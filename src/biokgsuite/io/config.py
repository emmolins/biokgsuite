"""Configuration loading for biokgsuite.

The repo keeps configuration in two files under ``config/``:

    config/kgs.yaml     knowledge-graph definitions
    config/params.yaml  analysis parameters

``load_config`` merges them into the single dict shape the notebooks expect
(``knowledge_graphs``, ``analysis_params``, ``_base_dir``), so notebook code
is unchanged from the pre-split layout.
"""
import os
from pathlib import Path
import yaml

_ROOT_MARKERS = ('config/kgs.yaml', 'config.yaml', 'pyproject.toml')


def find_repo_root(start: Path = None) -> Path:
    """Walk upward from *start* until a repo-root marker is found."""
    start = Path(start or Path.cwd()).resolve()
    for d in [start, *start.parents]:
        if any((d / m).exists() for m in _ROOT_MARKERS):
            return d
    raise FileNotFoundError(
        f'No repo root found above {start}. Expected one of {_ROOT_MARKERS}.')


# Back-compat alias: notebooks call find_config(...) and pass the result
# straight to load_config(...).
def find_config(start: Path = None) -> Path:
    return find_repo_root(start)


def load_config(root=None, *, trial_subset=True):
    """Load and merge the config files. *root* may be a repo root or a path
    inside it; a legacy single ``config.yaml`` is still honoured.

    ``trial_subset`` applies ``BKG_KGS`` or ``.trial_kgs`` (see
    ``_apply_kg_subset``). It is for the dimension panel, notebooks 00 to 07,
    whose cost is the whole panel. The grounding experiment and the set
    builder name their graphs explicitly (``params.grounding.kgs``, ``--kgs``)
    and pass ``trial_subset=False``, so a laptop trial file never removes a
    graph they need.
    """
    root = find_repo_root(root if root is None else Path(root))
    if (root / 'config' / 'kgs.yaml').exists():
        cfg = yaml.safe_load((root / 'config' / 'kgs.yaml').read_text()) or {}
        cfg.update(yaml.safe_load((root / 'config' / 'params.yaml').read_text()) or {})
    else:                                       # legacy layout
        cfg = yaml.safe_load((root / 'config.yaml').read_text()) or {}
    cfg['_base_dir'] = root
    _normalise_entity_types(cfg)
    if trial_subset:
        _apply_kg_subset(cfg)
    return cfg


def _apply_kg_subset(cfg: dict) -> None:
    """Restrict the panel for trial runs: ``BKG_KGS`` or a ``.trial_kgs`` file.

    Loading the full panel takes hours because of MATRIX (81 M edges). With
    ``BKG_KGS=primekg`` in the environment, or a gitignored ``.trial_kgs``
    file at the repo root naming the graphs (``primekg,hetionet``), every
    notebook and tool sees a config that declares only those graphs, so
    nothing else is loaded and no notebook needs editing. The file form is for
    a laptop: it is never committed, so a clone on the cluster runs the full
    panel. The full declared panel is kept under ``_full_panel`` so the panel
    guard can record that the run was a deliberate subset rather than a silent
    one.
    """
    raw = os.environ.get('BKG_KGS', '').strip()
    source = 'BKG_KGS'
    trial_file = Path(cfg['_base_dir']) / '.trial_kgs'
    if not raw and trial_file.exists():
        raw = ' '.join(l.split('#')[0] for l in trial_file.read_text().splitlines()).replace(' ', '').strip(',')
        source = '.trial_kgs'
    graphs = cfg.get('knowledge_graphs') or {}
    if not raw or not graphs:
        return
    wanted = [k.strip() for k in raw.split(',') if k.strip()]
    unknown = [k for k in wanted if k not in graphs]
    if unknown:
        raise KeyError(f'{source} names graphs not in config/kgs.yaml: {unknown}; '
                       f'known: {list(graphs)}')
    cfg['_full_panel'] = list(graphs)
    cfg['knowledge_graphs'] = {k: graphs[k] for k in graphs if k in wanted}
    print(f'{source}={raw}: TRIAL RUN on {len(wanted)} of {len(graphs)} graphs; '
          f'checkpoints will record the reduced panel'
          + ('  (delete .trial_kgs for the full panel)' if source == '.trial_kgs' else ''))


def _normalise_entity_types(cfg: dict) -> None:
    """Give every canonical entity type a TUPLE of the KG labels it covers.

    A graph may name one canonical type with several labels: BioKG carries both
    `Disease` (MeSH) and `GeneticDisorder` (OMIM), which are the same kind of
    object. The config accepts either a single label or a list, and every
    consumer downstream sees a tuple, so the one-label case needs no special
    handling anywhere else.
    """
    for spec in (cfg.get('knowledge_graphs') or {}).values():
        types = spec.get('entity_types')
        if not isinstance(types, dict):
            continue
        spec['entity_types'] = {
            canon: (tuple(labels) if isinstance(labels, (list, tuple)) else (labels,))
            for canon, labels in types.items()}


def require(params: dict, key: str, section: str = None):
    """Fetch a required parameter, raising if it is absent.

    Replaces ``params.get(key, <default>)``. A silent fallback means a run can
    use a value the config never mentions, which is how ``top_k``,
    ``train_ratio``, ``neg_strategies`` and ``min_pairs_per_area`` were
    applied for the whole project without appearing in config.yaml.
    """
    d = params.get(section, {}) if section else params
    if key not in d:
        where = f"analysis_params['{section}']" if section else 'analysis_params'
        raise KeyError(
            f"'{key}' is not defined under {where} in config/params.yaml. "
            'Add it rather than relying on a default.')
    return d[key]


class PartialPanelError(RuntimeError):
    """Raised when fewer knowledge graphs loaded than `config/kgs.yaml` declares."""


def assert_full_panel(loaded, config, allow_partial=False, notebook=None):
    """Refuse to score a panel that is missing graphs, unless asked explicitly.

    Why this exists
    ---------------
    Every dimension notebook wraps its load loop in a broad `except` so that one
    unreadable graph does not take down the run. The cost is that a run which
    found two of six graphs completes, prints a dimension score, and writes a
    checkpoint that looks exactly like a full one. Nothing downstream records
    which graphs were actually scored, so a partial run is indistinguishable
    from a complete one once it reaches the summary notebook or a figure.

    That is not a hypothetical: every notebook in this repository currently
    carries stored output from a run that loaded hetionet and biokg only, while
    `config/kgs.yaml` declares six.

    Cross-graph quantities make this worse than a missing row. Scores that
    intersect over the panel (shared entity types, shared relation types) and
    any normalisation against a panel maximum change value when a graph drops
    out, so the remaining graphs' numbers are not the numbers they would have
    had in a full run.

    Parameters
    ----------
    loaded : iterable of str
        Names of the graphs that actually loaded.
    config : dict
        The parsed `kgs.yaml`, or its `knowledge_graphs` mapping.
    allow_partial : bool
        Set True to run a deliberate subset. The run is still reported, so the
        reduced panel is a stated choice rather than a silent outcome.
    notebook : str, optional
        Name used in the message.

    Returns
    -------
    dict with `expected`, `loaded` and `missing`, for the checkpoint.
    """
    graphs = config.get("knowledge_graphs", config) if isinstance(config, dict) else config
    expected = list(graphs)
    got = list(loaded)
    missing = [n for n in expected if n not in got]
    report = {"expected_kgs": expected, "loaded_kgs": got, "missing_kgs": missing}
    # A BKG_KGS trial run is a subset by construction: say so in the record
    # and in the message, so the checkpoint can never pass for a full run.
    full = config.get("_full_panel") if isinstance(config, dict) else None
    if full and len(full) > len(expected):
        report["full_panel"] = list(full)
        report["subset_via"] = "BKG_KGS"

    where = f"{notebook}: " if notebook else ""
    if not missing:
        if "subset_via" in report:
            print(f"{where}trial subset complete, {len(got)} of {len(full)} declared graphs "
                  f"loaded (BKG_KGS); not a full-panel run")
        else:
            print(f"{where}panel complete, {len(got)}/{len(expected)} graphs loaded")
        return report

    msg = (f"{where}{len(got)}/{len(expected)} graphs loaded; missing "
           f"{missing}. Scores that intersect or normalise across the panel "
           f"are not comparable with a full run. Pass allow_partial=True to "
           f"score this subset deliberately.")
    if not allow_partial:
        raise PartialPanelError(msg)
    print(f"WARNING  {msg}")
    return report
