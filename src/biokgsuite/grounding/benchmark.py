"""The grounding benchmark: positives, per-graph flags, distractors and queries.

Positives are the post-cutoff approval set (``data/gold_standards/post_cutoff_approvals/``),
the same 128 pairs notebook 07 scores, so the grounding experiment and the
generalisation dimension share one benchmark and one provenance trail. The
per-graph coverage and leak flags come from the builder's status table rather
than being recomputed here, so the coverage-controlled analysis uses exactly
the checks the paper describes for notebook 07.

Two families of queries are drawn from the same pairs:

* ``main``    every pair, in every arm; the analysis then restricts post hoc
              using the flags carried on each row.
* ``common``  only the pairs evaluable in every KG arm, so the paired
              comparison between graphs holds evidence availability constant
              by construction.

Evaluable here is stricter than for notebook 07. The builder's ``evaluable``
says the graph holds the drug at the degree floor and the disease, carries no
direct edge, and predates the approval. A grounding dossier renders only gene
and pathway neighbours, so a drug can be covered (degree from drug-drug or
side-effect edges) and still reach the model with no evidence line; the
builder records that as ``grounded``. For the grounding experiment
``evaluable_<kg>`` = builder evaluable and grounded, and the three flags are
carried separately (``covered_<kg>``, ``grounded_<kg>``, ``evaluable_<kg>``)
so an analysis can condition on whichever it means.

Distractors are approved drugs from Open Targets that are not indicated for
the target disease or any of its MONDO descendants (``tools/build_distractor_pool.py``).
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from ..gold.ontology_classifier import OntologyClassifier

SET_DIR = Path('data/gold_standards/post_cutoff_approvals')
BRIDGE_DIR = Path('data/gold_standards/bridges')

# Open Targets drugType, collapsed to the modalities the bias audit reports on.
_MODALITY = {
    'Small molecule': 'small molecule', 'Antibody': 'antibody', 'Protein': 'protein',
    'Enzyme': 'protein', 'Oligonucleotide': 'oligonucleotide', 'Oligosaccharide': 'other',
    'Vaccine component': 'vaccine', 'Cell': 'cell or gene therapy', 'Gene': 'cell or gene therapy',
}


def load_benchmark(base, *, kgs, strict=False) -> pd.DataFrame:
    """One row per pair with the columns every stage reads.

    ``pair_id`` is the 1-based row number of ``post_cutoff_approvals.csv``,
    which the builder writes in a fixed sort order, so the id is stable across
    rebuilds as long as the set itself does not change (the set's SHA-256 is
    in its PROVENANCE.txt). ``evaluable_<kg>`` is the builder's
    ``evaluable_strict`` when ``strict`` is set, otherwise ``evaluable``, and
    in either case also requires ``grounded``; ``common`` is their conjunction
    over ``kgs``. ``n_kg_covered`` and ``n_kg_grounded`` count graphs per pair.
    ``shown_<kg>`` is looser than ``grounded_<kg>``: the graph holds the drug
    with at least one renderable gene or pathway neighbour, whatever its degree
    and whether or not the disease is in the graph. It is what the runner's
    dossier renders for the positive on that arm, so it is the flag
    ``pos_has_evidence`` must agree with; ``grounded`` additionally requires
    coverage and is the flag the evaluable set is built on.
    """
    base = Path(base)
    gold = pd.read_csv(base / SET_DIR / 'post_cutoff_approvals.csv')
    status = pd.read_csv(base / SET_DIR / 'post_cutoff_pair_status.csv')
    missing = set(kgs) - set(status.kg)
    if missing:
        raise ValueError(f'post_cutoff_pair_status.csv has no rows for {sorted(missing)}; '
                         'run tools/build_post_cutoff_approvals.py --kgs all first')
    bench = gold[['drugbank_id', 'mondo_id', 'drug_name', 'disease_name', 'doid', 'mesh_id',
                  'approval_date', 'agency', 'agencies', 'n_agencies', 'combination_type']].copy()
    bench.insert(0, 'pair_id', np.arange(1, len(bench) + 1))
    bench['approval_year'] = bench.approval_date.str[:4].astype(int)
    # The regulator that decided first. `agencies` is sorted alphabetically, so its
    # first entry is not the first decision; `agency` is the earliest-dated one.
    bench['first_agency'] = bench['agency']

    flag = 'evaluable_strict' if strict else 'evaluable'
    if 'grounded' not in status.columns:
        raise ValueError('post_cutoff_pair_status.csv has no grounded column; rebuild it with '
                         'tools/build_post_cutoff_approvals.py after re-extracting the views')
    for kg in kgs:
        s = status[status.kg == kg].set_index(['drugbank_id', 'mondo_id'])
        # The builder leaves grounded empty for a graph whose view carries no
        # neighbour count; for a grounding arm that is a stale view, not a False.
        g = s['grounded']
        if g.isna().any():
            raise ValueError(f'grounded is not computed for {kg} ({int(g.isna().sum())} empty rows); '
                             f'rerun tools/build_post_cutoff_approvals.py --extract-views {kg}, then --kgs all')
        s = s.assign(grounded=g.map({True: True, False: False, 'True': True, 'False': False}).astype(bool),
                     shown=pd.to_numeric(s['drug_n_mech'], errors='coerce').fillna(0) > 0)
        key = pd.MultiIndex.from_frame(bench[['drugbank_id', 'mondo_id']])
        for col, name in (('covered', f'covered_{kg}'), ('leaked', f'leaked_{kg}'),
                          ('grounded', f'grounded_{kg}'), ('shown', f'shown_{kg}'), (flag, f'evaluable_{kg}')):
            bench[name] = s[col].reindex(key).fillna(False).astype(bool).to_numpy()
        bench[f'evaluable_{kg}'] &= bench[f'grounded_{kg}']
    bench['n_kg_covered'] = bench[[f'covered_{kg}' for kg in kgs]].sum(axis=1)
    bench['n_kg_grounded'] = bench[[f'grounded_{kg}' for kg in kgs]].sum(axis=1)
    bench['n_kg_shown'] = bench[[f'shown_{kg}' for kg in kgs]].sum(axis=1)
    bench['common'] = bench[[f'evaluable_{kg}' for kg in kgs]].all(axis=1)

    clf = OntologyClassifier(base)
    areas = [clf.classify_mondo(m, n) for m, n in zip(bench.mondo_id, bench.disease_name)]
    bench['area'] = [a or 'unclassified' for a, _ in areas]
    bench['area_method'] = [m for _, m in areas]

    pool_path = base / BRIDGE_DIR / 'approved_drug_pool.csv'
    if pool_path.exists():
        pool = pd.read_csv(pool_path).set_index('drugbank_id')
        bench['modality'] = bench.drugbank_id.map(pool.drug_type).map(_MODALITY).fillna('unknown')
    else:
        bench['modality'] = 'unknown'
    return bench


def load_distractor_pool(base):
    """The approved-drug universe and, per MONDO id, the drugs indicated for it.

    The indicated sets are closed downward over the MONDO ``is_a`` hierarchy at
    query time (see ``build_queries``), not here, so this stays a plain lookup.
    """
    base = Path(base)
    universe = (pd.read_csv(base / BRIDGE_DIR / 'approved_drug_pool.csv')[['drugbank_id', 'drug_name']]
                .rename(columns={'drugbank_id': 'drug_id'}).reset_index(drop=True))
    ind = pd.read_csv(base / BRIDGE_DIR / 'approved_indications.csv')
    indicated: dict[str, set[str]] = defaultdict(set)
    for db, m in zip(ind.drugbank_id, ind.mondo_id):
        indicated[m].add(db)
    return universe, indicated


def mondo_children(base) -> dict[str, set[str]]:
    hier = pd.read_csv(Path(base) / BRIDGE_DIR / 'mondo_hierarchy.csv')
    children: dict[str, set[str]] = defaultdict(set)
    for c, p in zip(hier.child_id, hier.parent_id):
        children[p].add(c)
    return children


def descendants(mondo_id: str, children: dict[str, set[str]], max_depth: int = 15) -> set[str]:
    seen, frontier = {mondo_id}, {mondo_id}
    for _ in range(max_depth):
        frontier = {c for n in frontier for c in children.get(n, ()) if c not in seen}
        if not frontier:
            break
        seen |= frontier
    return seen


def query_order(bench: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    """Round-robin over therapeutic areas after a seeded permutation within each.

    Oncology holds a third of the pairs; a plain shuffle would front-load it,
    and any run stopped early by a rate limit would then be an oncology run.
    """
    areas = sorted(bench.area.unique())
    pools = {a: list(rng.permutation(bench.index[bench.area == a].to_numpy())) for a in areas}
    picks = []
    while any(pools.values()):
        for a in areas:
            if pools[a]:
                picks.append(pools[a].pop())
    return bench.loc[picks].reset_index(drop=True)


def build_queries(bench: pd.DataFrame, pool, *, family: str, pool_size: int, seed: int,
                  children: dict[str, set[str]] | None = None) -> list[dict]:
    """One query per pair: the positive plus ``pool_size - 1`` distractors.

    Distractors exclude every drug approved for the target disease or a
    descendant of it in MONDO, and every other post-cutoff positive for that
    disease, so no distractor is a correct answer under either the historical
    or the post-cutoff record. Candidate order is fixed here (positive first);
    the runner shuffles it per arm and per shuffle index.
    """
    if family == 'main':
        sel = bench
    elif family == 'common':
        sel = bench[bench.common]
    else:
        raise ValueError(f"unknown family {family!r}; expected 'main' or 'common'")
    universe, indicated = pool
    rng = np.random.default_rng(seed)
    ordered = query_order(sel, rng)
    pos_by_disease = bench.groupby('mondo_id').drugbank_id.agg(set).to_dict()
    queries = []
    for row in ordered.itertuples(index=False):
        excl = set(pos_by_disease.get(row.mondo_id, set())) | {row.drugbank_id}
        for m in (descendants(row.mondo_id, children) if children else {row.mondo_id}):
            excl |= indicated.get(m, set())
        cands = universe[~universe.drug_id.isin(excl)]
        k = pool_size - 1
        if len(cands) < k:
            raise ValueError(f'only {len(cands)} distractors available for {row.disease_name}')
        dis = cands.iloc[rng.choice(len(cands), size=k, replace=False)]
        candidates = [{'drug_id': row.drugbank_id, 'drug_name': row.drug_name, 'is_pos': True}]
        candidates += [{'drug_id': d.drug_id, 'drug_name': d.drug_name, 'is_pos': False} for d in dis.itertuples()]
        queries.append({'pair_id': int(row.pair_id), 'drugbank_id': row.drugbank_id, 'mondo_id': row.mondo_id,
                        'disease_name': row.disease_name, 'drug_name': row.drug_name, 'area': row.area,
                        'candidates': candidates})
    return queries


def selection_frame(bench: pd.DataFrame) -> pd.DataFrame:
    """The kept-versus-dropped frame for the selection-bias audit.

    ``kept`` marks the pairs the ``common`` family keeps. The audit asks whether
    that selection favours some kinds of pair (modality, area, year, how many
    regulators approved it), which is the assumption behind reading the
    evidence-controlled result as a fair comparison between graphs.
    """
    df = bench.copy()
    df['kept'] = df.common
    df['Modality'] = df.modality
    df['Therapeutic area'] = df.area
    df['Approval year'] = df.approval_year.astype(str)
    df['Regulators'] = np.where(df.n_agencies > 1, 'two or more', 'one')
    df['First approver'] = df.first_agency
    return df
