"""Coverage of a knowledge graph against task-derived reference sets.

    entity coverage   = |KG ∩ Ref| / |Ref|   per entity type
    relation coverage = the same, over (head, tail) pairs

Every KG identifier is resolved into the reference namespace first
(:mod:`biokgsuite.io.resolve`), so a graph is measured on the biology it holds
rather than on the accessions it happens to use. Coverage is always reported
beside resolvability, because one coverage number otherwise conflates three
different things: the entity is absent, the entity is present but unmappable,
and the entities are present but unlinked.

Which types are scored, and which are merely reported, is a decision the
calling notebook owns and passes in as ``eligible``. Nothing here decides it.

Used by ``biokgsuite_dimensions/01_coverage``.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, NamedTuple

import pandas as pd

NA = float('nan')

_BUILD_HINT = ('Run:  bash tools/fetch_opentargets.sh  &&  '
               'python tools/build_reference_sets.py')


# --------------------------------------------------------------- reference sets
@dataclass
class ReferenceSets:
    """Task-derived reference sets, keyed by canonical entity type / relation.

    Asymmetric by design: every disease is in scope, because repurposing exists
    to serve diseases with no approved drug, while only approved drugs are, because
    only an existing drug can be repurposed. See data/REFERENCE_SETS.md.
    """
    entities: dict = field(default_factory=dict)     # 'Drug' -> {chembl ids}
    relations: dict = field(default_factory=dict)    # 'drug_gene' -> {(a, b)}
    strata: dict = field(default_factory=dict)       # 'rare' -> {mondo ids}

    def summary(self) -> pd.DataFrame:
        rows = [{'Kind': 'entity', 'Set': k, 'N': len(v)}
                for k, v in self.entities.items()]
        rows += [{'Kind': 'disease stratum', 'Set': k, 'N': len(v)}
                 for k, v in self.strata.items()]
        rows += [{'Kind': 'relation', 'Set': k, 'N': len(v)}
                 for k, v in self.relations.items()]
        return pd.DataFrame(rows)


def load_reference_sets(directory) -> ReferenceSets:
    """Load every reference set written by ``tools/build_reference_sets.py``."""
    d = Path(directory)

    def ids(name: str) -> set:
        p = d / name
        if not p.exists():
            raise FileNotFoundError(f'{p} missing. {_BUILD_HINT}')
        return {ln.strip() for ln in p.open() if ln.strip()}

    def pairs(name: str) -> set:
        p = d / name
        if not p.exists():
            raise FileNotFoundError(f'{p} missing. {_BUILD_HINT}')
        with p.open() as fh:
            rd = csv.reader(fh, delimiter='\t')
            next(rd, None)                      # header
            return {(a, b) for a, b in rd}

    return ReferenceSets(
        entities={
            'Drug':         ids('entity_drug.txt'),        # ChEMBL, approved only
            'Disease':      ids('entity_disease.txt'),     # MONDO, all diseases
            'Gene/Protein': ids('entity_gene.txt'),        # Ensembl, drug targets
            'Pathway':      ids('entity_pathway.txt'),     # Reactome
        },
        strata={
            'rare':      ids('entity_disease_rare.txt'),
            'drugged':   ids('entity_disease_drugged.txt'),
            'undrugged': ids('entity_disease_undrugged.txt'),
        },
        relations={r: pairs(f'relation_{r}.tsv') for r in (
            'drug_gene', 'drug_pathway', 'disease_gene',
            'disease_pathway', 'disease_disease', 'drug_disease')},
    )


def reactome_pathway_names(directory) -> dict:
    """Reactome id -> name, for the pathway vocabulary diagnostic. {} if absent."""
    p = Path(directory) / 'entity_pathway_names.tsv'
    if not p.exists():
        return {}
    with p.open() as fh:
        rd = csv.reader(fh, delimiter='\t')
        next(rd, None)
        return {pid: nm for pid, nm in rd if nm.strip()}


# -------------------------------------------------------------------- coverage
def coverage(kg_set: set, ref_set: set) -> dict:
    """|KG ∩ Ref| / |Ref|, keeping both matched and unmatched for diagnostics."""
    overlap = kg_set & ref_set
    return {'kg_n': len(kg_set), 'gold_n': len(ref_set), 'overlap_n': len(overlap),
            'coverage_pct': round(100 * len(overlap) / len(ref_set), 2) if ref_set else NA,
            'matched': overlap, 'unmatched': ref_set - kg_set}


def _na(ref_set: set, kg_n: int = 0) -> dict:
    """A cell that was computed but is not evaluable for this KG.

    ``overlap_n`` and ``coverage_pct`` are NaN rather than 0, which is what the
    downstream code tests to tell "not asked" from "asked and found nothing".
    ``unmatched`` is left empty for the same reason: nothing was compared, so
    nothing can be reported as missing.
    """
    return {'kg_n': kg_n, 'gold_n': len(ref_set), 'overlap_n': NA,
            'coverage_pct': NA, 'matched': set(), 'unmatched': set()}


def entity_coverage(entity_sets: dict, refs: ReferenceSets,
                    eligible: Callable[[str, str], bool]):
    """Per KG x entity type. Returns (summary frame, membership lookup).

    An ineligible type still gets a row, as N/A: a type that silently vanishes
    from the table reads as an oversight, one shown as N/A states a decision.
    ``membership`` is left unset for it, so the gap diagnostics skip it.
    """
    rows, membership = [], {}
    for kg, sets in entity_sets.items():
        for etype, kg_set in sets.items():
            ref = refs.entities.get(etype)
            if ref is None:
                continue
            if eligible(kg, etype):
                res = coverage(kg_set, ref)
                membership.setdefault(etype, {})[kg] = res
            else:
                res = _na(ref, kg_n=len(kg_set))
            rows.append({'KG': kg, 'Entity Type': etype, 'KG Entities': res['kg_n'],
                         'Gold Standard': res['gold_n'], 'Overlap': res['overlap_n'],
                         'Coverage (%)': res['coverage_pct']})
    return pd.DataFrame(rows), membership


# ------------------------------------------------------------------- relations
class RelationSpec(NamedTuple):
    """One scored edge type, named by the entity types at its endpoints."""
    name: str
    src: str
    dst: str
    symmetric: bool = False

    @property
    def endpoints(self) -> tuple:
        return (self.src, self.dst)


def resolve_pairs(edges: pd.DataFrame, type_map: dict, spec: RelationSpec, resolver):
    """Edges between two canonical entity types, both endpoints resolved.

    Selected by endpoint entity TYPE rather than by relation name, which is how
    the grounding dossier gathers neighbours and which avoids double-counting a
    graph that reuses one predicate across two relation classes. Returns None
    when the KG does not type both endpoints, distinguishing "cannot be asked"
    from "asked and found nothing".
    """
    st, dt = type_map.get(spec.src), type_map.get(spec.dst)
    if not st or not dt:
        return None
    # `st`/`dt` are tuples: one canonical type may be spelled by several KG
    # labels (BioKG names diseases both `Disease` and `GeneticDisorder`).
    x_src, y_dst = edges['x_type'].isin(st), edges['y_type'].isin(dt)
    x_dst, y_src = edges['x_type'].isin(dt), edges['y_type'].isin(st)
    fwd = edges[x_src & y_dst]
    rev = edges[x_dst & y_src]
    raw = list(zip(fwd['x_id'], fwd['y_id'])) + list(zip(rev['y_id'], rev['x_id']))
    if not raw:
        return set()

    # Resolve each distinct id once, then expand: mappings are many-to-many.
    smap = {v: resolver.resolve(spec.src, [v])[0] for v in {a for a, _ in raw}}
    dmap = {v: resolver.resolve(spec.dst, [v])[0] for v in {b for _, b in raw}}
    out = set()
    for a, b in raw:
        for x in smap.get(a, ()):
            for y in dmap.get(b, ()):
                out.add(tuple(sorted((x, y))) if spec.symmetric else (x, y))
    return out


def relation_coverage(kgs: dict, type_maps: dict, specs: Iterable[RelationSpec],
                      refs: ReferenceSets, resolver,
                      eligible: Callable[[str, str], bool]):
    """Per KG x relation. Returns (results by relation, summary frame)."""
    results, rows = {}, []
    for spec in specs:
        ref = refs.relations[spec.name]
        ref_cmp = {tuple(sorted(p)) for p in ref} if spec.symmetric else set(ref)
        per_kg = {}
        for kg, data in kgs.items():
            # A relation is only as evaluable as the entity types it joins.
            if not all(eligible(kg, e) for e in spec.endpoints):
                per_kg[kg] = _na(ref_cmp)
                note = 'not evaluable (endpoint type not scored for this KG)'
            else:
                pairs = resolve_pairs(data['kg'], type_maps[kg], spec, resolver)
                if pairs is None:
                    per_kg[kg] = _na(ref_cmp)
                    note = 'not evaluable (KG does not type both endpoints)'
                else:
                    ov = pairs & ref_cmp
                    per_kg[kg] = {
                        'kg_n': len(pairs), 'gold_n': len(ref_cmp), 'overlap_n': len(ov),
                        'coverage_pct': (round(100 * len(ov) / len(ref_cmp), 2)
                                         if ref_cmp else NA),
                        'matched': ov, 'unmatched': ref_cmp - ov}
                    note = ''
            r = per_kg[kg]
            # Every KG x relation gets a row, evaluable or not. A cell that
            # vanishes from the table reads as an oversight; one shown as N/A
            # with its reason states a decision. Same rule as entity_coverage.
            rows.append({'KG': kg, 'Relation': spec.name, 'KG pairs': r['kg_n'],
                         'Overlap': r['overlap_n'], 'Gold': r['gold_n'],
                         'Coverage (%)': r['coverage_pct'], 'Note': note})
        results[spec.name] = per_kg
    return results, pd.DataFrame(rows)


def gap_decomposition(results: dict, entity_sets: dict,
                      specs: Iterable[RelationSpec]) -> pd.DataFrame:
    """Why is a reference pair missing? Three mutually exclusive reasons.

    Both endpoints present means the graph simply does not assert the edge; one
    or both absent means it could not assert it. A single coverage number hides
    that difference, and only the first is a statement about curation.
    """
    endpoints = {s.name: s.endpoints for s in specs}
    rows = []
    for rel, per_kg in results.items():
        src_c, dst_c = endpoints[rel]
        for kg, res in per_kg.items():
            if res['overlap_n'] != res['overlap_n']:       # NaN: not evaluable
                continue
            src = entity_sets.get(kg, {}).get(src_c, set())
            dst = entity_sets.get(kg, {}).get(dst_c, set())
            both = one = neither = 0
            for a, b in res['unmatched']:
                n = (a in src) + (b in dst)
                both, one, neither = (both + (n == 2), one + (n == 1), neither + (n == 0))
            rows.append({'KG': kg, 'Relation': rel, 'Missing': len(res['unmatched']),
                         'Both present (edge absent)': both,
                         'One entity absent': one, 'Both entities absent': neither})
    return pd.DataFrame(rows)


# --------------------------------------------------- pathway vocabulary check
def pathway_diagnostic(kgs: dict, type_maps: dict, reactome_kgs: set,
                       membership: dict, relation_results: dict,
                       reactome_names: dict) -> pd.DataFrame:
    """Why pathway coverage is reported but never scored.

    Compares each graph's pathway NAMES against Reactome's alongside what its
    pathway coverage would have been. Agreement on names beside zero or N/A
    coverage is a vocabulary mismatch, not missing biology, and scoring it would
    record an accession choice as a quality deficit.
    """
    ref_names = {n.strip().lower() for n in reactome_names.values() if n.strip()}

    def pct(d, kg):
        r = d.get(kg)
        return r['coverage_pct'] if r else NA

    rows = []
    for kg, data in kgs.items():
        nodes = data['nodes']
        ptypes = type_maps[kg].get('Pathway')
        sel = (nodes['type'].isin(ptypes) if ptypes
               else pd.Series(False, index=nodes.index))
        ids = nodes.loc[sel, 'id'].astype(str)
        names = ({s.strip().lower() for s in nodes.loc[sel, 'name'].astype(str)
                  if s.strip() and s.strip().lower() != 'nan'}
                 if 'name' in nodes.columns else set())
        n_nodes = len(ids)
        n_rhsa = int(ids.str.fullmatch(r'(?:\w+::)?R-HSA-\d+', na=False).sum())
        n_match = len(names & ref_names)
        rows.append({
            'KG': kg, 'Reactome declared': kg in reactome_kgs,
            'Pathway nodes': n_nodes, 'R-HSA accessions': n_rhsa,
            '% R-HSA': round(100 * n_rhsa / n_nodes, 1) if n_nodes else NA,
            'Distinct names': len(names), 'Name matches': n_match,
            '% name match': round(100 * n_match / len(names), 1) if names and ref_names else NA,
            'Entity (%)': pct(membership.get('Pathway', {}), kg),
            'drug→pathway (%)': pct(relation_results.get('drug_pathway', {}), kg),
            'disease→pathway (%)': pct(relation_results.get('disease_pathway', {}), kg)})
    return pd.DataFrame(rows)


# ------------------------------------------------- did we read the whole graph?
def schema_coverage(kgs: dict, type_maps: dict) -> pd.DataFrame:
    """How much of each graph's own schema the config actually maps.

    Every metric in this suite is computed over nodes whose KG type label maps
    to a canonical type. A label the config does not mention is invisible: its
    entities are never counted and its edges are never checked, and nothing
    else in the pipeline says so. That is not a finding about the graph, it is
    a gap in `config/kgs.yaml`, and it is the one failure that makes a
    cross-graph comparison unfair rather than merely incomplete.

    Reported per graph so an unmapped type has to be declared or dismissed on
    the record instead of vanishing. See `unmapped_types` for the detail.
    """
    rows = []
    for kg, data in kgs.items():
        labels = {lab for labs in type_maps.get(kg, {}).values() for lab in labs}
        counts = data['nodes']['type'].value_counts()
        mapped = counts[counts.index.isin(labels)]
        unmapped = counts[~counts.index.isin(labels)]
        rows.append({
            'KG': kg,
            'Node types': int(len(counts)),
            'Mapped': int(len(mapped)),
            'Unmapped': int(len(unmapped)),
            'Nodes': int(counts.sum()),
            'Nodes mapped (%)': round(100 * mapped.sum() / counts.sum(), 1) if counts.sum() else NA,
            'Unmapped types': '; '.join(f'{t} x{n:,}' for t, n in unmapped.items()) or '-',
        })
    return pd.DataFrame(rows)


def unmapped_types(kgs: dict, type_maps: dict, n_examples: int = 3) -> pd.DataFrame:
    """Each unmapped node type with example identifiers, to classify it."""
    rows = []
    for kg, data in kgs.items():
        labels = {lab for labs in type_maps.get(kg, {}).values() for lab in labs}
        nodes = data['nodes']
        for label, n in nodes['type'].value_counts().items():
            if label in labels:
                continue
            ids = nodes.loc[nodes['type'] == label, 'id'].astype(str)
            rows.append({'KG': kg, 'Unmapped type': label, 'Nodes': int(n),
                         'Example ids': ', '.join(ids.head(n_examples))})
    return pd.DataFrame(rows)
