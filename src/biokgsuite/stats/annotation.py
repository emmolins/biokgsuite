"""Annotation accuracy: does a graph's own labelling survive an outside check?

Both metrics here ask one question. A graph asserts that some node is a Drug,
or that some edge joins two Drugs. Is that assertion corroborated by what the
identifier actually is?

The arbiter is :class:`~biokgsuite.io.resolve.IdResolver`, the same namespace
dispatcher nb01 uses, so the judgement comes from outside the graph being
judged. That matters: the previous implementation checked each graph against a
schema hand-written by reading that graph, which no well-formed graph can fail,
and one graph's schema was parsed out of its own relation names.

Two checks are applied, and only the first is scored:

    check 1 (scored)    does the identifier's SHAPE fit a scheme the declared
                        type accepts?  DB12494 is DB plus five digits, so a
                        DrugBank accession, so a valid drug identifier.  The
                        shapes are declared in config/id_schemes.yaml.

    check 2 (reported)  does the identifier resolve to a live entry OF that
                        type in a reference registry?  Reported for Drug,
                        Disease and Gene/Protein, never folded into the score,
                        for the reason at REGISTRY_CHECKED below.

    entity validity         of the nodes whose shape the check can adjudicate,
                            the fraction whose shape fits their declared type

    relational consistency  of the edges both of whose endpoints are
                            adjudicable, the fraction where both shapes fit the
                            types the edge declares, so a `drug_drug` edge
                            really does join two things shaped like drugs

These are Table 1's names for the two metrics and the names nb00 reads out of
the checkpoint; keep them in step with the manuscript.

An identifier matching no declared shape at all is neither right nor wrong: it
is excluded from both denominators and counted separately.  So is an
identifier whose only matching schemes are type-agnostic (UMLS, MeSH, EFO name
concepts of every kind, so a CUI under Anatomy is no evidence of mislabelling).
Those exclusions are load-bearing and are why the score must be read beside
`nodes_judged_pct` and `edges_checkable_pct`, which say how much of the graph
the check could see.

Used by ``biokgsuite_dimensions/02_annotation_accuracy``.
"""
from __future__ import annotations

from collections import Counter
from typing import Iterable

import pandas as pd

# Every canonical type check 1 can adjudicate, which is every type declared in
# config/id_schemes.yaml. Not all are resolver-backed: Symptom, Pharmacologic
# Class, Phenotype, Exposure and Complex have declared shapes but no registry,
# so they are shape-checked only and `resolve_types` skips them deliberately.
TYPES = ('Drug', 'Disease', 'Gene/Protein', 'Pathway',
         'Biological Process', 'Molecular Function', 'Cellular Component',
         'Anatomy', 'Side Effect', 'Symptom', 'Pharmacologic Class',
         'Phenotype', 'Exposure', 'Complex', 'Taxon', 'ATC class', 'Cell type')

# The four verdicts check 1 can return.
CORRECT      = 'correct'        # the shape fits a scheme the declared type accepts
WRONG_SCHEME = 'wrong_scheme'   # the shape fits some OTHER type's own scheme: a mislabelling
AGNOSTIC     = 'agnostic'       # the shape fits only schemes that name every kind of concept
UNKNOWN      = 'unknown'        # the shape fits no declared scheme at all

# Scored: a graph is charged only for shapes that belong to a type it did not
# declare. UNKNOWN is excluded from the denominator rather than counted against
# the graph, because an unrecognised shape is at least as likely to be a gap in
# config/id_schemes.yaml as an error by the curator -- and an arbiter that
# charges a graph for the arbiter's own blind spots is not an outside check.
# The cost of that choice is that `nodes_judged_pct` must be read with the score.
FAILURES = (WRONG_SCHEME,)

# ONLY check 1 is scored.
#
# Check 2, resolution against a registry, is computed for every type the
# resolver covers and reported beside the score, never folded into it. The
# reason is that its complement does not mean what a reader would assume. When
# BioKG confirms 27.8% of its proteins, the other 72.2% are not wrong: they are
# reviewed UniProt accessions that HGNC, being gene-centric and human, does not
# carry. Hetionet's unconfirmed genes are withdrawn Entrez ids and its
# unconfirmed drugs are peptides and vitamins outside ChEMBL. Scoring that would
# charge a graph for the scope of our reference databases and for time passing.
REGISTRY_CHECKED = ('Drug', 'Disease', 'Gene/Protein')


def resolve_types(resolver, ids: Iterable[str]) -> dict:
    """Check 2's machinery: id -> the canonical types it resolves as.

    Built once per graph over the distinct identifiers, so the reporting below
    is a dictionary lookup rather than repeated resolution. Nothing SCORED
    depends on this: check 1 reads the identifier's shape and nothing else.
    """
    out = {}
    for raw in {str(i) for i in ids}:
        # Only the types the resolver holds a registry for. Symptom and
        # Pharmacologic Class are shape-checked only, so asking the resolver
        # about them is a KeyError, not a negative result.
        types = {t for t in TYPES
                 if t in resolver.RESOLVERS and resolver.resolve(t, [raw])[0]}
        # Pathway validity is broader than pathway coverage: a WikiPathways
        # accession is a real pathway that reaches no Reactome id. nb01 asks
        # whether an identifier reaches the reference set and must keep the
        # narrower test, so the two predicates are deliberately separate.
        if resolver.is_valid_pathway(raw):
            types.add('Pathway')
        out[raw] = frozenset(types)
    return out


def verdict(raw: str, declared: str, schemes) -> str:
    """Check 1: does the identifier's shape fit this entity type?

    ``DB12494`` is DB plus five digits, so a DrugBank accession, so a valid drug
    identifier. A shape that fits no scheme this type declares is a mislabelling
    when it fits some *other* type's shape, and simply unrecognised otherwise.
    """
    scheme, _ambiguous = schemes.match(raw, declared)
    if scheme is not None:
        return CORRECT
    if schemes.types_evidencing(raw):
        return WRONG_SCHEME
    return AGNOSTIC if schemes.types_matching(raw) else UNKNOWN


def _canonical_of(type_map: dict) -> dict:
    """KG's own type label -> canonical type, from config/kgs.yaml.

    One canonical type may be spelled by several KG labels, so the mapping is
    inverted rather than read directly.
    """
    return {label: canon for canon, labels in type_map.items() if canon in TYPES
            for label in labels}


def _verdict_sets(schemes, ids) -> tuple:
    """canonical -> (ids verdicted CORRECT, ids the check can adjudicate at all).

    Precomputed per graph so the edge pass is vectorised set membership.

    Deliberately NOT memoised. It is called twice per graph, by
    `relational_consistency` and `checkability_breakdown`, and recomputing costs
    a few seconds. Caching it costs ~11 sets of every identifier in the graph,
    held at module scope for the rest of the kernel's life, which is the wrong
    trade on a machine that also holds the raw KG tables.
    """
    correct, adjudicable = {}, {}
    for t in TYPES:
        ok, adj = set(), set()
        for raw in ids:
            v = verdict(raw, t, schemes)
            if v == CORRECT:
                ok.add(raw); adj.add(raw)
            elif v in FAILURES:
                adj.add(raw)
        correct[t], adjudicable[t] = ok, adj
    return correct, adjudicable


def entity_validity(kgs: dict, type_maps: dict, resolvers: dict,
                    schemes) -> pd.DataFrame:
    """Per KG x declared entity type, both checks reported separately.

        entity validity = valid / (valid + wrong scheme)

    Identifiers matching no declared shape, or only a type-agnostic one, are
    excluded from that denominator and given their own columns. They are a gap on our side as readily as an
    error on the curator's, so `Nodes` is carried beside the rate to show how
    much of the type the accuracy actually speaks for.
    """
    rows = []
    for kg, data in kgs.items():
        canon_of = _canonical_of(type_maps[kg])
        lookup = resolvers[kg]
        nodes = data['nodes']
        for label, canon in canon_of.items():
            sub = nodes[nodes['type'] == label]
            if sub.empty:
                continue
            counts, actually = Counter(), Counter()
            n_registry = 0                       # check 2, reported not scored
            for raw in sub['id'].astype(str):
                v = verdict(raw, canon, schemes)
                counts[v] += 1
                if canon in lookup.get(raw, ()):        # check 2, reported only
                    n_registry += 1
                if v == WRONG_SCHEME:
                    actually.update(schemes.types_evidencing(raw))
            judged = counts[CORRECT] + counts[WRONG_SCHEME]
            # Registry confirmation is reported for the three types drug
            # repurposing turns on, and for no others. Everything else is
            # shape-checked and left at that: a registry for side effects,
            # pathways or GO terms would measure our reference coverage, not
            # the graph, and is not what this dimension is for.
            covered = canon in REGISTRY_CHECKED
            rows.append({
                'KG': kg, 'Declared type': label, 'Canonical': canon,
                'Nodes': len(sub),
                'Valid': counts[CORRECT],
                'Wrong scheme': counts[WRONG_SCHEME],
                'Type-agnostic scheme': counts[AGNOSTIC],
                'Unrecognised shape': counts[UNKNOWN],
                'Validity (%)': round(100 * counts[CORRECT] / judged, 2) if judged else float('nan'),
                # Reported, never scored. The complement is NOT an error rate:
                # it is the share our reference registries do not carry.
                'Registry confirmed (%)': (round(100 * n_registry / len(sub), 1)
                                           if covered else float('nan')),
                'Actually': '; '.join(f'{t} x{n:,}' for t, n in actually.most_common(3)),
            })
    return pd.DataFrame(rows)


def pool_by_canonical(entity_acc: pd.DataFrame) -> pd.DataFrame:
    """Collapse per-declared-label rows to one row per KG x canonical type.

    `entity_validity` reports per DECLARED LABEL, which is the right
    granularity for diagnosis: it shows which of a graph's own labels is
    misfiled. Scoring and plotting need one row per CANONICAL type, because a
    graph that spells one type with two labels (BioKG files diseases as both
    `Disease` and `GeneticDisorder`) would otherwise weight that type twice in
    a macro-average, and cannot be pivoted at all.

    Validity is recomputed from the pooled counts rather than averaged from the
    per-label rates, so the pooled figure is the rate over all that type's
    nodes. Registry confirmation is node-weighted, and stays NaN for the types
    where it is not reported.
    """
    keys = ['KG', 'Canonical']
    g = (entity_acc.groupby(keys, as_index=False)
         .agg({'Nodes': 'sum', 'Valid': 'sum', 'Wrong scheme': 'sum',
               'Type-agnostic scheme': 'sum', 'Unrecognised shape': 'sum'}))
    judged = g['Valid'] + g['Wrong scheme']
    g['Validity (%)'] = (100 * g['Valid'] / judged.where(judged > 0)).round(2)

    w = entity_acc.assign(_c=entity_acc['Registry confirmed (%)'] * entity_acc['Nodes'])
    agg = w.groupby(keys, as_index=False).agg(
        _c=('_c', 'sum'), _n=('Nodes', 'sum'),
        _any=('Registry confirmed (%)', lambda s: s.notna().any()))
    g = g.merge(agg, on=keys, how='left')
    g['Registry confirmed (%)'] = (g['_c'] / g['_n'].where(g['_n'] > 0)).where(g['_any']).round(1)
    return g.drop(columns=['_c', '_n', '_any'])


def identifier_ambiguity(resolvers: dict) -> pd.DataFrame:
    """How often one identifier resolves into more than one canonical type.

    The assumption the accuracy metrics rest on. An identifier that resolves as
    both a gene and a disease cannot witness against a mislabelling, because
    whichever type the graph declares will be corroborated. A shape list that
    admits a bare integer as seven entity types at once fails for exactly that
    reason.
    """
    rows = []
    for kg, lookup in resolvers.items():
        resolved = [v for v in lookup.values() if v]
        multi = [v for v in resolved if len(v) > 1]
        pairs = Counter(tuple(sorted(v)) for v in multi)
        rows.append({
            'KG': kg, 'Identifiers': len(lookup), 'Resolvable': len(resolved),
            'Ambiguous': len(multi),
            'Ambiguous (%)': round(100 * len(multi) / len(resolved), 2) if resolved else float('nan'),
            'Most common overlap': ' + '.join(pairs.most_common(1)[0][0]) if pairs else '',
        })
    return pd.DataFrame(rows)


# Registry each bare number is read against, by the type the graph declares.
# A bare number carries no namespace, so the declared type is the only thing
# that says what it is; this table is that decision, made once and reported.
BARE_READ_AS = {'Disease': 'MONDO', 'Gene/Protein': 'NCBI Entrez',
                'Biological Process': 'GO', 'Molecular Function': 'GO',
                'Cellular Component': 'GO', 'Anatomy': 'UBERON', 'Phenotype': 'HPO',
                'Drug': 'PubChem CID', 'Taxon': 'NCBI Taxonomy', 'Side Effect': 'MedDRA'}


def bare_number_types(kgs: dict, type_maps: dict) -> pd.DataFrame:
    """Every declared type a graph keys by bare numbers, and what they are read as.

    Flagged rather than silently accepted. Unprefixed integers are the one
    identifier shape that fits several registries at once (a bare 8150 is a
    valid MONDO, Entrez, GO, UBERON and HPO number), so the metrics can only
    read them by trusting the declared type. That trust is checked by
    ``tools/audit_identifiers.py`` (registry membership and name agreement)
    and is stated here so a reader of the tables knows it was made.
    """
    rows = []
    for kg, data in kgs.items():
        canon_of = _canonical_of(type_maps[kg])
        nodes = data['nodes']
        for label, canon in canon_of.items():
            ids = nodes.loc[nodes['type'] == label, 'id'].astype(str)
            if ids.empty:
                continue
            stripped = ids.str.split('::').str[-1].str.strip()
            bare = stripped.str.fullmatch(r'\d+')
            grouped = stripped.str.fullmatch(r'\d+(_\d+)+')
            n_bare, n_grp = int(bare.sum()), int(grouped.sum())
            if not n_bare and not n_grp:
                continue
            rows.append({'KG': kg, 'Declared type': label, 'Canonical': canon,
                         'Nodes': len(ids), 'Bare numbers': n_bare, 'Grouped numbers': n_grp,
                         'Bare (%)': round(100 * (n_bare + n_grp) / len(ids), 1),
                         'Read as': BARE_READ_AS.get(canon, 'shape only, no registry')})
    return pd.DataFrame(rows)


def _endpoint_masks(ids: pd.Series, canon: pd.Series, ok_ids: dict, adj_ids: dict):
    """(conforms, adjudicable) for one endpoint column, as boolean Series.

    An endpoint whose KG type is unmapped has NaN canonical type, matches no
    branch, and stays False in both: not adjudicable, so its edge is excluded.
    """
    ok = pd.Series(False, index=ids.index)
    adj = pd.Series(False, index=ids.index)
    for t in TYPES:
        m = canon == t
        if not m.any():
            continue
        sel = ids[m]
        ok[m] = sel.isin(ok_ids[t])
        adj[m] = sel.isin(adj_ids[t])
    return ok, adj


def _resolved_as(lookup: dict, raw) -> str:
    """What an identifier actually resolved to, for the mismatch column."""
    return '/'.join(sorted(lookup.get(raw, ()))) or 'unresolved'


def relational_consistency(kgs: dict, type_maps: dict, resolvers: dict,
                           schemes) -> pd.DataFrame:
    """Per KG x relation: does the edge join the entity types it claims to?

    Both endpoints face the same two checks the nodes face, so a relation named
    for a drug-drug interaction is confirmed to join two things that really are
    drugs. Nothing is asserted about what a relation *ought* to join: the edge's
    own ``x_type``/``y_type`` are the claim and the identifiers are the witness,
    so every graph is held to one rule and none can score well by having had a
    permissive schema written for it.
    """
    rows = []
    for kg, data in kgs.items():
        canon_of = _canonical_of(type_maps[kg])
        lookup = resolvers[kg]
        ok_ids, adj_ids = _verdict_sets(schemes, lookup.keys())

        e = data['kg']
        cx, cy = e['x_type'].map(canon_of), e['y_type'].map(canon_of)
        ok_x, adj_x = _endpoint_masks(e['x_id'], cx, ok_ids, adj_ids)
        ok_y, adj_y = _endpoint_masks(e['y_id'], cy, ok_ids, adj_ids)

        checkable = adj_x & adj_y
        conforming = checkable & ok_x & ok_y

        for rel, idx in e.groupby('relation', sort=False).indices.items():
            chk, conf = checkable.values[idx], conforming.values[idx]
            n_check, n_conf = int(chk.sum()), int(conf.sum())
            # Only the violating rows are materialised. Attaching these columns
            # to the full edge table instead would copy every id string in the
            # graph, which is the difference between running and being
            # OOM-killed on a 2M-edge graph.
            bad = idx[chk & ~conf]
            reasons = Counter()
            for i in bad:
                parts = []
                if not ok_x.values[i]:
                    parts.append(f'{cx.values[i]}<-{_resolved_as(lookup, e["x_id"].values[i])}')
                if not ok_y.values[i]:
                    parts.append(f'{cy.values[i]}<-{_resolved_as(lookup, e["y_id"].values[i])}')
                reasons[', '.join(parts)] += 1
            n_edges = len(idx)
            rows.append({
                'KG': kg, 'Relation': rel, 'Edges': n_edges, 'Checkable': n_check,
                'Conforming': n_conf, 'Violations': n_check - n_conf,
                'Consistency (%)': round(100 * n_conf / n_check, 2) if n_check else float('nan'),
                'Checkable (%)': round(100 * n_check / n_edges, 1) if n_edges else float('nan'),
                'Top mismatch': reasons.most_common(1)[0][0] if reasons else '',
            })
    return pd.DataFrame(rows)


def checkability_breakdown(kgs: dict, type_maps: dict, resolvers: dict,
                           schemes) -> pd.DataFrame:
    """Why an edge could not be checked, so the coverage figure is explained.

    Out-of-scope endpoints and unadjudicable identifiers are different failures
    with different remedies, and a single "checkable %" hides which dominates.
    """
    rows = []
    for kg, data in kgs.items():
        canon_of = _canonical_of(type_maps[kg])
        _, adj_ids = _verdict_sets(schemes, resolvers[kg].keys())

        e = data['kg']
        cx, cy = e['x_type'].map(canon_of), e['y_type'].map(canon_of)
        in_scope = cx.notna() & cy.notna()
        _, adj_x = _endpoint_masks(e['x_id'], cx, adj_ids, adj_ids)
        _, adj_y = _endpoint_masks(e['y_id'], cy, adj_ids, adj_ids)
        both = adj_x & adj_y
        n = len(e)

        # Which KG type labels cost the most edges. Counted per endpoint, so an
        # edge with two unmapped endpoints is charged to both.
        oos = ~in_scope
        cost = Counter()
        for col, canon in (('x_type', cx), ('y_type', cy)):
            cost.update(e.loc[oos & canon.isna(), col].value_counts().to_dict())

        rows.append({
            'KG': kg, 'Edges': n,
            'Type out of scope': int(oos.sum()),
            'Out of scope (%)': round(100 * oos.sum() / n, 1) if n else float('nan'),
            'Endpoint unjudgeable': int((in_scope & ~both).sum()),
            'Unjudgeable (%)': round(100 * (in_scope & ~both).sum() / n, 1) if n else float('nan'),
            'Checkable': int(both.sum()),
            'Checkable (%)': round(100 * both.sum() / n, 1) if n else float('nan'),
            'Costliest out-of-scope types': '; '.join(f'{t} x{c:,}' for t, c in cost.most_common(4)),
        })
    return pd.DataFrame(rows)


def node_index_uniqueness(kgs: dict) -> pd.DataFrame:
    """Each internal index must name exactly one entity.

    A collision silently merges two entities in any matrix- or embedding-based
    model downstream, so this is an integrity check rather than a quality score.
    """
    rows = []
    for kg, data in kgs.items():
        nodes = data['nodes']
        n, uniq = len(nodes), nodes['idx'].nunique()
        rows.append({'KG': kg, 'Nodes': n, 'Unique indices': uniq,
                     'Collisions': n - uniq,
                     'Uniqueness (%)': round(100 * uniq / n, 2) if n else float('nan')})
    return pd.DataFrame(rows)


def _bare(series: pd.Series) -> pd.Series:
    return series.astype(str).str.split('::', n=1).str[-1]


def relation_exclusivity(kgs: dict, exclusive_pairs: dict) -> pd.DataFrame:
    """Contradictory relations sharing an entity pair.

    Diagnostic only. A graph that declares no mutually exclusive relations would
    score a vacuous 100%, which would reward schema minimalism over curation.
    """
    rows = []
    for kg, pairs in exclusive_pairs.items():
        if kg not in kgs:
            continue
        edges = kgs[kg]['kg']

        # Each relation's pair set is built once per graph and released with the
        # graph. Building it per conflict instead rebuilt the largest relations
        # several times over, which on a 2M-edge graph is the peak that matters.
        needed = {r for _, a, b in pairs for r in (a, b)}
        pair_sets = {}
        for rel in needed:
            sub = edges.loc[edges['relation'] == rel, ['x_id', 'y_id']]
            pair_sets[rel] = {(a, b) if a <= b else (b, a)
                              for a, b in zip(_bare(sub['x_id']), _bare(sub['y_id']))}

        for label, rel_a, rel_b in pairs:
            A, B = pair_sets[rel_a], pair_sets[rel_b]
            shared, union = len(A & B), len(A | B)
            rows.append({'KG': kg, 'Conflict': label, '|A|': len(A), '|B|': len(B),
                         'Shared pairs': shared,
                         'Conflict rate (%)': round(100 * shared / union, 3) if union else float('nan')})
        pair_sets.clear()
    return pd.DataFrame(rows)
