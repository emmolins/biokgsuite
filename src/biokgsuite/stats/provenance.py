"""Where does a graph say its edges came from?

Source diversity is only as good as the count that feeds it, and a number typed
in from a publication cannot be checked by anyone re-running this pipeline.
Worse, the publications count at different levels: Hetionet's 29 are primary
databases, MATRIX's first-level sources are five aggregator knowledge graphs,
and DRKG names Hetionet itself as one of its six. Comparing those numbers
directly compares bookkeeping conventions.

This module derives the count FROM THE DATA wherever a graph carries a
provenance field, and reports the documented figure BESIDE it rather than in
place of it. Two consequences worth stating in Methods:

  * The derived count is reproducible. Anyone with the release can recompute it.
  * Where derived and documented disagree, the disagreement is itself the
    finding, not something to reconcile away. DRKG documents "six existing
    databases ... and data collected from recent publications", so whether the
    literature is a seventh source is answered by whether the release tags it
    as one.

A graph with no provenance field cannot be measured this way at all, which is
exactly what edge traceability records. The two sub-metrics are therefore not
independent, and this module makes that visible instead of hiding it.

Used by ``biokgsuite_dimensions/03_trustworthiness``.
"""
from __future__ import annotations

import pandas as pd

NA = float('nan')

# How each graph exposes the source of an edge, if it does at all.
#   column          one or more columns naming the asserting source
#   relation_prefix the source is the prefix of the relation name (SRC::rel)
#   none            no per-edge provenance; only the publication can be cited
PROVENANCE = {
    'hetionet':   {'mode': 'column', 'cols': ('edge_source',), 'level': 'edge'},
    'biokg':      {'mode': 'none'},
    'primekg':    {'mode': 'column', 'cols': ('x_source', 'y_source'),
                   'level': 'node'},
    'drkg':       {'mode': 'relation_prefix', 'sep': '::', 'level': 'edge'},
    'openbilink': {'mode': 'column', 'cols': ('x_source', 'y_source'),
                   'level': 'edge'},
    'matrix':     {'mode': 'column', 'cols': ('edge_source',), 'level': 'edge'},
}

# Counts as stated by each resource's own documentation. Verified 2026-09-15;
# full citations and the counting rule are in data/SOURCE_COUNTS.md.
#
# COUNTING RULE. A source is a distinct upstream body of evidence that
# contributes edges, as enumerated by the resource itself. Literature-derived
# corpora count: they are independent bodies of evidence with their own failure
# modes, and excluding them would be inconsistent with EC-KG's 95, which
# reaches SemMedDB (literature-mined) through RTX-KG2. Ontologies and
# identifier-mapping resources that supply vocabulary rather than edges do not
# count, which is why OpenBioLink is 9 and not 16.
#
# These figures are NOT interchangeable with the derived counts above: a
# publication figure cannot be checked by re-running this pipeline, which is
# the whole reason `observed_sources` exists.
DOCUMENTED = {
    'hetionet': {
        'n': 29,
        'basis': '29 publicly available resources',
        'cite': 'Himmelstein et al. 2017, eLife 6:e26726; corroborated by '
                'Kaniewski et al. 2026, bioRxiv 2026.08.26.747253'},
    'primekg': {
        'n': 20,
        'basis': '20 high-quality resources',
        'cite': 'Chandak, Huang & Zitnik 2023, Sci Data 10:67; corroborated by '
                'Kaniewski et al. 2026'},
    'drkg': {
        'n': 7,
        # The README prose says "six existing databases ... and data collected
        # from recent publications", but its own statistics tables carry SEVEN
        # source columns (Drugbank, GNBR, Hetionet, STRING, IntAct, DGIdb,
        # Bibliography), and Bibliography contributes 9,464 entities and 84,756
        # triplets. The tables are the release; the sentence is prose about it.
        # Kaniewski et al. 2026 report 6 for DRKG, following the sentence.
        'basis': '6 databases + Bibliography, from the README statistics tables '
                 '(Bibliography: 9,464 entities, 84,756 triplets)',
        'cite': 'DRKG repository README, github.com/gnn4dr/DRKG (tables). '
                'NB Kaniewski et al. 2026 report 6, following the README prose'},
    'biokg': {
        'n': 13,
        'basis': '13 of the 18 databases in Table 1 are ticked in the '
                 '"In BioKG?" column',
        'cite': 'Walsh, Mohamed & Novacek 2020, CIKM ’20, Table 1'},
    'openbilink': {
        'n': 9,
        'basis': '9 distinct edge-contributing sources (STRING, GO, DisGeNet, '
                 'HPO, Bgee, STITCH, CTD, DrugCentral, SIDER); 16 names in the '
                 'table once ontology and id-mapping resources are included',
        'cite': 'OpenBioLink README source table, github.com/OpenBioLink/'
                'OpenBioLink; Breit et al. 2020, Bioinformatics 36:4097 '
                '(states no total)'},
    'matrix': {
        'n': 95,
        # EC-KG is the graph; MATRIX is Every Cure's platform. The paper
        # supersedes the docs page, which lists five first-level sources
        # including SPOKE and EmBiology: the paper integrates three (RTX-KG2,
        # ROBOKOP, PrimeKG) and names SPOKE as future work. 95 is the
        # de-duplicated union of their primary sources, not a sum.
        'basis': '95 aggregated knowledge sources across RTX-KG2, ROBOKOP and '
                 'PrimeKG; only 15/95 (16%) overlap between RTX-KG2 and ROBOKOP',
        'cite': 'Kaniewski et al. 2026, "Every Cure Knowledge Graph", bioRxiv '
                'doi:10.64898/2026.08.26.747253'},
}


def observed_sources(kg: pd.DataFrame, spec: dict) -> tuple:
    """(set of source names, fraction of edges carrying one) from the data."""
    mode = spec.get('mode', 'none')
    n = len(kg)
    if mode == 'none' or n == 0:
        return set(), NA

    if mode == 'column':
        cols = [c for c in spec['cols'] if c in kg.columns]
        if not cols:
            return set(), NA
        vals = pd.concat([kg[c] for c in cols], ignore_index=True)
        vals = vals.dropna().astype(str).str.strip()
        vals = vals[vals != '']
        # An edge counts as attributed if ANY of its source columns is filled.
        attributed = kg[cols].notna().any(axis=1).sum()
        return set(vals), attributed / n

    if mode == 'relation_prefix':
        sep = spec.get('sep', '::')
        rel = kg['relation'].astype(str)
        has = rel.str.contains(sep, regex=False)
        pref = rel[has].str.split(sep, n=1).str[0].str.strip()
        return set(pref), has.sum() / n

    raise ValueError(f'unknown provenance mode: {mode!r}')


def source_inventory(kgs: dict, provenance: dict = None,
                     documented: dict = None) -> pd.DataFrame:
    """Per KG: sources derived from the data, beside the documented count.

    `Derived` is NaN where the release carries no provenance field, which is a
    statement about the graph and not a zero.
    """
    provenance = PROVENANCE if provenance is None else provenance
    documented = DOCUMENTED if documented is None else documented
    rows = []
    for kg_name, data in kgs.items():
        spec = provenance.get(kg_name, {'mode': 'none'})
        names, attributed = observed_sources(data['kg'], spec)
        doc = documented.get(kg_name, {})
        derived = len(names) if spec.get('mode', 'none') != 'none' else None
        rows.append({
            'KG': kg_name,
            'Provenance': spec.get('mode', 'none'),
            'Level': spec.get('level', '-' if spec.get('mode') == 'none' else 'edge'),
            'Derived': derived if derived is not None else NA,
            'Documented': doc.get('n', NA),
            'Agree': (NA if derived is None or 'n' not in doc
                      else derived == doc['n']),
            'Edges attributed (%)': (round(100 * attributed, 1)
                                     if attributed == attributed else NA),
            'Basis for documented': doc.get('basis', ''),
            'Derived names': '; '.join(sorted(names)[:12]) + (
                f' (+{len(names) - 12} more)' if len(names) > 12 else ''),
        })
    return pd.DataFrame(rows)


def source_names(kgs: dict, provenance: dict = None) -> pd.DataFrame:
    """Every source name a graph attributes an edge to, with its edge count."""
    provenance = PROVENANCE if provenance is None else provenance
    rows = []
    for kg_name, data in kgs.items():
        spec = provenance.get(kg_name, {'mode': 'none'})
        if spec.get('mode', 'none') == 'none':
            continue
        kg = data['kg']
        if spec['mode'] == 'relation_prefix':
            sep = spec.get('sep', '::')
            rel = kg['relation'].astype(str)
            s = rel[rel.str.contains(sep, regex=False)].str.split(sep, n=1).str[0]
        else:
            cols = [c for c in spec['cols'] if c in kg.columns]
            s = pd.concat([kg[c] for c in cols], ignore_index=True).dropna().astype(str)
        for name, n in s.value_counts().items():
            rows.append({'KG': kg_name, 'Source': str(name).strip(), 'Edges': int(n)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- traceability
# What a single edge's provenance supports. Edge traceability is the MEAN of
# this over a graph's edges, which is the expected traceability of a randomly
# chosen edge.
#
# Scoring the tier alone asks only what KIND of provenance a graph offers and
# ignores how much of the graph carries it; scoring the attributed fraction
# alone does the reverse, and would rank PrimeKG level with DRKG because
# PrimeKG populates x_source/y_source on every edge -- fields that name where
# each ENDPOINT came from, not who asserted the edge. The product keeps both:
# the level says what one attributed edge is worth, the fraction says how many
# there are.
#
# It also handles a graph whose provenance is uneven, which Hetionet is: 20 of
# its 24 metaedges carry a per-edge source on 100% of edges and four carry none,
# because those four draw on several databases at once and no single source
# could be named. Hetionet therefore scores 0.69, not 0 and not 1.
LEVEL_SCORE = {
    'edge': 1.0,    # names the database that asserted the relationship
    'node': 0.5,    # names only where each endpoint's identifier came from
    'none': 0.0,    # no provenance field at all
}


def edge_traceability(kg, spec: dict) -> tuple:
    """(score, attributed fraction, level) for one graph.

    score = LEVEL_SCORE[level] x fraction of edges carrying provenance.
    """
    mode = spec.get('mode', 'none')
    if mode == 'none':
        return 0.0, 0.0, 'none'
    # A declared mechanism with no level is a config error, not a zero. Defaulting
    # it to 'none' silently scored DRKG 0 for traceability despite its per-edge
    # relation prefixes, which is the failure this raise exists to prevent.
    level = spec.get('level')
    if level not in LEVEL_SCORE or level == 'none':
        raise ValueError(
            f"provenance spec declares mode={mode!r} but level={level!r}; set "
            f"'level' to 'edge' or 'node' in biokgsuite.stats.provenance.PROVENANCE")
    _names, attributed = observed_sources(kg, spec)
    if attributed != attributed:                      # NaN: field absent
        return 0.0, NA, level
    return round(LEVEL_SCORE[level] * attributed, 4), attributed, level


def traceability_table(kgs: dict, provenance: dict = None) -> pd.DataFrame:
    """Per graph: provenance level, share of edges carrying it, and the score."""
    provenance = PROVENANCE if provenance is None else provenance
    rows = []
    for kg_name, data in kgs.items():
        spec = provenance.get(kg_name, {'mode': 'none'})
        score, frac, level = edge_traceability(data['kg'], spec)
        rows.append({
            'KG': kg_name,
            'Provenance level': level,
            'Level score': LEVEL_SCORE[level],
            'Edges attributed (%)': (round(100 * frac, 1) if frac == frac else NA),
            'Edge traceability': score,
        })
    return pd.DataFrame(rows)
