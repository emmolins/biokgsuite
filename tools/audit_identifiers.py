#!/usr/bin/env python3
"""Audit every identifier namespace a graph uses, per declared entity type.

Why this exists
---------------
Three of the six graphs write some identifiers as bare numbers with no
namespace prefix: PrimeKG writes MONDO, Entrez, GO, UBERON and HPO ids that
way; Hetionet and DRKG write Entrez genes as ``Gene::1956``; MATRIX's loader
emits bare MONDO and Entrez numbers after bridging. A bare number is read
according to the type the graph declares for the node, and that reading is a
claim that has to be checked, not assumed: every metric downstream counts on
'8150' under biological_process being GO:0008150.

This script makes the check explicit and reproducible. For every graph,
declared type and identifier shape it reports the count, and for shapes the
resolver holds a registry for (MONDO, HGNC/Entrez, GO, UBERON, HPO, ChEMBL
via DrugBank and PubChem, Reactome) the share confirmed as a live identifier,
with named examples. For bare numbers it adds the decisive test: the share
whose graph name equals the registry name of the term the number is read as
(existence proves the number is valid; a matching name proves it denotes the
same thing; the remainder is labels renamed since the graph was built). For
PrimeKG it also tabulates the graph's own ``x_source`` column, which states
the namespace of every bare number.

Usage (from the repo root; heavy graphs need the memory their notebooks need):

    python tools/audit_identifiers.py                    # every graph in config
    python tools/audit_identifiers.py --kgs primekg,hetionet
    BKG_KGS=biokg python tools/audit_identifiers.py      # the trial-run switch works too

Writes results/tables/id_audit/identifier_namespaces.csv and prints a summary.
"""
import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from biokgsuite.io import load_config, load_kg              # noqa: E402
from biokgsuite.io.resolve import IdResolver, split_curie  # noqa: E402

OUT_DIR = ROOT / 'results' / 'tables' / 'id_audit'


def shape_of(resolver, raw):
    """Namespace of a prefixed id, else the shape of the unprefixed one."""
    s = str(raw).split('::')[-1].strip()
    ns, _ = split_curie(s)
    if resolver._GROUPED_MONDO.match(s):
        return 'grouped integers'
    label = resolver._label(ns, s)
    if label in ('Reactome', 'WikiPathways', 'PathwayCommons'):
        ns = None
    if ns:
        return f'{ns}:'
    label = resolver._label(None, s)
    return {'integer': 'bare integer', 'UniProt': 'UniProt accession',
            'DrugBank': 'DrugBank accession', 'MeSH': 'MeSH id', 'UMLS CUI': 'UMLS CUI',
            'Reactome': 'Reactome stable id', 'WikiPathways': 'WikiPathways accession',
            'PathwayCommons': 'PathwayCommons accession', 'bare': 'other'}[label]


def _obo_names(path, prefix):
    """term id -> name, from an OBO file."""
    names, cur = {}, None
    for line in open(path, encoding='utf-8', errors='replace'):
        if line.startswith('id: ' + prefix):
            cur = line[4:].strip()
        elif line.startswith('[Term]'):
            cur = None
        elif cur and line.startswith('name: '):
            names[cur] = line[6:].strip()
    return names


def registry_names(crosswalks):
    """canonical type -> (id prefix, {id: name}); the names a bare number
    should carry if it is being read as the right term."""
    cw = Path(crosswalks)
    out = {}
    for canon, fname, prefix in (('Biological Process', 'go.obo', 'GO:'),
                                 ('Molecular Function', 'go.obo', 'GO:'),
                                 ('Cellular Component', 'go.obo', 'GO:'),
                                 ('Anatomy', 'uberon.obo', 'UBERON:'),
                                 ('Phenotype', 'hp.obo', 'HP:'),
                                 ('Disease', 'mondo.obo', 'MONDO:')):
        if (cw / fname).exists():
            out[canon] = (prefix, _obo_names(cw / fname, prefix))
    hg = cw / 'hgnc_complete_set.txt'
    if hg.exists():
        sym = {r['entrez_id']: r['symbol'] for r in csv.DictReader(open(hg), delimiter='\t') if r.get('entrez_id')}
        out['Gene/Protein'] = (None, sym)
    return out


def _norm(s):
    return re.sub(r'[^a-z0-9]+', ' ', str(s).lower()).replace('obsolete ', '').strip()


def name_agreement(canon, members, names_by_id, registry):
    """Share of bare-number nodes whose graph name equals the registry name of
    the term the number is read as. Existence in a registry shows the number
    is valid; a matching name shows it denotes the same thing. Disagreements
    are almost all labels renamed since the graph was built."""
    if canon not in registry:
        return None
    prefix, reg = registry[canon]
    same = total = 0
    for raw in members:
        s = str(raw).split('::')[-1].strip()
        if not s.isdigit():
            continue
        key = f'{prefix}{int(s):07d}' if prefix else s
        rn = reg.get(key)
        if rn is None:
            continue
        total += 1
        gn = names_by_id.get(raw)
        if gn is not None and _norm(rn) == _norm(gn):
            same += 1
    return round(100 * same / total, 1) if total else None


def primekg_sources(path):
    """PrimeKG's own statement of each node's namespace: (type, id) -> source."""
    csv.field_size_limit(10 ** 9)
    src = {}
    with open(path, newline='') as fh:
        for row in csv.DictReader(fh):
            for side in ('x', 'y'):
                src[(row[f'{side}_type'], row[f'{side}_id'])] = row[f'{side}_source']
    return src


def audit(kg_name, config, resolver, registry):
    kg, nodes = load_kg(kg_name, config)
    del kg
    spec = config['knowledge_graphs'][kg_name]
    canon_of = {label: canon for canon, labels in spec['entity_types'].items()
                for label in (labels if isinstance(labels, (list, tuple)) else (labels,))}
    sources = primekg_sources(ROOT / spec['path']) if kg_name == 'primekg' else None
    names = dict(zip(zip(nodes['type'], nodes['id'].astype(str)), nodes['name'].astype(str)))
    rows = []
    for label, canon in canon_of.items():
        ids = nodes.loc[nodes['type'] == label, 'id'].astype(str).unique()
        if not len(ids):
            continue
        by_shape = defaultdict(list)
        for raw in ids:
            by_shape[shape_of(resolver, raw)].append(raw)
        checkable = canon in resolver.RESOLVERS
        for shape, members in sorted(by_shape.items(), key=lambda kv: -len(kv[1])):
            confirmed, ok_ex, bad_ex = 0, [], []
            if checkable:
                for raw in members:
                    hit = resolver.resolve(canon, [raw])[0]
                    if hit:
                        confirmed += 1
                        if len(ok_ex) < 3:
                            ok_ex.append(f'{raw} ({names.get((label, raw), "?")[:40]})')
                    elif len(bad_ex) < 3:
                        bad_ex.append(f'{raw} ({names.get((label, raw), "?")[:40]})')
            agree = ''
            if shape == 'bare integer':
                v = name_agreement(canon, members, {m: names.get((label, m)) for m in members}, registry)
                agree = v if v is not None else ''
            graph_says = ''
            if sources is not None:
                said = Counter(sources.get((label, raw), '?') for raw in members)
                graph_says = '; '.join(f'{k} x{v:,}' for k, v in said.most_common(3))
            rows.append({
                'KG': kg_name, 'Declared type': label, 'Canonical': canon,
                'Shape': shape, 'Ids': len(members),
                'Registry': resolver.RESOLVERS.get(canon, '').replace('_', ' ') if checkable else 'none (shape only)',
                'Confirmed': confirmed if checkable else '',
                'Confirmed (%)': round(100 * confirmed / len(members), 1) if checkable else '',
                'Name agrees (%)': agree,
                'Graph says': graph_says,
                'Confirmed examples': '; '.join(ok_ex),
                'Unconfirmed examples': '; '.join(bad_ex),
            })
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--kgs', help='comma-separated graph names (default: every graph the config declares)')
    ap.add_argument('--out', default=str(OUT_DIR))
    a = ap.parse_args()

    config = load_config(ROOT)
    kgs = [k.strip() for k in a.kgs.split(',')] if a.kgs else list(config['knowledge_graphs'])
    unknown = [k for k in kgs if k not in config['knowledge_graphs']]
    if unknown:
        ap.error(f'unknown graph(s) {unknown}; known: {list(config["knowledge_graphs"])}')

    resolver = IdResolver(ROOT / 'data' / 'crosswalks', verbose=True, aux_registries=True)
    registry = registry_names(ROOT / 'data' / 'crosswalks')
    rows = []
    for kg in kgs:
        print(f'\n== {kg}')
        rows += audit(kg, config, resolver, registry)

    df = pd.DataFrame(rows)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / 'identifier_namespaces.csv', index=False)

    # The summary a reader wants: every bare-number shape, with what the graph
    # says it is and how much of it the registry confirms.
    bare = df[df['Shape'].isin(['bare integer', 'grouped integers'])]
    with pd.option_context('display.width', 200, 'display.max_colwidth', 60):
        print('\nBare-number identifiers, by graph and declared type:')
        print(bare[['KG', 'Declared type', 'Canonical', 'Shape', 'Ids', 'Registry',
                    'Confirmed (%)', 'Name agrees (%)', 'Graph says']].to_string(index=False))
    low = bare[(bare['Confirmed (%)'] != '') & (bare['Confirmed (%)'].astype(float) < 90)]
    if len(low):
        print('\nBelow 90% confirmed (read the examples in the CSV before trusting these):')
        print(low[['KG', 'Declared type', 'Shape', 'Ids', 'Confirmed (%)', 'Unconfirmed examples']]
              .to_string(index=False))
    print(f'\nwrote {out / "identifier_namespaces.csv"}  ({len(df)} rows)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
