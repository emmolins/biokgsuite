"""Extract the evidence neighbourhoods shown to the LLMs, one small graph per
query and knowledge graph, for the results dashboard.

For every prospective pair (family main, seed 0) and each grounding KG, the
script records the disease's rendered neighbours and each candidate drug's
rendered neighbours, using the same KGDossiers rules and cap as the run
(results/tables/G02_run), so the graph drawn on the dashboard is exactly the
evidence the prompt carried.

Usage (from the repo root; one KG at a time keeps memory low):
    python tools/extract_evidence_graphs.py --kg drkg
    python tools/extract_evidence_graphs.py --kg biokg
    python tools/extract_evidence_graphs.py --kg primekg
    python tools/extract_evidence_graphs.py --merge

Writes results/tables/G02_run/evidence_graphs/<kg>.json, then
results/tables/G02_run/evidence_graphs.json on --merge.
"""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from biokgsuite.grounding import (KGDossiers, build_queries, load_distractor_pool,  # noqa: E402
                                  load_entrez_to_symbol, mondo_children)
from biokgsuite.io import load_config  # noqa: E402

OUT_DIR = ROOT / 'results' / 'tables' / 'G02_run' / 'evidence_graphs'


def load_params():
    import yaml
    p = yaml.safe_load(open(ROOT / 'config' / 'params.yaml'))
    g = p['analysis_params']['grounding'] if 'grounding' in p.get('analysis_params', {}) else p['grounding']
    return g


def load_bench():
    import pickle
    ckpt = pickle.load(open(ROOT / 'results' / 'checkpoints' / 'G01_benchmark.pkl', 'rb'))
    return ckpt['benchmark'], ckpt['pool_size']


def neighbour_block(d, idx, role, exclude=()):
    """{phrase: [name, ...]} after the cap, as the prompt rendered it."""
    out = {}
    for phrase, nbrs in d._neighbours(idx, role, exclude).items():
        names = []
        for v in nbrs[:d.cap]:
            nm = str(d.name_map.get(v, v))
            kind = d.canon.get(d.idx_to_type.get(v), 'other')
            names.append([nm, kind])
        out[phrase] = names
    return out


def extract(kg):
    config = load_config(ROOT, trial_subset=False)
    g = load_params()
    bench, pool_size = load_bench()
    pool = load_distractor_pool(ROOT)
    children = mondo_children(ROOT)
    queries = build_queries(bench, pool, family='main', pool_size=pool_size, seed=0, children=children)
    entrez = load_entrez_to_symbol(ROOT)
    d = KGDossiers(kg, config, ROOT, cap=g['cap'], anonymize_genes=g['anonymize_genes'],
                   bridge_mode=g['bridge_mode'], entrez_to_symbol=entrez)
    out = {}
    for q in queries:
        x_idx = d._resolve(d.to_disease_ids(q['mondo_id']), d.disease_types)
        disease = neighbour_block(d, x_idx, 'disease') if x_idx is not None else None
        cands = []
        for c in q['candidates']:
            d_idx = d._resolve(d.to_drug_ids(c['drug_id']), d.drug_types)
            block = neighbour_block(d, d_idx, 'drug', {x_idx} if x_idx is not None else set()) if d_idx is not None else None
            cands.append({'name': c['drug_name'], 'pos': bool(c['is_pos']), 'in_kg': d_idx is not None, 'ev': block})
        out[str(q['pair_id'])] = {'disease': disease, 'candidates': cands}
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / f'{kg}.json'
    json.dump(out, open(path, 'w'), separators=(',', ':'))
    n_ev = sum(1 for v in out.values() for c in v['candidates'] if c['ev'])
    print(f'{kg}: {len(out)} queries, {n_ev} candidate dossiers with evidence -> {path}')


def merge():
    merged = {}
    for p in sorted(OUT_DIR.glob('*.json')):
        merged[p.stem] = json.load(open(p))
    path = OUT_DIR.parent / 'evidence_graphs.json'
    json.dump(merged, open(path, 'w'), separators=(',', ':'))
    print(f'merged {sorted(merged)} -> {path} ({path.stat().st_size / 1e6:.1f} MB)')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--kg')
    ap.add_argument('--merge', action='store_true')
    a = ap.parse_args()
    if a.kg:
        extract(a.kg)
    if a.merge:
        merge()
