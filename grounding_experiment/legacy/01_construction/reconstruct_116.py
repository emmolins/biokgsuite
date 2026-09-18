#!/usr/bin/env python3
"""
Faithful, validated reconstruction of the 253 -> 116 evaluable benchmark.

Rule (validated to 100% per-pair vs results .../coverage_annotation_v4.csv on
each KG):
  present_in_KG  = DrugBank_ID resolves to a Drug-type node
                 AND disease ID (KG scheme) resolves to a Disease-type node
  contaminated   = the drug->disease edge already exists in that KG
  evaluable_in_KG = present AND NOT contaminated
  kept           = evaluable in >=1 of PrimeKG / DRKG / BioKG

Resolution uses the pipeline's own loader (src.loading.load_kg) + config
entity_types, with alias expansion (suffix after '::', scheme-prefix strip,
leading-zero strip, and disease_id_separator splits).

Usage:  python scripts/reconstruct_116.py
Out:    results/tables/reconstruct_116_audit.csv
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pandas as pd
from collections import defaultdict
from src.loading import load_config, load_kg

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GS = os.path.join(ROOT, 'data', 'gold_standards')
cfg = load_config(ROOT)
pid = pd.read_csv(os.path.join(GS, 'legacy', 'pairs_with_ids_v2.csv'))
cov = pd.read_csv(os.path.join(GS, 'legacy', 'coverage_annotation_v4.csv'))   # 116 survivors, for validation
KGS = ['primekg', 'drkg', 'biokg']

def node_aliases(raw, sep=None):
    raw = str(raw); out = {raw, raw.lstrip('0')}
    tail = raw.split('::', 1)[-1] if '::' in raw else raw
    out.add(tail)
    for pre in ('MESH:', 'MeSH:', 'DOID:', 'MONDO:', 'UMLS:', 'OMIM:'):
        if tail.startswith(pre):
            out.add(tail[len(pre):]); out.add(tail[len(pre):].lstrip('0'))
    if sep and sep in raw:
        for p in raw.split(sep):
            out.add(p); out.add(p.lstrip('0'))
    return out

def cand_disease_keys(row, scheme):
    out = set()
    def add(x):
        if pd.isna(x): return
        s = str(x); out.add(s)
        t = s.split(':')[-1]; out.add(t); out.add(t.lstrip('0'))
    if scheme == 'mondo':   add(row.get('MONDO_ID'))
    elif scheme == 'mesh':  add(row.get('MeSH_ID'))
    elif scheme == 'doid_mesh': add(row.get('DOID_ID')); add(row.get('MeSH_ID'))
    return out

for kg in KGS:
    et = cfg['knowledge_graphs'][kg]['entity_types']
    sep = cfg['knowledge_graphs'][kg].get('disease_id_separator')
    scheme = cfg['knowledge_graphs'][kg]['disease_id_scheme']
    drug_t, dis_t = et['Drug'], et['Disease']
    kg_df, nodes = load_kg(kg, cfg)

    drug_set, dis_set = set(), set()
    for v in nodes[nodes.type == drug_t]['id']: drug_set |= node_aliases(v)
    for v in nodes[nodes.type == dis_t]['id']:  dis_set |= node_aliases(v, sep)

    # drug->disease edges for contamination: map drug-alias -> set(disease node-alias)
    dd = defaultdict(set)
    m = ((kg_df.x_type == drug_t) & (kg_df.y_type == dis_t))
    for x, y in zip(kg_df.loc[m, 'x_id'], kg_df.loc[m, 'y_id']):
        dd[str(x)].add(str(y)); dd[str(x).lstrip('0')].add(str(y))
    m2 = ((kg_df.x_type == dis_t) & (kg_df.y_type == drug_t))
    for x, y in zip(kg_df.loc[m2, 'x_id'], kg_df.loc[m2, 'y_id']):
        dd[str(y)].add(str(x)); dd[str(y).lstrip('0')].add(str(x))
    # disease-node alias -> set of canonical disease ids that contain it (for contamination match)
    del kg_df

    def present(row):
        drug = str(row['DrugBank_ID']) in drug_set
        dk = cand_disease_keys(row, scheme)
        dis = bool(dk & dis_set)
        return drug, dis, dk
    res = pid.apply(present, axis=1, result_type='expand')
    pid[f'{kg}_drug'], pid[f'{kg}_dis'], _dk = res[0], res[1], res[2]
    pid[f'{kg}_present'] = pid[f'{kg}_drug'] & pid[f'{kg}_dis']

    # contamination: candidate drug's disease-edge targets overlap candidate disease keys
    def contam(row, dk_series):
        d = str(row['DrugBank_ID'])
        if d not in dd: return False
        targets = dd[d]
        keys = cand_disease_keys(row, scheme)
        # expand edge targets to aliases for matching
        talias = set()
        for t in targets: talias |= node_aliases(t, sep)
        return bool(keys & talias)
    pid[f'{kg}_contam'] = pid.apply(lambda r: contam(r, None), axis=1)
    pid[f'{kg}_eval'] = pid[f'{kg}_present'] & ~pid[f'{kg}_contam']

    # validate present-coverage vs coverage_annotation_v4 (the 116 survivors)
    ck = cov[cov.kg == kg]
    print(f"[{kg}] 253: present={int(pid[f'{kg}_present'].sum())} contaminated={int(pid[f'{kg}_contam'].sum())} "
          f"eval={int(pid[f'{kg}_eval'].sum())}  | target covered among 116 = {int(ck.covered.sum())}")

pid['n_present'] = pid[[f'{k}_present' for k in KGS]].sum(axis=1)
pid['n_eval'] = pid[[f'{k}_eval' for k in KGS]].sum(axis=1)
pid['kept_presence'] = pid['n_present'] >= 1
pid['kept'] = pid['n_eval'] >= 1
print(f"\n=== UNION over 253 ===")
print(f"presence-only (>=1 KG): kept={int(pid.kept_presence.sum())}")
print(f"presence AND not-contaminated (>=1 KG): kept={int(pid.kept.sum())}   [target = 116]")

# compare kept set to the actual final 116
g = pd.read_csv(os.path.join(GS, 'legacy', 'gold_standard_v4_bigset.tsv'), sep='\t')
gk = set(g.drug_id.astype(str) + '|' + g.disease_id.astype(str))
pid['key'] = pid.DrugBank_ID.astype(str) + '|' + pid.DOID_ID.astype(str)
pid['in_final'] = pid['key'].isin(gk)
both = (pid.kept & pid.in_final).sum()
print(f"my-kept ∩ final-116 = {both}  | my-kept-only = {int((pid.kept & ~pid.in_final).sum())}  | final-only = {int((~pid.kept & pid.in_final).sum())}")
pid.to_csv(os.path.join(ROOT, 'results', 'tables', 'reconstruct_116_audit.csv'), index=False)
print("\nwrote results/tables/reconstruct_116_audit.csv")
