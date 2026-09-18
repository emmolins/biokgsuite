#!/usr/bin/env python3
"""
Faithful reconstruction of the 253 -> 234 -> 143 -> 116 funnel, enumerating the
pairs dropped at each stage.

  253 candidates (pairs_with_ids_v2.csv)
    -19  not strict / unresolved        -> 234 harness-ready (gold_standard_v2.tsv)
    -91  answer leakage (drug->disease edge in ANY of the 3 KGs, crosswalk-expanded)
                                         -> 143 leakage-free
    -27  drug or disease missing from all 3 KGs (coverage)
                                         -> 116 evaluable (gold_standard_v4_bigset.tsv)

Coverage (presence per KG) is taken from the validated resolver in
reconstruct_116.py (results/tables/reconstruct_116_audit.csv, n_present).
Leakage edges are streamed from the raw KG files; candidate diseases are expanded
through mondo.sssom.tsv + mesh_to_doid.csv so an edge keyed under an equivalent
ID is still caught.

Out: results/tables/build_116_dropped_step1_unresolved.csv
     results/tables/build_116_dropped_step2_leakage.csv
     results/tables/build_116_dropped_step3_coverage.csv
     results/tables/build_116_funnel.csv
"""
import os, pandas as pd
from collections import defaultdict
# This script sits two levels below the repository root (grounding/01_construction/),
# so ROOT needs three dirname() calls, not two. With two it resolved to
# grounding/ and every data path below it silently pointed at nothing.
ROOT = os.environ.get('BKS_ROOT') or os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__))))
GS = os.path.join(ROOT, 'data', 'gold_standards'); TAB = os.path.join(ROOT, 'results', 'tables')

p = pd.read_csv(os.path.join(GS, 'legacy', 'pairs_with_ids_v2.csv'))
v2 = pd.read_csv(os.path.join(GS, 'legacy', 'gold_standard_v2.tsv'), sep='\t')
g116 = pd.read_csv(os.path.join(GS, 'legacy', 'gold_standard_v4_bigset.tsv'), sep='\t')
pres = pd.read_csv(os.path.join(TAB, 'reconstruct_116_audit.csv'))[
    ['Pair_ID', 'n_present', 'primekg_drug', 'drkg_drug', 'biokg_drug']]
p = p.merge(pres, on='Pair_ID', how='left')

def doidnum(x):
    return None if pd.isna(x) else str(x).split(':')[-1].lstrip('0')
def nm(x):
    return str(x).strip().lower()
# multi-key membership: a candidate is "in" a gold set if drug_id matches AND
# (DOID matches OR disease name matches) — robust to disease-ID drift between files.
def keyset(df, drug='drug_id', doid='disease_id', name='disease_name'):
    ks = set()
    for _, r in df.iterrows():
        d = str(r[drug])
        ks.add(d + '|d|' + str(doidnum(r[doid])))
        ks.add(d + '|n|' + nm(r[name]))
    return ks
in234 = keyset(v2); in116 = keyset(g116)
def member(row, ks):
    d = str(row.DrugBank_ID)
    return (d + '|d|' + str(doidnum(row.DOID_ID)) in ks) or (d + '|n|' + nm(row.NewIndication) in ks)
p['dkey'] = p.DrugBank_ID.astype(str) + '|d|' + p.DOID_ID.map(doidnum).astype(str)  # for overlap report

# ---------- crosswalk maps ----------
def norm(x):
    x = str(x).split(' ')[0]  # strip sssom "{source=...}" suffixes
    return x
mondo2doid, mondo2mesh, doid2mondo, mesh2mondo = (defaultdict(set) for _ in range(4))
mesh2doid, doid2mesh = defaultdict(set), defaultdict(set)
with open(os.path.join(GS, 'bridges', 'mondo.sssom.tsv')) as f:
    for ln in f:
        if ln.startswith('#') or ln.startswith('subject_id'): continue
        c = ln.rstrip('\n').split('\t')
        if len(c) < 4: continue
        subj, obj = c[0], c[3]
        if subj.startswith('MONDO:'):
            m = subj.split(':')[-1].lstrip('0')
            if obj.startswith('DOID:'): d = obj.split(':')[-1].lstrip('0'); mondo2doid[m].add(d); doid2mondo[d].add(m)
            elif obj.upper().startswith('MESH'): me = obj.split(':')[-1]; mondo2mesh[m].add(me); mesh2mondo[me].add(m)
with open(os.path.join(GS, 'bridges', 'mesh_to_doid.csv')) as f:
    next(f)
    for ln in f:
        parts = ln.rstrip('\n').rsplit(',', 2)   # mesh_id may contain commas
        if len(parts) != 3: continue
        me, d, _src = parts
        me = norm(me).split(':')[-1]; d = str(d).split(':')[-1].lstrip('0')
        mesh2doid[me].add(d); doid2mesh[d].add(me)

def expand(row):
    M = None if pd.isna(row.MONDO_ID) else str(row.MONDO_ID).split(':')[-1].lstrip('0')
    D = None if pd.isna(row.DOID_ID) else str(row.DOID_ID).split(':')[-1].lstrip('0')
    Me = None if pd.isna(row.MeSH_ID) else str(row.MeSH_ID).split(':')[-1]
    mondos, doids, meshes = set(), set(), set()
    if M: mondos.add(M)
    if D: doids.add(D)
    if Me: meshes.add(Me)
    if M: doids |= mondo2doid[M]; meshes |= mondo2mesh[M]
    if D: mondos |= doid2mondo[D]; meshes |= doid2mesh[D]
    if Me: mondos |= mesh2mondo[Me]; doids |= mesh2doid[Me]
    return mondos, doids, meshes

# ---------- KG drug->disease edges ----------
pk_e, dr_e, bk_e = defaultdict(set), defaultdict(set), defaultdict(set)
with open(os.path.join(GS.replace('gold_standards', ''), 'biokg', 'biokg.links.tsv') if False else os.path.join(ROOT, 'data', 'biokg', 'biokg.links.tsv')) as f:
    for ln in f:
        c = ln.rstrip('\n').split('\t')
        if len(c) == 3 and c[1] == 'DRUG_DISEASE_ASSOCIATION':
            bk_e[c[0]].add(c[2])             # disease = MeSH D######
with open(os.path.join(ROOT, 'data', 'drkg', 'drkg.tsv')) as f:
    for ln in f:
        c = ln.rstrip('\n').split('\t')
        if len(c) != 3: continue
        a, b = c[0], c[2]
        if a.startswith('Compound::') and b.startswith('Disease::'):
            db = a.split('::', 1)[1]; dis = b.split('::', 1)[1]
            if dis.startswith('MESH:'): dr_e[db].add('M' + dis[5:])
            elif dis.startswith('DOID:'): dr_e[db].add('D' + dis[5:].lstrip('0'))
with open(os.path.join(ROOT, 'data', 'primekg', 'primekg.csv')) as f:
    next(f)
    for ln in f:
        c = ln.split(',')
        if len(c) < 11: continue
        if c[4] == 'drug' and c[9] == 'disease':
            for d in str(c[8]).split('_'): pk_e[c[3]].add(d.lstrip('0'))
        elif c[4] == 'disease' and c[9] == 'drug':
            for d in str(c[3]).split('_'): pk_e[c[8]].add(d.lstrip('0'))

# HetioNet (CtD/CpD): Compound::DB.. -> Disease::DOID:..   (disease = DOID)
het_e = defaultdict(set)
with open(os.path.join(ROOT, 'data', 'hetionet', 'edges.tsv')) as f:
    next(f)
    for ln in f:
        c = ln.rstrip('\n').split('\t')
        if len(c) == 3 and c[1] in ('CtD', 'CpD') and c[0].startswith('Compound::') and c[2].startswith('Disease::'):
            db = c[0].split('::', 1)[1]; doid = c[2].split('DOID:')[-1].lstrip('0')
            het_e[db].add(doid)
# OpenBioLink (DIS_DRUG): DOID:.. <- PUBCHEM.COMPOUND:CID  (map CID -> DrugBank)
cid2db = defaultdict(set)
with open(os.path.join(ROOT, 'data', 'openbilink', 'pubchem_to_drugbank.csv')) as f:
    next(f)
    for ln in f:
        a = ln.rstrip('\n').split(',')
        if len(a) >= 2: cid2db[a[1].strip()].add(a[0].strip())
ob_e = defaultdict(set)
with open(os.path.join(ROOT, 'data', 'openbilink', 'edges.csv')) as f:
    for ln in f:
        c = ln.rstrip('\n').split('\t')
        if len(c) >= 3 and c[1] == 'DIS_DRUG' and c[0].startswith('DOID:') and c[2].startswith('PUBCHEM'):
            doid = c[0].split('DOID:')[-1].lstrip('0'); cid = c[2].split(':')[-1]
            for db in cid2db.get(cid, ()): ob_e[db].add(doid)

def leaked(row):
    db = str(row.DrugBank_ID); mondos, doids, meshes = expand(row)
    if db in pk_e and (mondos & pk_e[db]): return True
    if db in bk_e and (meshes & bk_e[db]): return True
    if db in dr_e:
        keys = {'M' + m for m in meshes} | {'D' + d for d in doids}
        if keys & dr_e[db]: return True
    if db in het_e and (doids & het_e[db]): return True   # HetioNet
    if db in ob_e and (doids & ob_e[db]): return True     # OpenBioLink
    return False

p['leaked'] = p.apply(leaked, axis=1)
# Step 1 (exact): "unresolvable IDs" = missing DrugBank OR missing all disease IDs
# OR placeholder/TBD indication.  253 - 19 = 234 harness-ready.
def placeholder_ind(s):
    s = str(s).lower(); return any(k in s for k in ['tbd', 'population expansion', 'new indication ('])
p['unresolvable'] = p.DrugBank_ID.isna() | (p.MONDO_ID.isna() & p.DOID_ID.isna() & p.MeSH_ID.isna()) | p.NewIndication.map(placeholder_ind)
p['in234'] = ~p['unresolvable']
p['in116'] = p.apply(lambda r: member(r, in116), axis=1)
# Step 3 (-27): coverage. The slide phrases this as "true drug has no dossier in
# any of the 3 KGs"; in practice a pair is evaluable only if the drug AND disease
# both resolve in the same KG (this reproduces the per-KG 111/90/78 counts exactly).
# We keep co-presence in >=1 KG; drug-only presence is also reported for reference.
p['covered'] = p.n_present.fillna(0) >= 1
p['drug_anyKG'] = p[['primekg_drug', 'drkg_drug', 'biokg_drug']].fillna(False).any(axis=1)

# ---------- funnel ----------
s1_drop = p[~p.in234]                                   # step 1
s1 = p[p.in234]
s2_drop = s1[s1.leaked]                                 # step 2 (leakage)
s2 = s1[~s1.leaked]
s3_drop = s2[~s2.covered]                               # step 3 (coverage)
s3 = s2[s2.covered]

cols = ['Pair_ID', 'Drug', 'Active', 'NewIndication', 'TA_simplified', 'YearSet',
        'DrugBank_ID', 'MONDO_ID', 'DOID_ID', 'StrictRepurposing']
s1_drop[cols].to_csv(os.path.join(TAB, 'build_116_dropped_step1_unresolved.csv'), index=False)
s2_drop[cols + ['leaked']].to_csv(os.path.join(TAB, 'build_116_dropped_step2_leakage.csv'), index=False)
s3_drop[cols + ['n_present']].to_csv(os.path.join(TAB, 'build_116_dropped_step3_coverage.csv'), index=False)
p.to_csv(os.path.join(TAB, 'build_116_funnel.csv'), index=False)

print(f"253 candidates")
print(f"  step1 dropped (not in 234)     : {len(s1_drop):>3}   target -19   -> {len(s1)} harness-ready (target 234)")
print(f"  step2 dropped (answer leakage) : {len(s2_drop):>3}   target -91   -> {len(s2)} leakage-free (target 143)")
print(f"  step3 dropped (coverage)       : {len(s3_drop):>3}   target -27   -> {len(s3)} evaluable    (target 116)")
repro = set(s3.Pair_ID); truepid = set(p[p.in116].Pair_ID)
print(f"\nreproduced-116 ∩ true-116(mapped): {len(repro & truepid)}  | my-only: {len(repro-truepid)}  | true-only: {len(truepid-repro)}")
print("wrote build_116_dropped_step{1,2,3}*.csv + build_116_funnel.csv")
