#!/usr/bin/env python3
"""Size the candidate shared disease-gene retrieval benchmark.

Diagnostic only. Writes nothing. Answers one question: if notebook 06's
neighbourhood-retrieval queries were drawn from a shared, externally
adjudicated set instead of each KG's own disease-gene edges, how many
queries would survive, and how would they split by which KGs already
assert the pair?

Run from the repo root. Reads only data/ and data/gold_standards/.
"""
import collections, csv, gzip, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GS   = os.path.join(ROOT, 'data', 'gold_standards')
D    = os.path.join(ROOT, 'data')

# Hetionet: DOID diseases, NCBI gene ids. DaG/DuG/DdG are its disease-gene edges.
het = collections.defaultdict(set)
with open(os.path.join(D, 'hetionet', 'edges.tsv')) as fh:
    next(fh)
    for line in fh:
        p = line.split('\t')
        if p[1] in ('DaG', 'DuG', 'DdG'):
            het[p[0].split('::', 1)[1]].add(p[2].split('::', 1)[1])

# BioKG: MeSH descriptor diseases, UniProt proteins.
bio = collections.defaultdict(set)
with open(os.path.join(D, 'biokg', 'biokg.links.tsv')) as fh:
    for line in fh:
        p = line.rstrip('\n').split('\t')
        if len(p) >= 3 and p[1] == 'PROTEIN_DISEASE_ASSOCIATION':
            bio['MESH:' + p[2]].add(p[0])

# DOID -> MeSH. Many-to-many; a DOID with no mapping drops out of any shared
# set, which is crosswalk coverage loss and not a property of the graph.
d2m = collections.defaultdict(set)
for r in csv.DictReader(open(os.path.join(GS, 'bridges', 'mesh_to_doid.csv'))):
    d2m[r['doid'].strip()].add(r['mesh_id'].strip())
het_m = collections.defaultdict(set)
for dis, tgts in het.items():
    for m in d2m.get(dis, ()):
        het_m[m] |= tgts

sym2up = collections.defaultdict(set)
for r in csv.DictReader(open(os.path.join(GS, 'bridges', 'uniprot_genesproteins.csv'))):
    if r['gene_symbol']:
        sym2up[r['gene_symbol']].add(r['uniprot_id'])

ctd = set()
with gzip.open(os.path.join(GS, 'ctd', 'CTD_curated_genes_diseases.tsv.gz'), 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        p = line.rstrip('\n').split('\t')
        if len(p) >= 4:
            ctd.add((p[3], p[1], p[0]))

counts = collections.Counter()
diseases = set()
for dis, gid, sym in sorted(ctd):
    if dis not in het_m or dis not in bio:
        continue
    ups = sym2up.get(sym)
    if not ups:
        continue
    diseases.add(dis)
    h = gid in het_m[dis]
    b = bool(ups & bio[dis])
    counts['both' if (h and b) else 'neither' if not (h or b)
           else 'hetionet_only' if h else 'biokg_only'] += 1

total = sum(counts.values())
print(f'hetionet eligible diseases (>=5 targets) : '
      f'{sum(1 for t in het.values() if len(t) >= 5)}')
print(f'biokg    eligible diseases (>=5 targets) : '
      f'{sum(1 for t in bio.values() if len(t) >= 5)}')
print(f'CTD curated pairs resolvable in both     : {total} over {len(diseases)} diseases')
for k in ('both', 'neither', 'hetionet_only', 'biokg_only'):
    print(f'  {k:14s}: {counts[k]}')
