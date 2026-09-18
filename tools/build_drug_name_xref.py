#!/usr/bin/env python3
"""Build a drug name -> DrugBank accession crosswalk, offline.

WHY
---
data/gold_standards/bridges/drugbank_drugs.csv was rebuilt from Open Targets and
contains essentially NO monoclonal antibodies: of the 29 drugs in the
post-cutoff approval set, 16 fail to resolve and all 16 are biologics. The
cause is upstream, not a bug in that build: Open Targets simply carries no
DrugBank cross-reference for most biologics. Pembrolizumab's crossReferences
in ot_drug_molecule are EMA and DailyMed only, no drugbank, so any build keyed
on a DrugBank xref drops it correctly and silently.

The knowledge graphs themselves key drugs BY DrugBank accession and carry
names, so they can supply what Open Targets cannot, with no network and no
DrugBank account. That is also the right boundary for this benchmark: a drug
absent from every evaluated graph is unevaluable anyway, so a crosswalk scoped
to the graphs loses nothing the evaluation could have used.

Verified against the hand-curated gold standard: BioKG yields DB09037 for
pembrolizumab, DB09035 for nivolumab, DB12159 for dupilumab and DB13928 for
semaglutide, all matching the curated values exactly.

LIMITS, stated rather than discovered later
-------------------------------------------
Coverage is bounded by the graph snapshots. BioKG is a 2020 release, so drugs
approved since then (trastuzumab deruxtecan, glofitamab, the bispecifics) are
absent from it and will not resolve here. That is a property of the graphs, not
of this script, and it is exactly why those drugs cannot be scored on those
graphs either.

    python3 tools/build_drug_name_xref.py
"""
import csv
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from biokgsuite.io import find_config, load_config, load_kg   # noqa: E402

OUT = ROOT / 'data' / 'gold_standards' / 'bridges' / 'drug_name_to_drugbank.csv'
PROV = ROOT / 'data' / 'gold_standards' / 'bridges' / 'DRUG_NAME_XREF_PROVENANCE.txt'
LEGACY = ROOT / 'data' / 'gold_standards' / 'bridges' / 'drugbank_drugs.csv'
DB_RE = re.compile(r'\b(DB\d{5,})\b')


def norm(s):
    return re.sub(r'\s+', ' ', re.sub(r'[^a-z0-9 ]+', ' ', str(s).lower())).strip()


def main():
    cfg = load_config(find_config(str(ROOT)))
    entries = defaultdict(set)          # normalised name -> {(db_id, source)}
    raw_name = {}
    counts = {}

    for name in cfg['knowledge_graphs']:
        try:
            kg, nodes = load_kg(name, cfg)
        except Exception as exc:                        # noqa: BLE001
            print(f'{name:12s}: not loadable here ({type(exc).__name__}), skipped')
            continue
        dt = cfg['knowledge_graphs'][name]['entity_types'].get('Drug')
        if not dt:
            print(f'{name:12s}: no Drug entity type configured, skipped')
            continue
        drugs = nodes[nodes['type'].isin(dt)]
        n = 0
        for nid, nm in zip(drugs['id'].astype(str), drugs['name'].astype(str)):
            m = DB_RE.search(nid)
            if not m or not nm.strip():
                continue
            key = norm(nm)
            if not key:
                continue
            entries[key].add((m.group(1), name))
            raw_name.setdefault(key, nm)
            n += 1
        counts[name] = n
        print(f'{name:12s}: {n:,} drug nodes carrying a DrugBank accession')

    if LEGACY.exists():
        n = 0
        with open(LEGACY, encoding='utf-8') as fh:
            for r in csv.DictReader(fh):
                key = norm(r['name'])
                # Rows whose "name" is a bare CHEMBL id are not names.
                if not key or key.startswith('chembl'):
                    continue
                entries[key].add((r['drugbank_id'], 'drugbank_drugs.csv'))
                raw_name.setdefault(key, r['name'])
                n += 1
        counts['drugbank_drugs.csv'] = n
        print(f'{"legacy":12s}: {n:,} usable rows from drugbank_drugs.csv')

    rows, conflicts = [], 0
    for key in sorted(entries):
        ids = {db for db, _ in entries[key]}
        srcs = sorted({s for _, s in entries[key]})
        if len(ids) > 1:
            conflicts += 1
        for db in sorted(ids):
            rows.append({'name_normalised': key,
                         'name': raw_name[key],
                         'drugbank_id': db,
                         'sources': '|'.join(srcs),
                         'ambiguous': int(len(ids) > 1)})

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    unique_names = len(entries)
    n_mab = sum(1 for k in entries if k.endswith('mab'))
    PROV.write_text(
        '# Provenance for drug_name_to_drugbank.csv\n'
        f'# built:    {datetime.now(timezone.utc).isoformat(timespec="seconds")}\n'
        '# built_by: tools/build_drug_name_xref.py\n'
        '# sources:  each configured KG\'s own node table (drug-typed nodes whose\n'
        '#           id carries a DB accession), plus drugbank_drugs.csv with its\n'
        '#           CHEMBL-id-as-name rows dropped.\n'
        + ''.join(f'#   {k}: {v}\n' for k, v in counts.items()) +
        f'# unique names:      {unique_names}\n'
        f'# rows:              {len(rows)}\n'
        f'# antibody names:    {n_mab}\n'
        f'# ambiguous names:   {conflicts} (same name, >1 accession; flagged in\n'
        '#                    the ambiguous column, NOT silently resolved)\n'
        '# Coverage is bounded by the graph snapshots; drugs approved after a\n'
        '# graph was built are absent from it by construction.\n')

    print(f'\nunique names   : {unique_names:,}')
    print(f'rows           : {len(rows):,}')
    print(f'antibody names : {n_mab:,}  (was 0 in drugbank_drugs.csv)')
    print(f'ambiguous      : {conflicts:,}  flagged, not resolved')
    print(f'-> {OUT.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
