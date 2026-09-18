#!/usr/bin/env python3
"""Turn fetched approval events into a citable post-cutoff approval set.

Reads data/gold_standards/legacy/approvals_raw/ (written by tools/fetch_approvals.py)
and produces two files:

  approvals_auto.csv    records where BOTH the drug resolved to a DrugBank
                        accession AND the indication resolved to a DOID by an
                        exact name match. Safe to use without review.
  approvals_review.csv  everything else, with candidate matches ranked, for a
                        human to accept or reject.

DESIGN RULE
-----------
Nothing enters approvals_auto.csv on a fuzzy match. A gold standard whose DOIDs
were guessed by substring similarity is worse than a smaller one, because the
errors are invisible downstream and they bias exactly the metric the file
exists to support. Partial matches go to the review file with their evidence
and are not used until someone signs them off.

The FDA half needs indication TEXT, which drugsfda does not carry: an EFFICACY
supplement records that an indication changed, not what it changed to. Pass
--fetch-labels to pull indications_and_usage from openFDA drug/label.json per
application. That needs network, same as the fetch step.

    python3 tools/build_approval_gold_standard.py --fetch-labels
"""
import argparse
import csv
import json
import re
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GS = ROOT / 'data' / 'gold_standards' / 'legacy'   # superseded pipeline: reads and writes under legacy/
BRIDGES = ROOT / 'data' / 'gold_standards' / 'bridges'
RAW = GS / 'approvals_raw'
LABEL_API = 'https://api.fda.gov/drug/label.json'
UA = 'biokgsuite-bib/1.0 (academic benchmark; contact via repository)'


def norm(s):
    """Lowercase, strip punctuation and collapse whitespace."""
    return re.sub(r'\s+', ' ', re.sub(r'[^a-z0-9 ]+', ' ', str(s).lower())).strip()


def load_crosswalks():
    """Prefer drug_name_to_drugbank.csv, fall back to drugbank_drugs.csv.

    The preferred file is built by tools/build_drug_name_xref.py from the
    graphs' own node tables, which key drugs by DrugBank accession. It resolves
    25 of the 29 gold-standard drugs against drugbank_drugs.csv's 13, and every
    one of the 25 agrees with the hand-curated accession. Ambiguous names (one
    name, several accessions) are flagged in that file and skipped here rather
    than resolved by picking one.
    """
    drug = {}
    preferred = BRIDGES / 'drug_name_to_drugbank.csv'
    if preferred.exists():
        with open(preferred, encoding='utf-8') as fh:
            for r in csv.DictReader(fh):
                if r.get('ambiguous') == '1':
                    continue
                drug.setdefault(norm(r['name_normalised']), r['drugbank_id'])
    else:
        print('  NOTE: drug_name_to_drugbank.csv absent; falling back to '
              'drugbank_drugs.csv, which has no biologics. '
              'Run tools/build_drug_name_xref.py first.')
        with open(BRIDGES / 'drugbank_drugs.csv', encoding='utf-8') as fh:
            for r in csv.DictReader(fh):
                drug.setdefault(norm(r['name']), r['drugbank_id'])
    dis = {}
    with open(BRIDGES / 'do_diseases.csv', encoding='utf-8') as fh:
        for r in csv.DictReader(fh):
            dis.setdefault(norm(r['mondo_name']), r['doid'])
    return drug, dis


def fetch_label_indication(appno, cache):
    if appno in cache:
        return cache[appno]
    url = (f'{LABEL_API}?search=openfda.application_number:%22{appno}%22'
           f'&limit=1')
    try:
        req = urllib.request.Request(url, headers={'User-Agent': UA})
        with urllib.request.urlopen(req, timeout=60) as r:
            payload = json.loads(r.read())
        res = payload.get('results', [])
        text = ' '.join(res[0].get('indications_and_usage', [])) if res else ''
    except Exception:                                   # noqa: BLE001
        text = ''
    cache[appno] = text
    time.sleep(0.3)
    return text


def candidate_doids(text, dis, limit=5):
    """Exact match first; otherwise disease names contained in the text.

    Returned as CANDIDATES, never applied automatically. Containment is a hint
    for a human, not evidence: 'lung cancer' is contained in 'non-small cell
    lung cancer' and they are different DOIDs.
    """
    n = norm(text)
    if n in dis:
        return [(dis[n], n, 'exact')]
    hits = [(doid, name, 'contained') for name, doid in dis.items()
            if len(name) > 6 and name in n]
    hits.sort(key=lambda h: -len(h[1]))
    return hits[:limit]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fetch-labels', action='store_true')
    args = ap.parse_args()

    src = RAW / 'fda_indication_events.csv'
    if not src.exists():
        sys.exit(f'{src} not found. Run tools/fetch_approvals.py first.')
    rows = list(csv.DictReader(open(src, encoding='utf-8')))
    drug_x, dis_x = load_crosswalks()
    print(f'{len(rows):,} FDA events | {len(drug_x):,} DrugBank names | '
          f'{len(dis_x):,} disease names')

    # Preflight: drugbank_drugs.csv was rebuilt from Open Targets and contains
    # essentially NO monoclonal antibodies. Measured against the 29 drugs in the
    # existing hand-curated gold standard, 16 fail to resolve and all 16 are
    # biologics (pembrolizumab, nivolumab, durvalumab, trastuzumab deruxtecan,
    # the bispecifics, semaglutide). Modern indication approvals are dominated
    # by exactly those, so without a better crosswalk this step resolves the
    # wrong half of the data and does so quietly. Said out loud here rather
    # than emitting an empty approvals_auto.csv and looking like a no-op.
    n_mab = sum(1 for n in drug_x if n.endswith('mab'))
    if n_mab < 50:   # the graph-derived crosswalk carries 239
        print()
        print('  ' + '=' * 68)
        print('  WARNING: the DrugBank name crosswalk contains only '
              f'{n_mab} antibody names.')
        print('  data/gold_standards/bridges/drugbank_drugs.csv is small-molecule '
              'heavy and also')
        print('  carries ~2,400 rows whose "name" is a bare CHEMBL id. '
              'Biologic approvals')
        print('  will land in approvals_review.csv as "no drugbank match" '
              'rather than')
        print('  being resolved.')
        print('  Fix: download the DrugBank open vocabulary '
              '(drugbank_all_drugbank_vocabulary.csv,')
        print('  CC BY-NC 4.0, covers every accession with synonyms) and point '
              'load_crosswalks()')
        print('  at it. See data/POST_CUTOFF_APPROVALS.md.')
        print('  ' + '=' * 68)
        print()

    cache, auto, review = {}, [], []
    for i, r in enumerate(rows, 1):
        ingredients = [x for x in r['active_ingredients'].split('|') if x]
        db_hits = [(x, drug_x[norm(x)]) for x in ingredients if norm(x) in drug_x]
        text = (fetch_label_indication(r['application_number'], cache)
                if args.fetch_labels else '')
        cands = candidate_doids(text, dis_x) if text else []

        base = {k: r[k] for k in (
            'agency', 'application_number', 'submission_type',
            'submission_number', 'submission_class_code', 'approval_date',
            'sponsor', 'brand_names', 'active_ingredients', 'source_url')}

        if len(db_hits) == 1 and cands and cands[0][2] == 'exact':
            auto.append({**base,
                         'drug_name': db_hits[0][0],
                         'drugbank_id': db_hits[0][1],
                         'disease_name': cands[0][1],
                         'doid': cands[0][0],
                         'match_basis': 'exact drug + exact disease'})
        else:
            review.append({**base,
                           'drugbank_candidates': '|'.join(f'{n}={d}' for n, d in db_hits),
                           'n_drug_candidates': len(db_hits),
                           'indication_text': text[:400],
                           'doid_candidates': '|'.join(f'{d}={n}({how})'
                                                       for d, n, how in cands),
                           'reason': ('no drugbank match' if not db_hits else
                                      'ambiguous drug' if len(db_hits) > 1 else
                                      'no indication text' if not text else
                                      'no exact disease match')})
        if i % 200 == 0:
            print(f'  {i:,}/{len(rows):,}')

    stamp = datetime.now(timezone.utc).isoformat(timespec='seconds')
    for name, data in (('approvals_auto.csv', auto),
                       ('approvals_review.csv', review)):
        p = GS / name
        if data:
            with open(p, 'w', newline='', encoding='utf-8') as fh:
                w = csv.DictWriter(fh, fieldnames=list(data[0].keys()))
                w.writeheader()
                w.writerows(data)
        print(f'{name}: {len(data):,} rows')

    (GS / 'APPROVALS_BUILD.txt').write_text(
        f'# built:  {stamp}\n# by:     tools/build_approval_gold_standard.py\n'
        f'# source: {src.relative_to(ROOT)}\n'
        f'# auto:   {len(auto)} (exact drug AND exact disease match)\n'
        f'# review: {len(review)} (need a human decision)\n'
        f'# labels fetched: {args.fetch_labels}\n'
        f'# Nothing in approvals_auto.csv was matched fuzzily.\n')
    print('\nReview approvals_review.csv before using anything from it.')


if __name__ == '__main__':
    main()
