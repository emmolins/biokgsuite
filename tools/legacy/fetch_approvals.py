#!/usr/bin/env python3
"""Fetch regulatory approval records for the post-cutoff approval set.

WHY THIS EXISTS
---------------
The prospective axis of notebook 07 rested on two hand-curated sources with no
URLs, no document numbers, no retrieval dates and no version stamp (see
data/POST_CUTOFF_APPROVALS.md). One of them, fda_drugs_id_mapping.csv, is
missing outright. Hand curation is also why it cannot be refreshed: nobody can
re-derive it.

This script replaces the machine-derivable part with a fetch from official
sources, stamped the way tools/fetch_ctd.sh stamps CTD.

WHAT IT CAN AND CANNOT DO
-------------------------
  FDA   automated. openFDA drug/drugsfda.json. Gives NDA/BLA application
        number, supplement number, approval date, sponsor, brand name and
        active ingredients. An EFFICACY-class supplement is the regulatory
        event that corresponds to a NEW INDICATION, which is what a prospective
        evaluation needs; an ORIG approval is a new molecular entity.
  EMA   automated. Two published XLSX bulk files. The post-authorisation file
        carries variation opinions, which is where extensions of indication
        live. The medicines file carries initial authorisations.
  PMDA  NOT automated. Approvals are published as per-review-period PDFs with
        no stable bulk endpoint.
  TGA   NOT automated. tga.gov.au disallows automated fetching in robots.txt,
        and there is no bulk export of the new/extended-indications list.

For PMDA and TGA the existing 22 hand-curated records stand. They need citation
fields added by hand (PMDA approval-list reference, ARTG number). This script
does not invent them, and neither should anything else.

NETWORK
-------
api.fda.gov and ema.europa.eu are blocked from the agent sandboxes, so run this
yourself, the way you ran tools/fetch_ctd.sh:

    python3 tools/fetch_approvals.py --start 2023-01-01 --end 2026-12-31

Writes to data/gold_standards/legacy/approvals_raw/ plus a provenance file recording
each source URL, its bytes, its SHA-256 and the retrieval timestamp.
"""
import argparse
import hashlib
import json
import sys
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data' / 'gold_standards' / 'legacy' / 'approvals_raw'
PROV = ROOT / 'data' / 'gold_standards' / 'legacy' / 'approvals_raw' / 'APPROVALS_PROVENANCE.txt'

OPENFDA = 'https://api.fda.gov/drug/drugsfda.json'
EMA_FILES = {
    'ema_medicines.xlsx':
        'https://www.ema.europa.eu/en/documents/report/'
        'medicines-output-medicines-report_en.xlsx',
    'ema_post_authorisation.xlsx':
        'https://www.ema.europa.eu/en/documents/report/'
        'medicines-output-post_authorisation-report_en.xlsx',
}
UA = 'biokgsuite-bib/1.0 (academic benchmark; contact via repository)'


def _get(url, binary=False, tries=4):
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': UA})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read() if binary else json.loads(r.read())
        except Exception as exc:                       # noqa: BLE001
            if attempt == tries - 1:
                raise
            wait = 2 ** attempt
            print(f'  retry {attempt + 1}/{tries - 1} in {wait}s ({exc})')
            time.sleep(wait)
    return None


def fetch_fda(start, end, api_key=None):
    """Page through openFDA and flatten to one row per SUBMISSION.

    openFDA matches at APPLICATION level: a search on a submission field
    returns every application containing any matching submission, not the
    matching submissions. Filtering therefore has to happen here, per
    submission, or the output silently includes decades of unrelated history.
    """
    s, e = start.replace('-', ''), end.replace('-', '')
    query = f'submissions.submission_status_date:[{s}+TO+{e}]'
    rows, skip, limit = [], 0, 1000
    seen_apps = 0
    while True:
        url = f'{OPENFDA}?search={query}&limit={limit}&skip={skip}'
        if api_key:
            url += f'&api_key={api_key}'
        payload = _get(url)
        results = payload.get('results', [])
        if not results:
            break
        total = payload.get('meta', {}).get('results', {}).get('total', 0)
        seen_apps += len(results)
        for app in results:
            appno = app.get('application_number', '')
            sponsor = app.get('sponsor_name', '')
            products = app.get('products', []) or []
            brands = sorted({p.get('brand_name', '') for p in products if p.get('brand_name')})
            ingredients = sorted({
                i.get('name', '')
                for p in products for i in (p.get('active_ingredients') or [])
                if i.get('name')})
            for sub in app.get('submissions', []) or []:
                date = str(sub.get('submission_status_date', ''))
                if not (s <= date <= e):
                    continue
                if sub.get('submission_status') != 'AP':
                    continue
                cls = sub.get('submission_class_code', '')
                stype = sub.get('submission_type', '')
                # EFFICACY supplement = new or extended indication.
                # ORIG = new molecular entity. Everything else (labelling,
                # manufacturing) is not an indication event.
                if not (cls == 'EFFICACY' or stype == 'ORIG'):
                    continue
                rows.append({
                    'agency': 'FDA',
                    'application_number': appno,
                    'submission_type': stype,
                    'submission_number': sub.get('submission_number', ''),
                    'submission_class_code': cls,
                    'submission_class_description':
                        sub.get('submission_class_code_description', ''),
                    'review_priority': sub.get('review_priority', ''),
                    'approval_date': f'{date[:4]}-{date[4:6]}-{date[6:]}',
                    'sponsor': sponsor,
                    'brand_names': '|'.join(brands),
                    'active_ingredients': '|'.join(ingredients),
                    'source_url': (f'https://www.accessdata.fda.gov/scripts/cder/daf/'
                                   f'index.cfm?event=overview.process&ApplNo='
                                   f'{appno.lstrip("NDABLA")}'),
                })
        print(f'  applications {seen_apps:,} / {total:,}  '
              f'-> {len(rows):,} indication events')
        skip += limit
        if skip >= total or skip >= 25000:   # openFDA caps skip without a key
            if skip >= 25000 and skip < total:
                print('  NOTE: hit the 25,000 skip cap. Pass --api-key for the '
                      'full range, or narrow --start/--end and concatenate.')
            break
        time.sleep(0.4)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--start', default='2023-01-01')
    ap.add_argument('--end', default='2026-12-31')
    ap.add_argument('--api-key', default=None,
                    help='openFDA API key; raises rate and skip limits')
    ap.add_argument('--skip-ema', action='store_true')
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).isoformat(timespec='seconds')
    prov = [f'# Provenance for the prospective approval gold standard',
            f'# fetched:     {stamp}',
            f'# fetched_by:  tools/fetch_approvals.py',
            f'# date_range:  {args.start} .. {args.end}',
            '']

    print(f'FDA (openFDA drugsfda), {args.start} to {args.end}')
    try:
        fda = fetch_fda(args.start, args.end, args.api_key)
    except Exception as exc:                            # noqa: BLE001
        sys.exit(f'FDA fetch failed: {exc}\n'
                 f'If this is a network error, run it where api.fda.gov is '
                 f'reachable. Nothing was written.')
    import csv
    fda_path = OUT / 'fda_indication_events.csv'
    cols = list(fda[0].keys()) if fda else []
    with open(fda_path, 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(fda)
    digest = hashlib.sha256(fda_path.read_bytes()).hexdigest()
    n_eff = sum(1 for r in fda if r['submission_class_code'] == 'EFFICACY')
    prov += [f'== FDA',
             f'   endpoint:  {OPENFDA}',
             f'   query:     submissions.submission_status_date in range,',
             f'              submission_status=AP, class=EFFICACY or type=ORIG',
             f'   file:      {fda_path.relative_to(ROOT)}',
             f'   rows:      {len(fda)}  (EFFICACY supplements: {n_eff}, '
             f'ORIG: {len(fda) - n_eff})',
             f'   sha256:    {digest}',
             f'   note:      openFDA is updated continuously; cite the fetch',
             f'              date above and the per-row application number.',
             '']
    print(f'  wrote {len(fda):,} rows -> {fda_path.relative_to(ROOT)}')

    if not args.skip_ema:
        print('EMA (published bulk XLSX)')
        for fname, url in EMA_FILES.items():
            try:
                blob = _get(url, binary=True)
            except Exception as exc:                    # noqa: BLE001
                print(f'  {fname}: FAILED ({exc}); skipped')
                prov += [f'== EMA {fname}', f'   url:    {url}',
                         f'   status: FETCH FAILED', '']
                continue
            p = OUT / fname
            p.write_bytes(blob)
            prov += [f'== EMA {fname}',
                     f'   url:    {url}',
                     f'   bytes:  {len(blob)}',
                     f'   sha256: {hashlib.sha256(blob).hexdigest()}', '']
            print(f'  wrote {fname} ({len(blob):,} bytes)')

    prov += ['== PMDA', '   NOT FETCHED. Approvals are published as per-review-'
             'period PDFs', '   with no stable bulk endpoint. The 11 existing '
             'records are hand-', '   curated and still need an approval-list '
             'reference each.', '',
             '== TGA', '   NOT FETCHED. tga.gov.au disallows automated fetching '
             'in robots.txt', '   and publishes no bulk export of new/extended '
             'indications. The 11', '   existing records are hand-curated and '
             'still need an ARTG number each.', '']
    PROV.write_text('\n'.join(prov) + '\n')
    print(f'\nProvenance -> {PROV.relative_to(ROOT)}')
    print('Next: tools/build_approval_gold_standard.py maps these to DrugBank '
          'and DOID.')


if __name__ == '__main__':
    main()
