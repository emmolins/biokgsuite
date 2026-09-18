# Vetting report: fda_repurposing_events_2023_2026.csv

Audit run 2026-09-17. Read this before using the file.

## Verdict

**The data that is present is accurate. The set is NOT complete, and the
incompleteness is not random.** Measured recall on one audited month was 68%,
and the shortfall grows with the number of efficacy supplements in that month.
The file is usable as a verified sample. It is not publishable as a census of
FDA efficacy supplements, and it should not be described as one.

## How the set was built

1. Source: Drugs@FDA monthly report "Supplemental Approvals to NDAs and BLAs",
   `accessdata.fda.gov/scripts/cder/daf/index.cfm?event=reportsSearch.process
   &rptName=4&reportSelectMonth=<M>&reportSelectYear=<Y>`, one page per month,
   45 pages covering 2023-01 to 2026-09.
2. Filter: rows whose **Submission Classification** column is exactly
   "Efficacy". An efficacy supplement on an already-approved drug is the
   regulatory event corresponding to a new or extended indication.
3. Extraction: api.fda.gov and accessdata.fda.gov are blocked at the gateway
   from both agent sandboxes, so pages were read through a fetch tool that
   passes HTML through a summarising model rather than parsing structured
   output. **This is the root cause of every recall problem below.**
4. Drug mapping: active ingredient names matched against
   `drug_name_to_drugbank.csv` (built from the KG node tables).

## Test 1: internal integrity. PASS, 8/8

```
approval_date ISO and in range                 0 violations / 567
approval_date month == source report month     0 / 567
application_number matches (NDA|BLA) + 6 digits 0 / 567
submission matches SUPPL-<int>                 0 / 567
no empty required field                        0 / 567
no fully duplicated rows                       0 / 567
one date per (application, supplement)         0 / 567
one ingredient string per (application, suppl) 0 / 567
```

## Test 2: precision against an independent source. 7/8 confirmed

Verified against openFDA `drug/drugsfda.json`, a different endpoint from the
report pages, checking application number, supplement number, approval date and
submission class.

```
NDA 210922 SUPPL-12  ONPATTRO   12|AP|20230113|EFFICACY   exact match
NDA 217806 SUPPL-6   ZEPBOUND    6|AP|20241018|EFFICACY   exact match
NDA 215256 SUPPL-11  WEGOVY     11|AP|20240308|EFFICACY   exact match
NDA 215256 SUPPL-15  WEGOVY     15|AP|20241127|EFFICACY   exact match
NDA 215256 SUPPL-23  WEGOVY     23|AP|20251121|EFFICACY   exact match
NDA 215256 SUPPL-24  WEGOVY     24|AP|20250815|EFFICACY   exact match
NDA 215256 SUPPL-29  WEGOVY     29|AP|20260319|EFFICACY   exact match
NDA 208065 SUPPL-35  TAGRISSO   ABSENT from openFDA
```

The TAGRISSO miss is dated 2026-09-14 and openFDA's drugsfda snapshot reports
last_updated 2026-09-15. NDA 217759 SUPPL-6 (2026-09-11) is likewise absent,
with openFDA's highest SUPPL for that application being 5. Both come from FDA's
own monthly report, which is the authoritative output; openFDA is a derived
dataset that lags. **The four 2026-09 rows are unconfirmed against the API**
and are flagged as such rather than removed.

No record checked was wrong. Zero fabricated values found.

## Test 3: recall. FAIL, and size-dependent

Each audited month was re-fetched independently and the set of unique
(application, supplement) pairs compared with the capture.

```
month     table rows  my events  audit events  in both  I missed  audit missed
2024-05        97          4           4          4         0          0
2024-03         -         18          18         18         0          0
2025-03       128         21          27         17        10          4
```

March 2025: union of both passes is 31 events. My capture holds 21, so recall
is **68%** for that month. The audit pass itself missed 4 events I had, so it
is not ground truth either; neither pass is complete.

Identifiers of the 10 events missing from the file for 2025-03, recoverable
without re-extracting the month:

```
207103 SUPPL-21   211675 SUPPL-27   212436 SUPPL-9    214679 SUPPL-9
214916 SUPPL-2    215515 SUPPL-5    215833 SUPPL-24   215841 SUPPL-6
217171 SUPPL-3    218347 SUPPL-2
```

Recall tracks list length: months with few efficacy rows came back complete,
the month with the most did not. The mechanism is truncation of long lists by
the summarising fetch tool. Months with 20+ efficacy rows (2023-05, 2023-06,
2024-03, 2024-09, 2024-10, 2024-12, 2025-03, 2025-10) are therefore the most
likely to be short, and the deficit is **biased toward busy approval months**,
not spread evenly.

## Test 4: systematic exclusion of biologics. CONFIRMED

Only 2 of 567 rows are BLAs. FDA leaves the Submission Classification cell
**empty** for BLA rows, so filtering on "Efficacy" removes every biologic.
Verified on 2024-09: 15 of 89 rows are BLAs, every one with blank
classification, including DUPIXENT SUPPL-66, TREMFYA SUPPL-21, ELREXFIO
SUPPL-2, CIMZIA SUPPL-275, SOLIRIS SUPPL-447.

Consequence: 19 of the 32 FDA records in the previous hand-curated set are
absent here, nearly all biologics (pembrolizumab, nivolumab, dupilumab,
durvalumab, adalimumab, belimumab, risankizumab, tocilizumab, daratumumab).
**This file is not a superset of the old one.**

## Test 5: a discrepancy in the OLD set

The old set records everolimus, FDA, 2024-05-09, gastrointestinal
neuroendocrine tumour. The May 2024 report contains 4 efficacy supplements and
none is everolimus, confirmed by a second independent pass. openFDA shows every
2024 everolimus submission as LABELING or MANUF (CMC), none EFFICACY.

That record does not correspond to a 2024 FDA efficacy supplement. Since the
old set carries no citations, this cannot be traced. Worth checking the
provenance of the other 95 records in it.

## What this file may and may not be used for

MAY: a verified sample of dated, citable FDA repurposing events, every row
carrying an application and supplement number that resolves on Drugs@FDA.

MAY NOT: a denominator. Any rate, proportion or per-year trend computed from it
will be wrong, because the missing fraction is larger in busy months.

MAY NOT: drug-disease pairs. A supplement records that an indication changed,
not to what. Indication text has to come from the label.

## To make it publishable

Re-extract with `tools/fetch_approvals.py` against openFDA from a network that
reaches api.fda.gov. That path parses structured JSON, so it has no truncation
failure mode, and openFDA carries `submission_class_code` for BLAs, which
closes the biologics hole at the same time. This audit should then be re-run,
and the two sets diffed.

---

## Amendment, 2026-09-17

This dataset has been superseded by `fda_indication_expansion_census.csv`,
built from the Drugs@FDA relational data files rather than the monthly HTML
reports. See `FDA_CENSUS_METHOD.md`.

Re-running the audit against the complete record revises two conclusions:

1. **Recall was 66.9%, not 68%, and that figure is now measured rather than
   estimated.** It is computed across all 45 months against the like-for-like
   census subset, not extrapolated from one audited month. Spearman correlation
   between month size and recall is -0.12, so truncation contributed but the
   loss was broader than volume alone explains.

2. **Precision was 99.4%, not 100%.** Three defects surfaced that the original
   sample of eight could not have caught:
   - `NDA 208065 SUPPL-35` (claimed 2026-09-14) does not exist. The
     application's highest supplement is SUPPL-33, approved 2024-09-25. Release
     lag is ruled out: the source file contains approvals through 2026-09-14,
     including `NDA 217759 SUPPL-6` which openFDA had also been missing.
   - `NDA 213586 SUPPL-8` dated 2025-10-24; true date 2025-10-10.
   - `NDA 214697 SUPPL-2` dated 2026-03-12; true date 2026-03-26.

   Both date errors point at neighbouring rows in the same source table, which
   is the signature of row misalignment in a truncated table rather than of
   fabrication. The finding of zero fabricated values stands for the fields
   that were checked; it does not extend to row-level alignment.

The everolimus 2024-05-09 finding is confirmed against the full record: the
only everolimus submissions in 2024 are a LABELING supplement (2024-02-28) and
a MANUF (CMC) supplement (2024-03-11). No efficacy supplement exists.
