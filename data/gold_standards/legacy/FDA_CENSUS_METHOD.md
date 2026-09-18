# FDA indication-expansion census, 2023-01-01 to 2026-09-14

Method record for `fda_indication_expansion_census.csv` and its denominator
`fda_indication_expansion_census_frame.csv`. Written 2026-09-17.

## Why this replaces the previous dataset

`fda_repurposing_events_2023_2026.csv` was assembled by fetching the Drugs@FDA
monthly approvals pages and keeping rows classified `EFFICACY`. Vetting
(`FDA_REPURPOSING_VETTING.md`) established that the rows it contained were
accurate but that the set was incomplete in a non-random way: overall recall of
66.9% against the full record, and biologics excluded almost entirely because
BLA supplements carry a blank Submission Classification in those pages. A set
with those properties can support a claim about the drugs it contains. It
cannot serve as a denominator.

## Source

Drugs@FDA Data Files, the relational export of the Drugs@FDA database.

| | |
|---|---|
| Publisher | U.S. Food and Drug Administration, CDER |
| Landing page | https://www.fda.gov/drugs/drug-approvals-and-databases/drugsfda-data-files |
| File | `dafdata20260916.zip`, 12 tab-delimited tables |
| Data last updated | 2026-09-16 (latest approval present: 2026-09-14) |
| SHA-256 of tables | `e1739f79478a01bb925e83e5f5657b3f08f700d75310b5ddbac69a85c81c9eef` |
| Retrieved | 2026-09-17 |

Complete by construction. There is no pagination, no rendering layer and no
search ranking between the database and the file, which is precisely what the
previous route could not guarantee.

## Definition

A repurposing event is a molecule that already holds FDA approval gaining
approval for an indication it did not previously hold.

The operational rule is built entirely from FDA's own submission taxonomy
(`ActionTypes_Lookup.SupplCategoryLevel1Code` / `SupplCategoryLevel2Code`), so
the definition is the regulator's rather than ours. The same underlying event
reaches the database through three different codes, hence three tiers.

| Tier | Rule | n | Adjudication |
|---|---|---|---|
| A | `SUPPL` tagged `Efficacy-New Indication` | 303 | not required |
| B | `SUPPL` tagged `Efficacy-Accelerated Approval`, excluding confirmatory studies | 25 | required |
| C | `ORIG` with class `TYPE 9`, `TYPE 9- BLA`, `TYPE 10` | 13 | required |
| | **total** | **341** | 38 flagged |

Tier B exists because FDA tags accelerated approvals by the approval *pathway*
rather than by the change made. Inspection shows these are new indications
(KEYTRUDA, OPDIVO, YERVOY, BRAFTOVI and similar). They are disproportionately
oncology, so an `INDICATION`-only rule would bias the set by therapeutic area.

Tier C exists because a new indication can be filed as a separate application
rather than a supplement, which FDA classes as "new indication submitted as a
distinct NDA". ZEPBOUND against MOUNJARO and YEZTUGO against SUNLENCA are both
in this tier. A supplement-only rule misses this route entirely. A minority of
these are reformulations (OPDIVO QVANTIG is subcutaneous nivolumab), which is
why the tier is flagged rather than assumed.

Excluded, with the reason recorded rather than left implicit:

| Category | Reason |
|---|---|
| `Efficacy-New Patient Population` (147) | same disease, wider age band |
| `Efficacy-New Dosing Regimen` (63) | same disease |
| `Efficacy-Labeling Change With Clinical Data` (274) | new evidence, existing claim |
| `Efficacy-New Route Of Administration` (3), `Rx To OTC Switch` (3) | same disease |
| `Efficacy-Accelerated Approval Confirmatory Study` | converts an existing approval to full; indication already held |
| all `LABELING`, `MANUF (CMC)`, `REMS`, `BIOEQUIV` | not efficacy changes |

## Funnel

Every filter step with its row count, so the path from the raw file to the
final set is auditable line by line.

```
193,752  all submissions in the source
 20,514  status date in [2023-01-01, 2026-09-14]
 20,003  SubmissionStatus == AP
  5,372  ApplType in (NDA, BLA)              <- the denominator, written to _frame.csv
    341  matches the tiered indication rule  <- the census
```

## Verification

The previous scraped dataset is high precision and low recall, which makes it a
useful probe: every row it contains should reappear in the census, and any that
does not is a census defect requiring explanation. Run via
`tools/validate_fda_census.py`.

| Check | Result |
|---|---|
| A. Containment | 514/515 scraped rows present in the frame |
| B. Date agreement | 512/514 shared rows agree |
| C. The 10 identifiers the scrape lost in 2025-03 | 10/10 recovered |
| D. The 5 biologics supplements the scrape could not see (2024-09) | 5/5 recovered |
| E. Biologics overall | 1,212 BLA rows in the frame against 2 in the scrape |
| F. Everolimus 2024-05-09 | no FDA counterpart, confirming the entry was spurious |

The three residual discrepancies were run down individually and all three are
defects in the scrape, not the census:

- `NDA 208065 SUPPL-35` (claimed 2026-09-14). Application 208065 has no
  supplement beyond SUPPL-33, approved 2024-09-25. The file does contain approvals dated 2026-09-14, and it contains
  `NDA 217759 SUPPL-6` dated 2026-09-11 which openFDA had also been missing, so
  this is not release lag. The row is a phantom.
- `NDA 213586 SUPPL-8`, scrape date 2025-10-24, true date 2025-10-10.
- `NDA 214697 SUPPL-2`, scrape date 2026-03-12, true date 2026-03-26.

This amends the earlier vetting: scrape precision is 99.4%, not 100%. Both date
errors point at neighbouring rows in the same source table, consistent with row
misalignment when a long table is truncated.

Scrape recall measured against the like-for-like census subset (approved
`EFFICACY`-class supplements) is **66.9%**, matching the 68% estimated from a
single audited month. Spearman correlation between month size and scrape recall
is -0.12, so truncation was a contributor but the loss was broad rather than
purely volume-driven.

## Contents

341 events across 229 applications and 205 distinct active-ingredient strings.
145 carry an orphan-designation flag. 331 of 341 carry the URL of the FDA
approval letter for that exact submission, so any row can be verified against
its own primary document without re-running anything.

## DrugBank mapping

`fda_indication_expansion_census_drugbank.csv`, built by
`tools/map_census_drugbank.py`, kept separate so the census does not change
when the crosswalk does.

| Match level | n |
|---|---|
| exact | 172 |
| heuristic: salt stripped | 67 |
| heuristic: biologic suffix stripped | 46 |
| heuristic: suffix and prefix stripped | 5 |
| unmatched | 44 |
| incomplete combination | 6 |
| no ingredient listed | 1 |

290 of 341 events (85.0%) map fully. The residual 35 distinct ingredients are
not name-matching failures: dostarlimab, amivantamab, epcoritamab, tislelizumab,
retifanlimab, inclisiran and the rest are simply absent from the DrugBank
release backing the crosswalk. This is a real coverage limit and it skews
against recent biologics, which should be stated wherever the mapped subset is
used.

## Limitations to state in the paper

1. **Frame boundary.** Drugs@FDA covers CDER-regulated drugs and therapeutic
   biologics. CBER products (vaccines, blood products, most cell and gene
   therapies) are not in this database. The census is complete inside that
   frame; the frame is stated, not implied.
2. **Tiers B and C are flagged, not adjudicated.** 38 events need one-by-one
   confirmation against their approval letters before the set is used as a
   strict denominator.
3. **No disease attached yet.** The census records that an indication was
   added, not which one. Attaching the disease requires the approval letter or
   a label version diff and is the next stage.
4. **One malformed row** in `ApplicationDocs.txt` (line 78694, 9 fields where 8
   are expected) was skipped by the parser. It affects at most one document URL
   and no event record.
5. **Right-censoring.** Approvals in the days immediately before the release
   date may not yet be posted. The latest approval present is 2026-09-14.

## Reproduction

```bash
python3 tools/build_fda_census.py profile --zip <dafdata dir or zip>
python3 tools/build_fda_census.py build   --zip <dafdata dir or zip> \
    --start 2023-01-01 --out data/gold_standards/fda_indication_expansion_census.csv
python3 tools/validate_fda_census.py \
    --frame  data/gold_standards/fda_indication_expansion_census_frame.csv \
    --census data/gold_standards/fda_indication_expansion_census.csv \
    --scrape data/gold_standards/fda_repurposing_events_2023_2026.csv
python3 tools/map_census_drugbank.py \
    --census data/gold_standards/fda_indication_expansion_census.csv \
    --xref   data/gold_standards/drug_name_to_drugbank.csv
```

Counts reproduce exactly against the same dated release. A later release will
add rows at the right edge and may correct earlier ones, which is why the
release date and hash are recorded above.
