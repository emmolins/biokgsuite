# History of the prospective gold standard before the scripted build

Moved here from `data/PROSPECTIVE_GOLD_STANDARD.md` on 2026-09-18. These are the
2026-09-16 trace that found the hand-assembled sets uncitable and the 2026-09-17
rebuild notes. Paths in them refer to the flat `data/gold_standards/` layout that
existed at the time; the files now sit in this folder (see `README.md`). The
current set and its construction are described in `data/PROSPECTIVE_GOLD_STANDARD.md`.


# Prospective gold standard: provenance and version history

Traced 2026-09-16. This is the post-cutoff approval set behind notebook 07's
prospective generalisation axis.

## Verdict first

**This set is not currently citable.** There are no URLs, no document numbers,
no retrieval dates and no version stamp anywhere in it. Every other data source
in this repo carries provenance: CTD has a release date and SHA-256, MONDO has
a release tag, Hetionet and BioKG have release files. This one has a docstring
saying "accessed 2025".

A BiB reviewer can verify none of it, and neither can you in six months.

## Where it comes from

Two separate sources, which is the first thing to know.

### 1. `data/gold_standards/legacy/fda_drugs_id_mapping.csv`  — MISSING

The main FDA set. Hand-curated, covering FDA novel drug approvals June 2021 to
April 2026. Read by notebook 07 cell 31.

Columns: `No`, `Combo_Flag`, `DrugBank_ID`, `Active Ingredient`,
`Approval Date`, `Indication`, `Drug Name (Brand)`.

The notebook applies three exclusions to it: combination products flagged in
`Combo_Flag`, entry No. 203 (Kygevvi, DB02594, ambiguous nucleoside mapping),
and nine non-disease indications (imaging agents, contrast agents, a cosmetic
and a procedure).

This file is absent from the repo. It is the one genuine reconstruction blocker
left: the underlying FDA "Novel Drug Approvals" lists are public and republished
annually, but the DrugBank accession mapping, the combination flags and the
curated indication strings are hand work that is not recoverable from them.

### 2. `src/biokgsuite/gold/regulatory_approvals.py`

A Python module with the approval tuples written inline. Self-contained, reads
no files. 96 records:

| agency | records | unique (drug, disease) contribution |
|--------|---------|-------------------------------------|
| FDA    | 32      | used only as a dedup filter, never as data |
| EMA    | 42      | European public assessment reports |
| PMDA   | 11      | Japan, New Drug Approval Information |
| TGA    | 11      | Australia, ARTG |

Date range 2023-01-23 to 2025-02-20. 29 distinct drugs, 34 distinct diseases,
47 unique (drug, disease) pairs, 32 of them approved by more than one agency.
All 96 records carry a DOID.

Note the FDA list here is a **different and smaller** set than the CSV above
(32 records, 2023-2025, versus the CSV's June 2021 to April 2026). The notebook
never uses it as data. It exists only so `get_non_fda_pairs()` can drop
EMA/PMDA/TGA pairs that the FDA set already covers.

Per the module docstring, approval dates are the date of positive opinion for
EMA and the date of formal approval for PMDA and TGA. Drug IDs are DrugBank
accessions cross-referenced via DrugBank 6.0; disease IDs are Disease Ontology
DOIDs.

## What each field would need to be citable

Nothing below is currently recorded for any of the 96 records.

- **EMA**: the EPAR product page URL and the CHMP positive-opinion date from the
  meeting highlights, or the EPAR procedure number. EMA publishes both.
- **PMDA**: the approval list PDF (published per review period) and the item
  number within it.
- **TGA**: the ARTG entry number, which is a stable public identifier.
- **FDA**: the Drugs@FDA application number (NDA/BLA) and supplement number,
  which is what makes a *new indication* citable rather than just the drug.
- **All**: the date each record was retrieved, and a version stamp for the set
  as a whole.

## A bug this trace uncovered

`get_non_fda_pairs()` removes every EMA/PMDA/TGA pair that also appears in the
module's 32-record FDA list, so the two sources are not double counted when
merged. That is correct when the FDA CSV is present.

When the CSV is **absent**, as now, the filter still ran, and it was discarding
most of the only data left:

```
non-FDA records available      64
after the FDA-overlap filter   21     (what the notebook was using)
unique pairs available         43
unique pairs after the filter  15
```

15 pairs, minus the ones that fail to resolve into each graph, never reached
`MIN_NEW_PAIRS=10`. That is why both graphs scored N/A on the prospective axis,
and it was a plumbing artefact on top of the missing file rather than a fact
about the graphs. Notebook 07 now uses all 43 unique non-FDA pairs when the FDA
set is unavailable, and says so.

## Recommended next steps

1. Recover `fda_drugs_id_mapping.csv` from wherever it was curated. This is the
   highest-value item: it restores the prospective axis to full strength.
2. Add source columns to `regulatory_approvals.py`: a URL or document number and
   a retrieval date per record. The module is already a list of tuples, so this
   is two extra fields, and it converts an uncitable asset into a citable one.
3. Stamp the set with a version and a build date, as CTD and MONDO are, and
   record that stamp in notebook 07's checkpoint alongside `fda_available`.
4. State the date-of-positive-opinion versus date-of-approval distinction in the
   methods. Mixing EMA opinion dates with PMDA/TGA/FDA approval dates shifts
   some records by weeks relative to a KG snapshot cutoff, which matters for a
   prospective evaluation and nowhere else.


---

# Rebuild attempt, 2026-09-17

Asked to rebuild the 2023-2026 approval set for FDA + EMA + PMDA + TGA. Result:
the machine-derivable half is now scripted, two agencies cannot be automated,
and one crosswalk turns out to be unfit for the job. No records were written
from model memory, and none should be: an invented approval date or DOID that
looks plausible is worse than a missing file, because it is invisible
downstream and it biases the one metric this file exists to support.

## Network

api.fda.gov, ema.europa.eu, pmda.go.jp and tga.gov.au are all blocked from both
agent sandboxes (gateway answers 403 to CONNECT). So these scripts are written
to be run by you, on an unrestricted network, exactly as tools/fetch_ctd.sh was.

openFDA itself is live and current: its drugsfda dataset reported
last_updated 2026-09-16.

## What is now automated

`tools/fetch_approvals.py`

  FDA   openFDA `drug/drugsfda.json`. Filters to submission_status=AP and to
        submission_class_code=EFFICACY (a new or extended indication) or
        submission_type=ORIG (a new molecular entity). Yields the NDA/BLA
        application number, the supplement number, the approval date, sponsor,
        brand names and active ingredients, plus a Drugs@FDA URL per row. That
        is a citable regulatory event, which is what the old set lacked.

        One subtlety encoded in the script: openFDA matches at APPLICATION
        level, so a search on a submission field returns every application
        containing any matching submission. Filtering has to happen per
        submission client-side or the output quietly includes decades of
        unrelated history.

  EMA   two published bulk XLSX files, URLs confirmed from EMA's own download
        page:
          medicines-output-medicines-report_en.xlsx            (authorisations)
          medicines-output-post_authorisation-report_en.xlsx   (variations)
        The post-authorisation file is the one that matters here: an extension
        of indication is a Type II variation, not a new authorisation.

Both are SHA-256 stamped with a retrieval timestamp into
`data/gold_standards/legacy/approvals_raw/APPROVALS_PROVENANCE.txt`.

`tools/build_approval_gold_standard.py` maps those events to DrugBank and DOID
and splits the output in two: `approvals_auto.csv` for records where BOTH the
drug and the disease matched EXACTLY, and `approvals_review.csv` for everything
else, with ranked candidates and the reason. Nothing reaches the auto file on a
fuzzy match. Containment is offered as a hint for a human and never applied:
"lung cancer" is contained in "non-small cell lung cancer" and they are
different DOIDs.

## What cannot be automated

  PMDA  approvals are published as per-review-period PDFs with no stable bulk
        endpoint.
  TGA   tga.gov.au disallows automated fetching in robots.txt, and there is no
        bulk export of the new/extended-indications list. The right page exists
        ("Prescription medicines: new or extended uses...") but it is HTML and
        off-limits to a crawler.

The 22 existing PMDA and TGA records stand. They need a citation field each,
added by hand: a PMDA approval-list reference and an ARTG number. Both are short
lookups. I have not invented either.

## The crosswalk problem, which blocks this independently of network

`data/gold_standards/bridges/drugbank_drugs.csv` (rebuilt earlier from Open Targets
`ot_drug_molecule`) contains essentially **no monoclonal antibodies**. Measured
directly against the 29 drugs in the current gold standard:

```
resolvable by name via drugbank_drugs.csv : 13
NOT resolvable                            : 16
```

All 16 failures are biologics: adalimumab, tocilizumab, belimumab, nivolumab,
pembrolizumab, daratumumab, durvalumab, dupilumab, sacituzumab govitecan,
semaglutide, risankizumab, trastuzumab deruxtecan, glofitamab, teclistamab,
epcoritamab, talquetamab.

Modern indication approvals are dominated by exactly these, so the drug-mapping
step would resolve the wrong half of the data. The file also carries roughly
2,400 rows whose "name" is a bare CHEMBL identifier.

`tools/build_approval_gold_standard.py` now detects this and says so loudly
rather than emitting an empty auto file that looks like a no-op.

Two fixes, either works:

  a. DrugBank's open vocabulary (`drugbank_all_drugbank_vocabulary.csv`,
     CC BY-NC 4.0) covers every accession with synonyms, biologics included.
     Needs a free academic account.
  b. DrugCentral's identifiers dump, which needs no account and which this
     repo already knows how to fetch: `tools/build_drugbank_xref.py` pulls it
     for a different purpose (id-to-id cross references for MATRIX). Its
     structures and synonyms tables carry names for biologics too, so the same
     download can serve name-to-DrugBank.

Point `load_crosswalks()` at whichever you take. This is
also worth doing for its own sake, since `drugbank_drugs.csv` is a listed gold
standard used beyond this notebook, and REFERENCE_SETS.md currently describes
it without noting the biologic gap.

## Order of work

1. Download the DrugBank open vocabulary. Without it the mapping resolves
   small molecules only, and this is the cheapest of the three steps.
2. Run `tools/fetch_approvals.py --start 2023-01-01 --end 2026-12-31` on a
   machine with network access. Add `--api-key` if you have an openFDA key;
   without one the skip cap is 25,000 and the range may need splitting.
3. Run `tools/build_approval_gold_standard.py --fetch-labels`, then work
   through `approvals_review.csv`. The indication text comes from the label
   endpoint, since drugsfda records that an indication changed but not to what.
4. Add ARTG numbers and PMDA references to the 22 existing records by hand.
5. Version-stamp the result and record the stamp in notebook 07's checkpoint
   next to `fda_available`.
