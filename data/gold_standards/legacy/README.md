# data/gold_standards/legacy/

Every gold-standard attempt that preceded `post_cutoff_approvals/` (128 pairs, built by
`tools/build_post_cutoff_approvals.py` from the vendored `source/` release). Frozen on
2026-09-18; kept so that earlier numbers in notes, checkpoints and drafts can be
traced, not for use. Nothing in `src/`, `config/`, the notebooks or the paper
reads these files, except that notebook 07 falls back to `fda_drugs_id_mapping.csv`
when the post-cutoff set is absent. The tools that produced them
(`clean_116.py`, `build_v5_evidence.py`, `build_fda_census.py`,
`validate_fda_census.py`, `map_census_drugbank.py`, `corroborate_*.py`,
`parse_pmda.py`, `fetch_approvals.py`, `build_approval_gold_standard.py`) live in
`tools/legacy/` and are themselves superseded. Paths written inside the `.md`
records below refer to the flat layout that existed before the move.

The full account of why each generation was replaced is `HISTORY.md` in this
folder (the 2026-09-16 trace and the 2026-09-17 rebuild notes); the short version is that none of these
sets could be rebuilt from a cited source by a script, and the ones that could
mixed line extensions and population expansions with new indications.

## Generation 1: the 253 candidate pairs and the v2 to v4 sets (2026-09-16 and earlier)

| file | what it was |
|---|---|
| `pairs_with_ids_v2.csv` | 253 candidate post-2023 approval pairs with DrugBank, DOID, MONDO and MeSH ids, produced with Biomni; `PROVENANCE_253_candidate_pairs.md` records how |
| `pairs_with_ids_2023_2026.csv` | earlier 170-row candidate table with the same kind of columns (DrugBank, PubChem, DOID, MONDO, MeSH) |
| `usable_pairs_v2_audited.csv` | an evaluation-flag file over the 253; stale, and the source of a counting confusion recorded in `PROVENANCE_253_candidate_pairs.md` |
| `gold_standard_v2.tsv` | 234 harness-ready pairs (253 minus rows the loaders could not key) |
| `gold_standard_v3_common3.tsv` | 76 pairs present in all three graphs then compared, for the coverage-controlled arm |
| `gold_standard_v4_bigset.tsv` | 116 evaluable pairs after the KG-coverage and contamination audit |
| `coverage_annotation_v4.csv` | per-graph coverage of the 116 (drug in graph, disease in graph, mechanism edges) |
| `expanded_prospective_gold_standard.tsv` | 49-row hand-assembled expansion attempt with agencies and dates; never used downstream |
| `fda_drugs_id_mapping.csv` | 214 FDA approvals with a hand-set `Combo_Flag`; hand-curated, no download reproduces it |

Why superseded: candidate generation was not scripted from a regulator source,
so the 253 could not be regenerated; the audit that took 253 to 116 had no
build script.

## Generation 2: cleaning the 116 (2026-09-17)

| file | what it was |
|---|---|
| `gold_standard_v5_dedup.tsv` | 99 pairs after collapsing rows that recorded the same approval from two agencies |
| `gold_standard_v5_verified.tsv` | 94 pairs after removing five that the primary register contradicted (`EMA_PMDA_VERIFICATION.md`) |
| `gold_standard_v5_indication.tsv` | 90 pairs after removing formulation changes, vaccine strain updates and accelerated-to-full conversions |
| `gold_standard_v5_strict.tsv` | 80 pairs after also removing age-band and severity expansions |
| `gold_standard_v5_strict_evidence.csv` | the 80 with a regulator identifier and URL per row (`tools/build_v5_evidence.py`) |
| `cleaning_manifest.csv` | every row removed by `tools/clean_116.py`, with the step and reason |
| `EMA_PMDA_VERIFICATION.md` | the register checks behind the five factual removals |
| `RECHECK_CHECKLIST.md` | second-reviewer batches drafted 2026-09-18 for the 115-pair post-cutoff set (v0.5.0); the `prospective_recheck.csv` it asks for was never produced, the fourth-pass review inside the source package and the second-model identifier proposals in `source/` did that work instead |
| `pmda_approvals.csv` | 1,889 rows parsed from PMDA's List of Approved Drugs (`tools/parse_pmda.py`), used for the PMDA verification |

Why superseded: the cleaning was sound but the input (the 116) still could not
be rebuilt from a source, so the whole chain rested on an unreproducible table.

## Generation 3: the FDA census (2026-09-17)

| file | what it was |
|---|---|
| `approvals_raw/fda_efficacy_supplements_2023_2026.psv` | Drugs@FDA monthly supplemental-approval reports, 2023 onwards |
| `fda_repurposing_events_2023_2026.csv` | 567 events scraped from Drugs@FDA; vetted in `FDA_REPURPOSING_VETTING.md`, which found it incomplete |
| `fda_indication_expansion_census.csv`, `.meta.json` | 341 FDA efficacy supplements 2023 onwards, one row per supplement, from the Drugs@FDA relational files |
| `fda_indication_expansion_census_frame.csv` | the 5,372-row denominator the census was drawn from |
| `fda_indication_expansion_census_drugbank.csv`, `.mapping.json` | the census with DrugBank ids attached (`tools/map_census_drugbank.py`) |
| `fda_census_validation_2026-09-17.log` | output of `tools/validate_fda_census.py` |
| `FDA_CENSUS_METHOD.md` | the method record |

Why superseded: FDA only, and Drugs@FDA does not say which indication a
supplement added, so the census could not classify rows without the approval
letters. That work moved into the `repurposing-dataset` package, which fetches
letters and EPAR text for five regulators.

## Generation 4: earlier releases of the source package (2026-09-17 and 2026-09-18)

| file | what it was |
|---|---|
| `approvals_raw/repurposing-dataset-v0.2.0/`, `repurposing-dataset-0.2.0+20260917.zip` | first scripted release: PMDA, TGA and Health Canada only; produced the 272-pair set |
| `approvals_raw/repurposing-dataset-v0.5.0/`, `repurposing-dataset-0.5.0+20260918.zip` | ingredient-level classification of all five regulators after the independent re-read; produced the 115-pair set |

The current release, v0.6.0 (128 pairs), is in `source/`. Each release's
CHANGELOG and PROVENANCE inside its zip describe what changed.
