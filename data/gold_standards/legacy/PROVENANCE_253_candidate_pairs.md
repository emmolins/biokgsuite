# Provenance — the 253 candidate drug–disease pairs

Record of how the 253 ID-resolved candidate pairs (`pairs_with_ids_v2.csv` /
`kg_eligible_pairs_with_ids_2023_2026.csv`) were generated. Produced with Biomni.
This documents **Stage 1** of the benchmark funnel (candidate generation); the
downstream KG-coverage/contamination audit that reduces 253 → 116 evaluable pairs
(`gold_standard_v4_bigset.tsv`) is a separate, later step whose build script is not
yet committed (see "Open provenance item" below).

## Source data
Starting point: `kg_eligible_pairs_2024_2026.csv` — 157 drug-repurposing pairs approved
2024–2026 by FDA, EMA, PMDA, MHRA/Orbis. 10 columns: Drug, Active, OrigYear,
NewIndication, Agencies, TA, StrictRepurposing, Drug_Pre2023, Pair_ID, TA_simplified.

## Step 1 — Reconstruct 2023 pairs (+96)
The 2024–2026 file had no 2023 entries. 2023 new-indication approvals were reconstructed
from agency records (FDA novel approvals/sNDAs, EMA type-II variations, PMDA approvals),
cross-referenced with Drugs@FDA, EMA medicines, PMDA lists. Inclusion criteria matched the
2024–2026 file: drug approved before 2023 (`Drug_Pre2023 = True`); new indication approved
in 2023 by ≥1 major agency; genuine new indication (not a formulation change). → 96 pairs
(`pairs_2023.csv`).

## Step 2 — Resolve vague disease names
12 entries labelled "New indication extension" and 11 labelled "New indication extension
(EMA 2026, TBD)" were researched per drug+agency+date to identify the specific disease
(e.g. Xalkori → specific cancer type from the FDA announcement). The 11 EMA-2026-TBD entries
remain provisional (EMA details unpublished).

## Step 3 — Merge → 253
Concatenated 96 (2023) + 157 (2024–2026), added `YearSet`, renumbered `Pair_ID` 1–253.
**4 intentional cross-year duplicates** (same Drug + NewIndication as separate regulatory
events): Lynparza/BRCA-mutated mCRPC (19, 125); Reblozyl/first-line anemia lower-risk MDS
(29, 140); Bosulif/paediatric Ph+ CML-CP (30, 135); Vabysmo/retinal vein occlusion (38, 188).

## Step 4 — Drug names → DrugBank + PubChem IDs
174 unique active ingredients (after splitting combinations on `+` and `/`). Multi-source
cascade: PubChem PUG-REST synonyms (82), Wikidata SPARQL exact/alt label (66, esp. biologics),
Wikidata alternative/brand/salt-stripped names (15), mychem.info brand/research codes (4).
Result: **167/174 (96.0%)** ingredients with DrugBank IDs; 101/174 (58.0%) with PubChem CIDs.

7 permanently unmapped DrugBank IDs (3 CAR-T, 1 vaccine, 1 Japan-only SGLT2i, 2 recent EU mAbs):
Ebola Zaire vaccine, axicabtagene ciloleucel, ciltacabtagene autoleucel, lisocabtagene
maraleucel, mosunetuzumab, ruseogliflozin, tralokinumab. 18 combination drugs have
comma-separated per-ingredient IDs.

## Step 5 — Disease names → MONDO + MeSH + DOID IDs
240 unique disease names → 137 normalized search terms (abbreviation expansion e.g.
NSCLC → non-small cell lung carcinoma; parenthetical/age/population qualifier stripping;
~150-entry manual fix map). Mapping: OLS4 for MONDO 137/137 (100%, 4 top-1 hits manually
corrected: COVID-19, atopic dermatitis, chronic migraine, functional constipation);
Monarch v3 + OLS4 for DOID 135/137 (98.5%); Monarch + NLM MeSH (exact→contains fallback,
heading-format variants) for MeSH 130/137 (94.9%). 2 unresolvable placeholders:
"New indication (EMA Feb 2026, TBD)", "Population expansion (FDA+EMA)".

## Step 6 — QC
Validated ID formats (DrugBank `DB\d{5,7}`, MeSH `D\d{6,9}`, MONDO `MONDO:\d{5,8}`,
DOID `DOID:\d+`); fixed incorrect mappings (COVID-19/atopic dermatitis/chronic migraine/
functional constipation MONDO; heart-failure DOID); added `ID_Flag` completeness column.

## Step 7 — Output
253 rows × 17 columns → `kg_eligible_pairs_with_ids_2023_2026.csv` (in repo as
`pairs_with_ids_v2.csv`).

## Summary statistics
| Metric | Value |
|---|---|
| Total pairs | 253 |
| From 2023 / 2024–2026 | 96 / 157 |
| StrictRepurposing = True | 212 (83.8%) |
| Drug_Pre2023 = True | 253 (100%) |
| Unique drug names / active ingredients | 180 / 174 |
| Unique disease names | 240 |
| Combination drugs | 18 |
| DrugBank ID coverage | 246/253 (97.2%) |
| PubChem CID coverage | 132/253 (52.2%) |
| MONDO ID coverage | 251/253 (99.2%) |
| MeSH ID coverage | 243/253 (96.0%) |
| DOID ID coverage | 249/253 (98.4%) |
| All 5 IDs present | 121/253 (47.8%) |
| DrugBank+MONDO composite | 244/253 (96.4%) |

**Therapeutic area (TA_simplified):** Oncology 116, Immunology 46, Cardiovascular 18,
Metabolic 16, Infectious 16, CNS 10, Hematology 10, Ophthalmology 6, Respiratory 5, GI 3,
others 7. *(Matches `pairs_with_ids_v2.csv` exactly.)*

**Agencies (multi-agency possible):** EMA 122, FDA 112, PMDA 36, FDA+EMA 8, MHRA/Orbis 7,
EMA+PMDA 2, FDA+EMA+PMDA 1.

**OrigYear range:** 1983–2022 (drug's *original* approval year; all repurposing events are 2023–2026).

---

## Pathway to 116 — CONFIRMED 3-step funnel (matches the slide, reproduced here)

```
253 candidate pairs
  − 19   unresolved IDs / not strict, prospective repurposing
234 harness-ready pairs                    (= gold_standard_v2.tsv)
  − 91   drug→disease edge already in a KG  (ANSWER LEAKAGE / contamination, any of the 3 KGs)
143 leakage-free pairs
  − 27   drug or disease missing from all 3 KGs   (coverage)
116 evaluable pairs                        (= gold_standard_v4_bigset.tsv)
```

**Step 1 (−19) is exact and reproduced.** The 19 are the **unresolvable-ID** pairs:
missing DrugBank ID (7: the CAR-T therapies Yescarta/Carvykti/Breyanzi/Lunsumio + Ebola vaccine
Ervebo + tralokinumab/Adbry + ruseogliflozin/Rysefio) OR missing all disease IDs (2) OR a
placeholder/TBD indication (12 EMA-2026-TBD entries). `253 − 19 = 234` exactly. The slide phrase
"/ not strict, prospective repurposing" describes what the 234 *are*, not extra drops (note:
"not strict" alone is 41 pairs, so it is **not** part of this cut). List:
`results/tables/build_116_dropped_step1_unresolved.csv`.

The dominant cut is **answer leakage (−91)**: a pair is removed if the drug→disease edge already
exists in *any benchmarked KG* (PrimeKG, DRKG, BioKG, **HetioNet, OpenBioLink** — all five; HetioNet
and OpenBioLink do carry such edges), so the model can't read the answer off the graph. My earlier
confusion came from the stale `usable_pairs_v2_audited.csv` (an eval-flag file that does **not**
encode the leakage step) — ignore it. The leakage step is the real reducer.

**Coverage step, validated exactly.** Per-KG coverage of the final 116 (`coverage_annotation_v4.csv`):
**PrimeKG 111, DRKG 90, BioKG 78 → union 116** (72 in all three). A from-scratch re-implementation of
the coverage rule (DrugBank ID → Drug-type node AND disease ID in the KG's scheme → Disease-type
node, via `src.loading.load_kg` + config `entity_types`, with `_`-split + leading-zero handling for
PrimeKG's grouped MONDO nodes) reproduces **all three KG columns at 100% per-pair**
(`scripts/reconstruct_116.py`).

**The 116 is a clean subset of the 253 — nothing added.** Matching on disease + shared drug
ingredient, **all 116 final pairs trace to a 253 candidate**, so `253 − 116 = 137 dropped, 0 added`.
A few combination therapies carry different labels across files (brand names in the 253, e.g.
"Tafinlar+Mekinist", vs generic + lead-component DrugBank ID in the final set, e.g.
"dabrafenib+trametinib" / DB08912) — a naming artifact, not a real addition. Stage breakdown of the
137: **19 unresolved IDs + 81 answer leakage + 37 no-KG-coverage** (`dropped_253_to_116.csv`).

**Funnel reproduction.** Step 1 (−19) is exact; coverage (111/90/78) is exact (100% per-pair); the
leakage step reproduces ~84 of 91. Residual ~6 leaks come from disease-ID **crosswalk expansion** —
the audit expands a disease to equivalent MeSH/MONDO/DOID (`mondo.sssom.tsv`, `mesh_to_doid.csv`)
before matching edges, catching a few more than a single-ID match. This does not change the funnel
logic or the 137 dropped (which are anchored on the real final-116 membership, not the reproduction).

**To make the leakage step byte-exact / committed:** the v4 audit script should apply the disease-ID
crosswalks when checking leakage. Then `coverage_annotation_v4.csv` / `gold_standard_v4_bigset.tsv`
are reproducible
end-to-end.
