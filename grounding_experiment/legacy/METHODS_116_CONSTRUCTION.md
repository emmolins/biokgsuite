# Construction of the prospective repurposing benchmark (253 to 116)

> **Scope note, 2026-09-17.** This document covers the 116-pair evaluation set
> used by the LLM grounding experiment. Notebook 07's prospective axis no longer
> uses it: it reads the regulator-derived set described in
> `data/POST_CUTOFF_APPROVALS.md` (128 pairs, `tools/build_post_cutoff_approvals.py`),
> whose funnel reproduces from its inputs.

Materials and Methods record. Written 2026-09-17. Every count below was
re-derived by re-running the construction, not copied from a previous figure.
Where a published number could not be reproduced, that is stated rather than
smoothed over.

---

## 1. Candidate set (253 pairs)

Unit of analysis is a **drug-disease pair**: a drug already approved before the
evaluation window that subsequently gained approval for a new indication.

`data/gold_standards/legacy/pairs_with_ids_v2.csv`, 253 rows, 17 columns. This is a
single candidate set covering 2023 to 2026 under one set of inclusion criteria,
not a merger of separately defined datasets.

**Inclusion criteria**, applied uniformly across the window:

1. the drug held marketing approval before 2023 (`Drug_Pre2023 = True`), so the
   approval under study is an indication expansion rather than a first approval;
2. the new indication was approved within 2023 to 2026 by at least one of FDA,
   EMA, PMDA or MHRA/Orbis;
3. the approval granted a genuine new indication, not a formulation, dosage form
   or presentation change.

**Sources consulted**, the same set throughout: agency approval records (FDA
novel approvals and efficacy supplements, EMA type-II variations, PMDA approval
lists), cross-referenced against Drugs@FDA, the EMA medicines register and PMDA
listings.

Assembly was incremental. The 2024 to 2026 rows were compiled first and the 2023
rows were added in a later pass, against the same criteria and the same sources.
This is a fact about build order, not a difference in definition, and the set is
treated as one group throughout. `Pair_ID` runs 1 to 253 across the whole set.

Four cross-year duplicates are retained deliberately, because each is a separate
regulatory event for the same drug and indication: Lynparza (19, 125), Reblozyl
(29, 140), Bosulif (30, 135), Vabysmo (38, 188).

**Assembly was LLM-assisted (Biomni)**, as recorded in
`PROVENANCE_253_candidate_pairs.md`. This applies to the candidate set as a
whole and must be disclosed. Section 7 gives the independent corroboration
against the primary regulatory record that addresses it.

`YearSet` labels which pass a row entered in, `2023` or `2024-2026`. It is a
build-order label, not an approval date, and carries no methodological meaning.
`Agencies` distribution: FDA only 93, EMA only 93, PMDA only 25, remainder
multi-agency. 121 of the 253 are FDA-attributed.

## 2. Drug identifier resolution

174 unique active ingredients after splitting combinations on `+` and `/`.
Resolved by a documented cascade: PubChem PUG-REST synonyms (82), Wikidata
SPARQL exact or alternate label (66, mainly biologics), Wikidata alternate,
brand and salt-stripped names (15), mychem.info brand and research codes (4).

Result: **167/174 (96.0%)** ingredients carry a DrugBank ID; 101/174 (58.0%)
carry a PubChem CID. Seven are permanently unresolved (three CAR-T products, one
vaccine, one Japan-only SGLT2 inhibitor, two recent EU monoclonals). Eighteen
combination products carry comma-separated per-ingredient IDs.

## 3. Disease identifier resolution and crosswalks

Each pair carries `MONDO_ID`, `MeSH_ID` and `DOID_ID` where resolvable. Vague
labels were resolved per drug, agency and date; 11 entries labelled
"New indication extension (EMA 2026, TBD)" remain provisional because EMA had
not published details.

For matching, an identifier is expanded across vocabularies so that an edge
keyed under an equivalent ID is still detected:

| Crosswalk | File | Version |
|---|---|---|
| MONDO to DOID and MeSH | `mondo.sssom.tsv` | MONDO releases/2026-08-04 |
| MeSH to DOID | `mesh_to_doid.csv` | built 2026-09-14, `tools/build_crosswalks.py` |
| MONDO is_a hierarchy | `mondo_hierarchy.csv` | MONDO releases/2026-08-04, source SHA-256 `fce4bccd97c4eb66…` |

## 4. Filter 1: harness-ready (253 to 234, minus 19)

A candidate is dropped if **any** of the following holds:

1. no DrugBank ID, or
2. no disease identifier in any of MONDO, DOID or MeSH, or
3. the indication text is a placeholder, matching `tbd`,
   `population expansion`, or `new indication (`.

**Reproduces exactly: 19 dropped, 234 retained.** Dropped rows are enumerated in
`results/tables/build_116_dropped_step1_unresolved.csv`.

## 5. Filter 2: answer leakage

A pair is dropped if a direct drug-to-disease edge already exists in any
evaluated knowledge graph, since the model could then retrieve the answer rather
than infer it. Disease identifiers are crosswalk-expanded first (section 3).

Edge definitions, per graph:

| Graph | Edge tested | Disease keyed by |
|---|---|---|
| PrimeKG | `drug`-`disease` rows, either orientation | MONDO |
| BioKG | `DRUG_DISEASE_ASSOCIATION` | MeSH |
| DRKG | `Compound::` to `Disease::` | MeSH or DOID |
| Hetionet | `CtD`, `CpD` | DOID |
| OpenBioLink | `DIS_DRUG`, PubChem CID mapped to DrugBank | DOID |

Leakage hits among the 234, by graph: PrimeKG 41, BioKG 49, DRKG 44,
Hetionet 2, OpenBioLink 14. **Union of all five is 84, identical to the union of
PrimeKG, BioKG and DRKG alone**, so Hetionet and OpenBioLink contribute no
unique exclusions and describing the filter as three graphs is accurate.

## 6. Filter 3: knowledge-graph coverage

A pair is evaluable only if the drug and the disease both resolve within the
same graph, for at least one graph. Presence is taken from the validated
resolver in `reconstruct_116.py`
(`results/tables/reconstruct_116_audit.csv`, field `n_present`).

Coverage of the evaluable set, per graph: PrimeKG 111, DRKG 90, BioKG 78, with
72 covered by all three. The 72 are the common support on which evidence quality
is comparable across graphs.

## 7. Independent corroboration of the candidate set

Because assembly was LLM-assisted, every FDA-attributed pair was checked
against the FDA indication-expansion census built from the Drugs@FDA relational
data files (release 2026-09-16; see `FDA_CENSUS_METHOD.md`). Script:
`tools/corroborate_253.py`.

| Outcome | n |
|---|---|
| Directly corroborated in the census | 98 |
| Real FDA efficacy event in the window, different action type | 18 |
| Real FDA action in the window, not classed efficacy | 3 |
| CBER-regulated, outside the Drugs@FDA frame | 2 |
| **Unsupported by any FDA record** | **0** |

No fabricated pair was found anywhere in the set. 119 of the 121 FDA-attributed
pairs have a documented FDA regulatory action behind them; the remaining two are
CBER-regulated products, Yescarta (axicabtagene ciloleucel, a CAR-T cell
therapy) and Asceniv (an immune globulin), which are absent
from Drugs@FDA altogether, so they are outside the reference frame rather than
unsupported. Corroboration does not vary by which pass a row entered in: of the
55 FDA pairs labelled 2023 and the 66 labelled 2024-2026, neither group contains
an unsupported pair and each contains exactly one out-of-frame product.

The 18 in row two are genuine events of a different kind, chiefly
paediatric or patient-population expansions (Evkeeza, Rozlytrek, Zoryve, Adbry,
Cresemba) and one route change (Entyvio subcutaneous); they are the same disease
in a new population, which is why the `StrictRepurposing` flag matters.

EMA, PMDA and MHRA/Orbis pairs have not yet been corroborated against an
authoritative register.

## 8. Evaluation set and query construction

`data/gold_standards/legacy/gold_standard_v4_bigset.tsv`, 116 rows.

One anonymised ranking query per pair: the approved drug plus 7 distractors
drawn from Open Targets (Platform release 26.06), across 74 distinct diseases.
Four arms: no-KG plus one per graph. Three seeds by two presentation orders.

116 queries x 4 arms x 3 seeds x 2 orders = **2,784 calls per model**.

**Field caveat.** The column named `approval_date` in this file holds the
**drug's original approval year**, ranging 1983 to 2022, not the date of the new
indication. The new-indication window is carried by `year_set`. The column name
invites misreading and should be renamed or defined explicitly in the legend.

## 9. Reproduction status

Re-run on 2026-09-17 via `scripts/build_116.py` against the graph snapshots in
the `biokgsuite` working tree.

| Step | Published | Re-derived | Status |
|---|---|---|---|
| 253 to 234 | -19 | -19 | reproduces exactly |
| 234 to 143 (leakage) | -91 | -84 | does not reproduce |
| 143 to 116 (coverage) | -27 | -43 | does not reproduce |
| final evaluable set | 116 | 107 | see below |

Of the 116, **110 map back onto a row of the 234** on DrugBank ID plus either
DOID or disease name; 6 do not, owing to identifier drift between files. Against
those 110, the re-derived rule (not leaked AND covered) yields **107 true
positives, 0 false positives and 3 false negatives**. The three are Praluent
(alirocumab, heterozygous familial hypercholesterolaemia), Keytruda
(pembrolizumab, neoadjuvant and adjuvant HNSCC) and Noxafil (posaconazole,
paediatric invasive fungal infection), each flagged as leaked by the re-derived
test but present in the shipped 116.

**The published intermediate split cannot be correct as drawn.** 44 of the 234
pairs have no coverage in any graph, and a pair absent from every graph cannot
carry a leakage edge, so the coverage filter must remove approximately 44, not
27. An exhaustive search over all 31 graph subsets and both filter orderings
found no combination producing -91 then -27.

### Consequence for the manuscript

The 116 evaluable set itself is a stored artifact and is unaffected. Nothing in
the grounding experiment depends on the intermediate counts. Two defensible
options:

1. **Report the funnel as 253 to 234 to 116**, stating the leakage and coverage
   filters jointly as the criterion for evaluability, and giving 118 as the
   number removed. Every number in that version reproduces.
2. **Re-derive the 116 from the documented rule** and report the resulting
   split (-84, then -43, giving 107), noting the 3 pairs and 6 unmapped rows
   above as differences from the shipped set. This requires re-running the
   grounding experiment on the revised set and is not recommended unless a
   reviewer demands strict end-to-end reproducibility.

Option 1 is recommended. It removes an unverifiable claim without changing any
result.

## 10. Provenance gap to close

`data/PROVENANCE.md` in this repository pins releases and SHA-256 digests for
Hetionet and BioKG but lists PrimeKG, DRKG and OpenBioLink as "not yet
downloaded". The 116 was built against snapshots of all five that live in the
`biokgsuite` working tree. Those three snapshots need release identifiers,
retrieval timestamps and digests recorded before submission, otherwise sections
5 and 6 cannot be independently repeated.

## 11. Scripts

| Script | Role |
|---|---|
| `scripts/build_116.py` (in `biokgsuite`) | re-derives the funnel; the copy at `grounding/01_construction/build_116.py` has a broken `ROOT` and resolves to `grounding/` rather than the repository root |
| `reconstruct_116.py` | per-graph presence audit, produces `reconstruct_116_audit.csv` |
| `tools/corroborate_253.py` | checks FDA-attributed pairs against the Drugs@FDA census |
| `tools/build_fda_census.py` | builds the FDA census used in section 7 |
