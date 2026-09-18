# Verification of the EMA and PMDA candidate pairs

Written 2026-09-17. Closes the gap identified earlier: the FDA-attributed pairs
had been corroborated against Drugs@FDA, but 132 of the 253 candidates are
attributed only to EMA, PMDA or MHRA/Orbis and had never been checked against
any register. 78 of those sit inside the evaluated 116.

## Sources

| Agency | Source | Retrieved | Covers |
|---|---|---|---|
| EMA | Medicines register, `medicines-output-medicines-report_en.xlsx`, 2,339 human medicines | 2026-09-17 | current authorised indication text and MeSH therapeutic areas |
| PMDA | List of Approved Drugs, `000281190.pdf`, 1,889 records parsed | 2026-09-17 | Apr 2004 to Feb 2026, with New Approval / Partial Change flag and indication notes |

**One source was inspected and rejected.** EMA's "post-authorisation procedures"
table is the natural analogue of an FDA efficacy supplement, but it is a 140-row
snapshot of procedures currently in flight, with no historical rows and no
procedure scope text. It cannot corroborate anything approved in 2023 to 2025.

## What each check can and cannot establish

The EMA register reflects the **current** label, not variation history, so it
establishes that a drug is EU-authorised for a disease but not the year the
indication was added. EMA keeps that history only in per-product "procedural
steps" documents. PMDA is stronger: each record carries an approval date and a
Partial Change flag, so both the event and its date are checkable. Neither check
was allowed to borrow the other's strength; grades are reported separately.

## Result

| | all non-FDA pairs (172 agency-checks) | inside the 116 (78) |
|---|---|---|
| indication matched | 133 | 63 |
| partial match | 19 | 6 |
| drug and year only | 7 | 4 |
| drug outside window | 6 | 3 |
| drug found, disease not on label | 5 | 2 |
| drug not in register | 2 | 0 |

Every non-match was then adjudicated by hand against the primary record. Most
resolved as matcher artifacts rather than data faults: clinical abbreviations
("Ph+ CML-CP", "GVHD", "XLH") do not appear in spelled-out MeSH terms, EMA
renames products in place ("Onivyde pegylated liposomal (previously Onivyde)"),
and combination entries ("Opdivo+Yervoy") needed splitting. Matching the
ontology label carried by each pair, rather than its free-text indication,
resolved most of them.

## Adjudication: five removals, not ten

An earlier pass flagged ten. Re-checking each one adversarially rather than
confirming it reduced that to five. Two flags were wrong and three were
misclassified. The corrections are recorded here rather than quietly dropped.

### Confirmed removals (5)

Each verified against the primary record, and each checked for an in-window
record that would rescue it.

| Pair | Drug | Claimed | Evidence |
|---|---|---|---|
| 95 | Tremfya | Palmoplantar pustulosis (PMDA) | PMDA partial change **21 Nov 2018**. Tremfya's only in-window PMDA records are Crohn's disease (24 Jun 2025) and ulcerative colitis (27 Mar 2025, 19 Feb 2026). No in-window record for this indication |
| 96 | Lenvatinib | Unresectable thymic carcinoma (PMDA) | PMDA partial change, Lenvima, **23 Mar 2021**. No lenvatinib record of any kind after 2022 |
| 252 | Fibrinogen HT | Acquired hypofibrinogenemia in cardiovascular surgery | PMDA **28 Mar 2022**, and for obstetrical haemorrhage, not cardiovascular surgery. No fibrinogen record after 2022 |
| 158 | Omalizumab | EGPA | Xolair's Japanese records are 2009, 2013, chronic idiopathic urticaria (2017), seasonal allergic rhinitis (2019). Japanese EGPA approvals belong to mepolizumab (FY2018) and benralizumab (27 Dec 2024). No omalizumab EGPA record |
| 238 | Privigen | Measles post-exposure prophylaxis | Full 1,019-character EMA indication text read end to end: PID, hypogammaglobulinaemia, post-HSCT, congenital AIDS, ITP, Guillain-Barre, Kawasaki, CIDP. No measles. Pair is EMA-only |

Three of these five break the time-split premise directly: a 2018, 2021 or 2022
approval cannot postdate a graph snapshot.

### Retracted (2)

Both were flagged because the agency I checked does not support the claim. Both
are multi-agency rows, and the other agency does support it. **A gap at one
agency does not invalidate a pair attributed to two.**

| Pair | Drug | Why the flag was wrong |
|---|---|---|
| 119 | Welireg | Pair is **FDA, PMDA**. FDA NDA215383 SUPPL-12, **14 May 2025**, belzutifan for pheochromocytoma/paraganglioma. The PMDA finding stands (Japan's June 2025 approval is a new active ingredient for RCC and VHL tumours) but only the PMDA attribution is wrong, not the pair |
| 149 | Uplizna | Pair is **FDA, EMA**. FDA BLA761142 SUPPL-3, **3 Apr 2025**, Efficacy-New Indication, inebilizumab for IgG4-related disease. EMA has not authorised it, so again only the EMA attribution is wrong |

### Reclassified as a policy decision, not an error (3)

Cosentyx IV formulation (75), Spikevax strain update (78) and Wegovy oral tablet
(213) are genuinely not new indications. But they already carry
`StrictRepurposing = False` in the candidate file, so the dataset had already
identified them. They are not isolated errors; they are three instances of a
larger finding below, and picking them off individually would leave the rest.

## The larger finding: the strictness filter was never applied

**15 of the 110 candidate rows that map into the 116 carry
`StrictRepurposing = False`.**

| | |
|---|---|
| Keytruda 11 | MSI-H/dMMR solid tumours, full approval conversion |
| Veklury 24 | COVID-19 in severe renal impairment |
| Zoryve 51 | psoriasis in children 6 to 11 |
| Leqvio 54 | high LDL-C, expanded risk population |
| Praluent 63 | heterozygous familial hypercholesterolaemia |
| Evrysdi 71 | SMA paediatric expansion |
| Cosentyx 75 | IV formulation |
| Spikevax 78 | strain update |
| Vemlidy 79 | chronic hepatitis B, paediatric |
| HyQvia 81 | primary immunodeficiency, paediatric |
| Upadacitinib 163 | paediatric atopic dermatitis |
| Tzield 206 | delay T1D onset in children |
| Asceniv 207 | primary humoral immunodeficiency, paediatric |
| Palynziq 211 | PKU expanded to 12 years and older |
| Wegovy 213 | oral tablet formulation |

The construction figure labels the 234 as "strict repurposing, time-split", and
the 253 to 234 step as removing pairs that are "not strict repurposing". The
reconstructed filter (section 4 of the Methods record) tests only for missing
identifiers and placeholder text. It never reads `StrictRepurposing`. These 15
are the consequence.

Most are paediatric or population expansions: the same disease in a wider group.
Under a new-disease definition of repurposing they do not qualify. This is a
definitional decision to make once, deliberately, for all 15, not a defect to
patch three rows at a time.

## A false negative in the FDA census, found on the way

Belzutifan's pheochromocytoma approval of 14 May 2025 is tagged by FDA as
`Efficacy-Labeling Change With Clinical Data`, not `Efficacy-New Indication`.
That category holds 274 rows in the window and the census rule excludes it. It
therefore contains at least some genuine new indications, so the census as built
under-counts. `FDA_CENSUS_METHOD.md` records that exclusion with the reason
"new evidence, existing claim", which is now known to be true only in general,
not always.

## Scripts

| Script | Role |
|---|---|
| `tools/corroborate_ema.py` | EMA medicines register check, ontology-label matching |
| `tools/parse_pmda.py` | parses the PMDA PDF into 1,889 dated records |
| `tools/corroborate_pmda.py` | PMDA check with partial-change and indication matching |

Outputs: `results/tables/corroborate_ema.csv`, `corroborate_pmda.csv`,
`corroborate_ema_pmda.csv` (consolidated, with adjudications),
`data/gold_standards/pmda_approvals.csv`.
