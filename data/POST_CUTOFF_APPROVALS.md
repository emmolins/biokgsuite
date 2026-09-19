# Post-cutoff approval set

Current set: `data/gold_standards/post_cutoff_approvals/post_cutoff_approvals.csv`, 128 drug-disease
pairs, built 2026-09-18 by `tools/build_post_cutoff_approvals.py` from
`repurposing-dataset` v0.6.0 (section 1 below). The earlier hand-assembled sets
and the audit that found them uncitable are in `data/gold_standards/legacy/`
(`README.md` and `HISTORY.md` there).

## 1. Source

`data/gold_standards/source/repurposing-dataset-v0.6.0/candidates_reviewed.csv`
(SHA-256 `bbb9fb8748247517…`, 1,288 rows). One row per regulator decision that
added an indication to an authorised medicine, 2023-01-01 to 2026-09-17, from
five regulators, each mapped to the closest thing it publishes to an FDA
efficacy supplement of sub-type "New Indication":

| regulator | what a row is | rows |
|---|---|---|
| FDA (CDER) | efficacy supplement, sub-type New Indication (calendar-year reports + Drugs@FDA), plus 2026 BLAs not yet sub-typed | 347 |
| EMA | CHMP positive opinion under an "extension of indication" heading (meeting highlights) | 333 |
| Health Canada | SNDS whose reason text mentions an indication | 246 |
| PMDA | approval whose English notes say "new indication" | 173 |
| TGA | registration of type EOI (extension of indications) | 189 |

Every row carries the regulator's document URL. Indication text comes from the
regulator's own record: the Health Canada NOC reason, the PMDA notes, the TGA
indication, the FDA approval letter ("provides for" sentence; 335 of 347
letters) and the EPAR procedural-steps entry for the CHMP opinion (204 of 333
cards matched). Drug and disease identifiers were resolved with the NCATS Name
Resolver 1.5.2 and Node Normalizer 2.4.1 (responses cached, so the mapping is
frozen with the data) and the MONDO SSSOM file at release 2026-09-01, then
reviewed in four passes: 2,309 recorded decisions with a rationale each, all
`status = proposed` pending a second reviewer. The full package (raw regulator
files, parsers, tests, byte-identical rebuild, `PROVENANCE.md` with the
retrieval timeline) is the `repurposing-dataset` zip; only the reviewed table
is vendored here.

**What the later review passes changed.** The first pass had
taken each regulator's own "new indication" category at face value for Health
Canada, PMDA and TGA. Read at the level of the active ingredient, most of those
rows are line extensions (a new line, stage, regimen or combination for a
disease the drug already treats), population extensions, or indications the
ingredient already held in another product or that another regulator had
approved before 2023. The second pass classified every candidate row against
the ingredient's prior indications, with the earlier approval named in each
rationale, and harmonised disease terms across regulators (one MONDO term per
disease, whichever regulator reported it). The third pass was an independent
re-read of every surviving pair, cross-checked against the FY2022 rows of the
PMDA list shipped with the package (which is how the September 2022 Japanese
approval of aflibercept for retinopathy of prematurity was found); it removed
ten more pairs, merged terms that had split one indication into two pairs
(gastric/GEJ adenocarcinoma, CLL/SLL), and left an explicit confirmation on
every kept row. The fourth pass resolved the 73 rows whose identifier had
been pending a live lookup: a second, independent model proposed identifier
and class for each, every proposal was checked against the MONDO label index
and the Node Normalizer (one proposed DrugBank accession turned out to be
vortioxetine, not dextromethorphan), and the missing Translator responses
were fetched live and cached; 13 pairs entered, 18 of the 73 rows were
vaccines, diagnostics or sealants, and 30 failed the ingredient rule. The set
went 272 → 125 → 115 → 128 pairs, and its composition changed: 68 pairs are
FDA-approved and 41 EMA-approved (57 and 25 of them first seen there), 63 were
approved by two or more of the five regulators, and no pair is a setting or
population change of an existing indication.
What remains is closer to what "repurposing" means for a knowledge graph: an
ingredient linked to a disease it had no approval for before 2023-01-01, in
any of the five jurisdictions.

**Two limits of the source that a reader should know.** (i) "Before 2023,
anywhere" means the five jurisdictions; approvals elsewhere are not consulted,
and China's NMPA in particular approved several tislelizumab indications
years before the EU did (four tislelizumab pairs are in the set on their EU
dates). (ii) `approval_date` is the earliest regulator row in the source, and
the source does not hold every approval: the FDA feed is CDER efficacy
supplements, so original NDAs/BLAs of a new product containing an old
ingredient (Syfovre, Iwilfin, Opfolda, Zoryve foam, Zepbound) and CBER
products (Breyanzi) are absent, some CDER supplements are missing from the
calendar-year reports (Nucala COPD, Vtama atopic dermatitis, the Tafinlar
low-grade glioma supplement of March 2023), one letter is not on Drugs@FDA
(TNKase stroke), and 129 CHMP cards have no matched EPAR text (Olumiant JIA,
Lynparza endometrial cancer among them). For those pairs the recorded date is
months later than the true first approval. Every such date is still inside
2023 to 2026, so check 1 is unaffected for the five dated graphs, but a graph
with a content date inside the window would need the dates completed first.

## 2. Construction funnel (the flow diagram, `07h_post_cutoff_funnel`)

Global stages, independent of any graph. Every count is re-derived by the
builder from the files on disk; `post_cutoff_funnel.csv` holds the table and
`post_cutoff_dropped.csv` names every row that left and why.

| stage | n | dropped | why |
|---|---|---|---|
| S0 regulator rows | 1,288 | | |
| S1 reviewed new-indication rows, both identifiers resolved | 282 | 1,006 | 309 line extensions; 203 population extensions; 114 indications the ingredient already held (another product, or approved before 2023 by another regulator); 161 no indication text yet (unmatched EMA cards, FDA supplements without a letter, bare "New indication" NOCs); 49 no disease target (diagnostic, procedural, symptomatic, sealant); 39 unclear; 131 other out-of-scope classes (wording change, formulation, veterinary, commitment fulfilment, biosimilar alignment, withdrawn, first approval of a new substance) |
| S2 KG-mappable rows | 256 | 26 | 9 vaccines (product-level ids; no graph models them as compounds), 8 drugs without a DrugBank accession, 6 tumour-agnostic approvals whose target is the root "cancer" term, 2 diagnostic imaging agents, 1 fixed combination |
| S3 unique drug-disease pairs | 128 | 128 | 268 (row, disease) events collapse to one (DrugBank id, MONDO id) pair each, keeping the earliest decision date and listing every regulator; 63 pairs were approved by more than one regulator |

Then, per graph, three checks applied to the 128:

| check | test | Hetionet | BioKG |
|---|---|---|---|
| C2 covered | drug node present (DrugBank id, degree ≥ 5) and disease node present (MONDO, DOID or MeSH, whichever the graph keys by; every identifier the MONDO term maps to is tried) | 2 | 66 |
| C1 prospective | approval after the graph's content date (`config/kg_snapshots.yaml`) and no direct drug-disease edge of any relation to any node the disease maps to | 2 | 58 |
| C3 repurposing | new MONDO term more than 2 `is_a` steps from every disease the graph already links the drug to | 2 | 57 |
| common support | evaluable in both graphs | 1 | |

`post_cutoff_pair_status.csv` has the pair × graph detail (which identifier
resolved the disease, the drug's degree, the leaked edge, the number of prior
indications the graph holds for the drug, the nearest prior MONDO term). It
also carries `grounded`, with `drug_n_mech` and `disease_n_mech`: whether the
drug node has at least one gene or pathway neighbour that a grounding dossier
would render, counted by the dossier code (`biokgsuite.grounding.KGDossiers`).
A drug can clear the degree floor on drug-drug or side-effect edges alone and
still have nothing a dossier can show. The flag is on record here and as the
`G` row of the funnel table; notebook 07 does not use it (its scorer needs
the node, not evidence text), and the grounding experiment folds it into its
own evaluability (`grounding_experiment/01_benchmark/`).
Hetionet's coverage is low because it holds 1,552 compounds, almost no
biologics, and 137 diseases; 106 of the 128 pairs name a drug it does not
have. BioKG's 8 leaks are the Japanese "public knowledge-based" approvals of
long-standing off-label uses (cyclophosphamide for GVHD, carboplatin for
endometrial cancer, methylprednisolone for Kawasaki disease, mycophenolate for
nephrotic syndrome, levofolinate for gastric cancer), iobenguane for
neuroblastoma, and two palmoplantar-pustulosis approvals; BioKG's 66,867
`DRUG_DISEASE_ASSOCIATION` edges are associations rather than indications, so
a "leak" there is a known association, not a prior approval, which is exactly
what check 1 is meant to remove. PrimeKG, DRKG, OpenBioLink and MATRIX have
not been checked yet (their files are on the cluster; run
`--extract-views <kg>` there, then `--kgs all --figure`).

## 3. The four checks, as implemented

**Check 1, prospective relative to each graph.** Two conditions, both per
graph. (a) The approval date is after the graph's content date in
`config/kg_snapshots.yaml`: Hetionet 2017-06-08 and BioKG 2020-06-09 are read
from the release records in this repo; DRKG, OpenBioLink and PrimeKG carry
literature dates flagged `verified: false`; MATRIX has none and must be set to
the release actually downloaded before it is scored. Every 2023 to 2026
approval clears the five dated graphs, so today the date test removes nothing,
but it is what stops a MATRIX release from 2025 being scored on 2024
approvals. (b) No direct edge of any relation between the two nodes
(`views/<kg>_edges.csv` keeps every drug-disease relation; only the indication
relation counts as a prior indication for check 3). `approval_date` is the
earliest regulator decision; for EMA-first pairs that is the CHMP
positive-opinion date, which precedes the Commission decision by about two
months (`date_type` column).

**Check 2, coverage reported, not hidden.** A pair is evaluable in a graph only
if both nodes exist there and the drug node has degree at least 5 (the floor
the notebook already used). Coverage is written to the funnel and the
checkpoint as a result in its own right, and `prospective_common_support: true`
in `config/params.yaml` scores every loaded graph on the pairs evaluable in all
of them, so a small small-molecule graph is not rewarded for the pairs it
cannot be asked about.

**Check 3, new disease versus line extension.** For each pair the drug's prior
indications are taken from the graph itself (its indication edges), mapped to
MONDO through DOID (`do_diseases.csv`), MeSH (`mesh_to_doid.csv`, SSSOM exact
matches) or directly, and the new term is compared on the MONDO `is_a`
hierarchy (`mondo_hierarchy.csv`, release 2026-08-04): within
`prospective_max_hops = 2` steps (same term, parent or child, grandparent,
grandchild or sibling) it is related, and the pair is a line extension from
that graph's point of view. The graph's own edges are the reference because
they are what the model saw. Reported by default; `prospective_strict: true`
makes it a filter. The rule can only fire when the graph holds prior
indications for the drug (BioKG: 30 of the 66 covered pairs; the rest have a
drug node with no indication edge, which `n_prior_indications` makes visible).
`setting_change_noted` in the gold file flags pairs whose review note records
a change of treatment setting (neoadjuvant, first-line, paediatric, …) for a
disease the drug may already treat; 12 pairs, reported, not filtered. A
graph-independent version of this check needs a dated prior-indication source
and is listed under next steps.

**Check 4, evaluation design.** Unchanged scorer (Adamic-Adar, no retraining)
and type-constrained negatives across all graphs; AUROC, MRR and Hits@K;
bootstrap intervals clustered by disease and now also by drug (trametinib,
tislelizumab and upadacitinib contribute 4 pairs each), both in the checkpoint, the wider to be reported;
`min_new_pairs` (10) is now a parameter rather than a constant in the notebook;
the checkpoint records `gold_standard_version`, `prospective_strict`,
`prospective_common_support` and the per-graph filter counts.

## 4. Files

| file | role |
|---|---|
| `tools/build_post_cutoff_approvals.py` | S0 to S3, checks 1 to 3, figure; docstring gives the usage |
| `config/kg_snapshots.yaml` | content date per graph, with source and `verified` flag |
| `data/gold_standards/post_cutoff_approvals/post_cutoff_approvals.csv` | the 128 pairs in the loader schema plus `mondo_id`, `mesh_id`, `agencies`, `n_agencies`, `date_type`, `evidence_ref` (source row ids), `evidence_url`, `setting_change_noted` |
| `data/gold_standards/post_cutoff_approvals/post_cutoff_pair_status.csv` | pair × graph checks |
| `data/gold_standards/post_cutoff_approvals/post_cutoff_funnel.csv`, `post_cutoff_dropped.csv` | funnel counts; every dropped row with its reason |
| `data/gold_standards/post_cutoff_approvals/PROVENANCE.txt` | source version and SHA-256, bridges, parameters, graphs checked |
| `results/checkpoints/07_kg_views/` | compact per-graph views (drug and disease nodes with degree; every drug-disease edge with its relation) written by `--extract-views` |
| `results/figures/07_generalization/07h_post_cutoff_funnel.{pdf,png}` | the funnel with the per-graph checks |
| `tools/make_post_cutoff_flow_figure.py`, `results/figures/07_generalization/07i_post_cutoff_flow.{pdf,png}` | the full flow: rows per regulator, indication-text sources, each cleaning step with its count, the pairs with their regulator composition; every number re-derived from the shipped tables |
| `tools/make_post_cutoff_flow_pptx.py`, `tools/post_cutoff_flow_pptx.js`, `results/figures/07_generalization/07i_post_cutoff_flow.pptx` | the same flow as an editable PowerPoint figure in Nature artwork style (183 mm double column, Arial, 5 to 7 pt text, 8 pt bold lowercase panel labels, 0.5 pt lines; native shapes and tables, no images); `07i_post_cutoff_flow_pptx_preview.jpg` is a render of it |
| `src/biokgsuite/gold/post_cutoff_approvals.py` | loader; now keeps a row that has a MONDO id but no DOID (19 rows), invisible to the DOID-keyed graphs rather than dropped |
| `tests/test_build_post_cutoff_approvals.py` | funnel accounting, distance rule, per-graph checks on synthetic inputs |

## 5. Next steps, in order of value

1. Second reviewer on `review/review_decisions.csv` in the source package
   (`status = accepted`). The second-pass class decisions rest on the
   reviewer's knowledge of each ingredient's approval history; every
   `ingredient_prior_indication` and `line_extension` rationale names the prior
   approval (regulator and month), so each can be checked against Drugs@FDA or
   the EPAR. Also the 19 pairs without a DOID.
2. Extract views for PrimeKG, DRKG, OpenBioLink and MATRIX on the cluster and
   fill `config/kg_snapshots.yaml` from the release records, then re-run
   `--kgs all --figure`.
3. Fill the remaining source-text gaps in the package (129 EMA cards without a
   matched variation, 10 FDA supplements without a letter, the bare "New
   indication" NOCs); these are the 161 S1 drops that are not yet classified
   either way. No identifier is pending any longer.
4. A graph-independent line-extension flag from a dated indication source
   (Open Targets or ChEMBL drug indications with first-approval year), which
   would replace the reviewer's memory of approval histories in step 1 with a
   citable table.

## 6. Version log

| built | source | pairs | note |
|---|---|---|---|
| 2026-09-17 | repurposing-dataset v0.2.0 | 272 | first reproducible set; PMDA, TGA and Health Canada only (FDA and EMA rows had no indication text) |
| 2026-09-17 | v0.3.0 (not shipped) | 325 | FDA letters and EMA EPAR text added; every FDA and matched EMA row classified; superseded the same day when the ingredient-level review began |
| 2026-09-18 | v0.4.0 (not shipped) | 125 | ingredient-level classification of all regulators; disease terms harmonised; tumour-agnostic targets dropped at S2 |
| 2026-09-18 | v0.5.0 | 115 | independent re-read of every pair with a PMDA FY2022 cross-check; ten pairs removed, split terms merged; every kept row confirmed; every mapped identifier tried in the per-graph checks |
| 2026-09-18 | v0.6.0 | 128 | the 73 rows with a pending identifier resolved (second-model proposals checked against the Node Normalizer; live Translator responses cached); 13 pairs added |
