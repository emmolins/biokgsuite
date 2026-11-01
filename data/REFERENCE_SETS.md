# Reference sets — complete specification

Every reference set (the coverage dimension's gold standards) in `data/reference_sets/`, what it means, where it
comes from, and whether it is independent of the graphs being scored.

**Design principle.** Coverage asks whether the entities and relations needed
for *drug repurposing* are present. Reference sets are therefore defined by the
task, not by the size of an ontology. The definitions are asymmetric on purpose:
all diseases are in scope (repurposing exists to serve undrugged diseases),
but only approved drugs are (you can only repurpose a drug that exists).

**Circularity.** A gold standard drawn from a database that a graph ingests
measures ingestion recency, not quality. Every source below is marked
independent or contaminated. The 2026-09 revision moves the drug, gene and
drug-target standards off DrugBank for exactly this reason.

---

## Entity gold standards

| Entity | Rule | Source | Release | Independent? | In score? |
|---|---|---|---|---|---|
| **Drug** | approved therapeutics, one row per active molecule: `maximumClinicalStage == 'APPROVAL'`, with at least one APPROVAL-stage indication on the molecule or its ChEMBL parent, not a vaccine component, collapsed to the ChEMBL parent | OT `drug_molecule` + `clinical_indication` | 26.06 | **Yes**, ChEMBL-derived | scored |
| **Disease** | every human disease, incl. undrugged: MONDO non-obsolete, non-grouping, descendant of `MONDO:0700096` (human disease or disorder) | MONDO `mondo.obo`, `is_a` only | v2026-08-04 | **Yes** | scored |
| **Gene/Protein** | the genes the task touches: targets of the approved drugs, plus genes with curated disease evidence (relation #3) | OT `drug_mechanism_of_action` + the seven `evidence_*` datasets | 26.06 | **Yes** | scored |
| **Pathway** | Reactome pathways containing any gene of that set | OT `target.pathways` (no separate Reactome download needed) | 26.06 | **Yes** | *reported, not scored* |

Sizes on 26.06 and MONDO v2026-08-04: 3,793 drugs (from 5,485 approved
molecule rows), 22,276 diseases, 6,210 genes (1,051 targets and 5,775 disease
genes), 2,217 pathways. `data/reference_sets/MANIFEST.json` carries the rule
and the count for every set, written by the builder.

Disease is stratified, and the strata are the interesting result:

| Stratum | n | from |
|---|---|---|
| all (the denominator) | 22,276 | MONDO human disease, non-obsolete, non-grouping |
| rare-tagged | 15,634 | MONDO `subset: rare` |
| with >=1 approved indication | 928 | OT `clinical_indication`, APPROVAL stage |
| **no known indication** | 21,348 | complement; the repurposing frontier |

So 95.8% of the human disease space has no approved drug, and 15,219 of the
undrugged diseases are rare-tagged.

### Why each drug rule is there

An APPROVAL-stage molecule with no indication record is, in Open Targets,
almost always a diagnostic or an excipient (technetium imaging agents, soybean
oil, polyethylene glycol, silver): 820 of the 5,485 rows. A repurposing
candidate is a therapeutic with a known use, so the rule asks for one. Vaccine
components (55) are excluded because the post-cutoff approval set excludes
vaccines, and one rule should decide what a drug is in both parts of the
paper. The parent collapse makes the unit the active molecule: ChEMBL gives
doxazosin and doxazosin mesylate separate ids, DrugBank gives them one, and
counting both in a denominator a DrugBank-keyed graph can only ever match once
would understate every such graph. The resolver expands a graph's drug to its
ChEMBL family for the same reason, on the other side of the intersection.

### Why the disease rule is `is_a` under human disease

MONDO v2026-08 carries a veterinary branch (`MONDO:0005583` non-human animal
disease: "wry nose, giraffe", "bullous pemphigoid, pig", 8,319 live terms), 491
OMIM susceptibility loci ("psoriasis 4, susceptibility to"), 41 disease
characteristics ("idiopathic", "iatrogenic") and 18 injuries, all non-obsolete
and non-grouping. None is a target for drug repurposing. Walking `is_a` from
`MONDO:0700096` keeps the 22,276 human diseases and drops the 8,871 others.
The walk uses the asserted `is_a` edges only: the reasoned hierarchy in
`data/gold_standards/bridges/mondo_hierarchy.csv` (built for notebook 07's
therapeutic-area classifier) includes logical-definition edges that place the
animal branch under human disease, and is not used here.

### Why the gene set is the union, and why it still saturates

Targets of approved drugs alone (1,051 genes) are held by every gene-centric
graph (PrimeKG 99.8%, BioKG 100%), so a third of entity coverage would be a
constant. Adding the genes with curated disease evidence gives 6,210, the set
the task actually touches. It still saturates for graphs that ingest the human
gene set (PrimeKG 99.7%, BioKG 96.0%); that is a property of the graphs, not
of the set, and Methods should say that entity coverage of genes separates
graphs that carry genes from those that do not, while relation coverage
(disease-gene: PrimeKG 53%, BioKG 28%, Hetionet 2%) is where the signal is.

### What changed, and why
- **Drug** was DrugBank approved **+ investigational** (10,909), then every
  approved Open Targets molecule (5,485). Investigational drugs are not
  repurposing candidates, DrugBank is ingested by all six graphs, and the
  5,485 counted salts twice and diagnostics at all. Now 3,793 approved
  therapeutics, ChEMBL-derived, one per active molecule.
- **Disease** was DOID (11,932), then all of MONDO (31,147). DOID is reachable
  from only 35% of Orphanet and 52% of OMIM/UMLS, which silently penalised
  graphs that do not speak DOID and capped rare-disease measurement; all of
  MONDO included 8,871 terms that are not human diseases. Now 22,276.
- **Gene** was the whole UniProt human proteome, then the 1,076 targets of
  approved drugs. Neither is the set the task touches. Now 6,210.
- **Pathway** was all of Reactome, then the pathways of the targets. Now the
  pathways of the gene set, 2,217.

---

## Relation gold standards

Six relations, chosen because they are the edge types the grounding dossier
actually reads (`src/biokgsuite/grounding/dossiers.py`, `KGDossiers._triples`:
`allowed = ('gene', 'pathway') if role == 'drug' else ('gene', 'pathway', 'disease')`),
plus drug-disease which the link-prediction tasks in nb06/nb07 require.

| # | Relation | Gold standard | Source | Release | Independent? | In score? |
|---|---|---|---|---|---|---|
| 1 | **drug → gene** | drug and its molecular target | OT `drug_mechanism_of_action` | 26.06 | **Yes** — ChEMBL | scored |
| 2 | **drug → pathway** | *derived*: drug→gene ∘ gene→pathway | OT `drug_mechanism_of_action` + Reactome | 26.06 | **Yes** | *reported, not scored* |
| 3 | **disease → gene** | gene implicated in disease | union of 7 curated `evidence_*` datasets | 26.06 | **Yes** | scored |
| 4 | **disease → pathway** | *derived*: disease→gene ∘ gene→pathway | as #3 + Reactome | 26.06 | **Yes** | *reported, not scored* |
| 5 | **disease → disease** | ontological subtype (`is_a`) | MONDO `mondo.obo` | v2026-08-04 | **Yes** | scored |
| 6 | **drug → disease** | approved indication | OT `clinical_indication` | 26.06 | **Yes** | scored |

`drug → drug` (DDI) is **dropped**: the dossier never reads it, it is already
outside the dimension score, and its only comprehensive source (DrugBank) is
ingested by the graphs being scored.

---

## Why pathways are not scored

Pathway entities, drug→pathway and disease→pathway are built, computed, plotted
and written to the checkpoint exactly like everything else, but they are held
out of the dimension score. This is a deliberate measurement decision and it
belongs in Methods, not in a footnote.

### The argument

A pathway is not the same kind of object as the other three. A drug is a
molecule, a gene is a locus, a disease is an ontology concept. A pathway is a
curator's named grouping of genes, and where its boundaries fall is a judgement
call: Hetionet's pathway nodes have a median of 23 genes and a maximum of 1,956.
Two curators can partition the same biology differently and both be right.

The six graphs do not share a pathway vocabulary. Reactome, KEGG, SMPDB and
PathwayCommons all appear across them. Matching on accession therefore answers
"did this graph ingest Reactome?" rather than "does this graph represent pathway
biology?", and only the second question is a quality claim.

### The evidence

Hetionet is the clean case, and the notebook prints it:

| | |
|---|---|
| pathway nodes in the graph | 1,822 |
| carrying a Reactome `R-HSA-*` accession | 0 |
| pathway resolvability | 0.0000 (top unresolved: `PC7`, 1,528) |
| distinct pathway names | 1,808 |
| of those, exactly matching a name in Open Targets' Reactome subset | **873 (48.3%)** |

Nearly half of Hetionet's distinct pathway *names* are exact string matches to
Reactome pathway names, because PathwayCommons re-ingested Reactome and reissued
it under PC7 accessions. The content agrees and the identifiers do not. Scored on
accession, Hetionet would take a 0% on a dimension where it in fact holds
substantially the same pathway biology as the reference set. That is a
vocabulary artefact being recorded as a quality deficit.

48.3% is a floor, not the true agreement: the reference names come from Open
Targets' Reactome subset, which carries only the 2,298 Reactome pathways
containing a human drug target. A one-off check against PathwayCommons v7
itself, the release Hetionet drew from, put the figure near three quarters
(see below). The floor is enough to make the point, and it is the number the
notebook computes, so it is the number quoted here.

### What is scored instead

| | scored | computed and stored, not scored |
|---|---|---|
| Entities | Drug, Disease, Gene/Protein | Pathway |
| Relations | drug→gene, disease→gene, disease→disease, drug→disease | drug→pathway, disease→pathway |

Dropping the two pathway relations removes 94,720 of 167,493 reference pairs
(57% by count), which sounds severe and is not: on Hetionet both scored 0.00%
before the eligibility fix below and N/A after, so they contributed no
information to the mean either way. Hetionet's relation sub-score moves from
0.0599 to 0.0898 once the artefactual zeros are removed.

### Eligibility is applied consistently (fixed 2026-09-14)

Entity coverage always gated pathways on `pathway_ontology: reactome` in
`config/kgs.yaml`, so an undeclared KG showed N/A. The two pathway *relations*
were not gated and reported a hard 0%. A KG could therefore be marked N/A for
pathway entities and 0% for pathway relations in the same run. Both now use the
one gate, so an undeclared KG reads N/A throughout.

### What this costs, stated honestly

1. **The design justification weakens.** The six relations were chosen because
   they are what the grounding dossier reads (`src/biokgsuite/grounding/dossiers.py`,
   `allowed = ('gene', 'pathway')`). Scoring four while the method consumes five
   needs a sentence of explanation rather than a one-line principle.
2. **Mechanistic reasoning becomes unmeasured.** Pathways are the bridge when no
   direct disease→gene edge exists. The benchmark no longer says whether a graph
   supports indirect mechanistic inference, arguably the reason to use a KG at
   all rather than a lookup table.

Both are real. Neither is fixed by scoring a vocabulary check.

### The route back, if a reviewer asks for one

Match a KG pathway to a Reactome pathway when their **gene sets** overlap above
a Jaccard threshold, rather than when their accessions agree. That makes pathway
coverage vocabulary-independent and measures the thing the dimension claims to
measure. It costs one build step, a pathway→gene-set crosswalk per KG, and it is
the only option under which a graph storing genuine pathway structure in a
non-Reactome namespace can score above zero.

Name matching was considered and rejected as a *scoring* rule: most of
Hetionet's pathway names match Reactome *because* PathwayCommons re-ingested
Reactome, so as a metric it measures provenance lineage rather than independent
content agreement.

**One-off check against PathwayCommons v7 (2026-09-14).** Not part of the
pipeline, and not re-derivable from this repository: the PC7 GMT was downloaded
once, queried, and not retained, because nothing scored depends on it. Of
Hetionet's 1,822 pathways, 1,535 (84.2%) matched a PC7 pathway name exactly,
1,410 (77.4%) matched one by gene set at Jaccard >= 0.9, and the best match was
Reactome-sourced for 1,296 (71.1%); the remainder were mostly NCI-PID (221),
with small numbers from CTD, HumanCyc, KEGG and PANTHER.

This is background for the exclusion, not evidence behind a reported score.
Pathways do not enter the dimension score, so the decision rests on the
vocabulary argument above, which the notebook's own 48.3% already demonstrates.
Reproducing the PC7 figures would mean re-downloading
`https://download.baderlab.org/PathwayCommons/PC2/v7/`; they should be cited as
a one-off check, or dropped, rather than presented as a pipeline output.

Every figure in this section other than that check is regenerated by the
notebook, not asserted here.
The name comparison uses `data/reference_sets/entity_pathway_names.tsv`
(2,298 Reactome id/name pairs from `ot_target.pathways`, written by
`tools/build_reference_sets.py`), matched case-insensitively on stripped
strings. The diagnostic table is printed by the "Pathways: reported, not scored"
cell in `biokgsuite_dimensions/01_coverage/notebook.ipynb` and stored in the
checkpoint under `pathway_diagnostic`.

Everything needed to reverse this decision is retained. `entity_pathway.txt`,
`relation_drug_pathway.tsv` and `relation_disease_pathway.tsv` are still built
and still version-stamped, and the checkpoint carries `pathway_diagnostic`,
`scored_entity_types`, `scored_relations` and `reported_not_scored`. Re-including
pathways is a change to `SCORED_*` in the notebook, not a rebuild.

### gene → pathway is a building block, not a scored relation
The dossier reads drug→pathway and disease→pathway as *direct* edges, so
gene→pathway is not scored. But it is needed to CONSTRUCT golds #2 and #4.
Source: Reactome `UniProt2Reactome.txt` / `NCBI2Reactome.txt`
(https://reactome.org/download/current/), or the pathway annotations carried in
OT `target` — check the schema printed by `build_reference_sets.py` first; if
`target` carries them, no separate Reactome download is needed.

### Three decisions that must be stated in Methods

**direct, not indirect associations (#3).** Open Targets propagates evidence
through the disease ontology: an *indirect* association attributes a
descendant's evidence to the parent term. Using indirect would inflate gene
counts for general terms and systematically favour graphs holding
coarse-grained diseases. Use `..._direct`.

**curated evidence, not the overall association score (#3).** OT aggregates
genetics, somatics, expression, literature mining and more into a single score.
Literature-mined edges are noisy and expression edges are not causal, so the
overall score is not a defensible "this gene is implicated in this disease"
claim. The filter is applied at **download** time, not in the builder:
`tools/fetch_opentargets.sh` takes `GENE_DISEASE=curated|clinvar|genetic|all|none`
and only fetches the chosen `evidence_*` datasets;
`tools/build_reference_sets.py` then unions whatever `ot_evidence_*` directories
are present, with no `datatypeId` filter of its own. The current build used the
seven curated sources: `clingen`, `gene2phenotype`, `gene_burden`,
`genomics_england`, `orphanet`, `uniprot_literature`, `uniprot_variants`.
`data/reference_sets/MANIFEST.json` and the builder's per-source printout record
which ones actually contributed. Report an overall-score threshold as a
sensitivity check.

Note the datatype id in 26.06 is `genetic_literature`, not `genetic_association`;
earlier drafts of this file had the wrong name.

**#2 and #4 are derived, and the derivation is a modelling choice.** Their gold
is built by composing through genes. A KG may carry curated drug→pathway edges
that do not come from its drug's targets, in which case the gold and the graph
mean slightly different things. Disclose the composition rule. There is no
clean direct alternative: KEGG DRUG pathway maps are licence-restricted for
bulk use.

**#5 is subtype, not comorbidity.** MONDO `is_a` gives ontological
parent/child. Comorbidity would need EHR or claims data and is a different
relation. Hetionet's `DrD` mixes resemblance semantics — say which you mean.

## Getting the newest version

From the repository root, on a network that can reach ftp.ebi.ac.uk:

```bash

# 1. Open Targets 26.06 — all six relation + entity sources (~1.2 GB;
#    CORE_ONLY=1 fetches just the ~110 MB core, without target-disease associations)
bash tools/fetch_opentargets.sh

# 2. MONDO + Disease Ontology, pinned by release tag (GitHub, works anywhere)
python tools/build_crosswalks.py

# 3. UniChem PubChem<->ChEMBL crosswalk (~30 MB)
#    REQUIRED if you score OpenBioLink: Open Targets' drug crossReferences has
#    no PubChem source, so without this every PubChem CID resolves to nothing
#    and OpenBioLink is scored as holding no approved drugs.
bash tools/fetch_unichem.sh

# 4. DDInter 2.0, only if you keep drug_drug — manual, no stable API
#    https://ddinter.scbdd.com/download/  ->  data/crosswalks/ddinter/

# 5. Build the reference sets from whatever is present
python tools/build_reference_sets.py
```

### Namespace coverage of the resolver
`IdResolver` dispatches on the identifier's namespace, so a missing crosswalk
shows up as a whole entity type resolving to nothing rather than as a partial
loss. nb01 raises this loudly (`RESOLVABILITY_FLOOR`, 5%) instead of letting it
pass as a coverage result.

| Namespace | Reaches | Via | Present? |
|---|---|---|---|
| ChEMBL, DrugBank | ChEMBL | OT `crossReferences` | yes |
| **PubChem CID** | ChEMBL | **UniChem `src1src22`** | **only after step 3** |
| UniProt, Entrez, HGNC, symbol | Ensembl | HGNC + OT `target.proteinIds` | yes |
| MONDO, DOID, MeSH, OMIM, Orphanet | MONDO | MONDO SSSOM | yes |
| Reactome | Reactome | native | yes |
| KEGG, SMPDB, PathwayCommons PC7 | — | none (by design) | pathways are not scored |

Steps 1 and 2 each write a release file and SHA-256 checksums, so every
reported number traces to a named, re-fetchable version. To move to a newer
Open Targets release, set `OT_RELEASE=26.09` and rerun; the checksums change
and the diff shows it.

### Pinned release identifiers
| Source | Release | Notes |
|---|---|---|
| Open Targets Platform | **26.06** (2026-06-24) | dataset names changed from earlier releases: there is no `known_drug` and no `drug_indication` |
| MONDO | **v2026-08-04** | `mondo.obo` is a release asset; `mondo.sssom.tsv` is **not**, pin it by tag path |
| Human Disease Ontology | **v2026-08-31** | only needed if you keep a DOID-keyed view |
| DDInter | 2.0 | manual download |
| ChEMBL | 36 | **not verified as current** — check `ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest` before citing |

### Citations
- Open Targets Platform. *Nucleic Acids Research* 2025;53(D1):D1467.
- Vasilevsky NA, et al. Mondo: Unifying diseases for the world, by the world. *medRxiv* 2022. doi:10.1101/2022.04.13.22273750
- Schriml LM, et al. The Human Disease Ontology 2022 update. *NAR* 2022;50(D1):D1255. doi:10.1093/nar/gkab1063
- Xiong G, et al. DDInter 2.0. *NAR* 2025;53(D1):D1356.
- Chambers J, et al. UniChem. *J Cheminform* 2013;5:3. doi:10.1186/1758-2946-5-3

---

## Column names: resolved 2026-09-14
The Open Targets parquet schemas are not published in the docs, so these were
confirmed against the schema `tools/build_reference_sets.py` prints on first run:

| Dataset | Field | Resolution |
|---|---|---|
| `drug_molecule` | ~~`isApproved`~~ | **does not exist in 26.06**; approval is `maximumClinicalStage == 'APPROVAL'` |
| `clinical_indication` | `maxClinicalStage` | confirmed; note the name differs from `drug_molecule`'s `maximumClinicalStage` |
| `drug_mechanism_of_action` | `targets` | confirmed |
| `evidence_*` | `datatypeId` | confirmed; the genetic value is `genetic_literature` |

`data/reference_sets/MANIFEST.json` records the rule actually applied for each
set, so it is the authority if this table and the code ever diverge.

---

## data/gold_standards/ — the sets notebooks 06 and 07 read

Everything above describes `data/reference_sets/`, built from Open Targets by
`tools/build_reference_sets.py`. Notebooks 01, 04 and 05 use it. Notebooks 06 and
07 read a separate `data/gold_standards/` directory, which predates it. Recorded
here so the two are not confused. Since 2026-09-18 it is split into subfolders;
`data/gold_standards/README.md` is the index and `legacy/README.md` explains
every superseded file.

    post_cutoff_approvals/   the 128-pair post-cutoff approval set (notebook 07); built by tools/build_post_cutoff_approvals.py
    source/        the vendored repurposing-dataset v0.6.0 table and zip the builder reads
    bridges/       DOID / MONDO / MeSH / DrugBank / UniProt bridges (config/kgs.yaml, notebooks 06 and 07, the builder)
    ctd/           CTD files for notebook 06
    derived/       multihop_benchmark_v1.csv, written by notebook 06 on first run
    legacy/        every earlier attempt (v2 to v5, the FDA census, older package releases), frozen

| File | What it is | How to get it |
|---|---|---|
| `ctd/CTD_chem_gene_ixns.tsv.gz` | CTD chemical–gene interactions, filtered to *Homo sapiens* | `bash tools/fetch_ctd.sh` |
| `ctd/CTD_curated_genes_diseases.tsv.gz` | CTD gene–disease, `DirectEvidence ∈ {marker/mechanism, therapeutic}` | `bash tools/fetch_ctd.sh` |
| `bridges/mesh_to_doid.csv` | `MESH:*` → `DOID:*` | `python tools/build_gold_standard_aux.py` (from `doid.obo`) |
| `bridges/do_diseases.csv` | `DOID:*` → MONDO id and label | same (from `mondo.sssom.tsv` + `mondo.obo`) |
| `bridges/drugbank_drugs.csv` | DrugBank id → drug name | same (from `ot_drug_molecule`; no DrugBank account needed). **Contains essentially no biologics** — see the note below. |
| `bridges/uniprot_genesproteins.csv` | UniProt accession → symbol, Ensembl | same (from `ot_target`) |
| `bridges/mondo.sssom.tsv` | MONDO cross-map | same (copied from `data/crosswalks/`) |
| `bridges/drug_name_to_drugbank.csv` | drug name → DrugBank id, biologics included | `python tools/build_drug_name_xref.py` |
| `bridges/mondo_hierarchy.csv` | MONDO `is_a` edges for the builder's distance rule | `python tools/build_mondo_hierarchy.py` |
| `bridges/drugbank_xref.csv` | DrugBank → UNII / RxCUI / ATC / ChEMBL (not present) | `python tools/build_drugbank_xref.py` |
| `prospective/post_cutoff_approvals.csv` | the 128 post-2023 (drug, disease) pairs; see `data/POST_CUTOFF_APPROVALS.md` | `python tools/build_post_cutoff_approvals.py --kgs all --figure` |
| `legacy/fda_drugs_id_mapping.csv` | FDA approvals with a hand-set `Combo_Flag`; notebook 07 falls back to it only when the post-cutoff set is absent | **hand-curated; no download reproduces it** |
| `derived/multihop_benchmark_v1.csv` | 500 (drug, gene, disease) triples | generated by notebook 06 on first run |

### Versioning CTD

CTD releases monthly and puts no version in the filename, so the file alone is
not citable. `tools/fetch_ctd.sh` writes `CTD_PROVENANCE.txt` holding, per file,
the URL, fetch timestamp, byte count, SHA-256, and the report's own `#` header
block. **Cite the "Report created" date from that header, not the day it was
downloaded.** Record it here and in the manuscript's data-availability section
once fetched.

CTD asks to be cited as Davis AP, Wiegers TC, Johnson RJ, Sciaky D, Wiegers J,
Mattingly CJ. *Comparative Toxicogenomics Database (CTD): update 2023.* Nucleic
Acids Research 2023;51(D1):D1257–D1262. Confirm the current citation at
https://ctdbase.org/about/publications/ and the terms at
https://ctdbase.org/about/legal.jsp before redistributing anything.

### An open decision

The multi-hop standard could instead be derived from the Open Targets reference
sets already here: joining `relation_drug_gene` with `relation_disease_gene`,
restricted to approved indications, yields **396 (drug, gene, disease) triples**
over 249 drugs, 84 genes and 79 diseases. That keeps one arbiter across every
dimension and needs no download, at the cost of a narrower gene set (the targets of
approved drugs). CTD was kept for breadth; if a reviewer asks why this
dimension uses a second source, that is the trade that was made.


## `drugbank_drugs.csv` has no biologics (found 2026-09-17)

Measured against the 29 drugs in the post-cutoff approval set, this file
resolves 13 and misses 16, and all 16 misses are biologics: adalimumab,
tocilizumab, belimumab, nivolumab, pembrolizumab, daratumumab, durvalumab,
dupilumab, sacituzumab govitecan, semaglutide, risankizumab, trastuzumab
deruxtecan, glofitamab, teclistamab, epcoritamab, talquetamab. It also carries
roughly 2,400 rows whose `name` is a bare CHEMBL identifier.

The cause is upstream and the build is not at fault. Open Targets carries no
DrugBank cross-reference for most biologics: `ot_drug_molecule` gives
pembrolizumab (CHEMBL3137343, drugType Antibody) a `crossReferences` list of
EMA and DailyMed only. Any build keyed on a DrugBank xref drops it correctly,
and silently.

**Use `drug_name_to_drugbank.csv` instead** for name lookups. Built offline by
`tools/build_drug_name_xref.py` from the graphs' own node tables, which key
drugs by DrugBank accession, plus the usable rows of this file:

```
unique names   13,969        antibody names  239   (this file: 0)
rows           13,996        ambiguous        26   flagged, not resolved
```

Validated against the hand-curated gold standard: 25 of 29 drugs resolve
(against 13 here) and all 25 agree with the curated accession exactly, zero
disagreements. The 4 that do not resolve are bispecific antibodies approved
after BioKG's 2020 snapshot, so they are absent from every graph in the panel
and are unevaluable regardless.

Coverage is therefore bounded by the graph snapshots, which is the correct
boundary for this benchmark: a drug no evaluated graph contains cannot be
scored on it either.
