# data/gold_standards/

Inputs for notebooks 06 (multi-hop reasoning) and 07 (prospective generalisation).

    post_cutoff_approvals/   the 128-pair post-cutoff approval set, the set used in the paper
    source/        the vendored repurposing-dataset v0.6.0 release it is built from
    bridges/       DOID / MONDO / MeSH / DrugBank / UniProt bridges
    ctd/           CTD files for notebook 06
    derived/       multihop_benchmark_v1.csv, written by notebook 06 on first run
    legacy/        every superseded attempt, frozen; see legacy/README.md and legacy/HISTORY.md

## post_cutoff_approvals/

`post_cutoff_approvals.csv` holds the 128 (drug, disease) pairs with MONDO,
MeSH and DrugBank ids, regulators, approval date and evidence links.
`post_cutoff_pair_status.csv` is the pair x graph coverage and leak check,
`post_cutoff_funnel.csv` the stage counts, `post_cutoff_dropped.csv` every
removed row with its reason, and `PROVENANCE.txt` the source SHA-256,
parameters and graphs checked. All five come from one command:

    python tools/build_post_cutoff_approvals.py --kgs all --figure

The gold standard and dropped tables are byte-identical on every rebuild; the
status and funnel tables are identical when the same graph views exist in
`results/checkpoints/07_kg_views/` (the shipped files used all six). Only the
`built:` line of `PROVENANCE.txt` changes. Construction, checks and version
history: `data/POST_CUTOFF_APPROVALS.md`.

## source/

`repurposing-dataset-v0.6.0/candidates_reviewed.csv` is the 1,288-row table of
post-2023 approvals from FDA, EMA, Health Canada, PMDA and TGA, classified at the
active-ingredient level (SHA-256 `bbb9fb87...55e2`). The zip beside it is the
full package; `NORMALISE=1 OFFLINE=1 ./run_all.sh` rebuilds the table
byte-identically from cached responses. `review_proposals_second_model_2026-09-18.csv`
holds the second-model identifier proposals for 73 rows, each checked against the
Node Normalizer before use. To upgrade: unpack the new release here, point
`SOURCE_DEFAULT` in `tools/build_post_cutoff_approvals.py` at it, rebuild, log the
version in `data/POST_CUTOFF_APPROVALS.md`, move the old release to
`legacy/approvals_raw/`.

## bridges/

| file | what it is | produced by |
|---|---|---|
| `do_diseases.csv` | DOID to MONDO id and label | `tools/build_gold_standard_aux.py` |
| `mesh_to_doid.csv` | MeSH to DOID | same |
| `mondo.sssom.tsv` | MONDO cross-map (copy of `data/crosswalks/`) | same |
| `uniprot_genesproteins.csv` | UniProt to symbol and Ensembl | same |
| `drugbank_drugs.csv` | DrugBank id to name, small molecules only | same |
| `drug_name_to_drugbank.csv` | drug name to DrugBank id, biologics included | `tools/build_drug_name_xref.py` |
| `mondo_hierarchy.csv` | MONDO `is_a` edges for the distance rule | `tools/build_mondo_hierarchy.py` |
| `approved_drug_pool.csv`, `approved_indications.csv` | approved drugs with DrugBank ids, and their approved indications as MONDO ids, from the vendored Open Targets release; the grounding experiment's distractor pool | `tools/build_distractor_pool.py` |
| `entrez_to_symbol.tsv` | Entrez gene id to HGNC symbol, for rendering DRKG and OpenBioLink gene nodes in dossiers; from the vendored HGNC set | `tools/build_entrez_symbol_map.py` |

Each builder writes a `*_PROVENANCE.txt` beside its output. Read by
`config/kgs.yaml`, `src/biokgsuite/gold/ontology_classifier.py`, the builder
(`BRIDGE_DIR`), notebooks 06 and 07, and `src/biokgsuite/grounding/`. `drugbank_xref.csv`, named in
`config/kgs.yaml`, has never been built.

## ctd/ and derived/

`ctd/` holds the two CTD files and `CTD_PROVENANCE.txt` from `bash tools/fetch_ctd.sh`
(needs an unrestricted network); cite the "Report created" date in the
provenance file. `derived/` is empty in a fresh checkout.

## legacy/

Nothing in `src/`, `config/` or the notebooks reads these, except that notebook
07 falls back to `legacy/fda_drugs_id_mapping.csv` when the post-cutoff set is
absent. Reorganised 2026-09-18 by `reorganise_gold_standards.py` (tag `v0.9-legacy`);
nothing was deleted.
