# tools/

Scripts run from the repo root. Grouped by what they produce; each script's
docstring gives its usage.

## Fetching inputs (network)

| script | writes |
|---|---|
| `fetch_kgs.py` | `data/raw/<kg>/<date>/` and `data/<kg>/` for PrimeKG, DRKG, OpenBioLink, MATRIX, at pinned releases with checksums |
| `convert_hetionet_json.py` | Hetionet's flat node and edge tables from its hetnet JSON |
| `build_crosswalks.py` | `data/crosswalks/`: MONDO, DOID, HGNC, GO, UBERON, WikiPathways at pinned releases |
| `fetch_opentargets.sh`, `fetch_unichem.sh`, `download_mondo_sssom.sh` | the Open Targets datasets, the PubChem to ChEMBL map, the MONDO SSSOM |
| `fetch_ctd.sh` | `data/gold_standards/ctd/` (needs an unrestricted network; CTD is outside the analysis allowlist) |

## Building derived inputs (offline)

| script | writes |
|---|---|
| `build_reference_sets.py` | `data/reference_sets/`, the coverage dimension's reference sets |
| `build_gold_standard_aux.py` | `data/gold_standards/bridges/{do_diseases, mesh_to_doid, drugbank_drugs, uniprot_genesproteins}.csv`, `mondo.sssom.tsv` |
| `build_mondo_hierarchy.py`, `build_drug_name_xref.py`, `build_entrez_symbol_map.py` | the MONDO `is_a` edges, the drug-name bridge, the Entrez-to-symbol map, each with a provenance note |
| `build_drugbank_xref.py` | optional: `bridges/drugbank_xref.csv` for MATRIX resolution; not built in this release |
| `build_post_cutoff_approvals.py` | `data/gold_standards/post_cutoff_approvals/`, the 128-pair set, its status, funnel and dropped tables and PROVENANCE.txt; `--extract-views` writes the per-graph views under `results/checkpoints/07_kg_views/` |
| `build_distractor_pool.py` | `bridges/approved_drug_pool.csv` and `approved_indications.csv` from the vendored Open Targets release |
| `make_post_cutoff_flow_figure.py`, `make_post_cutoff_flow_pptx.py` (+ `post_cutoff_flow_pptx.js`) | the construction flow chart, as PDF/PNG and as an editable PowerPoint |

## Records and housekeeping

| script | does |
|---|---|
| `provenance.py` | regenerates `data/PROVENANCE.md` from the files on disk |
| `results_manifest.py` | regenerates `results/MANIFEST.md` from the files on disk |
| `regenerate_figs.py` | redraws every dimension figure from the checkpoints, without re-running the notebooks |
| `run_notebook.py` | executes a notebook headlessly where jupyter is not installed |
| `reversal_sensitivity.py`, `concordance_panels.py` | Supplementary Table S4d and the five-graph concordance rows, from the `00_summary` tables |
| `task_contrasts.py` | Supplementary Table S8d, pairwise graph contrasts between the two tasks, from the `06_task_performance` tables |
| `audit_identifiers.py` | per graph and declared type, every identifier shape and the share each registry confirms; the check behind reading a bare number as the type the graph declares. Writes `results/tables/id_audit/` |
| `measure_shared_retrieval_set.py` | diagnostic: sizes a shared disease-gene retrieval benchmark; writes nothing |
| `cleanup.sh` | removes `_to_delete/`, caches and OS files (`--yes` to delete) |

The tools of the superseded 116-pair benchmark are preserved at the tag
`v0.9-legacy`; they read and write under `data/gold_standards/legacy/`, which
stays on `main` as the record of that earlier set.

- `build_dashboard.py`: builds `docs/index.html`, the interactive results page served by GitHub Pages, from the tables under `results/tables` and the template `dashboard_template.html`.
