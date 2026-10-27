# Traceability: manuscript display items and statistics

Where every figure, table and reported number in the BiB manuscript and its
supplement comes from. Paths are relative to the repository root; notebooks are
`biokgsuite_dimensions/NN_*/notebook.ipynb` and `grounding_experiment/G0N_*/notebook.ipynb`.
`results/MANIFEST.md` lists every file with its checksum.

## Main text

| Item | Source files | Produced by |
|---|---|---|
| Table 1 (metrics by dimension) | `config/params.yaml`; values in Supplementary Table S3a | notebooks 01 to 07 |
| Table 2 (graph properties) | `data/PROVENANCE.md`, `config/kg_snapshots.yaml`, `results/tables/04_topology/raw_statistics.csv` | `tools/provenance.py`, notebook 04 |
| Fig. 1 | `paper/figures/Fig1_framework.pptx` | drawn by hand |
| Fig. 2 (profiles) | `results/figures/00_summary/quality_profiles.*`; data `results/tables/00_summary/dimension_scores.csv` | notebook 00 |
| Fig. 3 (metric heatmap) | `results/figures/00_summary/metric_heatmap.*`; data `results/tables/0N_*/dimension_scores.csv` | notebook 00 |
| Fig. 4a | `results/tables/06_task_performance/link_prediction.csv`, `neighbourhood_retrieval.csv` | `paper/figures/make_fig4.py` |
| Fig. 4b | `results/tables/04_topology/hop_recovery.csv` (reachability) and `link_prediction.csv` | `paper/figures/make_fig4.py` |
| Fig. 4c | `results/tables/07_generalization/sparse_tiers.csv` | `paper/figures/make_fig4.py` |
| Fig. 5 (query design) | `paper/figures/Fig5_evaluation_set.pptx`; counts from `results/tables/G01_benchmark/benchmark.csv`, `graph_status.csv` | notebook G01 |
| Fig. 6 (grounding) | `results/figures/G03_analysis/G03a_grounding_gain.*`, `G03b_rank_distribution.*`; data `results/tables/G03_analysis/arms_main.csv`, `contrasts_main.csv` | notebook G03 |
| Fig. 7 (capacity) | `results/figures/G04_scaling/G04a_capacity.*`; data `results/tables/G04_scaling/capacity.csv`, `scaling.csv` | notebook G04 |
| Fig. 8 (prompt) | `paper/figures/Fig8_example_prompt.*`; query from `results/tables/G01_benchmark/queries_seed0.csv` | `src/biokgsuite/grounding/prompt.py` |

Key statistics in the Results, by paragraph:

| Statistic | File |
|---|---|
| Dimension scores (MATRIX 0.920, PrimeKG 0.716, ...) and coverage 0.68 / 0.48 / 0.47 | `results/tables/00_summary/dimension_scores.csv` |
| Reversal rates 49% (five graphs), 40% and 32% (six), per dimension pair | `results/tables/00_summary/dimension_dissociations.csv` |
| Gap-restricted reversal rates (42% at 0.05 and 0.10, five graphs) | `results/tables/00_summary/reversal_sensitivity.csv` (`tools/reversal_sensitivity.py`) |
| Kendall's W, mean rho, permutation p | `results/tables/00_summary/concordance_specifications.csv`; permutation null (20,000 draws, seed 0) in notebook 00 |
| Rank span per graph | `results/tables/00_summary/rank_dispersion.csv` |
| Link prediction AUROC 0.998 / 0.754 etc. | `results/tables/06_task_performance/link_prediction.csv` (Adamic-Adar, type-constrained) |
| Retrieval Recall@0.5% (0.823, 0.553, 0.550, 0.380, 0.173, 0.149) | `results/tables/06_task_performance/neighbourhood_retrieval.csv` |
| Significant reversals between the two tasks | `results/tables/06_task_performance/task_contrasts.csv` (`tools/task_contrasts.py`) |
| Graph reachability 0.44 vs 0.55 to 0.98 | `results/tables/04_topology/hop_recovery.csv`, `dimension_scores.csv` |
| OpenBioLink edge retention by relation (3 to 5% of drug-gene) | `src/biokgsuite/io/loading.py` (`_openbilink_hq_mask`); counts in `data/openbilink/OPENBILINK_RELEASE.txt` |
| Sparse-tier AUROC (0.558 to 0.875), drops of 0.11 to 0.37 | `results/tables/07_generalization/sparse_tiers.csv` |
| Cross-domain per-area AUROC and CV | `results/tables/07_generalization/cross_domain.csv`, `domain_equity.csv` |
| Prospective AUROC (0.776, 0.868, 0.677, 0.831) | `results/tables/07_generalization/prospective.csv`, `prospective_coverage.csv` |
| Grounding: pairs with no profile (30 / 41 / 39; 9 / 9 / 6 held) | `results/tables/G03_analysis/diagnostics.csv`, `absence_penalty.csv` |
| MRR contrasts, BH over 9 contrasts per reading | `results/tables/G03_analysis/contrasts_main.csv`, `contrasts_common.csv`, `significance_tests.csv` |
| Lift with and without evidence (-0.21 to +0.02 vs +0.38 to +0.54) | `results/tables/G03_analysis/absence_penalty.csv` |
| KG and LLM as additive levers | `results/tables/G03_analysis/lever_decomposition.csv` |
| Llama capacity ladder | `results/tables/G04_scaling/capacity.csv` |
| Gemini below-chance omission rate | `results/tables/G03_analysis/fabrication_rate.csv` |

## Supplement

| Item | Source files | Produced by |
|---|---|---|
| Table S1 | written from the literature; no computation | |
| Table S2 | `data/REFERENCE_SETS.md`, `results/tables/01_coverage/reference_sets.csv`, `data/gold_standards/derived/multihop_benchmark_v1.meta.json` | `tools/build_reference_sets.py`, notebook 06 |
| Table S3a | `results/tables/0N_*/dimension_scores.csv`, summarised in `00_summary/dimension_scores.csv` | notebooks 01 to 07 |
| Table S3b to e | `results/tables/04_topology/graph_cohesion.csv`, `clustering_coefficient.csv`, `community_purity.csv`, `gold_pair_mapping.csv`, `hop_recovery.csv` | notebook 04 |
| Table S3f | `results/tables/07_generalization/cross_domain.csv` | notebook 07 |
| Table S4a, b, c | `results/tables/00_summary/dimension_dissociations.csv`, `rank_dispersion.csv`, `concordance_specifications.csv` | notebook 00 |
| Table S4d | `results/tables/00_summary/reversal_sensitivity.csv` | `tools/reversal_sensitivity.py` |
| Table S5 | `results/tables/G01_benchmark/selection_audit.csv` | notebook G01 |
| Table S6 | `results/tables/G03_analysis/*.csv`, `G04_scaling/capacity.csv` | notebooks G03, G04 |
| Table S7 | `results/tables/06_task_performance/link_prediction_embeddings.csv`, `embedding_scores/*.npz` | `supplementary/embeddings/` |
| Table S8a to c | `results/tables/06_task_performance/link_prediction.csv` (all scorers and strategies), `neighbourhood_retrieval.csv`, `multihop_reasoning.csv`, `05_stability/stability.csv` | notebooks 05, 06 |
| Table S8d | `results/tables/06_task_performance/task_contrasts.csv` | `tools/task_contrasts.py` |
| Table S9 | `results/tables/id_audit/identifier_namespaces.csv`, `01_coverage/resolvability.csv` | `tools/audit_identifiers.py`, notebook 01 |
| Fig. S1 | `results/figures/04_topology/hop_recovery.*` | notebook 04 |
| Fig. S2 | `results/figures/07_generalization/S2_sparse_tiers_and_domain.*` | notebook 07 |
| Fig. S3 | `results/figures/07_generalization/07i_post_cutoff_flow.*`, `G01_benchmark/G01d_benchmark_overview.*` | `tools/build_post_cutoff_approvals.py --figure`, notebook G01 |
| Fig. S4 | `results/figures/G03_analysis/G03a_grounding_gain.*` | notebook G03 |
| Fig. S5 | `results/figures/G03_analysis/G03c_evidence_use.*` | notebook G03 |
| Fig. S6 | `results/figures/07_generalization/07d_prospective_auroc.*`, `07e_prospective_hits.*` | notebook 07 |
| Fig. S7 | `results/figures/05_stability/stability.*` | notebook 05 |
| Fig. S8 | `results/figures/06_task_performance/link_prediction_heatmap.*` | notebook 06 |

## Audit, 2026-10-27

All 746 non-trivial numbers in the two documents (decimals with two or more
places, integers of 10 or more) were matched against the result tables,
reference-set records and configuration in this tree. Every one was found,
allowing for rounding, except the Data Availability commit hash, a MIM
identifier quoted inside the example prompt, and the multi-hop benchmark counts
(41 triples; 226 resolvable in all graphs; 170 removed), which match
`data/gold_standards/derived/multihop_benchmark_v1.meta.json` in the working
copy that produced the results rather than the previously committed version.
