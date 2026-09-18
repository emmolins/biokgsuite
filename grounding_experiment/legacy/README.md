# grounding_experiment/legacy/

The grounding experiment as it stood before 2026-09-18, frozen. It ran on the
hand-assembled 116-pair set (`data/gold_standards/legacy/gold_standard_v4_bigset.tsv`)
and its 76-pair common-3 subset, both superseded by the post-cutoff approval
set because neither could be rebuilt from a cited source (see
`data/gold_standards/legacy/README.md`). Its runs are therefore not comparable
with the current stages and are not used anywhere.

| folder | what it was |
|---|---|
| `01_construction/` | scripts that reconstructed the 253 to 116 funnel and its per-KG coverage table |
| `02_bias_audit/` | kept-versus-dropped audit of the 116 against the 253 candidates |
| `03_run/` | the ranking loop, prompt, clients and KG dossier code, all inside one notebook |
| `04_analysis/` | `build_09_analysis.py`, which generated the analysis notebook; the query-order replay it needed is gone because run rows now carry `pair_id` |
| `05_scaling/`, `06_diagnostic/` | Llama ladder and hard-pair diagnostics on the 116 runs |
| `lib/` | the four modules the stages imported (they also imported a `figstyle` and a flat `src/` layout that no longer exist, so they do not run as-is) |
| `runs/` | every run CSV of the 116-pair and common-3 experiments (GPT-4.1-mini, Gemini, Llama 3.x ladder); `README_old.md` lists the stages |
| `previous_run/` | the derived tables of that experiment's analysis notebook |
| `results_116/` | the construction figure and the corroboration and audit tables the 116-pair funnel produced |
| `METHODS_116_CONSTRUCTION.md` | the methods record of the 253 to 116 construction |

Everything reusable was ported to `src/biokgsuite/grounding/` and
`src/biokgsuite/viz/grounding.py` with the changes recorded in the module
docstrings; the stages were rebuilt as `../01_benchmark` to `../05_diagnostic`.
