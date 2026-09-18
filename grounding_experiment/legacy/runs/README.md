# LLM run outputs

`09_big_*.csv`, one file per model x seed, produced by `grounding/03_run/`.
Copied from the pre-review run (`biokgsuite/results/tables/09_llm_runs/`) on
2026-09-14 so that stages 04/05/06 can run without repeating the sweeps.

These are INPUTS to `04_analysis`, `05_scaling` and `06_diagnostic`.
Re-running `03_run/` overwrites them. If you re-run the sweeps for the
submission, delete these first so there is no mix of old and new.
