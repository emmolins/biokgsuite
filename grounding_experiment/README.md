# Grounding experiment

Does knowledge-graph evidence help a language model rank the true post-cutoff
repurposed drug above distractors, and does the answer depend on which graph
supplied the evidence? Five stages, run in order; only `02_run` calls models.

| Stage | Notebook | Produces |
|---|---|---|
| `01_benchmark/` | the 128 post-cutoff pairs, their per-graph status, the candidate pools, the selection audit | `G01_benchmark.pkl`, Fig G01 |
| `02_run/` | the ranking loop over models, arms, seeds and shuffles (`mock` or `real`) | `results/tables/G02_run/runs/` |
| `03_analysis/` | paired bootstrap, graph contrasts, coverage control, FDR | tables and Fig G03 (main text) |
| `04_scaling/` | the Llama capacity ladder, pooled versus covered lift | Fig G04 (supplement) |
| `05_diagnostic/` | which pairs are hard and why | Fig G05 (supplement) |

The benchmark is the same post-cutoff approval set notebook 07 scores
(`data/POST_CUTOFF_APPROVALS.md`), so the two parts of the paper share one
provenance trail. Every parameter is in `config/params.yaml` under
`analysis_params.grounding`; the notebooks read nothing else. Shared code is
the package: `src/biokgsuite/grounding/` (benchmark, prompt, clients, dossiers,
runner, analysis, bias) and `src/biokgsuite/viz/grounding.py` (figures), with
tests under `tests/test_grounding_*.py`.

The notebooks follow `biokgsuite_dimensions/CONVENTIONS.md`. Each writes
through `Outputs('G0N_<stage>', BASE)`, so figures, tables and checkpoints land
under `results/{figures,tables}/G0N_<stage>/` and `results/checkpoints/`.

## Run

    # once: the distractor pool (Open Targets) and the gene-symbol map (HGNC), both from vendored releases
    python tools/build_distractor_pool.py
    python tools/build_entrez_symbol_map.py

    # stages, from the repo root; 02 in mock mode first, then real
    jupyter nbconvert --to notebook --execute grounding_experiment/01_benchmark/notebook.ipynb
    BKG_MODE=mock jupyter nbconvert --to notebook --execute grounding_experiment/02_run/notebook.ipynb
    BKG_MODE=real BKG_MODELS=gpt-4.1-mini BKG_SEEDS=0 jupyter nbconvert --to notebook --execute grounding_experiment/02_run/notebook.ipynb
    ...
    jupyter nbconvert --to notebook --execute grounding_experiment/03_analysis/notebook.ipynb

`02_run` is idempotent per (family, model, seed) file; `BKG_MODELS`, `BKG_SEEDS`,
`BKG_FAMILIES` and `BKG_FORCE=1` scope one launch. Hosted models read their
key from the environment or from `.env` at the repo root (gitignored).

## legacy/

The 2026-09 experiment on the hand-assembled 116-pair set: its six stage
notebooks, `lib/`, and every run file. Kept for the record; nothing current
reads it. See `legacy/README.md`.
