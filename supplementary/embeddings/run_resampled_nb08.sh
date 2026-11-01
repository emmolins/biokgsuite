#!/usr/bin/env bash
# Run the resampled (multi-rerun) version of nb08 end-to-end on your machine.
#
# What this does:
#   Re-executes notebook 08 with the new resampling section enabled. Each of
#   the 5 reruns draws a different random 10% of drug-disease indications,
#   retrains TransE and RotatE from scratch, re-scores Gemma against the new
#   test set, and saves per-rerun results.
#
#   The single-run analysis (cells 2-29) runs first and produces bootstrap
#   CIs. The resampled section (cells 30-44) then produces empirical CIs
#   across the 5 reruns, and the comparison table at the end quantifies how
#   much the bootstrap underestimates variance.
#
# Prereqs:
#   1. cd into the biokgsuite repo root
#   2. conda activate biokgsuite
#   3. pip install torch transformers sentence-transformers  (for Gemma)
#   4. Accept EmbeddingGemma license + export HF_TOKEN  (see scripts/run_gemma_benchmark.sh)
#
# Expected wall-clock (TransE/RotatE training time × 5 reruns):
#   Hetionet     ~5 min × 5 = ~25 min          per model
#   PrimeKG      ~35 min × 5 = ~3 hr           per model
#   DRKG         ~10 min × 5 = ~50 min         per model
#   OpenBioLink  ~8 min × 5 = ~40 min          per model
#   BioKG        ~3 min × 5 = ~15 min          per model
#   MATRIX       ~30 min × 5 = ~2.5 hr         per model (with PER_KG_OVERRIDES; full hyperparams: ~10 hr × 5)
#
# Per model totals: ~7-8 hr for small KGs + MATRIX combined.
# With 2 models (TransE + RotatE): ~15 hr on CPU.
# Gemma reruns are cheap (~30 min total, encoding is cached).
#
# On a GPU you can expect ~3-6 hr total.
#
# Resumable: per-rerun cache files at results/cache/embedding_<kg>_resampled.json
# and results/cache/embedding_<kg>_gemma_resampled.json. Crash → restart → continues.

set -euo pipefail

cd "$(dirname "$0")/.."

if [[ -z "${HF_TOKEN:-}" ]]; then
    echo "WARNING: HF_TOKEN not set. Gemma reruns will be skipped." >&2
    echo "  Set HF_TOKEN to include them: export HF_TOKEN=hf_..." >&2
fi
export HUGGING_FACE_HUB_TOKEN="${HF_TOKEN:-}"

echo
echo "Executing nb08 (existing single-run + new resampled section)..."
echo "This is long. Progress is written incrementally to:"
echo "  results/cache/embedding_<kg>_resampled.json        (TransE/RotatE)"
echo "  results/cache/embedding_<kg>_gemma_resampled.json  (Gemma)"
echo "Safe to Ctrl-C and re-run later."
echo

jupyter nbconvert --to notebook --execute --inplace \
        --ExecutePreprocessor.timeout=172800 \
        supplementary/embeddings/notebook.ipynb

echo
echo "Done."
echo
echo "Resampled outputs:"
echo "  - results/tables/08_embedding_comparison_resampled.csv         (long-form, one row per rerun)"
echo "  - results/figures/08_resampled_auroc.{pdf,png}              (AUROC by KG x model)"
echo "  - results/figures/08_lift_over_gemma.{pdf,png}              (lift over Gemma prior)
#   - results/figures/08_heuristic_vs_embedding_scatter.{pdf,png}"
echo
echo "The 'bootstrap vs rerun' table in the notebook quantifies how much"
echo "the single-run bootstrap CIs were too tight. Cite the rerun-empirical"
echo "CIs as the headline stability numbers in your writeup."
