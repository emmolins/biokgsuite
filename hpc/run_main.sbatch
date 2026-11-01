#!/usr/bin/env bash
#SBATCH --job-name=biokgsuite-main
#SBATCH --partition=def
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --output=logs/main-%j.out
#SBATCH --error=logs/main-%j.err
#
# Runs the main benchmark notebooks 01-07 in order, then 00 (which aggregates
# their checkpoints), on one CPU compute node. No GPU is requested: embedding
# training and Gemma encoding live in supplementary/embeddings (hpc/run_embeddings.sbatch), not
# here. MATRIX is the memory driver: evaluated on its full released
# graph (config/kgs.yaml keep_categories: all) it peaks above 150 GB across
# nb01+nb02, hence 256 GB.
#
# Submit:   sbatch hpc/run_main.sbatch
# Monitor:  squeue -u $USER
#           tail -f logs/main-<jobid>.out
# Resumable: each notebook caches to results/checkpoints/ and results/cache/,
#            so a re-submit skips work that already finished.

set -euo pipefail

REPO_DIR="${REPO_DIR:-$HOME/biokgsuite}"
cd "$REPO_DIR"
mkdir -p logs
[ ! -f .trial_kgs ] || { echo "FATAL: .trial_kgs is present; the cluster runs the full panel. Remove it."; exit 1; }

# Refuse to run against anything that is not THIS repo. The default used to
# point at $HOME/biokgsuite, a different checkout with a different notebook
# layout, so a submit would silently execute superseded code and produce
# results that looked valid. Fail loudly instead.
for marker in config/kgs.yaml biokgsuite_dimensions src/biokgsuite; do
    if [ ! -e "$REPO_DIR/$marker" ]; then
        echo "FATAL: $REPO_DIR does not look like biokgsuite (missing $marker)."
        echo "       Set REPO_DIR to the checkout you mean to run."
        exit 1
    fi
done

echo "================================================================="
echo "HOST:  $(hostname)  DATE: $(date)  JOBID: ${SLURM_JOB_ID:-(no slurm)}"
echo "REPO:  $REPO_DIR"
echo "================================================================="

# ── Environment ─────────────────────────────────────────────────────
module purge
module load python/3.12.11
module load uv/latest
source .venv/bin/activate
echo "  python: $(which python) ($(python --version))"

# 00-07 need networkx (topology + community detection) and pytest (tests/).
# The embedding setup does not install them; add on demand so a fresh venv works.
python -c "import networkx" 2>/dev/null || uv pip install "networkx>=3.2,<4"
python -c "import pytest"   2>/dev/null || uv pip install "pytest>=8.0"

export PYTHONUNBUFFERED=1     # stream prints to the Slurm log immediately
export TQDM_MININTERVAL=30    # one progress line per ~30s (readable in logs)
export MPLBACKEND=Agg         # headless matplotlib

# ── Run 01-07 in order, then 00 (aggregation) ───────────────────────
# set -e stops the job if any notebook errors, so 00 never aggregates a
# partial set of checkpoints.
# Notebooks live one per directory as <dir>/notebook.ipynb; nbconvert is run
# from inside that directory so each notebook's root-finding bootstrap and its
# relative paths behave exactly as they do interactively.
# NOTEBOOKS overrides the list, for a rerun from the notebook that failed:
#   sbatch --export=ALL,NOTEBOOKS=06_task_performance+07_generalization+00_summary hpc/run_main.sbatch
# Each notebook reads its inputs from results/checkpoints, so the earlier
# ones need not be rerun as long as their checkpoints are from this panel.
NOTEBOOKS="${NOTEBOOKS:-01_coverage 02_annotation_accuracy 03_trustworthiness 04_topology 05_stability 06_task_performance 07_generalization 00_summary}"
NOTEBOOKS="${NOTEBOOKS//+/ }"      # + separates names on the sbatch line
for nb in $NOTEBOOKS; do
    target="biokgsuite_dimensions/${nb}/notebook.ipynb"
    if [ ! -f "$target" ]; then
        echo "FATAL: $target not found."; exit 1
    fi
    echo
    echo "================================================================="
    echo ">>> Executing ${target}   ($(date +%H:%M:%S))"
    echo "================================================================="
    ( cd "biokgsuite_dimensions/${nb}" && \
      jupyter nbconvert \
        --to notebook --execute --inplace \
        --log-level=INFO \
        --ExecutePreprocessor.timeout=86400 \
        notebook.ipynb )
done

echo
echo ">>> Metric unit tests"
python -m pytest tests -q

echo
echo "================================================================="
echo "Done. Outputs in results/figures, results/tables, results/checkpoints,"
echo "and the executed notebooks (00-07) with inline cell outputs."
echo
echo "Pull results back to your Mac:"
echo "  rsync -avz ashburn-login:biokgsuite/results/ ~/Desktop/biokgsuite/results/"
echo "  rsync -avz ashburn-login:biokgsuite/biokgsuite_dimensions/ ~/Desktop/biokgsuite/biokgsuite_dimensions/"
echo "================================================================="
