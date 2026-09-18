#!/usr/bin/env bash
#SBATCH --job-name=biokgsuite-stab-emb
#SBATCH --partition=def
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --gres=gpu:H100:1
#SBATCH --output=logs/stab-emb-%j.out
#SBATCH --error=logs/stab-emb-%j.err
#
# Supplementary: does notebook 05's stability result survive a learned readout?
#
# Notebook 05 measures stability with neighbourhood heuristics, which are a
# deterministic function of the graph. Embedding training is stochastic, so this
# job first establishes a RETRAINING-VARIANCE FLOOR (several seeds on the
# unperturbed graph, correlated with each other) and only then reads the
# perturbed runs against it. A perturbed number above the floor means the model
# cannot reproduce itself, not that the graph is robust.
#
# CALIBRATE BEFORE YOU TRUST ANYTHING. At dim 8 / 1 epoch the TransE floor sits
# at about r = -0.02, i.e. chance. Raise --epochs until the floor is high on one
# KG, then fix that value for the whole panel. Step 1 below does exactly this;
# do not skip to step 2 and report the output.
#
# Submit:   sbatch hpc/run_stability_embeddings.sbatch
# Monitor:  squeue -u $USER ; tail -f logs/stab-emb-<jobid>.out

set -euo pipefail

REPO_DIR="${REPO_DIR:-$HOME/biokgsuite-bib}"
cd "$REPO_DIR"
mkdir -p logs

for marker in config/kgs.yaml src/biokgsuite/models/embedding.py; do
    if [ ! -e "$REPO_DIR/$marker" ]; then
        echo "FATAL: $REPO_DIR does not look like biokgsuite-bib (missing $marker)."
        exit 1
    fi
done

echo "================================================================="
echo "HOST:  $(hostname)  DATE: $(date)  JOBID: ${SLURM_JOB_ID:-(no slurm)}"
echo "REPO:  $REPO_DIR"
echo "================================================================="
nvidia-smi -L 2>&1 | head -3 || echo "  (no GPU on this node?)"

module purge
module load python/3.12.11
module load uv/latest
module load cuda 2>/dev/null || true
source .venv/bin/activate
echo "  python: $(which python)"
echo "  torch:  $(python -c 'import torch; print(torch.__version__, "cuda=", torch.cuda.is_available())' 2>&1)"

export PYTHONUNBUFFERED=1
export TQDM_MININTERVAL=30
export MPLBACKEND=Agg

RUN=supplementary/embeddings/run_stability_embeddings.py
KGS="${KGS:-hetionet biokg}"
EPOCHS="${EPOCHS:-50}"
DIM="${DIM:-32}"

# ── Step 1: floor calibration on one KG ─────────────────────────────
# Reports the retraining floor at the chosen epoch count. Read it before
# believing step 2: if the floor is low, every perturbed number is noise.
CAL_KG="$(echo $KGS | awk '{print $1}')"
echo
echo ">>> Step 1: floor calibration on ${CAL_KG} at ${EPOCHS} epochs, dim ${DIM}"
python "$RUN" --kg "$CAL_KG" --epochs "$EPOCHS" --dim "$DIM"

# ── Step 2: the panel ───────────────────────────────────────────────
for kg in $KGS; do
    echo
    echo "================================================================="
    echo ">>> ${kg}   ($(date +%H:%M:%S))"
    echo "================================================================="
    python "$RUN" --kg "$kg" --epochs "$EPOCHS" --dim "$DIM"
done

echo
echo "Done -> results/tables/supp_stability_embeddings.csv"
echo "Check the 'retraining floor' rows FIRST. Perturbed rows only carry meaning"
echo "where they fall below their model's floor."
