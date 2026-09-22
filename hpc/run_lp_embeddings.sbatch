#!/usr/bin/env bash
#SBATCH --job-name=biokg-lp-emb
#SBATCH --partition=def
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=96G
#SBATCH --gres=gpu:H100:1
#SBATCH --output=logs/lp-emb-%j.out
#SBATCH --error=logs/lp-emb-%j.err
#
# Supplementary Table S7 only: notebook 06's link-prediction test pairs rescored
# with TransE and RotatE (supplementary/embeddings/run_link_prediction_embeddings.py).
# The full embeddings notebook (run_embeddings.sbatch) also trains the 60
# resampled models and the Gemma priors, which did not fit in 24 h; this job
# produces the one table the paper cites.
#
# The script writes results/tables/06_task_performance/link_prediction_embeddings.csv
# after every graph and skips graphs already in it, so a job that hits the wall
# clock keeps what it finished and a resubmit carries on.
#
# Submit (all graphs still missing from the table):
#   sbatch hpc/run_lp_embeddings.sbatch
# One graph, or a few, per job (MATRIX on its own is wise):
#   sbatch --export=ALL,KGS=matrix hpc/run_lp_embeddings.sbatch
#   sbatch --export=ALL,KGS=primekg+drkg+openbilink hpc/run_lp_embeddings.sbatch
# Redo graphs already in the table (for example ones scored at another dimension):
#   sbatch --export=ALL,KGS=hetionet+biokg,FORCE=1 hpc/run_lp_embeddings.sbatch

set -euo pipefail
REPO_DIR="${REPO_DIR:-$HOME/biokgsuite-bib}"
cd "$REPO_DIR"; mkdir -p logs
KGS="${KGS:-}"; KGS="${KGS//+/,}"          # + not comma: --export splits on commas
FORCE="${FORCE:-}"
echo "HOST $(hostname)  DATE $(date)  JOB ${SLURM_JOB_ID:-none}  KGS ${KGS:-all}"

module purge
module load python/3.12.11
module load uv/latest
module load cuda 2>/dev/null || true
source .venv/bin/activate
echo "  torch: $(python -c 'import torch; print(torch.__version__, "cuda=", torch.cuda.is_available())' 2>&1)"
export PYTHONUNBUFFERED=1 TQDM_MININTERVAL=30 MPLBACKEND=Agg

python supplementary/embeddings/run_link_prediction_embeddings.py ${KGS:+--kgs "$KGS"} ${FORCE:+--force}
echo "DONE"
