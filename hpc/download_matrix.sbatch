#!/usr/bin/env bash
#SBATCH --job-name=matrix-dl
#SBATCH --partition=def
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/matrix-dl-%j.out
#SBATCH --error=logs/matrix-dl-%j.err
#
# Download Every Cure's MATRIX KG (parquet shards) and convert to the
# TSV format the BioKGSuite loader expects.
#
# Runs as a Slurm job on a compute node — much faster than the login
# node because of better NFS access + dedicated network bandwidth.
# Expected wall time: 5-15 min.
#
# Submit:    sbatch scripts/hpc/download_matrix.sbatch
# Monitor:   tail -f logs/matrix-dl-<jobid>.out

set -euo pipefail

REPO_DIR="${REPO_DIR:-$HOME/biokgsuite}"
cd "$REPO_DIR"
mkdir -p logs

echo "================================================================="
echo "HOST:  $(hostname)"
echo "DATE:  $(date)"
echo "REPO:  $REPO_DIR"
echo "================================================================="

module purge
module load python/3.12.11
module load uv/latest 2>/dev/null || true

source .venv/bin/activate
echo "  python: $(which python) ($(python --version))"

# Faster HF downloads
export HF_HUB_ENABLE_HF_TRANSFER=1
export HF_HOME="$HOME/.cache/huggingface"
export PYTHONUNBUFFERED=1

python -u scripts/hpc/download_matrix.py
