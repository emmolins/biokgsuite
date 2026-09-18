#!/usr/bin/env bash
#SBATCH --job-name=biokgsuite-jupyter
#SBATCH --partition=def
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --output=logs/jupyter-%j.out
#SBATCH --error=logs/jupyter-%j.err
#
# Launches a Jupyter Lab server on a CPU compute node so you can run the main
# notebooks (00-07) interactively and see the figures render inline. No GPU is
# requested (embedding training lives in nb08). MATRIX is the memory driver,
# hence 120 GB.
#
# Usage:
#   1. sbatch scripts/hpc/run_jupyter.sbatch
#   2. tail -f logs/jupyter-<jobid>.out      # wait for the "JUPYTER READY" block
#   3. On your Mac, run the printed `ssh -N -L ...` tunnel command (new terminal)
#   4. Open the printed http://127.0.0.1:<port>/lab?token=... URL in your browser
#   5. Run the notebooks cell by cell; figures appear inline
#
# When done: scancel <jobid>   (or just let the 8 hr wall time expire)

set -euo pipefail

REPO_DIR="${REPO_DIR:-$HOME/biokgsuite}"
cd "$REPO_DIR"
mkdir -p logs

module purge
module load python/3.12.11
module load uv/latest
source .venv/bin/activate

# Main notebooks need networkx + an interactive Jupyter front end.
python -c "import networkx" 2>/dev/null || uv pip install "networkx>=3.2,<4"
python -c "import jupyterlab" 2>/dev/null || uv pip install jupyterlab

PORT=$(shuf -i 8888-8999 -n 1)
NODE=$(hostname)
LOGIN="ashburn-login"   # change to your login-node ssh alias if different

cat <<EOF

=================== JUPYTER READY ===================
Compute node : ${NODE}
Port         : ${PORT}

1) On your Mac, in a NEW terminal, open the tunnel:

     ssh -N -L ${PORT}:${NODE}:${PORT} ${LOGIN}

2) Then open the URL printed just below (the 127.0.0.1:${PORT} one,
   including its ?token=...) in your browser.

Leave this job running while you work. scancel ${SLURM_JOB_ID:-<jobid>} to stop.
=====================================================

EOF

jupyter lab --no-browser --ip=0.0.0.0 --port="${PORT}"
