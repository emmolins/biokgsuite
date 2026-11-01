#!/usr/bin/env bash
#SBATCH --job-name=biokg-ladder
#SBATCH --partition=def
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --gres=gpu:H100:1
#SBATCH --output=logs/ladder-%j.out
#SBATCH --error=logs/ladder-%j.err
#
# Same-family Llama 3.x capability-scaling sweep, ONE (model, seed) per job, single H100.
# Use for the rungs a single GPU can serve: llama3.2:1b, llama3.2:3b, llama3.1:8b, llama3.3:70b.
# 405B needs the multi-GPU script run_llm_405b.sbatch. Seeds are those in config/params.yaml
# (grounding.seeds); one job per (model, seed) writes one run file, so jobs are independent.
#
# Submit one job per (model, seed) so they run in PARALLEL and finish fast:
#   for M in llama3.2:1b llama3.2:3b llama3.1:8b; do
#     for S in 0 1 2; do
#       sbatch --export=ALL,MODEL=$M,SEED=$S hpc/run_llm_ladder.sbatch
#     done
#   done
#
# The 70B rung doubles as the third headline model, so it needs both
# families and more time than the small rungs:
#   sbatch --time=16:00:00 --export=ALL,MODEL=llama3.3:70b,SEED=$S,FAMILIES=main+common hpc/run_llm_ladder.sbatch
#
# Monitor:   tail -f logs/ladder-<jobid>.out
# Resumable: just resubmit a job if it dies (idempotent output name per model+seed).

set -euo pipefail
MODEL="${MODEL:?set MODEL, e.g. --export=ALL,MODEL=llama3.1:8b,SEED=1}"
SEED="${SEED:-1}"
FAMILIES="${FAMILIES:-main}"; FAMILIES="${FAMILIES//+/,}"   # main = all 128 pairs; main+common adds the paired family (+ not comma: --export splits on commas)
SAFE="$(echo "$MODEL" | tr 'A-Z' 'a-z' | sed -E 's/[^a-z0-9]+/_/g; s/^_//; s/_$//')"   # same slug as runner.model_slug
OUT="results/tables/G02_run/runs/main__${SAFE}__seed${SEED}.csv"   # written by grounding_experiment/02_run

REPO_DIR="${REPO_DIR:-$HOME/biokgsuite-bib}"
cd "$REPO_DIR"; mkdir -p logs
echo "HOST $(hostname)  DATE $(date)  JOB ${SLURM_JOB_ID:-none}  MODEL $MODEL  SEED $SEED  -> $OUT"

module purge 2>/dev/null || true
module load python/3.12.11 2>/dev/null || module load python 2>/dev/null || true
module load uv/latest 2>/dev/null || true
source .venv/bin/activate
python -c "import sys; print('python', sys.version.split()[0])" || { echo "FATAL: venv python unusable on $(hostname)"; exit 1; }
export PYTHONUNBUFFERED=1 TQDM_MININTERVAL=15

# ── Ollama: install to ~/.local if missing (no root needed) ─────────────────
export OLLAMA_MODELS="${OLLAMA_MODELS:-$HOME/.ollama/models}"
# One server per job. Slurm can place several ladder jobs on one node, and
# every job binding 11434 would silently share the first job's server (the
# others' `ollama serve` fails on the port), so one instance ends up swapping
# models between jobs and every job crawls. The port is derived from the job id.
OLLAMA_PORT=$(( 20000 + ${SLURM_JOB_ID:-0} % 20000 ))
export OLLAMA_HOST=127.0.0.1:${OLLAMA_PORT}
mkdir -p "$HOME/.local/bin" "$OLLAMA_MODELS"
if ! command -v ollama >/dev/null 2>&1 && [[ ! -x "$HOME/.local/bin/ollama" ]]; then
    echo "[ollama] downloading standalone binary..."
    curl -fL https://github.com/ollama/ollama/releases/latest/download/ollama-linux-amd64.tar.zst -o /tmp/ollama.tar.zst
    ( tar --zstd -xf /tmp/ollama.tar.zst -C "$HOME/.local" 2>/dev/null \
      || tar -I zstd -xf /tmp/ollama.tar.zst -C "$HOME/.local" )
fi
export PATH="$HOME/.local/bin:$PATH"
nvidia-smi || true

# ── start server, wait ──────────────────────────────────────────────────────
ollama serve > "logs/ollama-${SLURM_JOB_ID:-local}.log" 2>&1 &
OLLAMA_PID=$!
trap 'kill $OLLAMA_PID 2>/dev/null || true' EXIT
for i in $(seq 1 30); do curl -sf "http://${OLLAMA_HOST}/api/tags" >/dev/null && break; sleep 2; done
curl -sf "http://${OLLAMA_HOST}/api/tags" >/dev/null \
    || { echo "FATAL: ollama did not come up on ${OLLAMA_HOST}; see logs/ollama-${SLURM_JOB_ID:-local}.log"; exit 1; }
echo "[ollama] up on ${OLLAMA_HOST}"

# ── pull this rung's model (idempotent) ─────────────────────────────────────
ollama pull "$MODEL"

# ── run the sweep ───────────────────────────────────────────────────────────
BKG_MODE=real BKG_MODELS="$MODEL" BKG_SEEDS="$SEED" BKG_FAMILIES="$FAMILIES" \
    jupyter nbconvert --to notebook --execute grounding_experiment/02_run/notebook.ipynb \
    --output "$REPO_DIR/logs/02_run_${SAFE}_s${SEED}.ipynb" --ExecutePreprocessor.timeout=-1

echo "DONE -> $OUT"
