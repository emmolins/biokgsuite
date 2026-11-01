#!/usr/bin/env bash
#SBATCH --job-name=biokg-405b
#SBATCH --partition=def
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=256G
#SBATCH --gres=gpu:H100:4
#SBATCH --output=logs/405b-%j.out
#SBATCH --error=logs/405b-%j.err
#
# Top rung of the ladder: Llama-3.1-405B via Ollama, sharded across 4xH100.
# Llama-3.1-405B at q4 is ~230 GB of weights, so it needs >=4xH100 (320 GB) to
# hold weights + KV cache, and a one-time ~230 GB model pull. Ollama auto-shards
# across all visible GPUs (layer split), slower per token than a single-GPU 70B,
# hence the 24h wall and longer queue. Throughput-critical? swap Ollama for vLLM
# with tensor-parallel (TP=4); this script keeps the same harness for simplicity.
#
# PRE-FLIGHT (do once, interactively, before submitting):
#   - confirm OLLAMA_MODELS lives on a filesystem with >=250 GB free
#   - `df -h "$OLLAMA_MODELS"` and `squeue` for 4-GPU node availability
#
# Submit ONE seed at a time (each is long):
#   sbatch --export=ALL,SEED=1 hpc/run_llm_405b.sbatch
#   sbatch --export=ALL,SEED=2 hpc/run_llm_405b.sbatch
#   sbatch --export=ALL,SEED=3 hpc/run_llm_405b.sbatch

set -euo pipefail
MODEL="${MODEL:-llama3.1:405b}"
SEED="${SEED:-1}"
FAMILIES="${FAMILIES:-main}"; FAMILIES="${FAMILIES//+/,}"   # main = all 128 pairs; main+common adds the paired family (+ not comma: --export splits on commas)
SAFE="$(echo "$MODEL" | tr 'A-Z' 'a-z' | sed -E 's/[^a-z0-9]+/_/g; s/^_//; s/_$//')"   # same slug as runner.model_slug
OUT="results/tables/G02_run/runs/main__${SAFE}__seed${SEED}.csv"   # written by grounding_experiment/02_run

REPO_DIR="${REPO_DIR:-$HOME/biokgsuite}"
cd "$REPO_DIR"; mkdir -p logs
echo "HOST $(hostname)  DATE $(date)  JOB ${SLURM_JOB_ID:-none}  MODEL $MODEL  SEED $SEED  GPUS ${CUDA_VISIBLE_DEVICES:-?}  -> $OUT"

module purge 2>/dev/null || true
module load python/3.12.11 2>/dev/null || module load python 2>/dev/null || true
module load uv/latest 2>/dev/null || true
source .venv/bin/activate
python -c "import sys; print('python', sys.version.split()[0])" || { echo "FATAL: venv python unusable on $(hostname)"; exit 1; }
export PYTHONUNBUFFERED=1 TQDM_MININTERVAL=30

export OLLAMA_MODELS="${OLLAMA_MODELS:-$HOME/.ollama/models}"
# One server per job. Slurm can place several ladder jobs on one node, and
# every job binding 11434 would silently share the first job's server (the
# others' `ollama serve` fails on the port), so one instance ends up swapping
# models between jobs and every job crawls. The port is derived from the job id.
OLLAMA_PORT=$(( 20000 + ${SLURM_JOB_ID:-0} % 20000 ))
export OLLAMA_HOST=127.0.0.1:${OLLAMA_PORT}
export OLLAMA_SCHED_SPREAD=1          # spread layers across all 4 GPUs
mkdir -p "$HOME/.local/bin" "$OLLAMA_MODELS"
if ! command -v ollama >/dev/null 2>&1 && [[ ! -x "$HOME/.local/bin/ollama" ]]; then
    echo "[ollama] downloading standalone binary..."
    curl -fL https://github.com/ollama/ollama/releases/latest/download/ollama-linux-amd64.tar.zst -o /tmp/ollama.tar.zst
    ( tar --zstd -xf /tmp/ollama.tar.zst -C "$HOME/.local" 2>/dev/null \
      || tar -I zstd -xf /tmp/ollama.tar.zst -C "$HOME/.local" )
fi
export PATH="$HOME/.local/bin:$PATH"
nvidia-smi || true

ollama serve > "logs/ollama-405b-${SLURM_JOB_ID:-local}.log" 2>&1 &
OLLAMA_PID=$!
trap 'kill $OLLAMA_PID 2>/dev/null || true' EXIT
for i in $(seq 1 60); do curl -sf "http://${OLLAMA_HOST}/api/tags" >/dev/null && break; sleep 2; done
curl -sf "http://${OLLAMA_HOST}/api/tags" >/dev/null \
    || { echo "FATAL: ollama did not come up on ${OLLAMA_HOST}; see logs/ollama-${SLURM_JOB_ID:-local}.log"; exit 1; }
echo "[ollama] up on ${OLLAMA_HOST}"

echo "[ollama] pulling $MODEL (~230 GB first time, be patient)..."
ollama pull "$MODEL"

BKG_MODE=real BKG_MODELS="$MODEL" BKG_SEEDS="$SEED" BKG_FAMILIES="$FAMILIES" \
    jupyter nbconvert --to notebook --execute grounding_experiment/02_run/notebook.ipynb \
    --output "$REPO_DIR/logs/02_run_${SAFE}_s${SEED}.ipynb" --ExecutePreprocessor.timeout=-1

echo "DONE -> $OUT"
