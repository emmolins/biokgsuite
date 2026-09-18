#!/usr/bin/env bash
#SBATCH --job-name=biokg-llm-ollama
#SBATCH --partition=def
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --gres=gpu:H100:1
#SBATCH --output=logs/llm-%j.out
#SBATCH --error=logs/llm-%j.err
#
# Runs the LLM ranking sweep (eval_notebooks/run_llm_sweep.py) against a LOCAL
# llama-3.3-70b served by Ollama on the H100 — no API, no rate limits.
#
# Prereqs (one-time, see walkthrough): repo synced, .venv built (scripts/hpc/setup.sh),
# and the 3 KG data dirs present (data/primekg, data/drkg, data/biokg).
#
# Submit:   sbatch scripts/hpc/run_llm_ollama.sbatch
# Monitor:  tail -f logs/llm-<jobid>.out
# Resumable: FORCE_RERUN=True recomputes; just resubmit if it dies.

set -euo pipefail
REPO_DIR="${REPO_DIR:-$HOME/biokgsuite}"
cd "$REPO_DIR"; mkdir -p logs

echo "HOST $(hostname)  DATE $(date)  JOB ${SLURM_JOB_ID:-none}"
module purge 2>/dev/null || true
module load python/3.12.11 2>/dev/null || module load python 2>/dev/null || true
module load uv/latest 2>/dev/null || true
source .venv/bin/activate
# fail loudly only if the venv python is actually unusable on this node
python -c "import sys; print('python', sys.version.split()[0])" || { echo "FATAL: venv python unusable on $(hostname)"; exit 1; }
export PYTHONUNBUFFERED=1 TQDM_MININTERVAL=15

# ── Ollama: install to ~/.local if missing (no root needed) ─────────────────
export OLLAMA_MODELS="${OLLAMA_MODELS:-$HOME/.ollama/models}"   # ~40 GB lives here
export OLLAMA_HOST=127.0.0.1:11434
mkdir -p "$HOME/.local/bin" "$OLLAMA_MODELS"
if ! command -v ollama >/dev/null 2>&1 && [[ ! -x "$HOME/.local/bin/ollama" ]]; then
    echo "[ollama] downloading standalone binary (GitHub release, .tar.zst)..."
    curl -fL https://github.com/ollama/ollama/releases/latest/download/ollama-linux-amd64.tar.zst -o /tmp/ollama.tar.zst
    ( tar --zstd -xf /tmp/ollama.tar.zst -C "$HOME/.local" 2>/dev/null \
      || tar -I zstd -xf /tmp/ollama.tar.zst -C "$HOME/.local" )   # -> ~/.local/bin/ollama + ~/.local/lib/ollama
fi
export PATH="$HOME/.local/bin:$PATH"
nvidia-smi || true

# ── start the server, wait for it ──────────────────────────────────────────
ollama serve > logs/ollama-${SLURM_JOB_ID:-local}.log 2>&1 &
OLLAMA_PID=$!
trap 'kill $OLLAMA_PID 2>/dev/null || true' EXIT
for i in $(seq 1 30); do
    curl -sf http://127.0.0.1:11434/api/tags >/dev/null && break
    sleep 2
done
echo "[ollama] up"

# ── pull the model (idempotent; ~40 GB first time) ──────────────────────────
ollama pull llama3.3:70b

# ── run the sweep (bare tag -> harness routes to local Ollama) ───────────────
SEED="${SEED:-1}"   # one seed per job (~5h); submit once per seed via --export
python eval_notebooks/run_llm_sweep.py --model llama3.3:70b --gold gold_standard_v4_bigset.tsv --seeds "$SEED" --shuffles 2 --out "09_big_llama_s${SEED}.csv"

echo "DONE -> results/tables/09_llm_runs/09_ranking_llama3.3_70b.csv"
