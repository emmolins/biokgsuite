#!/usr/bin/env bash
#SBATCH --job-name=biokg-llm-76
#SBATCH --partition=def
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --gres=gpu:H100:1
#SBATCH --output=logs/llm76-%j.out
#SBATCH --error=logs/llm76-%j.err
#
# llama-3.3-70b on the 76-pair COMMON set (gold_standard_v3_common3.tsv),
# seeds 0,1,2 x 3 shuffles to match the gpt run exactly. Uses the CURRENT
# notebook (OT distractor pool) + the rank_call->ollama fix.

set -euo pipefail
REPO_DIR="${REPO_DIR:-$HOME/biokgsuite}"
cd "$REPO_DIR"; mkdir -p logs

echo "HOST $(hostname)  DATE $(date)  JOB ${SLURM_JOB_ID:-none}"
module purge
module load python/3.12.11
module load uv/latest 2>/dev/null || true
source .venv/bin/activate
export PYTHONUNBUFFERED=1 TQDM_MININTERVAL=15

export OLLAMA_MODELS="${OLLAMA_MODELS:-$HOME/.ollama/models}"
export OLLAMA_HOST=127.0.0.1:11434
mkdir -p "$HOME/.local/bin" "$OLLAMA_MODELS"
if ! command -v ollama >/dev/null 2>&1 && [[ ! -x "$HOME/.local/bin/ollama" ]]; then
    echo "[ollama] downloading standalone binary..."
    curl -fL https://github.com/ollama/ollama/releases/latest/download/ollama-linux-amd64.tar.zst -o /tmp/ollama.tar.zst
    ( tar --zstd -xf /tmp/ollama.tar.zst -C "$HOME/.local" 2>/dev/null \
      || tar -I zstd -xf /tmp/ollama.tar.zst -C "$HOME/.local" )
fi
export PATH="$HOME/.local/bin:$PATH"
nvidia-smi || true

ollama serve > logs/ollama76-${SLURM_JOB_ID:-local}.log 2>&1 &
OLLAMA_PID=$!
trap 'kill $OLLAMA_PID 2>/dev/null || true' EXIT
for i in $(seq 1 30); do
    curl -sf http://127.0.0.1:11434/api/tags >/dev/null && break
    sleep 2
done
echo "[ollama] up"
ollama pull llama3.3:70b

python eval_notebooks/run_llm_sweep.py --model llama3.3:70b \
    --gold gold_standard_v3_common3.tsv --seeds 0,1,2 --shuffles 3 \
    --out 09_big_llama_76.csv

echo "DONE -> results/tables/09_llm_runs/09_big_llama_76.csv"
