#!/usr/bin/env bash
#SBATCH --job-name=biokg-ground
#SBATCH --partition=def
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --dependency=singleton
#SBATCH --output=logs/ground-%j.out
#SBATCH --error=logs/ground-%j.err
#
# The grounding experiment (grounding_experiment/02_run) for ONE hosted model
# and ONE seed, both families (main: 128 pairs; common: the paired subset).
# No GPU: the work is API calls, the memory is the three graphs' dossier
# indexes. About 2,000 calls per job at 5 to 10 s each, so 3 to 6 hours.
#
# Submit one job per (model, seed). Give the jobs of one provider the same
# --job-name: with --dependency=singleton Slurm runs one job of that name at
# a time, so the three seeds of a provider run in sequence against its rate
# limit while the three providers run in parallel.
#
#   cd ~/biokgsuite-bib
#   for S in 0 1 2; do
#     sbatch --job-name=ground-gpt    --export=ALL,MODEL=gpt-4.1-mini,SEED=$S                    hpc/run_grounding_api.sbatch
#     sbatch --job-name=ground-gemini --export=ALL,MODEL=gemini:gemini-3.5-flash-lite,SEED=$S    hpc/run_grounding_api.sbatch
#   done
#   Llama 3.3 70B runs locally instead (hpc/run_llm_ladder.sbatch with
#   FAMILIES=main+common); Groq's hosted 70B is enterprise-only since 2026.
#
# Monitor:   squeue -u $USER;  tail -f logs/ground-<jobid>.out
# Resumable: a (family, model, seed) file that exists is loaded, not rerun, so
#            resubmit a job that died or was rate-limited (the log says
#            "incomplete (rate-limited) arms") and it completes the file.
# Keys:      .env at the repo root (OPENAI_API_KEY, GEMINI_API_KEY).

set -euo pipefail
MODEL="${MODEL:?set MODEL, e.g. --export=ALL,MODEL=gpt-4.1-mini,SEED=0}"
SEED="${SEED:?set SEED, e.g. --export=ALL,MODEL=gpt-4.1-mini,SEED=0}"
FAMILIES="${FAMILIES:-main,common}"
SAFE="$(echo "$MODEL" | tr 'A-Z' 'a-z' | sed -E 's/[^a-z0-9]+/_/g; s/^_//; s/_$//')"   # same slug as runner.model_slug

REPO_DIR="${REPO_DIR:-$HOME/biokgsuite-bib}"
cd "$REPO_DIR"; mkdir -p logs
for marker in config/kgs.yaml grounding_experiment/02_run/notebook.ipynb src/biokgsuite; do
    [ -e "$marker" ] || { echo "FATAL: $REPO_DIR is not biokgsuite-bib (missing $marker)"; exit 1; }
done
[ -f .env ] || { echo "FATAL: no .env at $REPO_DIR (copy it from the laptop; see hpc/README.md)"; exit 1; }
[ -f results/checkpoints/G01_benchmark.pkl ] || { echo "FATAL: run grounding_experiment/01_benchmark first (hpc/run_grounding_analysis.sbatch STAGE=benchmark)"; exit 1; }
[ ! -f .trial_kgs ] || { echo "FATAL: .trial_kgs is present; the cluster runs the full panel. Remove it."; exit 1; }

echo "HOST $(hostname)  DATE $(date)  JOB ${SLURM_JOB_ID:-none}  MODEL $MODEL  SEED $SEED  FAMILIES $FAMILIES"
module purge
module load python/3.12.11
module load uv/latest
source .venv/bin/activate
export PYTHONUNBUFFERED=1 TQDM_MININTERVAL=30 MPLBACKEND=Agg

# Each job writes its own executed copy under logs/; the notebook in the tree
# is never edited in place here, because several jobs run it at once.
BKG_MODE=real BKG_MODELS="$MODEL" BKG_SEEDS="$SEED" BKG_FAMILIES="$FAMILIES" \
    jupyter nbconvert --to notebook --execute grounding_experiment/02_run/notebook.ipynb \
        --output "$REPO_DIR/logs/02_run_${SAFE}_s${SEED}.ipynb" --ExecutePreprocessor.timeout=-1

echo "DONE  $(date)"
for fam in ${FAMILIES//,/ }; do ls -la "results/tables/G02_run/runs/${fam}__${SAFE}__seed${SEED}.csv"; done
