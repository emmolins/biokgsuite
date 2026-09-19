#!/usr/bin/env bash
#SBATCH --job-name=biokg-ground-analysis
#SBATCH --partition=def
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=logs/ground-analysis-%j.out
#SBATCH --error=logs/ground-analysis-%j.err
#
# The grounding stages that are not model calls, executed in place so the
# committed notebooks carry their outputs.
#
#   STAGE=benchmark   01_benchmark only: the pair status, pools and selection
#                     audit, and the G01 checkpoint every run job needs.
#                     Submit this first and wait for it.
#   STAGE=analysis    (default) 02_run once more with no scoping, so every run
#                     file in runs/ is loaded and the G02 checkpoint lists them
#                     all, then 03_analysis, 04_scaling, 05_diagnostic.
#                     Submit after every run_grounding_api / run_llm_* job
#                     has finished.
#
#   sbatch --export=ALL,STAGE=benchmark hpc/run_grounding_analysis.sbatch
#   sbatch hpc/run_grounding_analysis.sbatch

set -euo pipefail
STAGE="${STAGE:-analysis}"
REPO_DIR="${REPO_DIR:-$HOME/biokgsuite-bib}"
cd "$REPO_DIR"; mkdir -p logs
for marker in config/kgs.yaml grounding_experiment src/biokgsuite; do
    [ -e "$marker" ] || { echo "FATAL: $REPO_DIR is not biokgsuite-bib (missing $marker)"; exit 1; }
done
[ ! -f .trial_kgs ] || { echo "FATAL: .trial_kgs is present; the cluster runs the full panel. Remove it."; exit 1; }
echo "HOST $(hostname)  DATE $(date)  JOB ${SLURM_JOB_ID:-none}  STAGE $STAGE"

module purge
module load python/3.12.11
module load uv/latest
source .venv/bin/activate
export PYTHONUNBUFFERED=1 TQDM_MININTERVAL=30 MPLBACKEND=Agg

run() {   # run <dir>: execute grounding_experiment/<dir>/notebook.ipynb in place
    echo; echo "=== grounding_experiment/$1  ($(date +%H:%M:%S))"
    ( cd "grounding_experiment/$1" && jupyter nbconvert --to notebook --execute --inplace \
        --log-level=INFO --ExecutePreprocessor.timeout=-1 notebook.ipynb )
}

case "$STAGE" in
  benchmark)
    run 01_benchmark ;;
  analysis)
    n=$(ls results/tables/G02_run/runs/*.csv 2>/dev/null | wc -l)
    [ "$n" -gt 0 ] || { echo "FATAL: no run files under results/tables/G02_run/runs/"; exit 1; }
    echo "$n run files present"
    BKG_MODE=real run 02_run          # every file exists: loaded, not recomputed; writes the full checkpoint
    run 03_analysis
    run 04_scaling
    run 05_diagnostic ;;
  *) echo "FATAL: STAGE must be benchmark or analysis"; exit 1 ;;
esac

echo; echo ">>> tests"; python -m pytest tests -q -k grounding
echo "DONE  $(date)"
