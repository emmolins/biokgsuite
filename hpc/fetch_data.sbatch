#!/usr/bin/env bash
#SBATCH --job-name=biokg-fetch
#SBATCH --partition=def
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=logs/fetch-%j.out
#SBATCH --error=logs/fetch-%j.err
#
# Fetch the graphs on the cluster instead of pushing them from a laptop:
# tools/fetch_kgs.py downloads each one at its pinned release, verifies the
# publisher checksum where one exists, derives the files config/kgs.yaml
# reads and records a SNAPSHOT under data/raw/<kg>/<date>/. Then every
# derived file is hashed and compared with the record committed in
# data/<kg>/<KG>_RELEASE.txt, which is the laptop's copy, so the cluster is
# known to hold byte-identical inputs before anything runs on them.
#
#   sbatch hpc/fetch_data.sbatch                                  # primekg openbilink matrix
#   sbatch --export=ALL,KGS="primekg drkg" hpc/fetch_data.sbatch  # a subset
#
# Hetionet and BioKG are not fetched by tools/fetch_kgs.py (they were retrieved
# by hand, see data/PROVENANCE.md); copy those two from the laptop with rsync.
# MATRIX is tracked at Hugging Face main, so if Every Cure has pushed since the
# laptop's fetch the derived files will differ: the job then keeps any
# edges.tsv that was already on disk as edges.tsv.laptop and says so, rather
# than silently running on a different graph.

set -euo pipefail
KGS="${KGS:-primekg openbilink matrix}"
REPO_DIR="${REPO_DIR:-$HOME/biokgsuite}"
cd "$REPO_DIR"; mkdir -p logs
for marker in config/kgs.yaml tools/fetch_kgs.py src/biokgsuite; do
    [ -e "$marker" ] || { echo "FATAL: $REPO_DIR is not biokgsuite (missing $marker)"; exit 1; }
done
echo "HOST $(hostname)  DATE $(date)  JOB ${SLURM_JOB_ID:-none}  KGS $KGS"
df -h "$HOME" | tail -1

module purge
module load python/3.12.11
module load uv/latest
source .venv/bin/activate
export PYTHONUNBUFFERED=1

# A laptop copy of the MATRIX edge file is worth 11 GB of upload; keep it
# until the fetched one is shown to be identical.
if [[ " $KGS " == *" matrix "* && -f data/matrix/edges.tsv && ! -f data/matrix/edges.tsv.laptop ]]; then
    mv data/matrix/edges.tsv data/matrix/edges.tsv.laptop
    rm -f data/matrix/nodes.tsv           # a partial rsync leaves a truncated file under the final name
fi

for kg in $KGS; do
    echo; echo "=== fetch $kg  ($(date +%H:%M:%S))"
    python tools/fetch_kgs.py --kg "$kg"
done

echo; echo "=== verify derived files against the committed release records"
python - "$KGS" <<'PY'
import hashlib, subprocess, sys
from pathlib import Path
kgs = sys.argv[1].split()
bad = 0
for kg in kgs:
    rel = f'data/{kg}/{kg.upper()}_RELEASE.txt'
    committed = subprocess.run(['git', 'show', f'HEAD:{rel}'], capture_output=True, text=True).stdout
    expected, in_block = {}, False
    for line in committed.splitlines():
        if line.startswith('derived_files:'):
            in_block = True; continue
        if in_block and line.startswith('  data/'):
            parts = line.split()
            expected[parts[0]] = parts[parts.index('sha256') + 1]
        elif in_block and not line.startswith('  '):
            in_block = False
    for path, sha in expected.items():
        h = hashlib.sha256()
        with open(path, 'rb') as fh:
            for chunk in iter(lambda: fh.read(1 << 24), b''):
                h.update(chunk)
        got = h.hexdigest()[:len(sha)]
        ok = got == sha
        bad += not ok
        print(f'{"OK " if ok else "DIFF"} {path}  committed {sha}  fetched {got}')
sys.exit(1 if bad else 0)
PY
status=$?
if [ "$status" -eq 0 ]; then
    rm -f data/matrix/edges.tsv.laptop
    echo; echo "all derived files match the laptop's records; the laptop edge copy was not needed and was removed"
else
    echo; echo "MISMATCH: at least one fetched file differs from the committed record."
    echo "  For MATRIX that means Hugging Face main has moved since the laptop's fetch. Use the laptop's files:"
    echo "    mv data/matrix/edges.tsv.laptop data/matrix/edges.tsv   (if the edge file is the one that differs)"
    echo "  and rsync the other differing files from the laptop; do not run the panel on mixed inputs."
    exit 1
fi
echo "DONE  $(date)"
