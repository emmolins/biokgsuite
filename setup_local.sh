#!/usr/bin/env bash
# One-time local setup. No conda required.
#
#   bash setup_local.sh                # analysis environment (notebooks 00-07)
#   bash setup_local.sh --embeddings   # also installs torch, for supplementary/
#
# Creates .venv in this directory and installs from requirements.txt. Uses `uv`
# when it is on PATH (same tool hpc/setup.sh uses on the cluster) and falls back
# to python3 -m venv + pip otherwise. Idempotent: safe to re-run.
#
# Afterwards:
#   source .venv/bin/activate
#   python -m pytest tests -q
set -euo pipefail

cd "$(dirname "$0")"
WITH_EMB=0
[ "${1:-}" = "--embeddings" ] && WITH_EMB=1

# The repo ships a bare Python 3.9 Jupyter kernel at .venv; it has ipykernel and
# nothing else, so it cannot run the notebooks. Rebuild it rather than install
# into it, so the Python version matches environment.yml (3.11) where possible.
if [ -d .venv ] && [ ! -f .venv/.biokgsuite-setup ]; then
    echo "Existing .venv is not one this script made; moving it to .venv.bak"
    rm -rf .venv.bak && mv .venv .venv.bak
fi

if command -v uv >/dev/null 2>&1; then
    echo "[1/3] uv found; creating .venv on Python 3.11"
    uv venv --python 3.11 .venv 2>/dev/null || uv venv .venv
    PIP=(uv pip install --python .venv/bin/python)
else
    PY="$(command -v python3.11 || command -v python3.12 || command -v python3)"
    echo "[1/3] uv not found; using $PY ($("$PY" --version))"
    "$PY" -m venv .venv
    PIP=(.venv/bin/python -m pip install)
    .venv/bin/python -m pip install -q --upgrade pip
fi

echo "[2/3] installing requirements.txt"
"${PIP[@]}" -r requirements.txt

if [ "$WITH_EMB" = "1" ]; then
    echo "[2b/3] installing requirements-embeddings.txt (torch; this is large)"
    "${PIP[@]}" -r requirements-embeddings.txt
fi

touch .venv/.biokgsuite-setup

echo "[3/3] verifying"
.venv/bin/python - <<'PYCHECK'
import importlib, sys
sys.path.insert(0, 'src')
need = ['pandas', 'numpy', 'scipy', 'sklearn', 'networkx', 'matplotlib',
        'yaml', 'pyarrow', 'nbconvert', 'pytest']
missing = [m for m in need if not importlib.util.find_spec(m)]
if missing:
    raise SystemExit(f'MISSING: {missing}')
import biokgsuite.stats.stability, biokgsuite.io           # noqa: F401
print(f'  python {sys.version.split()[0]}  — imports OK')
PYCHECK

cat <<'DONE'

Ready. Next:
  source .venv/bin/activate
  python -m pytest tests -q

  # run one dimension notebook (05 takes 10-20 min)
  cd biokgsuite_dimensions/05_stability
  jupyter nbconvert --to notebook --execute --inplace \
      --ExecutePreprocessor.timeout=86400 notebook.ipynb
DONE
