#!/usr/bin/env bash
# One-time HPC setup. Run ONCE on the login node, from the repo root, before
# submitting any job. Idempotent: safe to rerun if it fails partway.
#
#   cd ~/biokgsuite && bash hpc/setup.sh
#
# Creates .venv (uv, Python 3.12) and installs requirements.txt plus the
# embedding stack (requirements-embeddings.txt, torch) that only
# supplementary/embeddings needs. Then lists the manual steps: the KG data
# (rsync from the laptop), the API keys (.env) and, for the embedding
# supplementary, the Hugging Face login.

set -euo pipefail

REPO_DIR="${REPO_DIR:-$HOME/biokgsuite}"
cd "$REPO_DIR"
for marker in config/kgs.yaml biokgsuite_dimensions src/biokgsuite requirements.txt; do
    [ -e "$marker" ] || { echo "FATAL: $REPO_DIR is not biokgsuite (missing $marker)"; exit 1; }
done

echo "== BioKGSuite HPC setup in $REPO_DIR"
module purge
module load python/3.12.11
module load uv/latest
echo "   python $(python --version 2>&1)   uv $(uv --version)"

if [[ ! -d .venv ]]; then
    uv venv --python 3.12
fi
source .venv/bin/activate
uv pip install -r requirements.txt
uv pip install -r requirements-embeddings.txt transformers sentence-transformers "huggingface_hub[cli]"
uv pip install requests            # tools/fetch_*.py; not needed by the notebooks
python -c "import pandas, numpy, scipy, sklearn, networkx, matplotlib, pyarrow, yaml, tqdm, nbconvert" \
    && echo "   notebook stack imports"
python -c "import torch; print('   torch', torch.__version__, 'cuda build', torch.version.cuda)"

mkdir -p logs data results/tables/G02_run/runs
echo
echo "== Done. Before submitting:"
echo "   1. KG data: from the laptop,"
echo "        rsync -avz --progress --exclude '_cache_matrix_*' --exclude 'raw' \\"
echo "              ~/Desktop/biokgsuite/data/ ashburn-login:biokgsuite/data/"
echo "   2. API keys for the grounding runs: from the laptop,"
echo "        scp ~/Desktop/biokgsuite/.env ashburn-login:biokgsuite/.env"
echo "   3. Embedding supplementary only: huggingface-cli login, and accept the"
echo "      EmbeddingGemma licence at https://huggingface.co/google/embeddinggemma-300m"
echo "   4. Check the graphs arrived:  du -sh data/*/"
echo "   Then follow hpc/README.md."
