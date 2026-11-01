#!/usr/bin/env bash
# Run EmbeddingGemma-300m as a word-priors baseline across the 5 small KGs.
# MATRIX runs separately via run_gemma_matrix.py (it needs subsampling /
# larger batches).
#
# Prereqs:
#   1. cd into the biokgsuite repo root
#   2. conda activate biokgsuite  (or your venv)
#   3. pip install torch transformers sentence-transformers
#   4. Accept the EmbeddingGemma license at:
#        https://huggingface.co/google/embeddinggemma-300m
#   5. export HF_TOKEN=<your_token>
#
# Expected runtime (CPU, M-series Mac or x86_64 server):
#   Hetionet   ~47K entities  →  ~15–30 min
#   PrimeKG    ~130K entities → ~45–75 min
#   DRKG       ~100K entities → ~30–60 min  (ID-only names, see caveat)
#   OpenBioL.  ~180K entities → ~60–90 min  (ID-only names, see caveat)
#   BioKG      ~100K entities → ~30–60 min  (ID-only names, see caveat)
#
# Total: roughly 3–5 hours on CPU. On a single A100, ~10–20× faster.
#
# Caveat for DRKG / OpenBioLink / BioKG:
#   These KGs' loaders return opaque ID strings (e.g. "2157", "DB00001",
#   "0001234") instead of human-readable names. EmbeddingGemma will embed
#   those literal strings, producing near-noise; they need external name
#   resolution to be meaningful. Run them anyway, the contrast (high signal on
#   Hetionet/PrimeKG, low signal on the others) confirms the experimental
#   hypothesis that the signal lives in word priors, not arbitrary IDs.

set -euo pipefail

if [[ -z "${HF_TOKEN:-}" ]]; then
    echo "ERROR: HF_TOKEN not set. Export your HuggingFace token first." >&2
    echo "  export HF_TOKEN=hf_..." >&2
    exit 1
fi
export HUGGING_FACE_HUB_TOKEN="$HF_TOKEN"

DIM=${GEMMA_DIM:-768}   # override with GEMMA_DIM=256 etc for Matryoshka

# Order: smallest first, so failures surface quickly.
for KG in hetionet biokg drkg primekg openbilink ; do
    echo
    echo "===================================================="
    echo "Running EmbeddingGemma on $KG  (dim=$DIM)"
    echo "===================================================="
    python scripts/run_emb_model.py "$KG" Gemma 0 "$DIM"
done

echo
echo "Done. Results cached to results/cache/embedding_<kg>.json"
echo "Re-execute supplementary/embeddings/notebook.ipynb to update figures."
