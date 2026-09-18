#!/usr/bin/env bash
# One-time download of MONDO's SSSOM crosswalk file.
#
# This file is ~30 MB and contains mappings between MONDO and external
# disease vocabularies (UMLS, OMIM, Orphanet, ICD9, NCIT, MESH, DOID, ...).
# load_matrix() consults it as the third disease-ID bridge so MATRIX disease
# nodes whose primary identifiers are UMLS/OMIM/etc. can still be canonicalised
# to bare MONDO numerics.
#
# Source: https://github.com/monarch-initiative/mondo (current release branch)
#
# Run from the repo root:
#     bash scripts/download_mondo_sssom.sh

set -euo pipefail
DEST_DIR="data/gold_standards/bridges"
DEST_FILE="$DEST_DIR/mondo.sssom.tsv"
URL="https://raw.githubusercontent.com/monarch-initiative/mondo/master/src/ontology/mappings/mondo.sssom.tsv"

mkdir -p "$DEST_DIR"
echo "Downloading $URL"
curl -L --fail --progress-bar -o "$DEST_FILE" "$URL"
size=$(wc -c < "$DEST_FILE" | tr -d ' ')
lines=$(wc -l < "$DEST_FILE" | tr -d ' ')
echo "Saved $DEST_FILE  (${size} bytes, ${lines} lines)"
