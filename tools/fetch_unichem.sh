#!/usr/bin/env bash
# PubChem CID -> ChEMBL crosswalk, from EBI UniChem.
#
# Why this is a separate script: Open Targets' drug crossReferences has no
# PubChem source, so without this mapping any graph that keys drugs by PubChem
# CID (OpenBioLink) resolves no drugs at all and is scored as holding none.
#
# Run this from a normal terminal on a machine with unrestricted internet.
# It will NOT work from the Claude session shells, which sit behind an egress
# allowlist that blocks ebi.ac.uk.
#
#   bash tools/fetch_unichem.sh
#
# UniChem source ids: 1 = ChEMBL, 22 = PubChem.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$ROOT/data/crosswalks"
URL="https://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1/src1src22.txt.gz"

mkdir -p "$OUT"
echo "Fetching UniChem ChEMBL<->PubChem mapping"
echo "  $URL"
curl -fSL --retry 3 -o "$OUT/src1src22.txt.gz" "$URL"

gunzip -c "$OUT/src1src22.txt.gz" > "$OUT/unichem_pubchem_chembl.tsv"
rm -f "$OUT/src1src22.txt.gz"

N=$(wc -l < "$OUT/unichem_pubchem_chembl.tsv" | tr -d ' ')
date -u +"retrieved_utc=%Y-%m-%dT%H:%M:%SZ" >  "$OUT/unichem_pubchem_chembl.release.txt"
echo "source_url=$URL"                        >> "$OUT/unichem_pubchem_chembl.release.txt"
echo "rows=$N"                                >> "$OUT/unichem_pubchem_chembl.release.txt"
if command -v shasum >/dev/null 2>&1; then
  shasum -a 256 "$OUT/unichem_pubchem_chembl.tsv" >> "$OUT/unichem_pubchem_chembl.release.txt"
fi

echo "Wrote $OUT/unichem_pubchem_chembl.tsv  ($N rows)"
echo "Release record: $OUT/unichem_pubchem_chembl.release.txt"
echo
echo "Then re-run:  python tools/provenance.py"
