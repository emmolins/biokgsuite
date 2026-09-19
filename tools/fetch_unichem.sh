#!/usr/bin/env bash
# PubChem CID -> ChEMBL and DrugBank -> ChEMBL crosswalks, from EBI UniChem.
#
# Why this is a separate script: Open Targets' drug crossReferences has no
# PubChem source, so without the first mapping any graph that keys drugs by
# PubChem CID (OpenBioLink) resolves no drugs at all and is scored as holding
# none. Its DrugBank source covers only about half of the approved molecules
# (varenicline, insulin glargine and enoxaparin are among the missing), so
# without the second mapping every DrugBank-keyed graph (PrimeKG, Hetionet,
# DRKG, BioKG) is capped at whatever the xref table happens to hold.
#
# Run this from a normal terminal on a machine with unrestricted internet.
# It will NOT work from the Claude session shells, which sit behind an egress
# allowlist that blocks ebi.ac.uk.
#
#   bash tools/fetch_unichem.sh
#
# UniChem source ids: 1 = ChEMBL, 2 = DrugBank, 22 = PubChem.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$ROOT/data/crosswalks"
BASE="https://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping/src_id1"
mkdir -p "$OUT"

fetch () {           # fetch <src file> <output stem>
  local url="$BASE/$1.gz" stem="$2"
  echo "Fetching $url"
  curl -fSL --retry 3 -o "$OUT/$1.gz" "$url"
  gunzip -c "$OUT/$1.gz" > "$OUT/$stem.tsv"
  rm -f "$OUT/$1.gz"
  local n; n=$(wc -l < "$OUT/$stem.tsv" | tr -d ' ')
  date -u +"retrieved_utc=%Y-%m-%dT%H:%M:%SZ" >  "$OUT/$stem.release.txt"
  echo "source_url=$url"                       >> "$OUT/$stem.release.txt"
  echo "rows=$n"                               >> "$OUT/$stem.release.txt"
  if command -v shasum >/dev/null 2>&1; then
    shasum -a 256 "$OUT/$stem.tsv"             >> "$OUT/$stem.release.txt"
  fi
  echo "Wrote $OUT/$stem.tsv  ($n rows); release record $OUT/$stem.release.txt"
}

fetch src1src22.txt unichem_pubchem_chembl
fetch src1src2.txt  unichem_drugbank_chembl

echo
echo "Then re-run:  python tools/provenance.py"
