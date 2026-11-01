#!/usr/bin/env bash
# Fetch the two CTD tables notebook 06's multi-hop benchmark is built from, and
# record enough provenance to cite a specific release.
#
#   bash tools/fetch_ctd.sh
#
# CTD publishes monthly and does NOT put a version in the filename, so the file
# alone is not citable. Every CTD report begins with '#' header lines carrying
# the field list and the report date; this script keeps those lines, the URL, the
# fetch date and a SHA-256 for each file in
# data/gold_standards/ctd/CTD_PROVENANCE.txt. Cite the report date from there, not
# the day you happened to download it.
#
# CTD terms of use: free for academic and research use with attribution. Check
# https://ctdbase.org/about/legal.jsp before redistributing these files, and do
# NOT commit them; data/ is inputs, not outputs.
#
# Needs a network that can reach ctdbase.org.
set -euo pipefail

cd "$(dirname "$0")/.."
GS=data/gold_standards/ctd
mkdir -p "$GS"

BASE_URL="${CTD_BASE_URL:-https://ctdbase.org/reports}"
FILES=(
  "CTD_chem_gene_ixns.tsv.gz"
  "CTD_curated_genes_diseases.tsv.gz"
)

PROV="$GS/CTD_PROVENANCE.txt"
PROV_TMP="$(mktemp)"
trap 'rm -f "$PROV_TMP"' EXIT
{
  echo "# CTD provenance for BioKGSuite"
  echo "# fetched:  $(date -u '+%Y-%m-%dT%H:%M:%SZ') UTC"
  echo "# base URL: $BASE_URL"
  echo "# NOTE: cite the 'Report created' date in the header block below,"
  echo "#       not the fetch date. CTD releases monthly."
  echo
} > "$PROV_TMP"

for f in "${FILES[@]}"; do
    echo ">>> $f"
    if ! curl -fSL --retry 3 --connect-timeout 20 -o "$GS/$f" "$BASE_URL/$f"; then
        echo
        echo "FAILED to fetch $f."
        echo "  If this is a network block, run this script from an unrestricted"
        echo "  network. If it is a 404, CTD has renamed the report: check"
        echo "  https://ctdbase.org/downloads/ and set CTD_BASE_URL or rename"
        echo "  the file to match what notebook 06 expects."
        exit 1
    fi
    {
        echo "== $f"
        echo "   url:    $BASE_URL/$f"
        echo "   bytes:  $(wc -c < "$GS/$f")"
        echo "   sha256: $( (sha256sum "$GS/$f" 2>/dev/null || shasum -a 256 "$GS/$f") | awk '{print $1}')"
        echo "   release: $( ( set +o pipefail
            gzip -dc "$GS/$f" 2>/dev/null | head -400 \
              | grep -m1 -i '^# Report created:' ) || true )"
        echo "   rows:    $( ( set +o pipefail
            gzip -dc "$GS/$f" 2>/dev/null | grep -vc '^#' ) || true )"
        echo "   header:"
        # `set -o pipefail` plus `head` closing the pipe early makes gzip die of
        # SIGPIPE and takes the whole script down with it, on a file that
        # downloaded perfectly. Read a bounded prefix in a subshell instead.
        ( set +o pipefail
          gzip -dc "$GS/$f" 2>/dev/null | head -400 | grep '^#' | head -20 \
            | sed 's/^/     /' ) || true
        echo
    } >> "$PROV_TMP"
done

# Only now, with both files fetched, does the provenance record become
# true. Writing it up front left a file claiming provenance for data that
# had not arrived.
mv "$PROV_TMP" "$PROV"
trap - EXIT

echo
echo "-> $GS/"
ls -la "$GS" | sed 's/^/   /'
echo
echo "Provenance written to $PROV. Record the 'Report created' date from it in"
echo "data/GOLD_STANDARDS.md and in the manuscript's data-availability section."
echo
echo "Citation (CTD asks to be cited as):"
echo "  Davis AP, Wiegers TC, Johnson RJ, Sciaky D, Wiegers J, Mattingly CJ."
echo "  Comparative Toxicogenomics Database (CTD): update 2023."
echo "  Nucleic Acids Research. 2023;51(D1):D1257-D1262."
echo "  Verify the current citation at https://ctdbase.org/about/publications/"
