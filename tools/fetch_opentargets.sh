#!/usr/bin/env bash
# Fetch the Open Targets Platform datasets the entity reference sets are built from.
#
#   bash tools/fetch_opentargets.sh [dest_dir]
#
# Needs a network that can reach ftp.ebi.ac.uk. Uses curl only.
# About 125 MB with the default GENE_DISEASE=curated; see the block below for
# the other gene-disease evidence options. Dataset names follow the 26.06
# index (there is no `known_drug` or `drug_indication` dataset in that release).
set -euo pipefail

REL="${OT_RELEASE:-26.06}"
BASE="https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/${REL}/output"
DEST="${1:-data/crosswalks}"

# Small core (~110 MB total)
CORE="drug_molecule clinical_indication drug_mechanism_of_action disease target reactome"
# --- disease->gene evidence (gold #3, and input to #4) -----------------------
# Pick ONE. The aggregated association table is 1.1 GB because it summarises
# ~9.3 GB of literature-mined evidence (evidence_europepmc, 125 files) that you
# do not want in a gold standard.
#
#   GENE_DISEASE=curated  (default, ~15 MB)
#       Expert-curated gene-disease assertions only: ClinGen, Orphanet,
#       Gene2Phenotype, Genomics England, UniProt, gene burden. Gene-level,
#       disease-specific, highest signal-to-noise.
#   GENE_DISEASE=clinvar  (~440 MB)  curated + ClinVar (evidence_eva).
#   GENE_DISEASE=genetic  (~720 MB)  + GWAS credible sets. The full
#       `genetic_association` datatype. GWAS is common-variant statistical
#       association, a different evidence class from a curated gene assertion.
#   GENE_DISEASE=all      (~1.1 GB)  the aggregated association table, which
#       includes literature mining. Not recommended as a gold standard.
CURATED="evidence_clingen evidence_orphanet evidence_gene2phenotype evidence_genomics_england evidence_uniprot_literature evidence_uniprot_variants evidence_gene_burden"
case "${GENE_DISEASE:-curated}" in
  curated) GD="$CURATED" ;;
  clinvar) GD="$CURATED evidence_eva" ;;
  genetic) GD="$CURATED evidence_eva evidence_gwas_credible_sets" ;;
  all)     GD="association_by_datatype_direct" ;;
  none)    GD="" ;;
  *) echo "GENE_DISEASE must be curated|clinvar|genetic|all|none"; exit 1 ;;
esac
if [ "${CORE_ONLY:-0}" = "1" ]; then DATASETS="$CORE"; else DATASETS="$CORE $GD"; fi

command -v curl >/dev/null || { echo "curl not found"; exit 1; }
if command -v sha256sum >/dev/null; then SHA=sha256sum; else SHA="shasum -a 256"; fi

mkdir -p "$DEST"
echo "Open Targets Platform ${REL}  ->  ${DEST}"
echo

total=0
failed=""
for d in $DATASETS; do
  printf '  %-34s ' "$d"
  dir="${DEST}/ot_${d}"
  mkdir -p "$dir"
  # Apache autoindex -> every .parquet href. Failures must not abort the run,
  # so the whole pipeline is guarded and the dataset is recorded and skipped.
  if ! listing=$(curl -fsSL --retry 2 "${BASE}/${d}/" 2>/dev/null); then
    echo "UNREACHABLE (check ${BASE}/${d}/)"; failed="$failed $d"; continue
  fi
  # Open Targets is inconsistent: most datasets use Spark's
  # part-00000-<uuid>-c000.snappy.parquet, but `disease` and
  # `clinical_indication` ship a single <dataset>.parquet. Match ANY .parquet
  # href, excluding absolute paths (the Parent Directory link).
  files=$(printf '%s' "$listing" | grep -o 'href="[^"/][^"]*\.parquet"' \
          | sed 's/^href="//; s/"$//' | sort -u || true)
  if [ -z "$files" ]; then
    echo "NO PARQUET FOUND"; failed="$failed $d"; continue
  fi
  want=$(printf '%s\n' $files | wc -w | tr -d ' ')
  n=0; bad=0
  for f in $files; do
    n=$((n+1))
    if [ ! -s "${dir}/${f}" ]; then
      printf '\r  %-34s %d/%s downloading...' "$d" "$n" "$want"
      curl -fsSL --retry 3 -o "${dir}/${f}" "${BASE}/${d}/${f}" || bad=$((bad+1))
    fi
  done
  printf '\r%-70s\r' " "
  sz=$(du -sh "$dir" 2>/dev/null | cut -f1)
  if [ "$bad" -gt 0 ]; then
    printf '  %-34s %s file(s), %s  (%d FAILED)\n' "$d" "$want" "$sz" "$bad"
    failed="$failed $d"
  else
    printf '  %-34s %s file(s), %s\n' "$d" "$want" "$sz"
  fi
  total=$((total+want))
done

echo
cat > "${DEST}/OPENTARGETS_RELEASE.txt" <<EOF
Open Targets Platform release: ${REL}
Source:    ${BASE}
Retrieved: $(date -u +%Y-%m-%dT%H:%M:%SZ)
Datasets:  ${DATASETS}
Files:     ${total}
Licence:   CC0 1.0
Citation:  Open Targets Platform. Nucleic Acids Research 2025;53(D1):D1467.
EOF
if ls "$DEST"/ot_*/*.parquet >/dev/null 2>&1; then
  ( cd "$DEST" && $SHA ot_*/*.parquet > OPENTARGETS_CHECKSUMS.txt )
  echo "Wrote ${DEST}/OPENTARGETS_RELEASE.txt  and  OPENTARGETS_CHECKSUMS.txt"
else
  echo "No parquet files downloaded — nothing to checksum."
fi
echo
if [ -n "$failed" ]; then
  echo
  echo "INCOMPLETE — these datasets did not download:$failed"
  echo "Re-run the same command; finished files are skipped."
fi
echo
echo "Next:  python tools/provenance.py && python tools/build_reference_sets.py"
echo "       (drug_drug gold, if kept: download DDInter 2.0 manually from"
echo "        https://ddinter.scbdd.com/download/ into ${DEST}/ddinter/)"
