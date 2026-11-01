"""Build a full DrugBank cross-reference table for MATRIX.

Output: ``data/gold_standards/bridges/drugbank_xref.csv`` with columns
    drugbank_id, namespace, external_id

Sources (in order of priority):
1. DrugCentral identifiers dump (free academic): pulls UNII / RxCUI / ATC /
   KEGG.DRUG / ChEBI / ChEMBL / CAS / PubChem mappings to DrugBank.
2. Existing ``data/openbilink/pubchem_to_drugbank.csv`` (PubChem only).
3. Optional: a user-provided ``drugbank_external_links.csv`` (download from
   drugbank.com with a free academic account at
   https://go.drugbank.com/releases/latest#external-links).

Usage:
    python scripts/build_drugbank_xref.py

Re-running is idempotent — safe to invoke after each DrugBank refresh.
"""
from __future__ import annotations
import argparse
import csv
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
GOLD = REPO_ROOT / "data" / "gold_standards" / "bridges"
OUT  = GOLD / "drugbank_xref.csv"

# DrugCentral's downloads page lives at https://drugcentral.org/Downloads
# The identifier crosswalk is `identifier.tsv.gz`, a TSV with columns:
#   identifier, id_type, struct_id, parent_match
# struct_id corresponds to a DrugCentral structure; we need the join to the
# `structure_external_id` table. Easier path: download `drug.target.interaction.tsv`
# which already includes DrugBank IDs and external_id pairs in one row. Or pull
# the full Postgres dump (drugcentral.dump.07292024.sql.gz), too big for here.
#
# For a self-contained free pull we use the public **Drugs@FDA / RxNorm
# crosswalk** maintained by NIH at https://rxnav.nlm.nih.gov, accessible via
# their REST API but requires per-drug calls. Skip that for batch builds.
#
# Pragmatic compromise: this script reads three local files if present:
#   1. data/openbilink/pubchem_to_drugbank.csv  (always present in this repo)
#   2. data/gold_standards/bridges/drugbank_external_links.csv  (user downloads from
#      drugbank.com, has all 16k DrugBank drugs with UNII/RxNorm/ATC/KEGG/
#      ChEBI/ChEMBL/CAS columns)
#   3. data/gold_standards/bridges/drugcentral_xref.csv  (optional, from DrugCentral)
# It then unions them into the canonical drugbank_xref.csv format.

# --- Column mappings for DrugBank's external_links.csv (the official format) -
# When a user downloads "external links" from drugbank.com, they get a CSV with
# columns like: 'DrugBank ID', 'Name', 'CAS Number', 'Drug Type', 'KEGG Drug ID',
# 'KEGG Compound ID', 'PubChem Compound ID', 'PubChem Substance ID',
# 'ChEBI ID', 'PharmGKB ID', 'HET ID', 'UniProt ID', 'GenBank ID', 'DPD',
# 'RxList Link', 'Pdrhealth Link', 'Wikipedia ID', 'Drugs.com Link',
# 'NDC ID', 'ATC Codes', 'AHFS Codes', 'PDB Entries', 'FDA label Link',
# 'MSDS Link', 'ChEMBL ID', 'Therapeutic Targets Database', 'UNII'
DRUGBANK_LINKS_COLUMNS = {
    "CAS Number":               "CAS",
    "KEGG Drug ID":             "KEGG.DRUG",
    "KEGG Compound ID":         "KEGG.COMPOUND",
    "PubChem Compound ID":      "PUBCHEM.COMPOUND",
    "PubChem Substance ID":     "PUBCHEM.SUBSTANCE",
    "ChEBI ID":                 "CHEBI",
    "ChEMBL ID":                "CHEMBL.COMPOUND",
    "ATC Codes":                "ATC",       # may contain pipe-delimited list
    "UNII":                     "UNII",
    "PharmGKB ID":              "PHARMGKB",
}


def _emit_links(writer, src_csv: Path) -> int:
    """Read DrugBank's external_links.csv and emit one row per (xref) cell."""
    n = 0
    with src_csv.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            db_id = row.get("DrugBank ID", "").strip()
            if not db_id:
                continue
            for col, namespace in DRUGBANK_LINKS_COLUMNS.items():
                raw = row.get(col, "").strip()
                if not raw:
                    continue
                # ATC and a few others may have multiple values pipe-delimited
                for ext in raw.split("|"):
                    ext = ext.strip()
                    if not ext:
                        continue
                    writer.writerow([db_id, namespace, ext])
                    n += 1
    return n


def _emit_pubchem(writer, src_csv: Path) -> int:
    """Read pubchem_to_drugbank.csv (drugbank_id, pubchem_cid) and emit rows."""
    n = 0
    with src_csv.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            db_id = row.get("drugbank_id", "").strip()
            cid   = row.get("pubchem_cid", "").strip()
            if not (db_id and cid):
                continue
            writer.writerow([db_id, "PUBCHEM.COMPOUND", cid])
            n += 1
    return n


def _emit_drugcentral(writer, src_csv: Path) -> int:
    """Read DrugCentral xref dump (drugbank_id, namespace, external_id)."""
    n = 0
    with src_csv.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            db_id = row.get("drugbank_id", "").strip()
            ns    = row.get("namespace", "").strip().upper()
            ext   = row.get("external_id", "").strip()
            if not (db_id and ns and ext):
                continue
            writer.writerow([db_id, ns, ext])
            n += 1
    return n


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--drugbank-links",
                    default=str(GOLD / "drugbank_external_links.csv"),
                    help="Path to DrugBank's external_links.csv "
                         "(download from go.drugbank.com/releases/latest#external-links)")
    ap.add_argument("--pubchem", default=str(REPO_ROOT / "data/openbilink/pubchem_to_drugbank.csv"))
    ap.add_argument("--drugcentral", default=str(GOLD / "drugcentral_xref.csv"),
                    help="Optional DrugCentral xref file (skipped if absent)")
    ap.add_argument("-o", "--output", default=str(OUT))
    args = ap.parse_args()

    GOLD.mkdir(parents=True, exist_ok=True)

    n_total = 0
    seen: set[tuple[str, str, str]] = set()

    with open(args.output, "w", newline="", encoding="utf-8") as fo:
        writer = csv.writer(fo)
        writer.writerow(["drugbank_id", "namespace", "external_id"])

        # Wrap the writer so duplicate rows are dropped when the same
        # (drugbank_id, namespace, external_id) appears in multiple sources.
        class DedupWriter:
            def __init__(self, real):
                self.real = real
            def writerow(self, row):
                key = (row[0], row[1], row[2])
                if key in seen:
                    return
                seen.add(key)
                self.real.writerow(row)

        dw = DedupWriter(writer)

        # Source 1: DrugBank external_links.csv (the fullest)
        if Path(args.drugbank_links).exists():
            n = _emit_links(dw, Path(args.drugbank_links))
            print(f"  [drugbank_external_links]  {n:,} rows")
            n_total += n
        else:
            print(f"  [drugbank_external_links]  not found at {args.drugbank_links}")
            print("     Download from https://go.drugbank.com/releases/latest#external-links "
                  "(free academic account) and place at that path to recover ~80% more drugs.")

        # Source 2: pubchem_to_drugbank.csv (always present)
        if Path(args.pubchem).exists():
            n = _emit_pubchem(dw, Path(args.pubchem))
            print(f"  [pubchem_to_drugbank]       {n:,} rows")
            n_total += n

        # Source 3: optional DrugCentral xref
        if Path(args.drugcentral).exists():
            n = _emit_drugcentral(dw, Path(args.drugcentral))
            print(f"  [drugcentral_xref]          {n:,} rows")
            n_total += n

    n_unique = len(seen)
    n_drugs  = len({k[0] for k in seen})
    print(f"\n  wrote {args.output}")
    print(f"    total rows seen:   {n_total:,}")
    print(f"    unique rows:       {n_unique:,}")
    print(f"    distinct DrugBank IDs covered: {n_drugs:,}")


if __name__ == "__main__":
    main()
