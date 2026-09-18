#!/usr/bin/env python3
"""Entrez gene id to HGNC symbol, for rendering numeric gene nodes in dossiers.

DRKG keys genes as ``Gene::3716`` and OpenBioLink as ``NCBIGENE:3716``; shown
to a language model as a number they carry no information, so the grounding
dossiers render them as symbols. The map is derived from the vendored HGNC
complete set (data/crosswalks/hgnc_complete_set.txt, release recorded in
data/PROVENANCE.md), approved symbols only, one symbol per Entrez id.

Writes data/gold_standards/bridges/entrez_to_symbol.tsv (two columns, no
header) and ENTREZ_TO_SYMBOL_PROVENANCE.txt. Run from the repo root:
    python tools/build_entrez_symbol_map.py
"""
from __future__ import annotations

import hashlib
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data/crosswalks/hgnc_complete_set.txt"
OUT = ROOT / "data/gold_standards/bridges/entrez_to_symbol.tsv"
PROV = ROOT / "data/gold_standards/bridges/ENTREZ_TO_SYMBOL_PROVENANCE.txt"


def main() -> None:
    hgnc = pd.read_csv(SRC, sep="\t", usecols=["symbol", "entrez_id", "status"], dtype=str, low_memory=False)
    approved = hgnc[(hgnc.status == "Approved") & hgnc.entrez_id.notna()].copy()
    approved["entrez_id"] = approved.entrez_id.str.strip()
    dup = approved.entrez_id.duplicated(keep=False).sum()
    table = approved.drop_duplicates("entrez_id").sort_values("entrez_id", key=lambda s: s.astype(int))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    table[["entrez_id", "symbol"]].to_csv(OUT, sep="\t", header=False, index=False)
    PROV.write_text("\n".join([
        f"# Provenance for {OUT.name}",
        f"# built:      {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
        "# built_by:   tools/build_entrez_symbol_map.py",
        f"# source:     data/crosswalks/hgnc_complete_set.txt sha256 {hashlib.sha256(SRC.read_bytes()).hexdigest()[:16]} (release in data/PROVENANCE.md)",
        "# rule:       status == Approved and entrez_id present; one symbol per Entrez id",
        f"# rows in:    {len(hgnc)}  approved with Entrez id: {len(approved)}  duplicate Entrez ids dropped: {int(dup)}",
        f"# rows out:   {len(table)}",
        f"# sha256 {OUT.name}: {hashlib.sha256(OUT.read_bytes()).hexdigest()[:16]}",
    ]) + "\n")
    print(PROV.read_text())


if __name__ == "__main__":
    main()
