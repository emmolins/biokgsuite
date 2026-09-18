#!/usr/bin/env python3
"""One-off migration (2026-09-18): reorganise data/gold_standards/ into

    prospective/   the 128-pair prospective gold standard and its provenance (built by tools/build_prospective_gold.py)
    source/        the vendored repurposing-dataset table and zip the builder reads, plus the second-model review proposals
    bridges/       ontology and identifier bridges notebooks 06/07 and the builder read (DOID/MONDO/MeSH, DrugBank names, UniProt)
    ctd/           CTD files for notebook 06
    legacy/        every superseded gold-standard attempt (v2 to v5, the FDA census, corroboration outputs, older vendored copies)
    derived/       tables a notebook builds and re-reads (multihop_benchmark_v1.csv)

and rewrite every path in src/, config/, tools/, tests/, the notebooks and the docs. Moves use os.rename and refuse to
overwrite; nothing is deleted. Run from the repo root:  python tools/reorganise_gold_standards.py [--dry-run]
"""
import argparse
import hashlib
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GS = ROOT / "data/gold_standards"

PROSPECTIVE = ["prospective_gold_standard.csv", "prospective_pair_status.csv", "prospective_funnel.csv", "prospective_dropped.csv"]
BRIDGES = ["do_diseases.csv", "mesh_to_doid.csv", "mondo_hierarchy.csv", "MONDO_HIERARCHY_PROVENANCE.txt", "mondo.sssom.tsv",
           "uniprot_genesproteins.csv", "drugbank_drugs.csv", "drug_name_to_drugbank.csv", "DRUG_NAME_XREF_PROVENANCE.txt"]
CTD = ["CTD_chem_gene_ixns.tsv.gz", "CTD_curated_genes_diseases.tsv.gz", "CTD_PROVENANCE.txt"]
LEGACY = ["gold_standard_v2.tsv", "gold_standard_v3_common3.tsv", "gold_standard_v4_bigset.tsv", "gold_standard_v5_dedup.tsv",
          "gold_standard_v5_indication.tsv", "gold_standard_v5_strict.tsv", "gold_standard_v5_strict_evidence.csv",
          "gold_standard_v5_verified.tsv", "pairs_with_ids_2023_2026.csv", "pairs_with_ids_v2.csv", "usable_pairs_v2_audited.csv",
          "expanded_prospective_gold_standard.tsv", "coverage_annotation_v4.csv", "cleaning_manifest.csv",
          "fda_drugs_id_mapping.csv", "fda_indication_expansion_census.csv", "fda_indication_expansion_census.csv.meta.json",
          "fda_indication_expansion_census_drugbank.csv", "fda_indication_expansion_census_drugbank.csv.mapping.json",
          "fda_indication_expansion_census_frame.csv", "fda_repurposing_events_2023_2026.csv", "fda_census_validation_2026-09-17.log",
          "pmda_approvals.csv", "EMA_PMDA_VERIFICATION.md", "FDA_CENSUS_METHOD.md", "FDA_REPURPOSING_VETTING.md",
          "PROVENANCE_253_candidate_pairs.md", "RECHECK_CHECKLIST.md"]
SOURCE_VERSION = "repurposing-dataset-v0.6.0"
SOURCE_ZIP = "repurposing-dataset-0.6.0+20260918.zip"

# (old relative path under gold_standards, new relative path)
MOVES = [(f, f"prospective/{f}") for f in PROSPECTIVE]
MOVES += [("PROSPECTIVE_GOLD_STANDARD_PROVENANCE.txt", "prospective/PROVENANCE.txt")]
MOVES += [(f, f"bridges/{f}") for f in BRIDGES]
MOVES += [(f, f"ctd/{f}") for f in CTD]
MOVES += [(f, f"legacy/{f}") for f in LEGACY]
MOVES += [(f"approvals_raw/{SOURCE_VERSION}", f"source/{SOURCE_VERSION}"), (f"approvals_raw/{SOURCE_ZIP}", f"source/{SOURCE_ZIP}"),
          ("pending_lookups_73_decisions.csv", "source/review_proposals_second_model_2026-09-18.csv"),
          ("approvals_raw/fda_efficacy_supplements_2023_2026.psv", "legacy/approvals_raw/fda_efficacy_supplements_2023_2026.psv"),
          ("approvals_raw/repurposing-dataset-v0.2.0", "legacy/approvals_raw/repurposing-dataset-v0.2.0"),
          ("approvals_raw/repurposing-dataset-0.2.0+20260917.zip", "legacy/approvals_raw/repurposing-dataset-0.2.0+20260917.zip"),
          ("approvals_raw/repurposing-dataset-v0.5.0", "legacy/approvals_raw/repurposing-dataset-v0.5.0"),
          ("approvals_raw/repurposing-dataset-0.5.0+20260918.zip", "legacy/approvals_raw/repurposing-dataset-0.5.0+20260918.zip")]
# an unpacked duplicate of the v0.6.0 zip; parked for the user to delete
DUPLICATE = ("approvals_raw/repurposing-dataset", ROOT / "_to_delete/repurposing-dataset_unpacked_duplicate_of_0.6.0_zip")

# textual path rewrites, longest first so prefixes do not shadow
REWRITES = [
    (f"gold_standards/approvals_raw/{SOURCE_VERSION}", f"gold_standards/source/{SOURCE_VERSION}"),
    (f"gold_standards/approvals_raw/{SOURCE_ZIP}", f"gold_standards/source/{SOURCE_ZIP}"),
    ("gold_standards/PROSPECTIVE_GOLD_STANDARD_PROVENANCE.txt", "gold_standards/prospective/PROVENANCE.txt"),
    ("gold_standards/pending_lookups_73_decisions.csv", "gold_standards/source/review_proposals_second_model_2026-09-18.csv"),
    ("gold_standards/approvals_raw/", "gold_standards/legacy/approvals_raw/"),
    ("gold_standards/multihop_benchmark_v1.csv", "gold_standards/derived/multihop_benchmark_v1.csv"),
]
REWRITES += [(f"gold_standards/{f}", f"gold_standards/prospective/{f}") for f in PROSPECTIVE]
REWRITES += [(f"gold_standards/{f}", f"gold_standards/bridges/{f}") for f in BRIDGES]
REWRITES += [(f"gold_standards/{f}", f"gold_standards/ctd/{f}") for f in CTD]
REWRITES += [(f"gold_standards/{f}", f"gold_standards/legacy/{f}") for f in LEGACY]
REWRITES.sort(key=lambda p: -len(p[0]))

# line-level edits where a directory is joined at run time
SPECIAL = {
    "src/biokgsuite/gold/ontology_classifier.py": [("gs = base_dir / 'data' / 'gold_standards'", "gs = base_dir / 'data' / 'gold_standards' / 'bridges'")],
    "tools/build_gold_standard_aux.py": [("GS = ROOT / 'data' / 'gold_standards'", "GS = ROOT / 'data' / 'gold_standards' / 'bridges'")],
    "tools/build_mondo_hierarchy.py": [("OUT = ROOT / 'data' / 'gold_standards' / 'mondo_hierarchy.csv'", "OUT = ROOT / 'data' / 'gold_standards' / 'bridges' / 'mondo_hierarchy.csv'"),
                                       ("PROV = ROOT / 'data' / 'gold_standards' / 'MONDO_HIERARCHY_PROVENANCE.txt'", "PROV = ROOT / 'data' / 'gold_standards' / 'bridges' / 'MONDO_HIERARCHY_PROVENANCE.txt'")],
    "tools/build_drug_name_xref.py": [("OUT = ROOT / 'data' / 'gold_standards' / 'drug_name_to_drugbank.csv'", "OUT = ROOT / 'data' / 'gold_standards' / 'bridges' / 'drug_name_to_drugbank.csv'"),
                                      ("PROV = ROOT / 'data' / 'gold_standards' / 'DRUG_NAME_XREF_PROVENANCE.txt'", "PROV = ROOT / 'data' / 'gold_standards' / 'bridges' / 'DRUG_NAME_XREF_PROVENANCE.txt'"),
                                      ("LEGACY = ROOT / 'data' / 'gold_standards' / 'drugbank_drugs.csv'", "LEGACY = ROOT / 'data' / 'gold_standards' / 'bridges' / 'drugbank_drugs.csv'")],
    "tools/build_drugbank_xref.py": [('GOLD = REPO_ROOT / "data" / "gold_standards"', 'GOLD = REPO_ROOT / "data" / "gold_standards" / "bridges"')],
    "tools/fetch_ctd.sh": [("GS=data/gold_standards\n", "GS=data/gold_standards/ctd\n")],
    "tools/download_mondo_sssom.sh": [('DEST_DIR="data/gold_standards"', 'DEST_DIR="data/gold_standards/bridges"')],
    "tools/measure_shared_retrieval_set.py": [("os.path.join(GS, 'mesh_to_doid.csv')", "os.path.join(GS, 'bridges', 'mesh_to_doid.csv')"),
                                              ("os.path.join(GS, 'uniprot_genesproteins.csv')", "os.path.join(GS, 'bridges', 'uniprot_genesproteins.csv')"),
                                              ("os.path.join(GS, 'CTD_curated_genes_diseases.tsv.gz')", "os.path.join(GS, 'ctd', 'CTD_curated_genes_diseases.tsv.gz')")],
    "tools/fetch_approvals.py": [("OUT = ROOT / 'data' / 'gold_standards' / 'approvals_raw'", "OUT = ROOT / 'data' / 'gold_standards' / 'legacy' / 'approvals_raw'"),
                                 ("PROV = ROOT / 'data' / 'gold_standards' / 'APPROVALS_PROVENANCE.txt'", "PROV = ROOT / 'data' / 'gold_standards' / 'legacy' / 'approvals_raw' / 'APPROVALS_PROVENANCE.txt'")],
    "tools/build_approval_gold_standard.py": [("GS = ROOT / 'data' / 'gold_standards'\n", "GS = ROOT / 'data' / 'gold_standards' / 'legacy'   # superseded pipeline: reads and writes under legacy/\nBRIDGES = ROOT / 'data' / 'gold_standards' / 'bridges'\n"),
                                              ("GS / 'drug_name_to_drugbank.csv'", "BRIDGES / 'drug_name_to_drugbank.csv'"),
                                              ("GS / 'drugbank_drugs.csv'", "BRIDGES / 'drugbank_drugs.csv'"),
                                              ("GS / 'do_diseases.csv'", "BRIDGES / 'do_diseases.csv'")],
    "tools/clean_116.py": [('ap.add_argument("--outdir", default="data/gold_standards")', 'ap.add_argument("--outdir", default="data/gold_standards/legacy")')],
    "tools/build_prospective_gold.py": [('GS_DIR = ROOT / "data/gold_standards"', 'GS_DIR = ROOT / "data/gold_standards/prospective"   # outputs\nBRIDGE_DIR = ROOT / "data/gold_standards/bridges"'),
                                        ("br = Bridges(GS_DIR)", "br = Bridges(BRIDGE_DIR)")],
    "biokgsuite_dimensions/06_task_performance/notebook.ipynb": [("'CTD_chem_gene_ixns.tsv.gz':", "'ctd/CTD_chem_gene_ixns.tsv.gz':"),
                                                                 ("'CTD_curated_genes_diseases.tsv.gz':", "'ctd/CTD_curated_genes_diseases.tsv.gz':"),
                                                                 ("'mesh_to_doid.csv':", "'bridges/mesh_to_doid.csv':"),
                                                                 ("'drugbank_drugs.csv':", "'bridges/drugbank_drugs.csv':"),
                                                                 ("_GS / 'CTD_chem_gene_ixns.tsv.gz'", "_GS / 'ctd' / 'CTD_chem_gene_ixns.tsv.gz'"),
                                                                 ("_GS / 'CTD_curated_genes_diseases.tsv.gz'", "_GS / 'ctd' / 'CTD_curated_genes_diseases.tsv.gz'"),
                                                                 ("MULTIHOP_BENCH = BASE / 'data' / 'gold_standards' / 'multihop_benchmark_v1.csv'", "MULTIHOP_BENCH = BASE / 'data' / 'gold_standards' / 'derived' / 'multihop_benchmark_v1.csv'")],
}
TEXT_DIRS = ["src", "config", "tools", "tests", "hpc", "supplementary", "biokgsuite_dimensions", "paper", "reference"]
TEXT_EXT = {".py", ".yaml", ".yml", ".sh", ".md", ".ipynb", ".toml", ".txt", ".cfg", ".ini", ".json"}
ROOT_DOCS = ["README.md", "NEXT_STEPS.md", "METHODS_BENCHMARK_CONSTRUCTION.md", "NOTEBOOK_SWEEP_2026-09-17.md", "EMBEDDING_CROSSCHECK.md",
             "data/REFERENCE_SETS.md", "data/MANIFEST.md", "data/PROVENANCE.md", "data/PROSPECTIVE_GOLD_STANDARD.md", "data/CROSSWALKS.md", "data/SOURCE_COUNTS.md"]


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()[:16]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    dry = a.dry_run
    # 0. fingerprints of the outputs we must reproduce after the move
    before = {f: sha(GS / f) for f in PROSPECTIVE if (GS / f).exists()}
    # 1. moves
    for old, new in MOVES:
        src, dst = GS / old, GS / new
        if not src.exists():
            print(f"skip (absent): {old}")
            continue
        if dst.exists():
            sys.exit(f"refusing to overwrite {dst}")
        print(f"move  {old}  ->  {new}")
        if not dry:
            dst.parent.mkdir(parents=True, exist_ok=True)
            os.rename(src, dst)
    src = GS / DUPLICATE[0]
    if src.exists():
        print(f"park  {DUPLICATE[0]}  ->  {DUPLICATE[1].relative_to(ROOT)}")
        if not dry:
            DUPLICATE[1].parent.mkdir(exist_ok=True)
            if DUPLICATE[1].exists():
                sys.exit(f"refusing to overwrite {DUPLICATE[1]}")
            os.rename(src, DUPLICATE[1])
    if not dry:
        (GS / "derived").mkdir(exist_ok=True)
        left = [p for p in (GS / "approvals_raw").iterdir()] if (GS / "approvals_raw").exists() else []
        if (GS / "approvals_raw").exists() and not left:
            os.rmdir(GS / "approvals_raw")
        elif left:
            print("approvals_raw still holds:", [p.name for p in left])
    # 2. special edits first (they contain the generic strings too)
    for rel, edits in SPECIAL.items():
        p = ROOT / rel
        if not p.exists():
            print(f"special: {rel} absent")
            continue
        t = p.read_text(encoding="utf-8")
        for old, new in edits:
            if old in t:
                t = t.replace(old, new)
                print(f"edit  {rel}: {old[:60]}")
            else:
                print(f"WARN  {rel}: pattern not found: {old[:70]}")
        if not dry:
            p.write_text(t, encoding="utf-8")
    # 3. generic rewrites
    files = []
    for d in TEXT_DIRS:
        for p in (ROOT / d).rglob("*"):
            if (p.is_file() and p.suffix in TEXT_EXT and ".venv" not in p.parts and "checkpoints" not in p.parts
                    and p.name != Path(__file__).name):
                files.append(p)
    files += [ROOT / d for d in ROOT_DOCS if (ROOT / d).exists()]
    n_files = 0
    for p in files:
        try:
            t = p.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        t2 = t
        for old, new in REWRITES:
            t2 = t2.replace(old, new)
        if t2 != t:
            n_files += 1
            print(f"paths {p.relative_to(ROOT)}")
            if not dry:
                p.write_text(t2, encoding="utf-8")
    print(f"{n_files} files rewritten")
    print("fingerprints before move:", before)


if __name__ == "__main__":
    main()
