#!/usr/bin/env python3
"""One-off rename (2026-09-18): "prospective gold standard" -> "post-cutoff approval set".

Renames the set's folder, files, builder, figure scripts, loader module, tests, figure files and
document, and rewrites every reference in src/, config/, tools/, tests/, the notebooks and the docs.
The evaluation axis in notebook 07 keeps its name (prospective evaluation, prospective_* parameters,
stats/prospective.py, figures 07d/07e): that is the test, not the data. Frozen material under
data/gold_standards/legacy/ and the vendored package under source/ are not touched.
Run from the repo root:  python tools/rename_post_cutoff.py [--dry-run]
"""
import argparse
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

MOVES = [
    ("data/gold_standards/prospective/prospective_gold_standard.csv", "data/gold_standards/prospective/post_cutoff_approvals.csv"),
    ("data/gold_standards/prospective/prospective_pair_status.csv", "data/gold_standards/prospective/post_cutoff_pair_status.csv"),
    ("data/gold_standards/prospective/prospective_funnel.csv", "data/gold_standards/prospective/post_cutoff_funnel.csv"),
    ("data/gold_standards/prospective/prospective_dropped.csv", "data/gold_standards/prospective/post_cutoff_dropped.csv"),
    ("data/gold_standards/prospective", "data/gold_standards/post_cutoff_approvals"),
    ("data/PROSPECTIVE_GOLD_STANDARD.md", "data/POST_CUTOFF_APPROVALS.md"),
    ("tools/build_prospective_gold.py", "tools/build_post_cutoff_approvals.py"),
    ("tools/make_prospective_flow_figure.py", "tools/make_post_cutoff_flow_figure.py"),
    ("tools/make_prospective_flow_pptx.py", "tools/make_post_cutoff_flow_pptx.py"),
    ("tools/prospective_flow_pptx.js", "tools/post_cutoff_flow_pptx.js"),
    ("src/biokgsuite/gold/prospective_gold.py", "src/biokgsuite/gold/post_cutoff_approvals.py"),
    ("tests/test_prospective_gold.py", "tests/test_post_cutoff_approvals.py"),
    ("tests/test_build_prospective_gold.py", "tests/test_build_post_cutoff_approvals.py"),
    ("results/figures/07_generalization/07h_prospective_funnel.pdf", "results/figures/07_generalization/07h_post_cutoff_funnel.pdf"),
    ("results/figures/07_generalization/07h_prospective_funnel.png", "results/figures/07_generalization/07h_post_cutoff_funnel.png"),
    ("results/figures/07_generalization/07i_prospective_flow.pdf", "results/figures/07_generalization/07i_post_cutoff_flow.pdf"),
    ("results/figures/07_generalization/07i_prospective_flow.png", "results/figures/07_generalization/07i_post_cutoff_flow.png"),
    ("results/figures/07_generalization/07i_prospective_flow.pptx", "results/figures/07_generalization/07i_post_cutoff_flow.pptx"),
    ("results/figures/07_generalization/07i_prospective_flow_pptx_preview.jpg", "results/figures/07_generalization/07i_post_cutoff_flow_pptx_preview.jpg"),
]

# ordered: longer / more specific first
PLAIN = [
    ("gold_standards/prospective/", "gold_standards/post_cutoff_approvals/"),
    ("prospective_gold_standard.csv", "post_cutoff_approvals.csv"),
    ("prospective_pair_status.csv", "post_cutoff_pair_status.csv"),
    ("prospective_funnel.csv", "post_cutoff_funnel.csv"),
    ("prospective_dropped.csv", "post_cutoff_dropped.csv"),
    ("PROSPECTIVE_GOLD_STANDARD_PROVENANCE.txt", "PROVENANCE.txt"),
    ("PROSPECTIVE_GOLD_STANDARD.md", "POST_CUTOFF_APPROVALS.md"),
    ("test_build_prospective_gold", "test_build_post_cutoff_approvals"),
    ("test_prospective_gold", "test_post_cutoff_approvals"),
    ("build_prospective_gold", "build_post_cutoff_approvals"),
    ("make_prospective_flow_figure", "make_post_cutoff_flow_figure"),
    ("make_prospective_flow_pptx", "make_post_cutoff_flow_pptx"),
    ("prospective_flow_pptx", "post_cutoff_flow_pptx"),
    ("load_prospective_gold_standard", "load_post_cutoff_approvals"),
    ("prospective_gold", "post_cutoff_approvals"),
    ("07h_prospective_funnel", "07h_post_cutoff_funnel"),
    ("07i_prospective_flow", "07i_post_cutoff_flow"),
    ("PROSPECTIVE_PATH", "POST_CUTOFF_PATH"),
    ("prospective gold-standard", "post-cutoff approval-set"),
    ("Prospective gold-standard", "Post-cutoff approval-set"),
    ("prospective gold standard", "post-cutoff approval set"),
    ("Prospective gold standard", "Post-cutoff approval set"),
    ("prospective set", "post-cutoff set"),
    ("Prospective set", "Post-cutoff set"),
    ("regulator-derived gold standard", "regulator-derived post-cutoff approval set"),
    ("prospective (post-snapshot) gold standard", "post-cutoff approval set (approvals dated after every graph's content date)"),
    ("## prospective/", "## post_cutoff_approvals/"),
]
REGEX = [(re.compile(r"gold_standards/prospective(?![A-Za-z0-9_/])"), "gold_standards/post_cutoff_approvals"),
         (re.compile(r"(?<![A-Za-z0-9_/])prospective/(?![A-Za-z0-9_])"), "post_cutoff_approvals/")]

TEXT_DIRS = ["src", "config", "tools", "tests", "hpc", "supplementary", "biokgsuite_dimensions", "paper", "reference", "data/gold_standards/post_cutoff_approvals"]
TEXT_EXT = {".py", ".yaml", ".yml", ".sh", ".md", ".ipynb", ".toml", ".txt", ".cfg", ".ini", ".json", ".js"}
ROOT_DOCS = ["README.md", "NEXT_STEPS.md", "METHODS_BENCHMARK_CONSTRUCTION.md", "NOTEBOOK_SWEEP_2026-09-17.md", "EMBEDDING_CROSSCHECK.md",
             "data/REFERENCE_SETS.md", "data/MANIFEST.md", "data/PROVENANCE.md", "data/POST_CUTOFF_APPROVALS.md", "data/CROSSWALKS.md",
             "data/SOURCE_COUNTS.md", "data/gold_standards/README.md"]
SKIP_PARTS = {".venv", "node_modules", "__pycache__", "checkpoints", "_to_delete", "legacy", "source"}


def rewrite(t: str) -> str:
    for old, new in PLAIN:
        t = t.replace(old, new)
    for rx, new in REGEX:
        t = rx.sub(new, t)
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    dry = ap.parse_args().dry_run
    for old, new in MOVES:
        src, dst = ROOT / old, ROOT / new
        if not src.exists():
            print(f"skip (absent): {old}")
            continue
        if dst.exists():
            sys.exit(f"refusing to overwrite {new}")
        print(f"move  {old} -> {new}")
        if not dry:
            os.rename(src, dst)
    files = []
    for d in TEXT_DIRS:
        base = ROOT / d
        if not base.exists():
            continue
        for p in base.rglob("*"):
            if p.is_file() and p.suffix in TEXT_EXT and not (SKIP_PARTS & set(p.relative_to(ROOT).parts)) and p.name not in (Path(__file__).name, "reorganise_gold_standards.py"):
                files.append(p)
    files += [ROOT / d for d in ROOT_DOCS if (ROOT / d).exists()]
    n = 0
    for p in files:
        try:
            t = p.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        t2 = rewrite(t)
        if t2 != t:
            n += 1
            print(f"edit  {p.relative_to(ROOT)}")
            if not dry:
                p.write_text(t2, encoding="utf-8")
    # legacy/README.md: only the pointers to the current set, not the historical file names
    lp = ROOT / "data/gold_standards/legacy/README.md"
    if lp.exists():
        t = lp.read_text(encoding="utf-8")
        t2 = t.replace("`prospective/`", "`post_cutoff_approvals/`").replace("`tools/build_prospective_gold.py`", "`tools/build_post_cutoff_approvals.py`")
        t2 = t2.replace("115-pair prospective set", "115-pair post-cutoff set").replace("the prospective set is absent", "the post-cutoff set is absent")
        if t2 != t:
            print("edit  data/gold_standards/legacy/README.md (pointers only)")
            if not dry:
                lp.write_text(t2, encoding="utf-8")
    print(f"{n} files rewritten")


if __name__ == "__main__":
    main()
