#!/usr/bin/env python3
"""Editable PowerPoint version of the post-cutoff approval-set flow chart (Nature artwork style).

Re-derives every count from the shipped tables (same code as make_post_cutoff_flow_figure.py), writes them to a
JSON file and calls tools/post_cutoff_flow_pptx.js (pptxgenjs) to draw native, editable shapes and tables.

    python tools/make_post_cutoff_flow_pptx.py   # -> results/figures/07_generalization/07i_post_cutoff_flow.pptx
"""
import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from build_post_cutoff_approvals import GS_DIR, SOURCE_DEFAULT  # noqa: E402
from make_post_cutoff_flow_figure import counts  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", type=Path, default=SOURCE_DEFAULT)
    ap.add_argument("--out", type=Path, default=ROOT / "results/figures/07_generalization/07i_post_cutoff_flow.pptx")
    a = ap.parse_args()
    src = pd.read_csv(a.source, dtype=str)
    dropped = pd.read_csv(GS_DIR / "post_cutoff_dropped.csv")
    gold = pd.read_csv(GS_DIR / "post_cutoff_approvals.csv")
    funnel = pd.read_csv(GS_DIR / "post_cutoff_funnel.csv")
    c = counts(src, dropped, gold)
    c["any"], c["rows_by_reg_s3"] = dict(c["any"]), dict(c["rows_by_reg_s3"])
    c["kg"] = {k: {r.stage: int(r.n) for r in funnel[funnel.scope == k].itertuples()} for k in funnel.scope.unique() if k not in ("global", "common_support")}
    c["common"] = {r.stage: int(r.n) for r in funnel[funnel.scope == "common_support"].itertuples()}
    c["version"] = (a.source.parent / "VERSION").read_text().strip() if (a.source.parent / "VERSION").exists() else "?"
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as f:
        json.dump(c, f)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["node", str(ROOT / "tools/post_cutoff_flow_pptx.js"), f.name, str(a.out)], check=True)


if __name__ == "__main__":
    main()
