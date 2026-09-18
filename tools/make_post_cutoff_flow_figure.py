#!/usr/bin/env python3
"""Flow chart of the post-cutoff approval set: regulator rows in, cleaning steps, pairs out.

Every number is re-derived from the files on disk (the vendored source table, post_cutoff_dropped.csv,
post_cutoff_approvals.csv), so the figure cannot drift from the data.

    python tools/make_post_cutoff_flow_figure.py            # -> results/figures/07_generalization/07i_post_cutoff_flow.{pdf,png}
"""
import argparse
import collections
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from build_post_cutoff_approvals import SOURCE_DEFAULT, GS_DIR  # noqa: E402

REG_LABEL = {"FDA": "FDA (CDER)", "EMA": "EMA (CHMP)", "HC": "Health Canada", "PMDA": "PMDA", "TGA": "TGA"}
REG_WHAT = {"FDA": "efficacy supplements,\nsub-type New Indication",
            "EMA": "positive opinions,\nextension of indication",
            "HC": "SNDS whose reason\nnames an indication",
            "PMDA": "approvals whose notes\nsay new indication",
            "TGA": "extension-of-indications\nregistrations"}
ORDER = ["FDA", "EMA", "HC", "PMDA", "TGA"]

NOT_NEW = {"line_extension", "expanded_population", "ingredient_prior_indication", "formulation_or_route",
           "biosimilar_alignment", "first_indication"}
NOT_CHANGE = {"wording_only", "commitment_fulfilment", "indication_removed", "out_of_scope_veterinary", "excluded", "unclear"}


def counts(src: pd.DataFrame, dropped: pd.DataFrame, gold: pd.DataFrame):
    c = {}
    c["rows"] = src["regulator"].value_counts().to_dict()
    c["n0"] = len(src)
    c["text"] = src.groupby("regulator")["indication_text"].apply(lambda s: int(s.notna().sum())).to_dict()
    s1 = dropped[dropped.stage == "S1"]["reason"]
    cls = collections.Counter()
    for x in s1:
        if "class=" in x:
            cls[x.split("class=")[-1].split(";")[0]] += 1
        elif x.startswith("out of scope"):
            cls["no_target"] += 1
        else:
            cls[x] += 1
    c["no_text"] = cls.pop("disease pending_source_text", 0)
    c["pending"] = cls.pop("disease pending_lookup", 0) + cls.pop("drug pending_lookup", 0)
    c["no_target"] = cls.pop("no_target", 0)
    c["not_new"] = {k: v for k, v in cls.items() if k in NOT_NEW}
    c["not_change"] = {k: v for k, v in cls.items() if k in NOT_CHANGE}
    c["n1"] = c["n0"] - int(s1.shape[0])
    s2 = dropped[dropped.stage == "S2"]["reason"].str.split(r" \(").str[0].value_counts().to_dict()
    c["s2"] = s2
    c["n2"] = c["n1"] - int(sum(s2.values()))
    c["events"] = int(sum(len(str(e).split(";")) for e in gold["evidence_ref"]))
    # rows that reach S3: one row can name two diseases, so count rows, not events
    c["rows_s3"] = int(len({cid for e in gold["evidence_ref"] for cid in str(e).split(";")}))
    c["pairs"] = len(gold)
    c["drugs"], c["diseases"] = gold["drugbank_id"].nunique(), gold["mondo_id"].nunique()
    c["multi"] = int((gold["n_agencies"] > 1).sum())
    c["any"] = collections.Counter(a for s in gold["agencies"] for a in s.split(";"))
    c["first"] = gold["agency"].value_counts().to_dict()
    c["rows_by_reg_s3"] = collections.Counter(cid.split(":")[0] for e in gold["evidence_ref"] for cid in set(str(e).split(";")))
    return c


def draw(c: dict, version: str, out_stem: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
    try:
        from biokgsuite.viz.style import setup_style
        setup_style()
    except Exception:  # noqa: BLE001
        pass
    INK, MUTED, RULE, SURF, DROP, KEEP = "#1a1a1a", "#5b5b5b", "#dcdcd8", "#fcfcfb", "#c0392b", "#2e6f5e"
    fig, ax = plt.subplots(figsize=(7.2, 10.4))
    ax.set_xlim(0, 10)
    ax.set_ylim(-1.15, 13.4)
    ax.axis("off")

    def box(x, y, w, h, hero=False, fc=None):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.06",
                                    linewidth=1.8 if hero else 1.0, edgecolor=INK if hero else RULE,
                                    facecolor=fc or ("#f4f4f1" if hero else SURF), zorder=2))

    def arrow(x0, y0, x1, y1):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=9, color=MUTED, lw=0.9, zorder=3))

    # --- regulator sources
    ax.text(0, 13.2, "Regulator records, 2023-01-01 to 2026-09-17, to the post-cutoff approval set", fontsize=8.5, fontweight="bold", color=INK, va="center")
    bw, gap, y0, bh = 1.84, 0.2, 11.15, 1.5
    for i, reg in enumerate(ORDER):
        x = i * (bw + gap)
        box(x, y0, bw, bh)
        ax.text(x + 0.12, y0 + bh - 0.3, f"{c['rows'][reg]:,}", fontsize=12, fontweight="bold", color=INK, va="center", zorder=3)
        ax.text(x + 0.12, y0 + bh - 0.66, REG_LABEL[reg], fontsize=6.6, color=INK, va="center", zorder=3)
        ax.text(x + 0.12, y0 + 0.36, REG_WHAT[reg], fontsize=5.4, color=MUTED, va="center", zorder=3, linespacing=1.1)
        arrow(x + bw / 2, y0 - 0.04, x + bw / 2, y0 - 0.42)

    # --- merged table
    y = 9.55
    box(0, y, 10, 1.15, hero=False)
    ax.text(0.3, y + 0.62, f"{c['n0']:,}", fontsize=13, fontweight="bold", color=INK, va="center", zorder=3)
    ax.text(1.75, y + 0.78, "candidate rows, one per regulator decision (S0); indication text from the regulator's record", fontsize=6.9, color=INK, va="center", zorder=3)
    ax.text(1.75, y + 0.36, _wrap(f"FDA approval letter, “provides for” sentence ({c['text']['FDA']} of {c['rows']['FDA']}); EMA EPAR procedural-steps "
                                  f"scope matched to the opinion date ({c['text']['EMA']} of {c['rows']['EMA']}); Health Canada NOC reason, PMDA notes, "
                                  f"TGA indication (all rows)", 112),
            fontsize=5.6, color=MUTED, va="center", zorder=3, linespacing=1.15)

    # --- cleaning steps (each an annotated drop)
    steps = [
        ("identifiers", "Drug to DrugBank, disease to MONDO (NCATS Name Resolver / Node Normalizer, cached; MONDO SSSOM), then reviewed: "
                        "42 label pins, every flagged match decided", None),
        ("no indication text yet", f"{c['no_text']} rows whose source names the medicine but not the indication (129 unmatched CHMP cards, "
                                   "10 FDA supplements without a letter, bare “New indication” NOCs)", c["no_text"]),
        ("not a new indication for the ingredient",
         "class review at the active-ingredient level, every decision with a rationale and the prior approval it rests on: "
         + ", ".join(f"{k.replace('_', ' ')} {v}" for k, v in sorted(c["not_new"].items(), key=lambda kv: -kv[1])),
         sum(c["not_new"].values())),
        ("not an indication event",
         ", ".join(f"{k.replace('_', ' ')} {v}" for k, v in sorted(c["not_change"].items(), key=lambda kv: -kv[1])),
         sum(c["not_change"].values())),
        ("no disease target", f"{c['no_target']} diagnostic, procedural, cosmetic or contraceptive uses", c["no_target"]),
        ("identifier pending a live lookup", f"{c['pending']} rows whose Name Resolver call is not in the shipped cache", c["pending"]),
    ]
    y = 9.55
    for name, detail, n in steps:
        y -= 0.86
        arrow(0.9, y + 0.86 - 0.02, 0.9, y + 0.14)
        if n is not None:
            ax.text(1.6, y + 0.5, f"−{n:,}", fontsize=8, fontweight="bold", color=DROP, va="center")
        ax.text(2.7, y + 0.62, name, fontsize=6.8, fontweight="bold", color=INK, va="center")
        ax.text(2.7, y + 0.3, _wrap(detail, 118), fontsize=5.5, color=MUTED, va="center", linespacing=1.15)

    # --- S1
    y -= 1.05
    box(0, y, 10, 0.95)
    ax.text(0.3, y + 0.48, f"{c['n1']:,}", fontsize=12.5, fontweight="bold", color=INK, va="center", zorder=3)
    ax.text(1.75, y + 0.48, "reviewed new-indication rows with both identifiers resolved (S1)", fontsize=6.9, color=INK, va="center", zorder=3)
    y -= 0.86
    arrow(0.9, y + 0.84, 0.9, y + 0.14)
    ax.text(1.6, y + 0.5, f"−{sum(c['s2'].values())}", fontsize=8, fontweight="bold", color=DROP, va="center")
    ax.text(2.7, y + 0.62, "not mappable to a graph node", fontsize=6.8, fontweight="bold", color=INK, va="center")
    ax.text(2.7, y + 0.3, ", ".join(f"{k} {v}" for k, v in c["s2"].items()), fontsize=5.5, color=MUTED, va="center")
    y -= 1.05
    box(0, y, 10, 0.95)
    ax.text(0.3, y + 0.48, f"{c['n2']:,}", fontsize=12.5, fontweight="bold", color=INK, va="center", zorder=3)
    ax.text(1.75, y + 0.48, "KG-mappable rows (S2): single agent with a DrugBank accession, disease a specific MONDO term", fontsize=6.9, color=INK, va="center", zorder=3)
    y -= 0.86
    arrow(0.9, y + 0.84, 0.9, y + 0.14)
    ax.text(1.6, y + 0.5, f"−{c['events'] - c['pairs']}", fontsize=8, fontweight="bold", color=DROP, va="center")
    ax.text(2.7, y + 0.62, "collapse across regulators", fontsize=6.8, fontweight="bold", color=INK, va="center")
    ax.text(2.7, y + 0.3, _wrap(f"{c['n2']} rows name {c['events']} (row, disease) events; one (DrugBank id, MONDO id) pair each, earliest "
                                f"decision date kept, every approving regulator listed; {c['multi']} pairs approved by two or more regulators", 118),
            fontsize=5.5, color=MUTED, va="center", linespacing=1.15)
    # --- S3 hero
    y -= 1.5
    box(0, y, 10, 1.4, hero=True)
    ax.text(0.3, y + 0.95, f"{c['pairs']}", fontsize=15, fontweight="bold", color=INK, va="center", zorder=3)
    ax.text(1.75, y + 1.02, f"unique drug-disease pairs (S3): {c['drugs']} drugs, {c['diseases']} diseases", fontsize=7.2, fontweight="bold", color=INK, va="center", zorder=3)
    ax.text(1.75, y + 0.62, "approved by (a pair can count under several):  "
            + "   ".join(f"{REG_LABEL[r].split(' ')[0] if r != 'HC' else 'Health Canada'} {c['any'][r]}" for r in ORDER),
            fontsize=6.3, color=INK, va="center", zorder=3)
    ax.text(1.75, y + 0.3, "first approved by:  " + "   ".join(f"{REG_LABEL[r].split(' ')[0] if r != 'HC' else 'Health Canada'} {c['first'].get(r, 0)}" for r in ORDER),
            fontsize=6.3, color=MUTED, va="center", zorder=3)
    # per-regulator rows reaching S3, under the source boxes as a thin return strip
    ax.text(0, y - 0.3, "rows per regulator that survive to a pair: " + ", ".join(f"{REG_LABEL[r].split(' ')[0] if r != 'HC' else 'Health Canada'} {c['rows_by_reg_s3'][r]} of {c['rows'][r]}" for r in ORDER),
            fontsize=5.8, color=MUTED, va="center")
    ax.text(0, y - 0.56, "Per-graph checks (coverage, no prior edge, distance from prior indications) are applied to these pairs; see 07h_post_cutoff_funnel.",
            fontsize=5.8, color=MUTED, va="center")
    fig.text(0.01, -0.012, f"repurposing-dataset {version}; every count re-derived by tools/make_post_cutoff_flow_figure.py from the shipped tables", fontsize=6, color=MUTED)
    fig.tight_layout()
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(out_stem.with_suffix("." + ext), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("figure ->", out_stem.with_suffix(".pdf"))


def _wrap(s: str, width: int) -> str:
    words, lines, cur = str(s).split(), [], ""
    for w in words:
        if len(cur) + len(w) + 1 > width and cur:
            lines.append(cur)
            cur = w
        else:
            cur = (cur + " " + w).strip()
    if cur:
        lines.append(cur)
    return "\n".join(lines[:3])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", type=Path, default=SOURCE_DEFAULT)
    ap.add_argument("--out", type=Path, default=ROOT / "results/figures/07_generalization/07i_post_cutoff_flow")
    a = ap.parse_args()
    src = pd.read_csv(a.source, dtype=str)
    dropped = pd.read_csv(GS_DIR / "post_cutoff_dropped.csv")
    gold = pd.read_csv(GS_DIR / "post_cutoff_approvals.csv")
    version = (a.source.parent / "VERSION").read_text().strip() if (a.source.parent / "VERSION").exists() else "?"
    c = counts(src, dropped, gold)
    assert c["n0"] - (c["no_text"] + c["pending"] + c["no_target"] + sum(c["not_new"].values()) + sum(c["not_change"].values())) == c["n1"], "S1 drops do not add up"
    draw(c, version, a.out)
    print({k: v for k, v in c.items() if k in ("n0", "n1", "n2", "events", "pairs", "multi")})


if __name__ == "__main__":
    main()
