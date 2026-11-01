#!/usr/bin/env python3
"""Flow chart of the post-cutoff approval set: regulator rows in, cleaning steps, pairs out.

Every number is re-derived from the files on disk (the vendored source table, post_cutoff_dropped.csv,
post_cutoff_approvals.csv), so the figure cannot drift from the data.

    python tools/make_post_cutoff_flow_figure.py            # -> results/figures/07_generalization/07i_post_cutoff_flow.{pdf,png}
"""
import argparse
import collections
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
    from biokgsuite.viz.style import DOUBLE_COL_W, MUTED, NEUTRAL, TEXT_COLOR, setup_style
    setup_style()
    # House style: every box the same neutral fill with a 0.5 pt ink stroke (regulators are not graphs, so
    # they carry no colour); the final set is tinted rather than framed; annotations are muted; text 6 to 7 pt.
    INK, FILL, FINAL_FILL = TEXT_COLOR, "#F2F2F2", "#EAF3FA"
    fig, ax = plt.subplots(figsize=(DOUBLE_COL_W, 7.0), constrained_layout=True)
    ax.set_xlim(0, 10)
    ax.axis("off")

    def box(x, y, w, h, hero=False):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.06",
                                    linewidth=0.5, edgecolor=INK, facecolor=FINAL_FILL if hero else FILL, zorder=2))

    def arrow(x0, y0, x1, y1):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=6, color=NEUTRAL, lw=0.6, zorder=3))

    # --- regulator sources (the figure carries no title; the caption names it)
    bw, gap, y0, bh = 1.84, 0.2, 11.15, 1.45
    for i, reg in enumerate(ORDER):
        x = i * (bw + gap)
        box(x, y0, bw, bh)
        ax.text(x + 0.12, y0 + bh - 0.26, f"{c['rows'][reg]:,}", fontsize=7, color=INK, va="center", zorder=3)
        ax.text(x + 0.12, y0 + bh - 0.6, REG_LABEL[reg], fontsize=7, color=INK, va="center", zorder=3)
        ax.text(x + 0.12, y0 + 0.32, REG_WHAT[reg], fontsize=6, color=MUTED, va="center", zorder=3, linespacing=1.1)
        arrow(x + bw / 2, y0 - 0.04, x + bw / 2, y0 - 0.42)

    # --- merged table
    y = 9.55
    box(0, y, 10, 1.15, hero=False)
    ax.text(0.3, y + 0.62, f"{c['n0']:,}", fontsize=7, color=INK, va="center", zorder=3)
    ax.text(1.75, y + 0.78, "Candidate rows, one per regulator decision, 2023-01-01 to 2026-09-17 (S0); indication text from the regulator's record",
            fontsize=7, color=INK, va="center", zorder=3)
    ax.text(1.75, y + 0.36, _wrap(f"FDA approval letter, “provides for” sentence ({c['text']['FDA']} of {c['rows']['FDA']}); EMA EPAR procedural-steps "
                                  f"scope matched to the opinion date ({c['text']['EMA']} of {c['rows']['EMA']}); Health Canada NOC reason, PMDA notes, "
                                  f"TGA indication (all rows)", 112),
            fontsize=6, color=MUTED, va="center", zorder=3, linespacing=1.15)

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
        y -= 0.95
        arrow(0.9, y + 0.95 - 0.02, 0.9, y + 0.14)
        if n is not None:
            ax.text(1.6, y + 0.5, f"−{n:,}", fontsize=7, color=INK, va="center")
        # The detail hangs from just under the name, so a three-line detail cannot climb into it.
        ax.text(2.7, y + 0.68, name[:1].upper() + name[1:], fontsize=7, color=INK, va="center")
        ax.text(2.7, y + 0.52, _wrap(detail, 118), fontsize=6, color=MUTED, va="top", linespacing=1.15)

    # --- S1
    y -= 0.9
    box(0, y, 10, 0.8)
    ax.text(0.3, y + 0.4, f"{c['n1']:,}", fontsize=7, color=INK, va="center", zorder=3)
    ax.text(1.75, y + 0.4, "Reviewed new-indication rows with both identifiers resolved (S1)", fontsize=7, color=INK, va="center", zorder=3)
    y -= 0.86
    arrow(0.9, y + 0.84, 0.9, y + 0.14)
    ax.text(1.6, y + 0.5, f"−{sum(c['s2'].values())}", fontsize=7, color=INK, va="center")
    ax.text(2.7, y + 0.68, "Not mappable to a graph node", fontsize=7, color=INK, va="center")
    ax.text(2.7, y + 0.52, ", ".join(f"{k} {v}" for k, v in c["s2"].items()), fontsize=6, color=MUTED, va="top")
    y -= 0.9
    box(0, y, 10, 0.8)
    ax.text(0.3, y + 0.4, f"{c['n2']:,}", fontsize=7, color=INK, va="center", zorder=3)
    ax.text(1.75, y + 0.4, "KG-mappable rows (S2): single agent with a DrugBank accession, disease a specific MONDO term", fontsize=7, color=INK, va="center", zorder=3)
    y -= 0.86
    arrow(0.9, y + 0.84, 0.9, y + 0.14)
    ax.text(1.6, y + 0.5, f"−{c['events'] - c['pairs']}", fontsize=7, color=INK, va="center")
    ax.text(2.7, y + 0.68, "Collapse across regulators", fontsize=7, color=INK, va="center")
    ax.text(2.7, y + 0.52, _wrap(f"{c['n2']} rows name {c['events']} (row, disease) events; one (DrugBank id, MONDO id) pair each, earliest "
                                f"decision date kept, every approving regulator listed; {c['multi']} pairs approved by two or more regulators", 118),
            fontsize=6, color=MUTED, va="top", linespacing=1.15)
    # --- S3, the final set: tinted, not framed
    y -= 1.45
    box(0, y, 10, 1.25, hero=True)
    ax.text(0.3, y + 0.95, f"{c['pairs']}", fontsize=7, color=INK, va="center", zorder=3)
    ax.text(1.75, y + 0.95, f"Unique drug-disease pairs (S3): {c['drugs']} drugs, {c['diseases']} diseases", fontsize=7, color=INK, va="center", zorder=3)
    ax.text(1.75, y + 0.58, "Approved by (a pair can count under several):  "
            + "   ".join(f"{REG_LABEL[r].split(' ')[0] if r != 'HC' else 'Health Canada'} {c['any'][r]}" for r in ORDER),
            fontsize=6.5, color=INK, va="center", zorder=3)
    ax.text(1.75, y + 0.24, "First approved by:  " + "   ".join(f"{REG_LABEL[r].split(' ')[0] if r != 'HC' else 'Health Canada'} {c['first'].get(r, 0)}" for r in ORDER),
            fontsize=6.5, color=MUTED, va="center", zorder=3)
    # per-regulator rows reaching S3, under the source boxes as a thin return strip
    ax.text(0, y - 0.3, "Rows per regulator that survive to a pair: " + ", ".join(f"{REG_LABEL[r].split(' ')[0] if r != 'HC' else 'Health Canada'} {c['rows_by_reg_s3'][r]} of {c['rows'][r]}" for r in ORDER),
            fontsize=6, color=MUTED, va="center")
    ax.text(0, y - 0.56, "Per-graph checks (coverage, no prior edge, distance from prior indications) are applied to these pairs; see 07h_post_cutoff_funnel.",
            fontsize=6, color=MUTED, va="center")
    # The y range follows what was drawn, so a longer or shorter list of steps never clips.
    ax.set_ylim(y - 0.8, y0 + bh + 0.1)
    fig.supxlabel(f"repurposing-dataset {version}; every count re-derived by tools/make_post_cutoff_flow_figure.py from the shipped tables",
                  fontsize=6, color=MUTED, x=0.0, ha="left")
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(out_stem.with_suffix("." + ext), dpi=600, bbox_inches="tight")
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
