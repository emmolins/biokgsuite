#!/usr/bin/env python3
"""
Decompose the nb09 KG effect into COVERAGE and QUALITY-GIVEN-COVERAGE.

The single-number "does KG X help the LLM?" conflates two things:

  coverage  -- how often the graph has any evidence for the queried pair
               (PrimeKG 111, DRKG 90, BioKG 78 of 116)
  quality   -- how good that evidence is WHEN IT EXISTS, which is only
               estimable on common support: the 72 pairs all three cover

A graph can look good because it speaks often, or because it speaks well.
Comparing arms on the full 116 confounds the two, because each arm is then
scored on a different effective subset. Restricting to the 72-pair common
core holds the queried pairs fixed, so the remaining difference between
arms is attributable to evidence quality.

Outcome is per-pair MRR lift over the no-KG arm, pooled over all model
families and seeds. Lift absorbs pair difficulty (a pair that is hard for
every arm contributes ~0 rather than dragging the level down).

Inputs   data/gold_standards/legacy/coverage_annotation_v4.csv   (Pair_ID x kg)
         data/gold_standards/legacy/gold_standard_v4_bigset.tsv  (the 116)
         results/tables/09_llm_runs/09_big_*.csv          (run outcomes)
Outputs  results/tables/09_coverage_vs_quality.csv        (per kg x stratum)
         results/tables/09_coverage_vs_quality_bypair.csv (per pair x kg)

Usage:  python scripts/coverage_vs_quality.py
"""
import glob
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
if str(BASE) not in sys.path:
    sys.path.insert(0, str(BASE))

GOLD = BASE / "data" / "gold_standards"
RUNS = BASE / "results" / "tables" / "09_llm_runs"
OUT = BASE / "results" / "tables"

KGS = ["primekg", "drkg", "biokg"]
N_BOOT = 2000
SEED = 42

# ---------------------------------------------------------------- coverage
cov = pd.read_csv(GOLD / "legacy" / "coverage_annotation_v4.csv")
gold = pd.read_csv(GOLD / "legacy" / "gold_standard_v4_bigset.tsv", sep="\t")

# covered = drug_in AND disease_in (the definition behind 111/90/78).
wide = (cov.pivot(index="Pair_ID", columns="kg", values="covered")
           .reindex(columns=KGS).fillna(0).astype(int))
meta = cov[["Pair_ID", "disease_name", "drug_id"]].drop_duplicates("Pair_ID") \
                                                  .set_index("Pair_ID")
wide = wide.join(meta)
wide["n_kg"] = wide[KGS].sum(axis=1)
COMMON = wide.index[wide["n_kg"] == 3]          # the 72

assert len(wide) == 116, f"expected 116 pairs, got {len(wide)}"

# ---------------------------------------------------------------- outcomes
frames = []
for f in sorted(glob.glob(str(RUNS / "09_big_*.csv"))):
    df = pd.read_csv(f)
    df["src"] = os.path.basename(f)
    frames.append(df)
runs = pd.concat(frames, ignore_index=True)
# Do NOT drop failed parses. The protocol scores them as misses and the
# run files already encode that (rank = pool_size + 1, RR = 1/9).
# Filtering them would remove 157 legitimate misses, unevenly across
# arms (no_kg 50, primekg 49, biokg 32, drkg 26), biasing the lift.
n_unparsed = int((~runs["parsed"].fillna(1).astype(bool)).sum())

# Arm label: 'no_kg' or the KG name.
runs["arm"] = np.where(runs["condition"] == "no_kg", "no_kg", runs["kg"])
runs = runs[runs["arm"].isin(["no_kg"] + KGS)]

# Collapse to one number per (pair-proxy, arm): mean RR over models x seeds
# x shuffles. Disease name is the join key; verified lossless because the 7
# diseases carrying two drugs agree on coverage across all three KGs.
cell = (runs.groupby(["disease", "arm"], as_index=False)["reciprocal_rank"]
            .agg(rr="mean", n_obs="size"))

pv = cell.pivot(index="disease", columns="arm", values="rr")
nb = cell.pivot(index="disease", columns="arm", values="n_obs")
pv = pv.dropna(subset=["no_kg"])

# Attach coverage by disease name.
# Prefix coverage columns: pv already has arm columns named after the KGs.
cov_by_dis = (wide.reset_index()
                  .groupby("disease_name")[KGS].max()      # max == any pair
                  .rename(columns={k: f"cov_{k}" for k in KGS}))
pv = pv.join(cov_by_dis, how="inner")

for kg in KGS:
    pv[f"lift_{kg}"] = pv[kg] - pv["no_kg"]

print(f"diseases joined: {len(pv)}  (runs cover {runs.disease.nunique()})")
print(f"observations pooled: {len(runs):,} rows from {len(frames)} run files, "
      f"{runs.model.nunique()} model families, seeds {sorted(runs.seed.unique())}")
print(f"failed parses retained as misses: {n_unparsed}")

# ---------------------------------------------------------------- bootstrap
rng = np.random.RandomState(SEED)


def boot_ci(vals, n_boot=N_BOOT):
    """Percentile 95% CI of the mean; (nan, nan) if fewer than 3 values."""
    v = np.asarray([x for x in vals if np.isfinite(x)], dtype=float)
    if v.size < 3:
        return np.nan, np.nan
    idx = rng.randint(0, v.size, size=(n_boot, v.size))
    means = v[idx].mean(axis=1)
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


rows = []
common_dis = set(wide.loc[COMMON, "disease_name"])

for kg in KGS:
    n_cov = int(wide[kg].sum())

    # --- stratum 1: all 116 (what a naive comparison reports) ---
    all_lift = pv[f"lift_{kg}"].dropna()
    lo_a, hi_a = boot_ci(all_lift)

    # --- stratum 2: pairs THIS graph covers ---
    own = pv[pv[f"cov_{kg}"] == 1][f"lift_{kg}"].dropna()
    lo_o, hi_o = boot_ci(own)

    # --- stratum 3: the 72-pair common core (quality given coverage) ---
    comm = pv[pv.index.isin(common_dis)][f"lift_{kg}"].dropna()
    lo_c, hi_c = boot_ci(comm)

    # --- stratum 4: pairs this graph does NOT cover (should be ~0) ---
    unc = pv[pv[f"cov_{kg}"] == 0][f"lift_{kg}"].dropna()
    lo_u, hi_u = boot_ci(unc)

    for stratum, s, lo, hi, note in [
        ("all_116", all_lift, lo_a, hi_a, "confounds coverage with quality"),
        ("covered_by_this_kg", own, lo_o, hi_o, "own support, not comparable"),
        ("common_core_72", comm, lo_c, hi_c, "QUALITY GIVEN COVERAGE"),
        ("not_covered", unc, lo_u, hi_u, "sanity: expect no lift"),
    ]:
        rows.append({
            "kg": kg,
            "coverage_n": n_cov,
            "coverage_pct": round(100 * n_cov / 116, 1),
            "stratum": stratum,
            "n_diseases": int(s.size),
            "mean_lift": round(float(s.mean()), 4) if s.size else np.nan,
            "ci_lo": round(lo, 4) if np.isfinite(lo) else np.nan,
            "ci_hi": round(hi, 4) if np.isfinite(hi) else np.nan,
            "note": note,
        })

res = pd.DataFrame(rows)
res.to_csv(OUT / "09_coverage_vs_quality.csv", index=False)

bypair = pv.reset_index()
bypair.to_csv(OUT / "09_coverage_vs_quality_bypair.csv", index=False)

# ---------------------------------------------------------------- report
print("\n" + "=" * 78)
print("COVERAGE  --  how often the graph has evidence (of 116 pairs)")
print("=" * 78)
for kg in KGS:
    n = int(wide[kg].sum())
    print(f"  {kg:<10s} {n:3d} / 116   ({100 * n / 116:4.1f}%)")
print(f"\n  covered by all three (common core): {len(COMMON)}")
print(f"  covered by exactly 2:               {int((wide.n_kg == 2).sum())}")
print(f"  covered by exactly 1:               {int((wide.n_kg == 1).sum())}")
print(f"  covered by none:                    {int((wide.n_kg == 0).sum())}")

print("\n" + "=" * 78)
print("QUALITY GIVEN COVERAGE  --  mean per-pair MRR lift over no-KG")
print("=" * 78)
print(f"{'kg':<10}{'stratum':<22}{'n':>5}{'lift':>9}{'95% CI':>20}")
print("-" * 78)
for _, r in res.iterrows():
    ci = (f"[{r.ci_lo:+.4f}, {r.ci_hi:+.4f}]"
          if pd.notna(r.ci_lo) else "n/a")
    print(f"{r.kg:<10}{r.stratum:<22}{r.n_diseases:>5}"
          f"{r.mean_lift:>+9.4f}{ci:>20}")

print("\n" + "-" * 78)
print("Ranking on the naive all-116 comparison vs the common core:")
naive = res[res.stratum == "all_116"].sort_values("mean_lift", ascending=False)
core = res[res.stratum == "common_core_72"].sort_values("mean_lift", ascending=False)
print("  all 116     :  " + "  >  ".join(naive.kg))
print("  common 72   :  " + "  >  ".join(core.kg))
if list(naive.kg) != list(core.kg):
    print("  -> ORDER CHANGES once coverage is held fixed.")
else:
    print("  -> order is unchanged; the coverage advantage is not driving it.")

# ------------------------------------------------------- exact decomposition
# lift(all) = P(covered) * lift(covered) + P(uncovered) * lift(uncovered)
# The first term is the coverage channel, the second is what the arm still
# earns on pairs it has no evidence for (distractor discrimination in a
# ranking task -- see note below).
print("\n" + "=" * 78)
print("DECOMPOSITION of the all-116 lift")
print("=" * 78)
print(f"{'kg':<10}{'P(cov)':>8}{'lift|cov':>10}{'covered':>10}"
      f"{'lift|unc':>10}{'uncov':>9}{'total':>9}")
print("-" * 78)
for kg in KGS:
    r = res[res.kg == kg].set_index("stratum")
    n_c, n_u = r.loc["covered_by_this_kg", "n_diseases"], r.loc["not_covered", "n_diseases"]
    p_c = n_c / (n_c + n_u)
    l_c = r.loc["covered_by_this_kg", "mean_lift"]
    l_u = r.loc["not_covered", "mean_lift"] if n_u else 0.0
    t_c, t_u = p_c * l_c, (1 - p_c) * l_u
    print(f"{kg:<10}{p_c:>8.3f}{l_c:>+10.4f}{t_c:>+10.4f}"
          f"{l_u:>+10.4f}{t_u:>+9.4f}{t_c + t_u:>+9.4f}")

print("\nPairwise gaps -- naive vs coverage-held-fixed:")
core = res[res.stratum == "common_core_72"].set_index("kg")["mean_lift"]
naive = res[res.stratum == "all_116"].set_index("kg")["mean_lift"]
for a, b in [("primekg", "drkg"), ("primekg", "biokg"), ("biokg", "drkg")]:
    gn, gc = naive[a] - naive[b], core[a] - core[b]
    share = 100 * (1 - gc / gn) if gn else float("nan")
    print(f"  {a:>8s} - {b:<9s} naive {gn:+.4f} -> common core {gc:+.4f}"
          f"   ({share:4.0f}% of the gap was coverage)")

print(f"\nWrote {OUT / '09_coverage_vs_quality.csv'}")
print(f"Wrote {OUT / '09_coverage_vs_quality_bypair.csv'}")
