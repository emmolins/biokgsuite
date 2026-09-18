#!/usr/bin/env python3
"""
Clean the 116 evaluable pairs into a defensible evaluation set.

Four steps, ordered so each is independently justifiable and each can be
stopped at. Every removed row is written out with its reason, so the cleaning
is auditable and reversible.

STEP 0  DEDUPLICATE (forced by the task definition)
  The evaluation asks "rank drug X for disease Y". Two rows with the same
  (drug_id, disease_id) are the same question, usually the same approval
  recorded by two agencies. They inflate n and are perfectly correlated in
  difficulty, which no clustered interval can repair. 116 rows -> 99 pairs.

STEP 1  FACTUAL REMOVALS (forced by evidence)
  Five pairs verified against the primary register as wrong: three approved
  before the window, one attributed to the wrong drug, one never authorised by
  the agency named. See EMA_PMDA_VERIFICATION.md.

STEP 2  NOT AN INDICATION CHANGE (forced by any definition)
  Formulation changes, vaccine strain updates and accelerated-to-full approval
  conversions. These add no drug-disease edge under any of the three candidate
  definitions of repurposing.

STEP 3  POPULATION AND SEVERITY EXPANSIONS (a definitional choice)
  The drug already treated that disease; the approval widened the age band or
  the severity group. These are invalid for a PROSPECTIVE test regardless of
  definition: the drug-disease pair is not new, so a model that ranks it highly
  is recalling an existing fact rather than predicting one. They survived the
  leakage filter only because the particular graph happened to lack that edge,
  which makes them a test of graph incompleteness, not of prediction.

USAGE
  python3 tools/clean_116.py --gold data/gold_standards/legacy/gold_standard_v4_bigset.tsv
"""
import argparse, json
from pathlib import Path
import pandas as pd

# Verified against the primary register; see EMA_PMDA_VERIFICATION.md
FACTUAL = {
    "95":  "PMDA approval 21 Nov 2018, outside the window",
    "96":  "PMDA approval 23 Mar 2021, outside the window",
    "252": "PMDA approval 28 Mar 2022, outside the window, and a different indication",
    "158": "EGPA in Japan belongs to benralizumab/mepolizumab; omalizumab has no such approval",
    "238": "measles prophylaxis is absent from the full EMA indication text; pair is EMA-only",
}
# Add no drug-disease edge under any definition of repurposing
NOT_INDICATION = {
    "11":  "accelerated-to-full approval conversion, same indication",
    "75":  "new IV formulation; the diseases were already authorised",
    "78":  "vaccine strain update, same disease",
    "213": "new oral formulation; obesity already authorised",
}
# Same disease, wider population or severity band
POPULATION = {
    "24":  "COVID-19 in renal impairment, severity subgroup",
    "51":  "psoriasis in children 6-11",
    "54":  "LDL-C, expanded risk population",
    "63":  "heterozygous familial hypercholesterolaemia, flagged non-strict",
    "71":  "SMA paediatric expansion",
    "79":  "chronic hepatitis B, paediatric",
    "81":  "primary immunodeficiency, paediatric",
    "163": "atopic dermatitis, paediatric",
    "206": "delay of T1D onset extended to children",
    "211": "PKU extended to 12 years and older",
}


def dn(x):
    return str(x).split(":")[-1].lstrip("0")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gold", required=True)
    ap.add_argument("--pairs", default="data/gold_standards/legacy/pairs_with_ids_v2.csv")
    ap.add_argument("--outdir", default="data/gold_standards/legacy")
    a = ap.parse_args()

    g = pd.read_csv(a.gold, sep="\t", dtype=str)
    p = pd.read_csv(a.pairs, dtype=str)
    g["_k"] = g.drug_id.astype(str) + "|" + g.disease_id.map(dn)
    p["_k"] = p.DrugBank_ID.astype(str) + "|" + p.DOID_ID.map(dn)
    g["Pair_ID"] = g._k.map(p.drop_duplicates("_k").set_index("_k").Pair_ID)

    removed = []
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)

    # STEP 0: deduplicate. Keep the earliest-window row, then the most
    # agency-attributed one, so the survivor is the best-evidenced instance.
    g["_rank"] = g.year_set.eq("2023").map({True: 0, False: 1})
    g["_nag"] = g.source_agencies.fillna("").str.count(r"[,+]")
    gs = g.sort_values(["_k", "_rank", "_nag"], ascending=[True, True, False])
    dup = gs[gs._k.duplicated(keep="first")]
    for _, r in dup.iterrows():
        removed.append({**r.to_dict(), "step": "0_duplicate",
                        "reason": "same (drug, disease) as a retained row"})
    s0 = gs[~gs._k.duplicated(keep="first")].copy()

    def drop(df, table, step):
        keep, gone = [], []
        for _, r in df.iterrows():
            pid = str(r.Pair_ID)
            if pid in table:
                gone.append({**r.to_dict(), "step": step, "reason": table[pid]})
            else:
                keep.append(r)
        removed.extend(gone)
        return pd.DataFrame(keep), len(gone)

    s1, n1 = drop(s0, FACTUAL, "1_factual")
    s2, n2 = drop(s1, NOT_INDICATION, "2_not_an_indication_change")
    s3, n3 = drop(s2, POPULATION, "3_population_expansion")

    cols = [c for c in g.columns if not c.startswith("_")]
    for name, df in (("v5_dedup", s0), ("v5_verified", s1),
                     ("v5_indication", s2), ("v5_strict", s3)):
        df[cols].to_csv(out / f"gold_standard_{name}.tsv", sep="\t", index=False)

    pd.DataFrame(removed)[[c for c in cols if c in pd.DataFrame(removed).columns]
                           + ["step", "reason"]].to_csv(
        out / "cleaning_manifest.csv", index=False)

    print(f"  {len(g):>4}  v4, as evaluated")
    print(f"  {len(s0):>4}  step 0: deduplicated        (-{len(g)-len(s0)} duplicate rows)")
    print(f"  {len(s1):>4}  step 1: factual removals    (-{n1})")
    print(f"  {len(s2):>4}  step 2: not an indication   (-{n2})")
    print(f"  {len(s3):>4}  step 3: population expansions (-{n3})   <- recommended set")
    print(f"\ndiseases: v4 {g.disease_id.nunique()} -> strict {s3.disease_id.nunique()}")
    print(f"drugs   : v4 {g.drug_id.nunique()} -> strict {s3.drug_id.nunique()}")
    print(f"\nwrote 4 sets + cleaning_manifest.csv ({len(removed)} removed rows) to {out}")
    (out / "cleaning_summary.json").write_text(json.dumps({
        "v4": len(g), "dedup": len(s0), "verified": len(s1),
        "indication": len(s2), "strict": len(s3),
        "n_removed": len(removed)}, indent=2))


if __name__ == "__main__":
    main()
