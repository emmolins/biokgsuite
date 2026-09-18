#!/usr/bin/env python3
"""
Corroborate the 253 candidate pairs against the FDA indication-expansion census.

PURPOSE, AND WHAT THIS IS NOT
-----------------------------
The 253 pairs are the benchmark's candidate set and the basis of the grounding
experiment. This script does NOT rebuild or replace them. It answers one
question the 253 could not answer on its own: for each pair attributed to FDA,
does an actual FDA indication-expansion event exist for that drug in that year?

That matters because PROVENANCE_253_candidate_pairs.md records that the 2023
portion (96 pairs) was reconstructed by an LLM agent from agency records rather
than read out of a regulatory database. Corroboration against the authoritative
record is the check a reviewer would ask for, and it costs no experiment reruns.

INTERPRETING A MISS
-------------------
A pair with no census match is NOT thereby shown to be false. Three innocent
explanations exist and are reported separately:
  * the event is real but outside the census tiers (new patient population,
    dosing regimen, labeling change with clinical data)
  * the product is CBER-regulated (vaccines, cell and gene therapy) and so is
    absent from Drugs@FDA entirely
  * the approval year in the 253 differs by one from the FDA action date
Only pairs with no FDA record of any kind for that drug in that year are
flagged for manual review.

USAGE
  python3 tools/corroborate_253.py \
    --pairs  data/gold_standards/legacy/pairs_with_ids_v2.csv \
    --census data/gold_standards/legacy/fda_indication_expansion_census_drugbank.csv \
    --frame  data/gold_standards/legacy/fda_indication_expansion_census_frame.csv
"""
import argparse, json, re
from pathlib import Path
import pandas as pd

CBER_HINTS = ("vaccine", "ciloleucel", "cel-t", "autoleucel", "maraleucel",
              "cabtagene", "leucel", "gene therapy")


def norm(s):
    return re.sub(r"\s+", " ", re.sub(r"[^a-z0-9 ]+", " ", str(s).lower())).strip()


def strip_bio(s):
    return re.sub(r"-[a-z]{4}$", "", str(s).strip(), flags=re.I)


def ingredient_keys(raw):
    """All normalised component names for a census or pair ingredient string."""
    out = set()
    for part in re.split(r"[|;,/+]| and ", str(raw)):
        p = norm(strip_bio(part))
        if p:
            out.add(p)
            toks = p.split()
            if len(toks) > 1:
                out.add(" ".join(toks[1:]))   # fam-trastuzumab deruxtecan
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--census", required=True)
    ap.add_argument("--frame", required=True)
    ap.add_argument("--out", default="results/tables/corroborate_253.csv")
    a = ap.parse_args()

    pairs = pd.read_csv(a.pairs, dtype=str)
    cen = pd.read_csv(a.census, dtype=str)
    frame = pd.read_csv(a.frame, dtype=str)

    cen["year"] = pd.to_datetime(cen["approval_date"], errors="coerce").dt.year
    frame["year"] = pd.to_datetime(frame["approval_date"], errors="coerce").dt.year

    # index: (ingredient key, year) -> rows, for the census and for the wider frame
    def build_index(df):
        idx = {}
        for i, r in df.iterrows():
            ks = ingredient_keys(r.get("ActiveIngredients", ""))
            ks |= ingredient_keys(r.get("DrugNames", ""))
            for k in ks:
                idx.setdefault((k, r["year"]), []).append(i)
        return idx

    cen_idx = build_index(cen)
    eff = frame[frame["SubmissionClassCode"].fillna("").eq("EFFICACY")]
    eff_idx = build_index(eff)
    any_idx = build_index(frame)

    # DrugBank index too, since the 253 carries DrugBank IDs
    db_idx = {}
    for i, r in cen.iterrows():
        for d in str(r.get("drugbank_ids", "") or "").split("|"):
            if d and d != "nan":
                db_idx.setdefault((d, r["year"]), []).append(i)

    rows = []
    for _, p in pairs.iterrows():
        agencies = str(p.get("Agencies", ""))
        is_fda = "FDA" in agencies.upper()
        # YearSet is a SET label, not a year: either "2023" or "2024-2026".
        # Treating it as a year silently skipped two thirds of the FDA pairs.
        ys = str(p.get("YearSet", "")).strip()
        if "-" in ys:
            lo, hi = ys.split("-")[:2]
            years = list(range(int(lo), int(hi) + 1))
        elif ys.isdigit():
            years = [int(ys)]
        else:
            years = []
        keys = ingredient_keys(p.get("Active", "")) | ingredient_keys(p.get("Drug", ""))
        dbids = [d.strip() for d in str(p.get("DrugBank_ID", "") or "").split(",")
                 if d.strip() and d.strip() != "nan"]

        def lookup(idx, yr):
            hits = []
            for k in keys:
                hits += idx.get((k, yr), [])
            return sorted(set(hits))

        status, evidence = "", ""
        def lookup_years(idx):
            hits = []
            for y in years:
                hits += lookup(idx, y)
            return sorted(set(hits))

        if not is_fda:
            status = "not_fda_attributed"
        elif not years:
            status = "no_year"
        else:
            hits = lookup_years(cen_idx)
            for d in dbids:
                for y in years:
                    hits += db_idx.get((d, y), [])
            hits = sorted(set(hits))
            if hits:
                status = "corroborated_in_census"
            else:
                e = lookup_years(eff_idx)
                if e:
                    status, hits = "efficacy_event_outside_tiers", e
                else:
                    anyh = lookup_years(any_idx)
                    if anyh:
                        status, hits = "fda_action_but_not_efficacy", anyh
                    elif any(h in norm(p.get("Active", "")) for h in CBER_HINTS):
                        status = "out_of_frame_cber"
                    else:
                        status = "NO_FDA_RECORD_review"
            if hits:
                src = cen if status.startswith("corroborated") else (
                    eff if status == "efficacy_event_outside_tiers" else frame)
                evidence = "; ".join(
                    f"{src.loc[i, 'approval_date']} {src.loc[i, 'application']} "
                    f"{src.loc[i, 'supplement']}" for i in hits[:3])

        rows.append({
            "Pair_ID": p.get("Pair_ID"), "Drug": p.get("Drug"),
            "Active": p.get("Active"), "DrugBank_ID": p.get("DrugBank_ID"),
            "YearSet": p.get("YearSet"), "Agencies": agencies,
            "NewIndication": p.get("NewIndication"),
            "StrictRepurposing": p.get("StrictRepurposing"),
            "corroboration": status, "fda_evidence": evidence,
        })

    res = pd.DataFrame(rows)
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out, index=False)

    print(f"{len(pairs)} candidate pairs\n")
    vc = res["corroboration"].value_counts()
    for k, v in vc.items():
        print(f"  {v:>4}  {k}")

    fda = res[res["corroboration"].ne("not_fda_attributed")]
    ok = fda["corroboration"].str.startswith("corroborated").sum()
    print(f"\nFDA-attributed pairs: {len(fda)}")
    print(f"  directly corroborated in the census: {ok} "
          f"({100*ok/max(len(fda),1):.1f}%)")
    flag = res[res["corroboration"].eq("NO_FDA_RECORD_review")]
    print(f"  flagged for manual review: {len(flag)}")
    if len(flag):
        print(flag[["Pair_ID", "Drug", "Active", "YearSet",
                    "NewIndication"]].to_string(index=False))
    print(f"\nwrote {out}")

    Path(str(out) + ".summary.json").write_text(json.dumps({
        "n_pairs": int(len(pairs)),
        "n_fda_attributed": int(len(fda)),
        "n_corroborated": int(ok),
        "counts": {k: int(v) for k, v in vc.items()},
    }, indent=2))


if __name__ == "__main__":
    main()
