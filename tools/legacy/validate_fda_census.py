#!/usr/bin/env python3
"""
Verify the Drugs@FDA census against the previously vetted HTML-scrape dataset.

The scrape (fda_repurposing_events_2023_2026.csv) was audited in
FDA_REPURPOSING_VETTING.md and found ACCURATE but INCOMPLETE:
  * 8/8 internal integrity checks passed, 0 violations
  * 7/8 records confirmed field-for-field against openFDA, 0 fabrications
  * recall ~68% on the audited busy month (2025-03), truncation-driven
  * BLA supplements excluded by construction (blank classification in reports)

That makes it a high-precision, low-recall probe, which is exactly what is
needed to test a census: every scrape row SHOULD reappear in the census. A
scrape row absent from the census is a census defect and must be explained.
The reverse direction is expected and is the whole point of rebuilding.

Checks run here:
  A. containment   - every scraped (application, supplement) present in the frame
  B. field agreement - approval dates agree where both sources have the row
  C. targeted recall - the 10 identifiers the scrape missed in 2025-03
  D. biologics      - the 2024-09 BLA supplements the scrape structurally dropped
  E. month profile  - per-month counts, census vs scrape, to expose any new bias
  F. known defect   - the everolimus 2024-05-09 row with no FDA counterpart

USAGE
  python3 tools/validate_fda_census.py \
      --frame data/gold_standards/legacy/fda_indication_expansion_census_frame.csv \
      --census data/gold_standards/legacy/fda_indication_expansion_census.csv \
      --scrape data/gold_standards/legacy/fda_repurposing_events_2023_2026.csv
"""
import argparse, sys
from pathlib import Path
import pandas as pd

# Recorded in FDA_REPURPOSING_VETTING.md, test 3: present in the source table for
# 2025-03 but absent from the scrape. A complete census must contain all ten.
MISSED_2025_03 = [
    ("207103", "21"), ("211675", "27"), ("212436", "9"), ("214679", "9"),
    ("214916", "2"), ("215515", "5"), ("215833", "24"), ("215841", "6"),
    ("217171", "3"), ("218347", "2"),
]

# Recorded in test 4: BLA supplements from 2024-09 carrying a blank Submission
# Classification in the monthly reports, and therefore invisible to the scrape.
BLA_2024_09 = [
    ("DUPIXENT", "66"), ("TREMFYA", "21"), ("ELREXFIO", "2"),
    ("CIMZIA", "275"), ("SOLIRIS", "447"),
]


def norm_appl(s: pd.Series) -> pd.Series:
    """Strip the NDA/BLA prefix and any zero padding so the two sources compare."""
    return (s.astype(str).str.upper().str.replace(r"^(NDA|BLA|ANDA)", "", regex=True)
             .str.strip().str.lstrip("0"))


def norm_suppl(s: pd.Series) -> pd.Series:
    return (s.astype(str).str.upper().str.replace(r"^SUPPL[-_ ]?", "", regex=True)
             .str.strip().str.lstrip("0"))


def norm_type(s: pd.Series) -> pd.Series:
    """ORIG and SUPPL share a numbering space: application 214070 has both an
    ORIG-1 and a SUPPL-1. Keying on number alone silently matches the original
    approval to a supplement and manufactures date disagreements."""
    t = s.astype(str).str.upper().str.strip()
    return t.where(t.isin(["ORIG", "SUPPL"]), "SUPPL")


def key(df, appl_col, suppl_col, type_col=None):
    if type_col is not None and type_col in df.columns:
        t = norm_type(df[type_col])
    else:
        t = df[suppl_col].astype(str).str.upper().str.contains("ORIG").map(
            {True: "ORIG", False: "SUPPL"})
    return norm_appl(df[appl_col]) + "|" + t + "|" + norm_suppl(df[suppl_col])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frame", required=True)
    ap.add_argument("--census", required=True)
    ap.add_argument("--scrape", required=True)
    a = ap.parse_args()

    frame = pd.read_csv(a.frame, dtype=str)
    census = pd.read_csv(a.census, dtype=str)
    scrape = pd.read_csv(a.scrape, dtype=str)

    print(f"frame  {len(frame):,} approved NDA/BLA supplements")
    print(f"census {len(census):,} indication-expansion events")
    print(f"scrape {len(scrape):,} rows (vetted, high precision, ~68% recall)\n")

    def pick(cands, pred, label):
        for c in cands:
            if c in scrape.columns:
                return c
        hit = next((c for c in scrape.columns if pred(c.lower())), None)
        if hit is None:
            raise SystemExit(f"cannot locate the {label} column in {list(scrape.columns)}")
        return hit

    sc_appl = pick(["application_number", "application"],
                   lambda c: "appl" in c, "application")
    sc_supp = pick(["submission", "supplement", "submission_number"],
                   lambda c: "suppl" in c or c == "submission", "supplement")
    sc_date = pick(["approval_date"], lambda c: "date" in c, "date")
    print(f"scrape columns used: application={sc_appl} supplement={sc_supp} date={sc_date}\n")

    frame["_k"] = key(frame, "ApplNo", "SubmissionNo", "SubmissionType")
    census["_k"] = key(census, "ApplNo", "SubmissionNo", "SubmissionType")
    scrape["_k"] = key(scrape, sc_appl, sc_supp)

    fail = 0

    # A. containment
    missing = sorted(set(scrape["_k"]) - set(frame["_k"]))
    print(f"[A] containment: {len(scrape['_k'].unique()) - len(missing)}"
          f"/{scrape['_k'].nunique()} scraped rows found in the census frame")
    if missing:
        fail += 1
        print(f"    {len(missing)} scraped rows NOT in the frame (investigate each):")
        for m in missing[:25]:
            print("      ", m)
    print()

    # B. date agreement
    fd = frame.drop_duplicates("_k").set_index("_k")["approval_date"]
    sd = scrape.drop_duplicates("_k").set_index("_k")[sc_date]
    both = fd.index.intersection(sd.index)
    sd_norm = pd.to_datetime(sd.loc[both], errors="coerce").dt.strftime("%Y-%m-%d")
    fd_norm = pd.to_datetime(fd.loc[both], errors="coerce").dt.strftime("%Y-%m-%d")
    disagree = both[(sd_norm.values != fd_norm.values)]
    print(f"[B] date agreement on {len(both):,} shared rows: "
          f"{len(both) - len(disagree):,} agree, {len(disagree):,} disagree")
    if len(disagree):
        fail += 1
        for d in list(disagree)[:25]:
            print(f"      {d}  census={fd.loc[d]}  scrape={sd.loc[d]}")
    print()

    # C. targeted recall of the ten the scrape missed
    want = {f"{a_.lstrip('0')}|SUPPL|{s_.lstrip('0')}" for a_, s_ in MISSED_2025_03}
    got = want & set(frame["_k"])
    print(f"[C] 2025-03 gap: {len(got)}/{len(want)} of the identifiers the scrape "
          f"missed are present in the census frame")
    if got != want:
        fail += 1
        print("    still missing:", sorted(want - got))
    print()

    # D. biologics
    print("[D] biologics recovered (2024-09 BLA supplements the scrape could not see):")
    names = frame["DrugNames"].fillna("") + " " + frame["ActiveIngredients"].fillna("")
    for brand, suppl in BLA_2024_09:
        hit = frame[names.str.upper().str.contains(brand, na=False)
                    & norm_suppl(frame["SubmissionNo"]).eq(suppl.lstrip("0"))]
        mark = "found" if len(hit) else "MISSING"
        extra = f" {hit.iloc[0]['approval_date']} {hit.iloc[0]['application']}" if len(hit) else ""
        print(f"    {brand:<10} SUPPL-{suppl:<4} {mark}{extra}")
        if not len(hit):
            fail += 1
    n_bla_frame = int(frame["ApplType"].eq("BLA").sum()) if "ApplType" in frame else -1
    n_bla_scrape = int(scrape[sc_appl].astype(str).str.upper().str.startswith("BLA").sum())
    print(f"    BLA rows: census frame {n_bla_frame:,} vs scrape {n_bla_scrape:,}\n")

    # E. month profile
    # The scrape could only ever see EFFICACY-class supplements, so comparing it
    # against the whole frame (which includes LABELING, CMC, REMS) understates its
    # recall and hides the pattern we are testing for. Compare like with like.
    eff = frame[frame["SubmissionClassCode"].fillna("").eq("EFFICACY")
                & frame["SubmissionType"].fillna("").str.upper().str.strip().eq("SUPPL")]
    fm = pd.to_datetime(eff["approval_date"], errors="coerce").dt.to_period("M")
    sm = pd.to_datetime(scrape[sc_date], errors="coerce").dt.to_period("M")
    prof = pd.concat([fm.value_counts().rename("efficacy_suppl"),
                      sm.value_counts().rename("scrape")], axis=1).fillna(0).astype(int)
    prof = prof.sort_index()
    prof["scrape_recall_pct"] = (100 * prof["scrape"] /
                                 prof["efficacy_suppl"].replace(0, pd.NA)).round(1)
    print("[E] per-month recall of the scrape against the like-for-like census "
          "subset (approved EFFICACY-class supplements). Correlation between "
          "month size and recall is the truncation signature.")
    print(prof.to_string())
    r = prof.dropna(subset=["scrape_recall_pct"])
    if len(r) > 2:
        rho = r["efficacy_suppl"].corr(r["scrape_recall_pct"], method="spearman")
        print(f"\n    Spearman(month size, scrape recall) = {rho:+.3f}  "
              f"(negative = the scrape lost more in busy months)")
        print(f"    overall scrape recall = "
              f"{100*prof['scrape'].sum()/prof['efficacy_suppl'].sum():.1f}%")
    print()

    # F. known defect
    ev = frame[names.str.upper().str.contains("EVEROLIMUS", na=False)]
    ev2024 = ev[ev["approval_date"].astype(str).str.startswith("2024")]
    print("[F] everolimus 2024 submissions in the census frame "
          f"({len(ev2024)} rows):")
    if len(ev2024):
        print(ev2024[["approval_date", "application", "supplement",
                      "SubmissionClassCode", "SupplCategoryL2"]].to_string(index=False))
    else:
        print("    none, which corroborates the vetting finding that the "
              "everolimus FDA 2024-05-09 entry has no FDA counterpart")
    print()

    print("VERDICT:", "PASS" if fail == 0 else f"{fail} CHECK(S) FAILED")
    return 0 if fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
