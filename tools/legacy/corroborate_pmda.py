#!/usr/bin/env python3
"""
Corroborate the PMDA-attributed candidate pairs against PMDA's own
"List of Approved Drugs".

Each PMDA record states whether the approval was a New Approval or a Partial
Change, and the Notes sentence names the indication. A partial change with a
"new indication" note is the Japanese regulatory event corresponding to
repurposing, so both the event type and the disease can be checked, not just the
drug and the year.

Grades mirror the EMA check:
  indication_match    drug matched, window matched, disease present in the Notes
  partial_match       drug and window matched, disease partially present
  drug_year_only      drug matched in window, Notes do not name the disease
  drug_outside_window drug appears in the list, but not in 2023-2026
  no_pmda_record      drug not found in the list at all

USAGE
  python3 tools/corroborate_pmda.py --pmda data/gold_standards/legacy/pmda_approvals.csv \
      --pairs data/gold_standards/legacy/pairs_with_ids_v2.csv
"""
import argparse, json, re
from pathlib import Path
import pandas as pd

MODIFIERS = {
    "1l", "2l", "3l", "first", "second", "third", "line", "adjuvant", "neoadjuvant",
    "maintenance", "advanced", "metastatic", "unresectable", "locally", "relapsed",
    "refractory", "recurrent", "previously", "treated", "untreated", "adult",
    "adults", "paediatric", "pediatric", "children", "monotherapy", "combination",
    "with", "and", "or", "the", "of", "in", "for", "patients", "patient",
    "treatment", "therapy", "expansion", "extension", "new", "indication",
    "years", "aged", "age", "stage", "including", "plus", "positive", "negative",
    "high", "low", "risk", "severe", "moderate", "mild", "chronic", "acute",
    "primary", "secondary", "prophylaxis", "prevention", "japan", "jpn",
}


def norm(s):
    return re.sub(r"\s+", " ", re.sub(r"[^a-z0-9 ]+", " ", str(s).lower())).strip()


def disease_core(text):
    t = re.sub(r"\([^)]*\)", " ", str(text))
    t = re.sub(r"\+\s*\w+", " ", t)
    return [w for w in norm(t).split() if w not in MODIFIERS and len(w) > 2]


def load_labels(path="data/gold_standards/bridges/do_diseases.csv"):
    import os
    if not os.path.exists(path):
        return {}
    d = pd.read_csv(path, dtype=str); lut = {}
    for _, r in d.iterrows():
        for k in (r.get("doid"), r.get("mondo_id")):
            if pd.notna(k) and pd.notna(r.get("mondo_name")):
                lut[str(k).strip()] = str(r["mondo_name"])
    return lut


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmda", required=True)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--gold", default="data/gold_standards/legacy/gold_standard_v4_bigset.tsv")
    ap.add_argument("--out", default="results/tables/corroborate_pmda.csv")
    ap.add_argument("--full", type=float, default=0.60)
    ap.add_argument("--partial", type=float, default=0.30)
    a = ap.parse_args()

    LAB = load_labels()
    d = pd.read_csv(a.pmda, dtype=str)
    d["year"] = pd.to_numeric(d["year"], errors="coerce")
    d["_hay"] = (d.block.fillna("") + " " + d.brand_guess.fillna("")).map(norm)

    pairs = pd.read_csv(a.pairs, dtype=str)
    g = pd.read_csv(a.gold, sep="\t", dtype=str)
    ks = set()
    for _, r in g.iterrows():
        ks.add(f"{r.drug_id}|d|{str(r.disease_id).split(':')[-1].lstrip('0')}")
        ks.add(f"{r.drug_id}|n|{str(r.disease_name).strip().lower()}")

    rows = []
    for _, p in pairs.iterrows():
        if "PMDA" not in str(p.get("Agencies", "")).upper():
            continue
        names = {norm(x) for f in (p.get("Drug"), p.get("Active"))
                 for x in re.split(r"[+/]", str(f or "")) if len(norm(x)) > 3}
        hit = d[d._hay.apply(lambda h: any(n in h for n in names))]
        inwin = hit[(hit.year >= 2023) & (hit.year <= 2026)]

        cores = [disease_core(p.get("NewIndication", ""))]
        for c in ("DOID_ID", "MONDO_ID"):
            lab = LAB.get(str(p.get(c) or "").strip())
            if lab:
                cores.append(disease_core(lab))
        cores = [c for c in cores if c]

        status, best, ev = "no_pmda_record", 0.0, ""
        if len(hit) and not len(inwin):
            status = "drug_outside_window"
            ev = f"years present: {sorted(set(hit.year.dropna().astype(int)))[:8]}"
        elif len(inwin):
            status = "drug_year_only"
            for _, r in inwin.iterrows():
                for core in cores:
                    frac = sum(1 for w in core if w in r._hay) / len(core)
                    if frac > best:
                        best = frac
                        ev = (f"{str(r.approval_date)[:10]} | "
                              f"{'PartialChange' if r.is_partial_change=='True' else 'NewApproval'}"
                              f" | {str(r.block)[:170]}")
            if best >= a.full:
                status = "indication_match"
            elif best >= a.partial:
                status = "partial_match"
            if not ev:
                r = inwin.iloc[0]
                ev = f"{str(r.approval_date)[:10]} | {str(r.block)[:170]}"

        dn = str(p.get("DOID_ID") or "").split(":")[-1].lstrip("0")
        rows.append({
            "Pair_ID": p.get("Pair_ID"), "Drug": p.get("Drug"), "Active": p.get("Active"),
            "Agencies": p.get("Agencies"), "YearSet": p.get("YearSet"),
            "NewIndication": p.get("NewIndication"), "pmda_status": status,
            "token_coverage": round(best, 3), "pmda_evidence": ev,
            "in_116": (f"{p.get('DrugBank_ID')}|d|{dn}" in ks or
                       f"{p.get('DrugBank_ID')}|n|{str(p.get('NewIndication')).strip().lower()}" in ks),
        })

    res = pd.DataFrame(rows)
    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out, index=False)
    print(f"PMDA-attributed pairs checked: {len(res)}")
    for k, v in res.pmda_status.value_counts().items():
        print(f"  {v:>4}  {k}")
    print("\nby membership of the evaluated 116:")
    print(pd.crosstab(res.pmda_status, res.in_116).to_string())
    Path(str(out) + ".summary.json").write_text(json.dumps(
        {"n": int(len(res)),
         "counts": {k: int(v) for k, v in res.pmda_status.value_counts().items()},
         "source": "PMDA List of Approved Drugs, April 2004 to February 2026",
         "caveat": "English list ends February 2026; later approvals are not covered"},
        indent=2))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
