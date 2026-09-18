#!/usr/bin/env python3
"""
Parse PMDA's "List of Approved Drugs" PDF into one record per approval.

WHY THIS SOURCE IS THE STRONGEST OF THE THREE
---------------------------------------------
Each PMDA record carries four things the other registers do not all provide:
an approval date, the brand and active ingredient, an explicit
"New Approval / Partial Change" flag, and a Notes sentence naming the
indication ("A drug with a new indication for the treatment of complement 3
glomerulopathy"). "Change" is the Japanese partial-change approval, which is
exactly the regulatory event corresponding to a new indication on an existing
drug, so repurposing is identifiable directly rather than inferred.

KNOWN LIMITS, STATED UP FRONT
-----------------------------
  * The English list ends February 2026, so the tail of a 2023-2026 window is
    not covered by this source at all.
  * It is a PDF, not a table. Records are reconstructed from laid-out text, so
    the parse is verified by re-deriving counts and by keeping the raw block for
    every record as evidence rather than trusting the extraction.

USAGE
  python3 tools/parse_pmda.py --pdf <pdf> --out data/gold_standards/legacy/pmda_approvals.csv
"""
import argparse, re, subprocess, sys
from pathlib import Path
import pandas as pd

# The document mixes two date styles: recent monthly tables spell the month out
# ("May 19, 2025") while the fiscal-year tables abbreviate it ("Jun. 26, 2023").
# Matching only the full names silently drops every fiscal-year section, which is
# the bulk of the file and all of 2023 and 2024.
MONTH = (r"(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Sept|Oct|Nov|Dec)"
         r"(?:uary|ruary|ch|il|e|y|ust|tember|ober|ember)?\.?")
HEADER = re.compile(rf"New Drugs Approved in (?:FY\s*)?({MONTH})?\s*(\d{{4}})")
# a record line carries a full date and then the record's index number
REC = re.compile(rf"^\s*(?P<cat>\S.*?)\s{{2,}}(?P<date>{MONTH}\s+\d{{1,2}},\s+\d{{4}})"
                 rf"\s+(?P<no>\d+)\s+(?P<rest>.*)$")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdf", required=True)
    ap.add_argument("--out", default="data/gold_standards/legacy/pmda_approvals.csv")
    a = ap.parse_args()

    txt = subprocess.run(["pdftotext", "-layout", a.pdf, "-"],
                         capture_output=True, text=True, check=True).stdout
    lines = txt.split("\n")

    recs, cur, month = [], None, (None, None)
    for ln in lines:
        h = HEADER.search(ln)
        if h:
            month = (h.group(1), h.group(2))
            continue
        m = REC.match(ln)
        if m:
            if cur:
                recs.append(cur)
            cur = {"review_category": m.group("cat").strip(),
                   "approval_date_raw": m.group("date"),
                   "record_no": m.group("no"),
                   "header_month": f"{month[0]} {month[1]}" if month[0] else None,
                   "block": m.group("rest").strip()}
        elif cur is not None and ln.strip():
            cur["block"] += " " + ln.strip()
    if cur:
        recs.append(cur)

    df = pd.DataFrame(recs)
    if not len(df):
        sys.exit("no records parsed; the PDF layout has changed")

    def parse_date(v):
        v = re.sub(r"\.", "", str(v))
        for fmt in ("%B %d, %Y", "%b %d, %Y"):
            try:
                return pd.to_datetime(v, format=fmt)
            except (ValueError, TypeError):
                continue
        return pd.NaT

    df["approval_date"] = df.approval_date_raw.map(parse_date)
    df["year"] = df.approval_date.dt.year
    # "Change" marks a partial-change approval: a new indication on a drug that
    # is already approved. "Approval" marks a wholly new product.
    df["is_partial_change"] = df.block.str.contains(r"\bChange\b", na=False)
    df["is_new_approval"] = df.block.str.contains(r"\bApproval\b", na=False)
    df["mentions_new_indication"] = df.block.str.contains(r"new indication", case=False, na=False)
    # brand name is the leading text of the record body, before the applicant
    df["brand_guess"] = (df.block.str.split(r"\(", regex=True).str[0]
                           .str.replace(r"\s{2,}", " ", regex=True).str.strip().str[:90])

    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)

    print(f"records parsed          : {len(df):,}")
    print(f"date parsed             : {int(df.approval_date.notna().sum()):,}")
    print(f"date range              : {df.approval_date.min()} .. {df.approval_date.max()}")
    print(f"partial-change records  : {int(df.is_partial_change.sum()):,}")
    print(f"mention 'new indication': {int(df.mentions_new_indication.sum()):,}")
    w = df[(df.year >= 2023) & (df.year <= 2026)]
    print(f"\nin 2023-2026            : {len(w):,} records, "
          f"{int(w.mentions_new_indication.sum()):,} mentioning a new indication")
    print(w.groupby("year").size().to_string())
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
