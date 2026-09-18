#!/usr/bin/env python3
"""
Build a complete census of FDA post-approval indication-expansion events
("repurposing events") from the official Drugs@FDA relational data files.

WHY THIS SOURCE
---------------
Earlier versions of this dataset were assembled by scraping the Drugs@FDA
monthly approvals HTML reports. Vetting (FDA_REPURPOSING_VETTING.md) showed
that route is accurate but NOT complete: recall was ~68% on a busy month
because long HTML tables get truncated, and BLA supplements carry a blank
Submission Classification in those reports and so were dropped wholesale.

drugsatfda.zip is the relational database behind Drugs@FDA. It is complete by
construction (no pagination, no rendering layer), it contains NDA and BLA
submissions in one table, and it ships FDA's own supplement action-type
taxonomy, which distinguishes a NEW INDICATION from other efficacy changes
(new dosing regimen, new patient population). That taxonomy is the operational
definition of a repurposing event used here, so the definition is FDA's, not
ours.

SCOPE BOUNDARY (state this in the paper)
----------------------------------------
Drugs@FDA covers CDER-regulated drugs and therapeutic biologics. Biologics
regulated by CBER (vaccines, blood products, most cell and gene therapies) are
NOT in this database and are therefore out of scope for this census. The census
is complete WITHIN that frame; the frame itself is stated, not implied.

USAGE
  python3 tools/build_fda_census.py profile --zip <path to drugsatfda.zip>
  python3 tools/build_fda_census.py build   --zip <path> --start 2023-01-01
"""
import argparse, hashlib, io, json, re, sys, zipfile
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

TABLES = [
    "ApplicationDocs", "ApplicationsDocsType_Lookup",
    "Applications", "Submissions", "SubmissionClass_Lookup", "Products",
    "Join_Submission_ActionTypes_Lookup", "ActionTypes_Lookup",
    "SubmissionPropertyType", "MarketingStatus", "MarketingStatus_Lookup",
]
# The FDA data-definitions page names this table in the singular
# ("Join_Submission_ActionType_Lookup") but ships it in the plural. Accept both
# so the script does not silently drop the action-type taxonomy.
TABLE_ALIASES = {
    "Join_Submission_ActionTypes_Lookup": ["Join_Submission_ActionType_Lookup"],
}


def sha256(path: Path) -> str:
    """Hash the source. For a directory, hash the sorted (name, content) stream so
    the digest is reproducible regardless of whether the user kept the zip."""
    h = hashlib.sha256()
    path = Path(path).expanduser()
    if path.is_dir():
        for f in sorted(p for p in path.rglob("*")
                        if p.is_file() and p.suffix.lower() == ".txt"):
            h.update(f.name.encode())
            with open(f, "rb") as fh:
                for chunk in iter(lambda: fh.read(1 << 20), b""):
                    h.update(chunk)
        return h.hexdigest()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def src_bytes(path: Path) -> int:
    path = Path(path).expanduser()
    if path.is_dir():
        return sum(f.stat().st_size for f in path.rglob("*")
                   if f.is_file() and f.suffix.lower() == ".txt")
    return path.stat().st_size


def _read_members(src: Path):
    """Yield (table stem, raw bytes) from either drugsatfda.zip or an unzipped dir."""
    src = Path(src).expanduser()
    if src.is_dir():
        found = {}
        for f in src.rglob("*"):
            if f.is_file() and f.suffix.lower() == ".txt":
                found.setdefault(f.stem, f)
        return {k: ("file", v) for k, v in found.items()}
    zf = zipfile.ZipFile(src)
    return {Path(n).stem: ("zip", (zf, n))
            for n in zf.namelist() if n.lower().endswith(".txt")}


def load_tables(zip_path: Path) -> dict:
    """Read the tab-delimited tables as all-string frames, from zip or directory."""
    out = {}
    members = _read_members(zip_path)
    missing = [t for t in TABLES
               if t not in members
               and not any(a in members for a in TABLE_ALIASES.get(t, []))]
    if missing:
        print(f"[warn] tables absent from source: {missing}", file=sys.stderr)
    if True:
        for t in TABLES:
            name = t
            if name not in members:
                name = next((a for a in TABLE_ALIASES.get(t, []) if a in members), None)
            if name is None:
                continue
            kind, ref = members[name]
            raw = ref.read_bytes() if kind == "file" else ref[0].read(ref[1])
            for enc in ("utf-8-sig", "latin-1"):
                try:
                    df = pd.read_csv(
                        io.BytesIO(raw), sep="\t", dtype=str, encoding=enc,
                        keep_default_na=False, na_values=[""], engine="python",
                        quoting=3, on_bad_lines="warn",
                    )
                    break
                except UnicodeDecodeError:
                    continue
            df.columns = [c.strip() for c in df.columns]
            for c in df.columns:
                df[c] = df[c].astype("string").str.strip()
            out[t] = df
    return out


def parse_dates(s: pd.Series) -> pd.Series:
    """Drugs@FDA ships dates in more than one shape across releases."""
    d = pd.to_datetime(s, errors="coerce", format="mixed", dayfirst=False)
    return d


def enrich(tabs: dict) -> pd.DataFrame:
    """One row per submission, joined to application, class and action type."""
    sub = tabs["Submissions"].copy()
    sub["SubmissionStatusDateParsed"] = parse_dates(sub["SubmissionStatusDate"])

    app = tabs["Applications"][["ApplNo", "ApplType", "SponsorName"]]
    sub = sub.merge(app, on="ApplNo", how="left", validate="many_to_one")

    cls = tabs["SubmissionClass_Lookup"].rename(
        columns={"SubmissionClassCodeID": "SubmissionClassCodeID"})
    sub = sub.merge(cls, on="SubmissionClassCodeID", how="left")

    # Action types are many-to-one per submission; keep them as a joined set.
    j = tabs.get("Join_Submission_ActionTypes_Lookup")
    at = tabs.get("ActionTypes_Lookup")
    if j is not None and at is not None:
        ja = j.merge(at, on="ActionTypes_LookupID", how="left")
        key = ["ApplNo", "SubmissionType", "SubmissionNo"]
        agg = (ja.groupby(key, dropna=False)
                 .agg(ActionTypeDescriptions=("ActionTypes_LookupDescription",
                                              lambda v: " | ".join(sorted({x for x in v if pd.notna(x)}))),
                      SupplCategoryL1=("SupplCategoryLevel1Code",
                                       lambda v: " | ".join(sorted({x for x in v if pd.notna(x)}))),
                      SupplCategoryL2=("SupplCategoryLevel2Code",
                                       lambda v: " | ".join(sorted({x for x in v if pd.notna(x)}))))
                 .reset_index())
        sub = sub.merge(agg, on=key, how="left")
    else:
        sub["ActionTypeDescriptions"] = pd.NA
        sub["SupplCategoryL1"] = pd.NA
        sub["SupplCategoryL2"] = pd.NA

    spt = tabs.get("SubmissionPropertyType")
    if spt is not None and "SubmissionPropertyTypeCode" in spt.columns:
        orph = (spt.assign(_o=spt["SubmissionPropertyTypeCode"].str.lower().eq("orphan"))
                   .groupby(["ApplNo", "SubmissionType", "SubmissionNo"], dropna=False)["_o"]
                   .any().reset_index().rename(columns={"_o": "Orphan"}))
        sub = sub.merge(orph, on=["ApplNo", "SubmissionType", "SubmissionNo"], how="left")
    else:
        sub["Orphan"] = pd.NA
    return sub


def window(sub: pd.DataFrame, start: str, end: str | None) -> pd.DataFrame:
    lo = pd.Timestamp(start)
    hi = pd.Timestamp(end) if end else pd.Timestamp(datetime.now().date())
    d = sub["SubmissionStatusDateParsed"]
    return sub[(d >= lo) & (d <= hi)]


def cmd_profile(args):
    zp = Path(args.zip).expanduser()
    tabs = load_tables(zp)
    print(f"# drugsatfda.zip  sha256={sha256(zp)}  bytes={src_bytes(zp)}")
    print(f"# profiled {datetime.now(timezone.utc).isoformat(timespec='seconds')}")
    for t, df in tabs.items():
        print(f"\n## {t}: {len(df):,} rows, columns={list(df.columns)}")

    sub = enrich(tabs)
    print(f"\n## Submissions joined: {len(sub):,} rows")
    print("\n### SubmissionStatusDate parse failures:",
          int(sub["SubmissionStatusDateParsed"].isna().sum()))
    print("\n### ApplType over all submissions")
    print(sub["ApplType"].value_counts(dropna=False).to_string())

    w = window(sub, args.start, args.end)
    print(f"\n## In window {args.start}..{args.end or 'today'}: {len(w):,} submissions")
    print("\n### ApplType in window")
    print(w["ApplType"].value_counts(dropna=False).to_string())
    print("\n### SubmissionStatus in window")
    print(w["SubmissionStatus"].value_counts(dropna=False).to_string())

    wa = w[w["SubmissionStatus"].eq("AP") & w["ApplType"].isin(["NDA", "BLA"])]
    print(f"\n## Approved NDA/BLA submissions in window: {len(wa):,}")
    print("\n### SubmissionClassCode (approved NDA/BLA in window)")
    print(wa["SubmissionClassCode"].value_counts(dropna=False).to_string())
    print("\n### Does class code survive for BLAs? (the monthly-report hole)")
    print(pd.crosstab(wa["ApplType"], wa["SubmissionClassCode"].fillna("<NULL>")).to_string())
    print("\n### SupplCategoryLevel1Code (approved NDA/BLA in window)")
    print(wa["SupplCategoryL1"].fillna("<NULL>").value_counts().to_string())
    print("\n### SupplCategoryLevel2Code (approved NDA/BLA in window)")
    print(wa["SupplCategoryL2"].fillna("<NULL>").value_counts().to_string())
    print("\n### ActionTypes_LookupDescription (approved NDA/BLA in window)")
    print(wa["ActionTypeDescriptions"].fillna("<NULL>").value_counts().head(40).to_string())
    print("\n### Full ActionTypes_Lookup vocabulary")
    print(tabs["ActionTypes_Lookup"].to_string())


# ---------------------------------------------------------------------------
# Inclusion rule
# ---------------------------------------------------------------------------
# A repurposing event is defined as: a molecule that already held FDA approval
# gains approval for an indication it did not previously hold. The rule below
# is built entirely out of FDA's own submission taxonomy, so it can be audited
# without trusting any judgement of ours. Three tiers, because the taxonomy
# routes the same underlying event through three different codes.
#
#   A  SUPPL tagged "Efficacy-New Indication" (SupplCategoryLevel2Code
#      = INDICATION). The canonical case, and unambiguous.
#
#   B  SUPPL tagged "Efficacy-Accelerated Approval". Inspection of these rows
#      (KEYTRUDA, OPDIVO, YERVOY, BRAFTOVI and similar) shows they are new
#      indications granted under subpart H/E, which FDA tags by the approval
#      PATHWAY rather than by the change. An INDICATION-only rule silently
#      drops them, and they are disproportionately oncology, so dropping them
#      would bias the set by therapeutic area. Flagged for adjudication rather
#      than assumed.
#
#   C  ORIG applications with SubmissionClassCode TYPE 9, TYPE 9- BLA or
#      TYPE 10, which FDA defines as "new indication submitted as a distinct
#      NDA". This is repurposing carried out through a separate application
#      rather than a supplement (ZEPBOUND vs MOUNJARO; YEZTUGO vs SUNLENCA).
#      A supplement-only rule misses this route entirely. Also flagged, since
#      a minority are reformulations rather than new indications.
#
# EXCLUDED, with the reason stated rather than left implicit:
#   Efficacy-New Patient Population  - same disease, wider age band
#   Efficacy-New Dosing Regimen      - same disease
#   Efficacy-New Route Of Administration, Rx To OTC Switch - same disease
#   Efficacy-Labeling Change With Clinical Data - new evidence, existing claim
#   Efficacy-Accelerated Approval Confirmatory Study - converts an existing
#       accelerated approval to full; the indication already existed
#   All LABELING, MANUF (CMC), REMS, BIOEQUIV categories
TIER_A = "Efficacy-New Indication"
TIER_B = "Efficacy-Accelerated Approval"
TIER_C_CLASSES = {"TYPE 9", "TYPE 9- BLA", "TYPE 10"}


def assign_tier(row) -> str:
    desc = (row.get("ActionTypeDescriptions") or "")
    stype = (row.get("SubmissionType") or "").upper()
    cls = (row.get("SubmissionClassCode") or "")
    if stype == "SUPPL" and TIER_A in desc:
        return "A_new_indication"
    if stype == "SUPPL" and TIER_B in desc and "Confirmatory" not in desc:
        return "B_accelerated_approval"
    if stype == "ORIG" and cls in TIER_C_CLASSES:
        return "C_new_indication_distinct_application"
    return ""


def attach_documents(ev: pd.DataFrame, tabs: dict) -> pd.DataFrame:
    """Attach the FDA approval letter and label URLs for each submission, so a
    reviewer can verify any row against its own primary document."""
    docs = tabs.get("ApplicationDocs")
    if docs is None:
        ev["approval_letter_url"] = pd.NA
        ev["label_url"] = pd.NA
        return ev
    d = docs.copy()
    d["SubmissionType"] = d["SubmissionType"].str.strip()
    d["SubmissionNo"] = d["SubmissionNo"].str.strip()
    key = ["ApplNo", "SubmissionType", "SubmissionNo"]
    for tid, col in (("1", "approval_letter_url"), ("2", "label_url")):
        sel = (d[d["ApplicationDocsTypeID"].eq(tid)]
               .sort_values("ApplicationDocsDate")
               .groupby(key, dropna=False)["ApplicationDocsURL"].last()
               .reset_index().rename(columns={"ApplicationDocsURL": col}))
        ev = ev.merge(sel, on=key, how="left")
    return ev


def cmd_build(args):
    zp = Path(args.zip).expanduser()
    digest = sha256(zp)
    tabs = load_tables(zp)
    sub = enrich(tabs)

    funnel = []
    def step(label, df):
        funnel.append((label, len(df)))
        return df

    step("all submissions in source", sub)
    w = step(f"status date in [{args.start}, {args.end or 'today'}]",
             window(sub, args.start, args.end))
    w = step("SubmissionStatus == AP (approved)", w[w["SubmissionStatus"].eq("AP")])
    w = step("ApplType in (NDA, BLA)", w[w["ApplType"].isin(["NDA", "BLA"])])

    prod = (tabs["Products"].groupby("ApplNo")
            .agg(DrugNames=("DrugName", lambda v: " | ".join(sorted({x for x in v if pd.notna(x)}))),
                 ActiveIngredients=("ActiveIngredient",
                                    lambda v: " | ".join(sorted({x for x in v if pd.notna(x)}))))
            .reset_index())

    w = w.merge(prod, on="ApplNo", how="left")
    w = w.assign(
        approval_date=w["SubmissionStatusDateParsed"].dt.strftime("%Y-%m-%d"),
        application=w["ApplType"] + w["ApplNo"],
        supplement=w["SubmissionType"].str.strip() + "-" + w["SubmissionNo"],
    )
    w["SubmissionType"] = w["SubmissionType"].str.strip()
    w = attach_documents(w, tabs)

    base_cols = ["approval_date", "application", "ApplType", "ApplNo", "SubmissionType",
                 "SubmissionNo", "supplement", "ReviewPriority", "SubmissionClassCode",
                 "SupplCategoryL1", "SupplCategoryL2", "ActionTypeDescriptions",
                 "Orphan", "SponsorName", "DrugNames", "ActiveIngredients",
                 "approval_letter_url", "label_url"]

    # The frame is the denominator: every approved NDA/BLA action in the window.
    frame_path = Path(str(args.out).replace(".csv", "_frame.csv")).expanduser()
    frame_path.parent.mkdir(parents=True, exist_ok=True)
    w[[c for c in base_cols if c in w.columns]].sort_values(
        ["approval_date", "application", "SubmissionNo"]).to_csv(frame_path, index=False)
    print(f"wrote {frame_path}  ({len(w):,} approved NDA/BLA actions = the denominator)")

    w["tier"] = w.apply(assign_tier, axis=1)
    ev = step("matches the tiered indication rule", w[w["tier"].ne("")]).copy()
    ev["adjudication_required"] = ev["tier"].ne("A_new_indication")

    cols = ["approval_date", "application", "tier", "adjudication_required"] + \
           [c for c in base_cols if c not in ("approval_date", "application")]
    ev = ev[[c for c in cols if c in ev.columns]].sort_values(
        ["approval_date", "application", "SubmissionNo"]).reset_index(drop=True)

    out = Path(args.out).expanduser()
    ev.to_csv(out, index=False)

    by_tier = ev["tier"].value_counts().to_dict()
    meta = {
        "built_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source": "Drugs@FDA Data Files (drugsatfda.zip, unzipped)",
        "source_url": "https://www.fda.gov/drugs/drug-approvals-and-databases/drugsfda-data-files",
        "source_sha256": digest,
        "source_bytes": src_bytes(zp),
        "window_start": args.start,
        "window_end": args.end or str(datetime.now().date()),
        "rule": {
            "A_new_indication": "SUPPL with action type 'Efficacy-New Indication'",
            "B_accelerated_approval": "SUPPL with action type 'Efficacy-Accelerated Approval', excluding confirmatory studies",
            "C_new_indication_distinct_application": "ORIG with SubmissionClassCode in TYPE 9 / TYPE 9- BLA / TYPE 10",
            "excluded": ["Efficacy-New Patient Population", "Efficacy-New Dosing Regimen",
                         "Efficacy-New Route Of Administration", "Efficacy-Rx To OTC Switch",
                         "Efficacy-Labeling Change With Clinical Data",
                         "Efficacy-Accelerated Approval Confirmatory Study",
                         "all LABELING / MANUF (CMC) / REMS / BIOEQUIV"],
        },
        "scope_note": ("CDER-regulated drugs and therapeutic biologics only; CBER "
                       "biologics (vaccines, blood products, most cell and gene "
                       "therapies) are not in Drugs@FDA and are out of frame."),
        "known_join_gap": ("Submissions whose ApplNo is absent from Applications.txt "
                           "carry a null ApplType; all observed cases are pre-1982 "
                           "antibiotic or legacy ANDA numbers, excluded by the "
                           "NDA/BLA filter regardless."),
        "funnel": [{"step": s, "rows": n} for s, n in funnel],
        "n_events": int(len(ev)),
        "n_by_tier": {k: int(v) for k, v in by_tier.items()},
        "n_applications": int(ev["application"].nunique()),
        "n_ingredients": int(ev["ActiveIngredients"].nunique()),
        "n_adjudication_required": int(ev["adjudication_required"].sum()),
    }
    Path(str(out) + ".meta.json").write_text(json.dumps(meta, indent=2))

    print("\nFUNNEL")
    for s, n in funnel:
        print(f"  {n:>9,}  {s}")
    print("\nBY TIER")
    for k, v in sorted(by_tier.items()):
        print(f"  {v:>9,}  {k}")
    print(f"\nwrote {out}  ({len(ev):,} events, {ev['application'].nunique():,} applications, "
          f"{ev['ActiveIngredients'].nunique():,} distinct ingredient strings)")
    print(f"wrote {out}.meta.json")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="cmd", required=True)

    pp = sp.add_parser("profile", help="describe the file and its vocabularies")
    pp.add_argument("--zip", required=True, help="drugsatfda.zip or the unzipped directory")
    pp.add_argument("--start", default="2023-01-01")
    pp.add_argument("--end", default=None)
    pp.set_defaults(func=cmd_profile)

    pb = sp.add_parser("build", help="emit the event census")
    pb.add_argument("--zip", required=True, help="drugsatfda.zip or the unzipped directory")
    pb.add_argument("--start", default="2023-01-01")
    pb.add_argument("--end", default=None)
    pb.add_argument("--out", default="data/gold_standards/legacy/fda_indication_expansion_census.csv")
    pb.set_defaults(func=cmd_build)

    a = p.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
