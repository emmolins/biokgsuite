#!/usr/bin/env python3
"""
Attach a vettable citation to every pair in the cleaned evaluation set.

The point is that a reviewer can check any row without rerunning anything and
without trusting us: each row carries the regulator's own identifier for the
approval and a URL that resolves to the regulator's page or document.

Evidence per agency:
  FDA   application and supplement number from the Drugs@FDA relational data
        files, plus the URL of the approval letter for that exact submission
  EMA   EMA product number (EMEA/H/C/......) and the EPAR page URL, from the
        EMA medicines register
  PMDA  the dated record in PMDA's List of Approved Drugs, with its record
        number, plus the URL of that published list

Verification strength is recorded per row rather than implied, because it is
not uniform:
  disease_level   the register names the disease, so drug AND disease are checked
  drug_window     the register confirms an approval for that drug in the window
                  but does not name the indication (Drugs@FDA does not record
                  which indication a supplement added)
  unverified      no register check was possible (MHRA/Orbis has no equivalent
                  machine-readable register here)
"""
import argparse, re
from pathlib import Path
import pandas as pd


def dn(x):
    return str(x).split(":")[-1].lstrip("0")


def norm(s):
    return re.sub(r"\s+", " ", re.sub(r"[^a-z0-9 ]+", " ", str(s).lower())).strip()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gold", default="data/gold_standards/legacy/gold_standard_v5_strict.tsv")
    ap.add_argument("--pairs", default="data/gold_standards/legacy/pairs_with_ids_v2.csv")
    ap.add_argument("--ema-xlsx", default="ema/medicines.xlsx")
    ap.add_argument("--pmda", default="data/gold_standards/legacy/pmda_approvals.csv")
    ap.add_argument("--fda-census", default="data/gold_standards/legacy/fda_indication_expansion_census.csv")
    ap.add_argument("--fda-corr", default="results/tables/corroborate_253.csv")
    ap.add_argument("--ema-corr", default="results/tables/corroborate_ema.csv")
    ap.add_argument("--pmda-corr", default="results/tables/corroborate_pmda.csv")
    ap.add_argument("--out", default="data/gold_standards/legacy/gold_standard_v5_strict_evidence.csv")
    a = ap.parse_args()

    g = pd.read_csv(a.gold, sep="\t", dtype=str)
    p = pd.read_csv(a.pairs, dtype=str)
    g["_k"] = g.drug_id.astype(str) + "|" + g.disease_id.map(dn)
    p["_k"] = p.DrugBank_ID.astype(str) + "|" + p.DOID_ID.map(dn)
    pm = p.drop_duplicates("_k").set_index("_k")
    for col in ("Pair_ID", "Active", "MONDO_ID", "MeSH_ID", "OrigYear", "TA_simplified"):
        g[col] = g._k.map(pm[col])

    # ---- EMA: product number and EPAR URL -------------------------------
    em = pd.read_excel(a.ema_xlsx, header=8, dtype=str).dropna(how="all")
    em.columns = [str(c).strip() for c in em.columns]
    em = em[em["Category"].fillna("").str.lower().eq("human")]
    eidx = {}
    for _, r in em.iterrows():
        nm = str(r.get("Name of medicine") or "")
        keys = {norm(nm), norm(re.split(r"\(previously", nm, flags=re.I)[0])}
        for c in ("International non-proprietary name (INN) / common name", "Active substance"):
            if pd.notna(r.get(c)):
                for part in re.split(r"[/,;]| and ", str(r[c])):
                    keys.add(norm(part))
        for k in keys:
            if len(k) > 2:
                eidx.setdefault(k, r)

    ema_status = pd.read_csv(a.ema_corr, dtype=str).set_index("Pair_ID")

    # ---- FDA: application, supplement, approval letter -------------------
    fc = pd.read_csv(a.fda_corr, dtype=str).set_index("Pair_ID")
    cen = pd.read_csv(a.fda_census, dtype=str)
    cen["_ref"] = cen.application.astype(str) + " " + cen.supplement.astype(str)
    letters = cen.drop_duplicates("_ref").set_index("_ref")

    # ---- PMDA -------------------------------------------------------------
    pc = pd.read_csv(a.pmda_corr, dtype=str).set_index("Pair_ID")

    PMDA_URL = "https://www.pmda.go.jp/files/000281190.pdf"
    rows = []
    for _, r in g.iterrows():
        pid = str(r.Pair_ID)
        ags = str(r.source_agencies or "")
        ev_src = ev_ref = ev_url = ev_date = status = note = ""

        if "FDA" in ags.upper() and pid in fc.index:
            f = fc.loc[pid]
            # Accept the weaker FDA corroboration grades too. The evidence is a
            # real FDA action either way; only the action TYPE differs, and that
            # is recorded in the note rather than used to discard the citation.
            if str(f.corroboration).startswith("corroborated") or \
               str(f.corroboration) in ("efficacy_event_outside_tiers",
                                        "fda_action_but_not_efficacy"):
                first = str(f.fda_evidence).split(";")[0].strip()
                m = re.match(r"(\d{4}-\d{2}-\d{2})\s+(\S+)\s+(\S+)", first)
                if m:
                    ev_date, appl, sup = m.groups()
                    ev_src, ev_ref = "Drugs@FDA", f"{appl} {sup}"
                    row = letters.loc[f"{appl} {sup}"] if f"{appl} {sup}" in letters.index else None
                    ev_url = str(row.approval_letter_url) if row is not None and pd.notna(
                        row.get("approval_letter_url")) else \
                        f"https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=overview.process&ApplNo={re.sub(r'[A-Z]','',appl)}"
                    status = "drug_window"
                    if not str(f.corroboration).startswith("corroborated"):
                        note = (f"FDA action of a different type "
                                f"({f.corroboration}); the approval is real but "
                                f"is not tagged Efficacy-New Indication")

        # Combination products carry comma-separated DrugBank IDs and so do not
        # map to a single candidate row. Look them up in the census by ingredient.
        if not ev_src and "FDA" in ags.upper():
            toks = [t for t in re.split(r"[+/]", str(r.drug_name)) if len(t) > 4]
            for t in toks:
                h = cen[cen.ActiveIngredients.fillna("").str.contains(
                    t.strip(), case=False, regex=False)]
                if len(h):
                    h = h.sort_values("approval_date")
                    row = h.iloc[0]
                    ev_src, ev_ref = "Drugs@FDA", row["_ref"]
                    ev_date = str(row.approval_date)
                    ev_url = (str(row.approval_letter_url)
                              if pd.notna(row.get("approval_letter_url")) else "")
                    status = "drug_window"
                    note = ("combination product; citation is the FDA action for "
                            f"the component '{t.strip()}'")
                    break

        if not ev_src and "EMA" in ags.upper():
            e = eidx.get(norm(r.drug_name))
            if e is None:
                e = eidx.get(norm(r.Active))
            if e is not None:
                ev_src = "EMA EPAR"
                ev_ref = str(e.get("EMA product number") or "")
                ev_url = str(e.get("Medicine URL") or "")
                ev_date = str(e.get("European Commission decision date") or "")
                st = str(ema_status.loc[pid].ema_status) if pid in ema_status.index else ""
                status = "disease_level" if st == "indication_match" else "drug_window"
                note = "EMA register reflects the current label, so the date is the "\
                       "latest decision, not necessarily this indication's"

        if not ev_src and "PMDA" in ags.upper():
            ev_src, ev_url, ev_ref = "PMDA List of Approved Drugs", PMDA_URL, ""
            if pid not in pc.index:
                # unmapped row (combination, or no candidate match): look the
                # drug up directly in the parsed PMDA records
                pm_rec = pd.read_csv(a.pmda, dtype=str)
                pm_rec["year"] = pd.to_numeric(pm_rec.year, errors="coerce")
                toks = [t for t in re.split(r"[+/]", str(r.drug_name)) if len(t) > 4]
                dcore = [w for w in norm(r.disease_name).split() if len(w) > 4]
                best = None
                for t in toks:
                    h = pm_rec[pm_rec.block.fillna("").str.contains(t.strip(), case=False, regex=False)
                               & pm_rec.year.between(2023, 2026)]
                    for _, hr in h.iterrows():
                        hay = norm(hr.block)
                        frac = (sum(1 for w in dcore if w in hay) / len(dcore)) if dcore else 0
                        if best is None or frac > best[0]:
                            best = (frac, hr)
                if best is not None:
                    frac, hr = best
                    ev_date = str(hr.approval_date)[:10]
                    ev_ref = f"record dated {ev_date}"
                    status = "disease_level" if frac >= 0.6 else "drug_window"
            if pid in pc.index:
                pe = str(pc.loc[pid].pmda_evidence or "")
                d = re.match(r"(\d{4}-\d{2}-\d{2})", pe)
                ev_date = d.group(1) if d else ""
                ev_ref = f"record dated {ev_date}" if ev_date else ""
                status = ("disease_level"
                          if str(pc.loc[pid].pmda_status) == "indication_match"
                          else "drug_window")

        if not ev_src:
            status, note = "unverified", "no machine-readable register checked for this agency"

        rows.append({
            "pair_id": pid, "drugbank_id": r.drug_id, "drug_name": r.drug_name,
            "active_ingredient": r.Active, "disease_name": r.disease_name,
            "doid": r.disease_id, "mondo_id": r.MONDO_ID, "mesh_id": r.MeSH_ID,
            "therapeutic_area": r.TA_simplified,
            "drug_first_approved": r.OrigYear, "approval_agencies": ags,
            "window": r.year_set, "evidence_source": ev_src,
            "evidence_ref": ev_ref, "evidence_date": ev_date,
            "evidence_url": ev_url, "verification": status, "note": note,
        })

    out = pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, index=False)
    print(f"{len(out)} pairs\n")
    print("evidence source:"); print(out.evidence_source.value_counts().to_string())
    print("\nverification strength:"); print(out.verification.value_counts().to_string())
    print(f"\nrows with a resolvable URL: {int((out.evidence_url.fillna('') != '').sum())}/{len(out)}")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
