#!/usr/bin/env python3
"""
Corroborate the EMA-attributed candidate pairs against EMA's register of
centrally authorised medicines.

WHICH EMA SOURCE, AND WHY NOT THE OTHER ONE
-------------------------------------------
EMA publishes two relevant tables. The "post-authorisation procedures" table
sounds like the right one, being the EU analogue of an FDA efficacy supplement,
but it is a 140-row snapshot of procedures currently in flight. It carries no
historical rows and no procedure scope text, so it cannot corroborate anything
approved in 2023 to 2025. It was inspected and rejected, not overlooked.

The medicines table is the usable one: 2,339 human medicines, each carrying the
full authorised `Therapeutic indication` text and a curated
`Therapeutic area (MeSH)` list.

WHAT THIS ESTABLISHES, AND WHAT IT DOES NOT
-------------------------------------------
ESTABLISHES: that the drug is centrally authorised in the EU and that the
claimed disease appears in its authorised indication. That is a genuine
drug-disease check, and it is STRONGER than the FDA corroboration, which could
only match drug plus year because Drugs@FDA does not record which indication a
supplement added.

DOES NOT ESTABLISH: the year the indication was added. This table reflects the
current label, not variation history. EMA keeps that history only in per-product
"procedural steps" documents, one per medicine. So a pair confirmed here is
confirmed as real; its placement inside the 2023 to 2026 window rests on the
original assembly and is reported as an open item rather than as verified.

Grades:
  indication_match       drug found AND disease present in indication or MeSH area
  partial_match          drug found, disease partially present, needs adjudication
  drug_found_no_match    drug is EU-authorised but the disease is not on its label
  drug_not_in_register   no centrally authorised medicine of that name

USAGE
  python3 tools/corroborate_ema.py --xlsx ema/medicines.xlsx \
      --pairs data/gold_standards/legacy/pairs_with_ids_v2.csv
"""
import argparse, json, re
from pathlib import Path
import pandas as pd

# Modifiers that describe HOW a drug is used rather than WHICH disease it treats.
# Left in, they inflate overlap and make every oncology pair look corroborated.
MODIFIERS = {
    "1l", "2l", "3l", "first", "second", "third", "line", "adjuvant",
    "neoadjuvant", "maintenance", "advanced", "metastatic", "unresectable",
    "locally", "relapsed", "refractory", "recurrent", "previously", "treated",
    "untreated", "adult", "adults", "paediatric", "pediatric", "children",
    "monotherapy", "combination", "with", "and", "or", "the", "of", "in", "for",
    "patients", "patient", "treatment", "therapy", "expansion", "extension",
    "new", "indication", "years", "aged", "age", "stage", "including", "plus",
    "positive", "negative", "high", "low", "risk", "severe", "moderate", "mild",
    "chronic", "acute", "primary", "secondary", "prophylaxis", "prevention",
}


def norm(s):
    return re.sub(r"\s+", " ", re.sub(r"[^a-z0-9 ]+", " ", str(s).lower())).strip()


def load_disease_labels(path="data/gold_standards/bridges/do_diseases.csv"):
    """DOID / MONDO to the ontology's spelled-out label.

    The candidate file writes indications the way a clinician would ("Ph+ CML-CP",
    "chronic GVHD", "XLH infants <1yr"). EMA writes them the way MeSH does
    ("Leukemia, Myelogenous, Chronic, BCR-ABL Positive"). Matching the raw strings
    fails on the abbreviation, not on the fact. The pairs already carry ontology
    IDs, so the ontology label is a spelled-out form of the same disease and is a
    fairer thing to match on than an ad hoc abbreviation dictionary.
    """
    import os
    if not os.path.exists(path):
        return {}
    d = pd.read_csv(path, dtype=str)
    lut = {}
    for _, r in d.iterrows():
        for key in (r.get("doid"), r.get("mondo_id")):
            if pd.notna(key) and pd.notna(r.get("mondo_name")):
                lut[str(key).strip()] = str(r["mondo_name"])
    return lut


def disease_core(text):
    """Content tokens naming the disease, with usage modifiers removed."""
    t = str(text)
    t = re.sub(r"\([^)]*\)", " ", t)          # drop parentheticals: (+trastuzumab+chemo)
    t = re.sub(r"\+\s*\w+", " ", t)           # drop combination partners
    toks = norm(t).split()
    return [w for w in toks if w not in MODIFIERS and len(w) > 2]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx", required=True)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--out", default="results/tables/corroborate_ema.csv")
    ap.add_argument("--gold", default="data/gold_standards/legacy/gold_standard_v4_bigset.tsv")
    ap.add_argument("--full", type=float, default=0.60,
                    help="fraction of disease tokens that must be present")
    ap.add_argument("--partial", type=float, default=0.30)
    a = ap.parse_args()

    LABELS = load_disease_labels()
    print(f"ontology labels loaded: {len(LABELS):,}")
    df = pd.read_excel(a.xlsx, header=8, dtype=str).dropna(how="all")
    df.columns = [str(c).strip() for c in df.columns]
    df = df[df["Category"].fillna("").str.lower().eq("human")]
    print(f"EMA human medicines in register: {len(df):,}")

    # index by brand name and by INN; a pair may carry either.
    # EMA renames products in place ("Onivyde pegylated liposomal (previously
    # Onivyde)", "Spikevax (previously COVID-19 Vaccine Moderna)"), so the bare
    # leading name is indexed too or those products look absent from the register.
    idx = {}
    for i, r in df.iterrows():
        raw_name = str(r.get("Name of medicine") or "")
        keys = {norm(raw_name), norm(re.split(r"\(previously", raw_name, flags=re.I)[0])}
        for col in ("International non-proprietary name (INN) / common name",
                    "Active substance"):
            v = r.get(col)
            if pd.notna(v):
                for part in re.split(r"[/,;]| and ", str(v)):
                    if len(norm(part)) > 2:
                        keys.add(norm(part))
        for k in keys:
            if k:
                idx.setdefault(k, []).append(i)

    pairs = pd.read_csv(a.pairs, dtype=str)

    # Which candidates actually survive to the evaluated set: a failure on a row
    # that filter 1 already removed does not affect any reported result.
    g = pd.read_csv(a.gold, sep="\t", dtype=str) if a.gold else None
    ks = set()
    if g is not None:
        for _, r in g.iterrows():
            dn = str(r.disease_id).split(":")[-1].lstrip("0")
            ks.add(f"{r.drug_id}|d|{dn}")
            ks.add(f"{r.drug_id}|n|{str(r.disease_name).strip().lower()}")

    def in116(p):
        dn = str(p.get("DOID_ID") or "").split(":")[-1].lstrip("0")
        return (f"{p.get('DrugBank_ID')}|d|{dn}" in ks or
                f"{p.get('DrugBank_ID')}|n|{str(p.get('NewIndication')).strip().lower()}" in ks)

    PLACEHOLDER = ("tbd", "population expansion", "new indication (")
    rows = []
    for _, p in pairs.iterrows():
        if "EMA" not in str(p.get("Agencies", "")).upper():
            continue
        # Combination entries are written "Opdivo+Yervoy"; each component is its
        # own authorised medicine, so look each up rather than the joined string.
        keys = set()
        for field in (p.get("Drug"), p.get("Active")):
            for part in re.split(r"[+/]", str(field or "")):
                k = norm(part)
                if len(k) > 2:
                    keys.add(k)
        hits = sorted({i for k in keys for i in idx.get(k, [])})
        if not hits:
            # EMA also renames to a longer form ("Onivyde pegylated liposomal"),
            # so fall back to prefix matching before declaring a drug absent.
            for k in keys:
                for rk, rows_ in idx.items():
                    if rk.startswith(k + " "):
                        hits.extend(rows_)
            hits = sorted(set(hits))

        # Try the free-text indication AND the ontology label for the same disease
        cores = [disease_core(p.get("NewIndication", ""))]
        for idcol in ("DOID_ID", "MONDO_ID"):
            lab = LABELS.get(str(p.get(idcol) or "").strip())
            if lab:
                cores.append(disease_core(lab))
        cores = [c for c in cores if c]
        status, best, ev, where = "drug_not_in_register", 0.0, "", ""
        if hits:
            status = "drug_found_no_match"
            for i in hits:
                r = df.loc[i]
                mesh = norm(r.get("Therapeutic area (MeSH)"))
                ind = norm(r.get("Therapeutic indication"))
                for label, hay in (("MeSH area", mesh), ("indication text", ind)):
                    if not hay:
                        continue
                    for ci, core in enumerate(cores):
                        frac = sum(1 for w in core if w in hay) / len(core)
                        if frac > best:
                            best = frac
                            where = label + ("" if ci == 0 else " via ontology label")
                        ev = (f"{r.get('Name of medicine')} "
                              f"[{r.get('Medicine status')}] "
                              f"{str(r.get('Therapeutic area (MeSH)'))[:110]}")
            if best >= a.full:
                status = "indication_match"
            elif best >= a.partial:
                status = "partial_match"

        rows.append({
            "Pair_ID": p.get("Pair_ID"), "Drug": p.get("Drug"),
            "Active": p.get("Active"), "Agencies": p.get("Agencies"),
            "YearSet": p.get("YearSet"), "NewIndication": p.get("NewIndication"),
            "ema_status": status, "token_coverage": round(best, 3),
            "matched_on": where, "ema_evidence": ev,
            "in_116": in116(p),
            "placeholder_text": any(t in str(p.get("NewIndication", "")).lower()
                                    for t in PLACEHOLDER),
        })

    res = pd.DataFrame(rows)
    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(out, index=False)

    print(f"\nEMA-attributed pairs checked: {len(res)}")
    for k, v in res.ema_status.value_counts().items():
        print(f"  {v:>4}  {k}")
    ok = int(res.ema_status.isin(["indication_match"]).sum())
    print(f"\nfully corroborated (drug + disease): {ok}/{len(res)} "
          f"({100*ok/max(len(res),1):.1f}%)")
    Path(str(out) + ".summary.json").write_text(json.dumps(
        {"n": int(len(res)), "n_indication_match": ok,
         "counts": {k: int(v) for k, v in res.ema_status.value_counts().items()},
         "thresholds": {"full": a.full, "partial": a.partial},
         "source": "EMA medicines register (medicines-output-medicines-report_en.xlsx)",
         "caveat": "confirms the drug-disease authorisation, not the year it was added"},
        indent=2))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
