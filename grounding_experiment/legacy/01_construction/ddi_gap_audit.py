"""Definitive drug-drug gap audit (v2): replicates what the loader does.

Step 1: stream nodes.tsv, build map  curie -> drugbank_id  (from primary id
        OR equivalent_identifiers). Same regex the loader uses.
Step 2: stream edges.tsv, look up subject and object. If BOTH map to a
        DrugBank ID and the predicate is in our DDI set, record the pair.
Step 3: also record pairs under ANY predicate (upper bound for Matrix's
        DDI ceiling under maximally permissive selection).
Step 4: compare against gold drugbank_ddi.csv.

This tells us whether the 33% gap is in the data (Matrix didn't ingest those
DDI edges) or in our pipeline (we filtered them out somewhere).
"""
from __future__ import annotations
import json
import re
import time
from collections import Counter
from pathlib import Path

import pandas as pd

NODES = Path("data/matrix/nodes.tsv")
EDGES = Path("data/matrix/edges.tsv")
GOLD  = Path("data/gold_standards/bridges/drugbank_ddi.csv")
OUT   = Path("results/ddi_gap_audit_v2.json")
PUBCHEM_BRIDGE = Path("data/openbilink/pubchem_to_drugbank.csv")
# Optional comprehensive multi-namespace DrugBank xref (UNII / RxCUI / ATC /
# KEGG.DRUG / ChEBI / ChEMBL / CAS / etc.) built by scripts/build_drugbank_xref.py
DRUGBANK_XREF = Path("data/gold_standards/bridges/drugbank_xref.csv")

DDI_PREDS    = {"interacts_with", "directly_physically_interacts_with"}
DRUGBANK_RE  = re.compile(r"\b(DB\d{5,})\b")
PUBCHEM_RE   = re.compile(r"PUBCHEM\.COMPOUND:(\d+)")
CURIE_RE     = re.compile(r"[A-Za-z][\w.]*:[\w.\-]+")


def main() -> None:
    print(f"[{time.strftime('%H:%M:%S')}] loading gold DDIs ...")
    df = pd.read_csv(GOLD)
    gold: set[tuple[str, str]] = set()
    for a, b in zip(df["drug1_id"], df["drug2_id"]):
        a, b = str(a).strip(), str(b).strip()
        gold.add((min(a, b), max(a, b)))
    print(f"  gold pairs: {len(gold):,}")

    pubchem_to_db: dict[str, str] = {}
    if PUBCHEM_BRIDGE.exists():
        pc = pd.read_csv(PUBCHEM_BRIDGE).dropna(subset=["drugbank_id", "pubchem_cid"])
        for r in pc[["pubchem_cid", "drugbank_id"]].itertuples(index=False):
            pubchem_to_db.setdefault(str(r.pubchem_cid).strip(),
                                     str(r.drugbank_id).strip())
        print(f"  loaded {len(pubchem_to_db):,} PubChem→DrugBank bridges")

    # Multi-namespace DrugBank xref (UNII / RxCUI / ATC / KEGG.DRUG / ChEBI / ...)
    db_xref: dict[str, dict[str, str]] = {}
    if DRUGBANK_XREF.exists():
        xr = pd.read_csv(DRUGBANK_XREF).dropna(
            subset=["drugbank_id", "namespace", "external_id"])
        for r in xr[["namespace", "external_id", "drugbank_id"]].itertuples(index=False):
            ns = str(r.namespace).strip().upper()
            db_xref.setdefault(ns, {}).setdefault(
                str(r.external_id).strip(), str(r.drugbank_id).strip())
        n_total = sum(len(v) for v in db_xref.values())
        print(f"  loaded {n_total:,} multi-namespace DrugBank xrefs "
              f"({sorted(db_xref.keys())})")
    else:
        print(f"  [no {DRUGBANK_XREF} — run scripts/build_drugbank_xref.py "
              f"to recover more drugs]")

    print(f"[{time.strftime('%H:%M:%S')}] streaming nodes.tsv to build CURIE→DB map ...")
    curie_to_db: dict[str, str] = {}
    n_nodes = 0
    n_with_db = 0
    for chunk in pd.read_csv(NODES, sep="\t", dtype=str,
                             usecols=["id", "category", "equivalent_identifiers"],
                             chunksize=500_000):
        chunk = chunk.fillna("")
        # Only consider drug-like nodes (replicates loader's category filter)
        is_drug = chunk["category"].str.contains(
            r"SmallMolecule|Drug|ChemicalEntity|MolecularMixture|Food",
            regex=True, na=False)
        sub = chunk[is_drug]
        n_nodes += len(sub)
        for nid, eq in zip(sub["id"], sub["equivalent_identifiers"]):
            haystack = nid + " " + eq
            m = DRUGBANK_RE.search(haystack)
            if m:
                curie_to_db[nid] = m.group(1)
                n_with_db += 1
                continue
            # PubChem→DrugBank fallback
            if pubchem_to_db:
                pc = PUBCHEM_RE.search(haystack)
                if pc and pc.group(1) in pubchem_to_db:
                    curie_to_db[nid] = pubchem_to_db[pc.group(1)]
                    n_with_db += 1
                    continue
            # Multi-namespace xref fallback (UNII / RxCUI / ATC / KEGG.DRUG / ...)
            if db_xref:
                hit = None
                for tok in CURIE_RE.findall(haystack):
                    ns, _, ext = tok.partition(":")
                    ns_dict = db_xref.get(ns.upper())
                    if ns_dict and ext in ns_dict:
                        hit = ns_dict[ext]
                        break
                if hit:
                    curie_to_db[nid] = hit
                    n_with_db += 1
        print(f"  ... {n_nodes:,} drug-like nodes scanned, "
              f"{n_with_db:,} mapped to DrugBank ({100*n_with_db/max(n_nodes,1):.1f}%)",
              flush=True)
    print(f"[{time.strftime('%H:%M:%S')}] done. "
          f"{len(curie_to_db):,} unique CURIEs map to DrugBank IDs")
    distinct_db = len(set(curie_to_db.values()))
    print(f"  → {distinct_db:,} distinct DrugBank drugs covered by Matrix")

    print(f"[{time.strftime('%H:%M:%S')}] streaming edges.tsv ...")
    raw_any: set[tuple[str, str]] = set()
    raw_ddi: set[tuple[str, str]] = set()
    pred_counts: Counter[str] = Counter()
    n_edges = 0
    n_dd = 0
    for i, chunk in enumerate(pd.read_csv(
        EDGES, sep="\t", dtype=str,
        usecols=["subject", "predicate", "object"],
        chunksize=2_000_000,
    ), start=1):
        chunk = chunk.fillna("")
        n_edges += len(chunk)
        sub = chunk["subject"].map(curie_to_db)
        obj = chunk["object"].map(curie_to_db)
        mask = sub.notna() & obj.notna()
        kept = chunk.loc[mask].assign(s_db=sub[mask].values, o_db=obj[mask].values)
        n_dd += len(kept)
        for s_db, p, o_db in zip(kept["s_db"], kept["predicate"], kept["o_db"]):
            pair = (min(s_db, o_db), max(s_db, o_db))
            raw_any.add(pair)
            pred = p.replace("biolink:", "")
            pred_counts[pred] += 1
            if pred in DDI_PREDS:
                raw_ddi.add(pair)
        print(f"  chunk {i}: {n_edges:,} total, {n_dd:,} drug-drug "
              f"(DB-mappable both sides)", flush=True)

    overlap_any = raw_any & gold
    overlap_ddi = raw_ddi & gold
    missing     = gold - raw_any

    payload = {
        "gold_pairs":                  len(gold),
        "matrix_db_drugs":             distinct_db,
        "matrix_drug_nodes_total":     n_nodes,
        "matrix_drug_nodes_with_db":   n_with_db,
        "matrix_dd_edges_with_db":     n_dd,
        "matrix_pairs_any_pred":       len(raw_any),
        "matrix_pairs_ddi_pred":       len(raw_ddi),
        "overlap_any_pred":            len(overlap_any),
        "overlap_ddi_pred":            len(overlap_ddi),
        "coverage_any_pred_pct":       round(100 * len(overlap_any) / len(gold), 2),
        "coverage_ddi_pred_pct":       round(100 * len(overlap_ddi) / len(gold), 2),
        "gold_pairs_truly_missing":    len(missing),
        "predicate_distribution":      pred_counts.most_common(30),
        "sample_truly_missing_pairs":  [list(p) for p in list(missing)[:20]],
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2))
    print(f"\n[{time.strftime('%H:%M:%S')}] wrote {OUT}\n")
    print("=" * 60)
    print("HEADLINE")
    print("=" * 60)
    print(f"  gold DDI pairs:                            {payload['gold_pairs']:>10,}")
    print(f"  Matrix DD pairs (ANY predicate):           {payload['matrix_pairs_any_pred']:>10,}")
    print(f"  Matrix DD pairs (DDI predicates only):     {payload['matrix_pairs_ddi_pred']:>10,}")
    print(f"  overlap (ANY pred):                        {payload['overlap_any_pred']:>10,} = "
          f"{payload['coverage_any_pred_pct']:5.1f}%   ← ceiling")
    print(f"  overlap (DDI pred):                        {payload['overlap_ddi_pred']:>10,} = "
          f"{payload['coverage_ddi_pred_pct']:5.1f}%   ← what loader does today")
    print(f"  gold pairs truly missing from Matrix raw:  {len(missing):>10,} "
          f"({100*len(missing)/len(gold):.1f}%)")
    print()
    print("INTERPRETATION:")
    if payload['coverage_any_pred_pct'] - payload['coverage_ddi_pred_pct'] > 5:
        print(f"  → Adding broader predicates would gain "
              f"{payload['coverage_any_pred_pct']-payload['coverage_ddi_pred_pct']:.1f} "
              f"percentage points")
    if 100 - payload['coverage_any_pred_pct'] > 5:
        print(f"  → {100 - payload['coverage_any_pred_pct']:.1f}% of gold DDIs are NOT "
              f"in Matrix's raw data (upstream curation gap)")


if __name__ == "__main__":
    main()
