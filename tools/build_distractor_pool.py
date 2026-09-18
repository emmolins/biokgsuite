#!/usr/bin/env python3
"""Build the distractor pool for the grounding experiment from Open Targets.

The ranking task hides each post-cutoff approval among approved drugs that are
NOT indicated for the target disease. The pool therefore needs two tables, both
derived from the Open Targets Platform release already vendored under
data/crosswalks/ (release and retrieval date in OPENTARGETS_RELEASE.txt):

    data/gold_standards/bridges/approved_drug_pool.csv
        one row per approved drug that carries a DrugBank id:
        drugbank_id, drug_name, chembl_id, drug_type, id_source
    data/gold_standards/bridges/approved_indications.csv
        one row per (drug, disease) approved indication with the disease as a
        MONDO id: drugbank_id, mondo_id, source_disease_id

DrugBank ids come from the Open Targets cross-references where present, and
otherwise from the name bridge data/gold_standards/bridges/drug_name_to_drugbank.csv
(unambiguous names only). Biologics rarely carry a DrugBank cross-reference in
Open Targets, so without the name bridge the pool would be small molecules
only, and a pool that never contains an antibody makes an antibody positive
easy to spot. Disease ids that are not MONDO (EFO, HP, Orphanet) are mapped
through the MONDO SSSOM exact matches; indications whose disease cannot be
expressed in MONDO are dropped and counted in the provenance note.

Run from the repo root:  python tools/build_distractor_pool.py
"""
from __future__ import annotations

import glob
import hashlib
import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
XW = ROOT / "data/crosswalks"
BRIDGES = ROOT / "data/gold_standards/bridges"
OUT_POOL = BRIDGES / "approved_drug_pool.csv"
OUT_IND = BRIDGES / "approved_indications.csv"
OUT_PROV = BRIDGES / "APPROVED_DRUG_POOL_PROVENANCE.txt"


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", str(s).lower()).strip()


def _drugbank_xrefs(xrefs) -> list[str]:
    out: list[str] = []
    if xrefs is None:
        return out
    for entry in xrefs:
        if str(entry.get("source", "")).lower() == "drugbank":
            ids = entry.get("ids")
            out.extend(str(i) for i in (list(ids) if ids is not None else []))
    return out


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()[:16]


def main() -> None:
    mol_files = sorted(glob.glob(str(XW / "ot_drug_molecule" / "*.parquet")))
    mol = pd.concat(pd.read_parquet(f, columns=["id", "name", "drugType", "crossReferences", "maximumClinicalStage"])
                    for f in mol_files)
    approved = mol[mol.maximumClinicalStage == "APPROVAL"].copy()
    approved["db_xref"] = approved.crossReferences.map(_drugbank_xrefs).map(lambda l: l[0] if l else None)

    bridge = pd.read_csv(BRIDGES / "drug_name_to_drugbank.csv")
    unambiguous = bridge[bridge.ambiguous == 0].drop_duplicates("name_normalised").set_index("name_normalised").drugbank_id
    approved["db_name"] = approved.name.map(_norm).map(unambiguous)
    conflicts = approved[approved.db_xref.notna() & approved.db_name.notna() & (approved.db_xref != approved.db_name)]
    approved["drugbank_id"] = approved.db_xref.fillna(approved.db_name)
    approved["id_source"] = approved.db_xref.map(lambda v: "opentargets_xref" if pd.notna(v) else None).fillna(
        approved.db_name.map(lambda v: "name_bridge" if pd.notna(v) else None))
    pool = (approved.dropna(subset=["drugbank_id"])
            .rename(columns={"id": "chembl_id", "name": "drug_name", "drugType": "drug_type"})
            [["drugbank_id", "drug_name", "chembl_id", "drug_type", "id_source"]])
    # Salts and formulations of one ingredient share a DrugBank id; one row per id,
    # the ChEMBL parent (shortest name) kept so the pool never lists a drug twice.
    pool = (pool.assign(_len=pool.drug_name.str.len()).sort_values(["drugbank_id", "_len"])
            .drop_duplicates("drugbank_id").drop(columns="_len"))
    pool["drug_name"] = pool.drug_name.str.lower()
    pool = pool.sort_values("drugbank_id").reset_index(drop=True)

    # Indications: Open Targets clinical_indication at APPROVAL stage, disease -> MONDO.
    ind = pd.read_parquet(XW / "ot_clinical_indication" / "clinical_indication.parquet",
                          columns=["drugId", "diseaseId", "maxClinicalStage"])
    ind = ind[ind.maxClinicalStage == "APPROVAL"]
    dis = pd.read_parquet(XW / "ot_disease" / "disease.parquet", columns=["id", "dbXRefs"]).set_index("id")
    sssom = pd.read_csv(BRIDGES / "mondo.sssom.tsv", sep="\t", comment="#", usecols=["subject_id", "object_id", "predicate_id"])
    sssom = sssom[sssom.predicate_id == "skos:exactMatch"]
    to_mondo: dict[str, str] = {}
    for s, o in zip(sssom.subject_id, sssom.object_id):
        if s.startswith("MONDO:") and not o.startswith("MONDO:"):
            to_mondo.setdefault(o.replace(":", "_", 1), s)

    def mondo_of(did: str) -> str | None:
        if did.startswith("MONDO_"):
            return did.replace("_", ":", 1)
        if did in to_mondo:
            return to_mondo[did]
        xrefs = dis.dbXRefs.get(did)
        if xrefs is not None:
            for x in xrefs:
                if str(x).startswith("MONDO"):
                    return str(x).replace("_", ":", 1)
        return None

    chembl_to_db = pool.set_index("chembl_id").drugbank_id
    ind = ind.assign(drugbank_id=ind.drugId.map(chembl_to_db), mondo_id=ind.diseaseId.map(mondo_of))
    n_no_drug = int(ind.drugbank_id.isna().sum())
    n_no_mondo = int(ind.drugbank_id.notna().sum() - ind.dropna(subset=["drugbank_id", "mondo_id"]).shape[0])
    indications = (ind.dropna(subset=["drugbank_id", "mondo_id"])
                   .rename(columns={"diseaseId": "source_disease_id"})
                   [["drugbank_id", "mondo_id", "source_disease_id"]]
                   .drop_duplicates().sort_values(["drugbank_id", "mondo_id"]).reset_index(drop=True))

    BRIDGES.mkdir(parents=True, exist_ok=True)
    pool.to_csv(OUT_POOL, index=False)
    indications.to_csv(OUT_IND, index=False)
    release = (XW / "OPENTARGETS_RELEASE.txt").read_text().strip()
    OUT_PROV.write_text("\n".join([
        f"# Provenance for {OUT_POOL.name} and {OUT_IND.name}",
        f"# built:        {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
        "# built_by:     tools/build_distractor_pool.py",
        "# source:       Open Targets Platform drug_molecule, clinical_indication, disease (data/crosswalks/)",
        *("#   " + line for line in release.splitlines()),
        f"# name_bridge:  drug_name_to_drugbank.csv sha256 {_sha(BRIDGES / 'drug_name_to_drugbank.csv')}",
        f"# sssom:        mondo.sssom.tsv sha256 {_sha(BRIDGES / 'mondo.sssom.tsv')}",
        f"# approved drugs in Open Targets:      {len(approved)}",
        f"# with a DrugBank id (pool):           {len(pool)} ({int((pool.id_source == 'opentargets_xref').sum())} by cross-reference, "
        f"{int((pool.id_source == 'name_bridge').sum())} by name bridge; {len(conflicts)} conflicts, cross-reference kept)",
        f"# approved indications in Open Targets: {len(ind)}",
        f"#   dropped, drug has no DrugBank id:   {n_no_drug}",
        f"#   dropped, disease not in MONDO:      {n_no_mondo}",
        f"# indications kept:                     {len(indications)} ({indications.drugbank_id.nunique()} drugs, {indications.mondo_id.nunique()} diseases)",
        f"# sha256 {OUT_POOL.name}: {_sha(OUT_POOL)}",
        f"# sha256 {OUT_IND.name}:  {_sha(OUT_IND)}",
    ]) + "\n")
    print(OUT_PROV.read_text())


if __name__ == "__main__":
    main()
