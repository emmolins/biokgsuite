#!/usr/bin/env python3
"""Build the post-cutoff approval set (approvals dated after every graph's content date) for notebook 07 from the
regulator-derived approval dataset, and run the four evaluability checks.

WHAT THIS REPLACES
------------------
The previous post-cutoff set (253 -> 116, LLM-assisted; see
METHODS_BENCHMARK_CONSTRUCTION.md) could not be rebuilt from its inputs. This
builder starts from `repurposing-dataset` v0.6.0 (data/gold_standards/
approvals_raw/), a table of 1,288 regulator decisions (FDA, EMA, Health Canada,
PMDA, TGA; 2023-01-01 to 2026-09-17) that is itself rebuilt byte-for-byte from
the regulators' own files, with every drug and disease identifier resolved
through the NCATS Translator services and reviewed. Nothing here is written
from memory: every row of the output cites a regulator document.

THE FUNNEL (what the flow diagram shows)
----------------------------------------
Global stages, independent of any graph:

  S0  regulator rows            every decision in the source table
  S1  reviewed new indication   in_scope, repurposing_class == new_indication,
                                drug and disease identifiers resolved (kg_ready)
  S2  KG-mappable rows          drop vaccines (product-level ids, no KG models
                                them as compounds), combination products (a
                                composite DrugBank id matches no node), drugs
                                without a DrugBank accession, diagnostic agents
                                and tumour-agnostic approvals (target = the root
                                "cancer" term)
  S3  unique drug-disease pairs one row per (DrugBank id, MONDO id), earliest
                                decision date across regulators, agencies listed

Per graph (check 1 to 3 from the design note), each applied to the S3 pairs:

  C2  covered        drug node present (DrugBank id, degree >= MIN_DEGREE) and
                     disease node present (MONDO / DOID / MeSH id, whichever the
                     graph keys by)
  C1  prospective    approval date after the graph's content date
                     (config/kg_snapshots.yaml) AND no direct drug-disease edge
                     of any relation in the graph (leakage test)
  C3  repurposing    the new MONDO term is more than MAX_HOPS is_a steps from
                     every disease the graph already links the drug to; the
                     complement ("related to a prior indication") is a line
                     extension from the graph's point of view

  evaluable        = covered and prospective
  evaluable_strict = evaluable and repurposing
  common support   = evaluable(_strict) in every graph of the panel

  G   grounded       the drug node has at least one neighbour a grounding
                     dossier would render (a gene or pathway, ADME relations
                     excluded where the graph labels them), counted by the
                     dossier code itself (biokgsuite.grounding.KGDossiers).
                     A drug can clear the degree floor on drug-drug or
                     side-effect edges alone and still reach a language model
                     with no evidence line. Computed only for the graphs in
                     analysis_params.grounding.kgs (the dossier index is the
                     whole graph's adjacency, which the other views do not
                     need); elsewhere the column is empty. Recorded per pair
                     and as a funnel row; folded into evaluability by the
                     grounding experiment only (load_benchmark), never by
                     notebook 07, whose link-prediction scorer needs the
                     node, not the evidence text.

Check 4 (evaluation design: one scorer, rank-based metrics, clustered
bootstrap) lives in the notebook; this script only produces the pair table and
the per-graph status the notebook reads.

OUTPUTS
-------
  data/gold_standards/post_cutoff_approvals/post_cutoff_approvals.csv       S3 pairs in the loader
        schema (post_cutoff_approvals.py) plus provenance columns
  data/gold_standards/post_cutoff_approvals/post_cutoff_pair_status.csv         pair x graph: covered,
        date_ok, leaked, related_to_prior, evaluable, evaluable_strict, grounded (with drug_n_mech,
        disease_n_mech)
  data/gold_standards/post_cutoff_approvals/post_cutoff_funnel.csv              stage counts (global
        and per graph) with the drop reasons; the figure is drawn from this
  data/gold_standards/post_cutoff_approvals/post_cutoff_dropped.csv             every source row that
        left the funnel, with the stage and reason
  data/gold_standards/post_cutoff_approvals/PROVENANCE.txt
  results/checkpoints/07_kg_views/<kg>_{nodes,edges}.csv  compact per-graph
        views (drug/disease nodes with degree and dossier neighbour count,
        drug-disease edges) so the checks rerun without reloading the graphs
  results/figures/07_generalization/07h_post_cutoff_funnel.{pdf,png}

USAGE
-----
  python tools/build_post_cutoff_approvals.py --extract-views hetionet biokg   # once per graph (loads the KG)
  python tools/build_post_cutoff_approvals.py --kgs hetionet biokg --figure    # checks from the views
  python tools/build_post_cutoff_approvals.py                                  # global stages only

Every count printed is re-derived from the files on disk; nothing is carried
over from a previous run.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from collections import defaultdict, deque
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

SOURCE_DEFAULT = ROOT / "data/gold_standards/source/repurposing-dataset-v0.6.0/candidates_reviewed.csv"
GS_DIR = ROOT / "data/gold_standards/post_cutoff_approvals"   # outputs
BRIDGE_DIR = ROOT / "data/gold_standards/bridges"
VIEW_DIR = ROOT / "results/checkpoints/07_kg_views"
FIG_DIR = ROOT / "results/figures/07_generalization"

MIN_DEGREE = 5        # same floor notebook 07 applies to a drug node before scoring it
MAX_HOPS = 2          # is_a distance at or below which a new indication counts as related to a prior one
DB_RE = re.compile(r"DB\d{5}")
# Wording in the review note that marks a change of treatment setting for a disease the drug may already
# treat (a line extension in the clinical sense). Reported as a column, not used as a filter: whether the
# graph already links the drug to that disease is what check 3 tests.
NAME_OVERRIDES = {"DB00028": "human immunoglobulin G"}  # regulator strings are process descriptions
GENERIC_DISEASES = {"MONDO:0004992", "MONDO:0005070", "MONDO:0000001", "MONDO:0700096"}  # cancer, neoplasm, disease, human disease
SETTING_RE = re.compile(r"neoadjuvant|adjuvant|perioperative|first-line|maintenance|earlier line|line-of-therapy|"
                        r"line qualifier|p(a)?ediatric|age (extension|expansion)|population|expansion", re.I)


SALT_RE = re.compile(r"\b(sodium|potassium|calcium|magnesium|hydrochloride|hydrobromide|mesylate|mesilate|maleate|fumarate|"
                     r"tartrate|citrate|acetate|succinate|sulfate|sulphate|phosphate|s-malate|malate|besylate|tosylate|"
                     r"hydrate|monohydrate|dihydrate|hemihydrate|dimethyl sulfoxide|dimesylate|disodium)\b", re.I)


def clean_name(s: str) -> str:
    """First ingredient of a regulator ingredient string, without salt or process qualifiers."""
    s = re.split(r";|/|\band\b", str(s), maxsplit=1)[0]
    s = re.sub(r"\((genetical recombination|activated|human|recombinant)\)", "", s, flags=re.I)
    s = re.sub(r"-[a-z]{4}\b", "", s, flags=re.I)  # FDA biologic suffixes (-nxki, -qvfc)
    s = SALT_RE.sub("", s)
    return re.sub(r"\s+", " ", s).strip(" ,").lower()


# --------------------------------------------------------------------------- helpers
def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def mondo_curie(x) -> str | None:
    """'0004975' / 4975 / 'MONDO:0004975' / 'MONDO_0004975' -> 'MONDO:0004975'."""
    s = str(x).strip()
    m = re.search(r"(\d+)$", s)
    if not m or (not s.upper().startswith("MONDO") and not s.isdigit()):
        return None
    return f"MONDO:{int(m.group(1)):07d}"


def norm_curie(x) -> str:
    """'Disease::DOID:14227' -> 'DOID:14227'; 'MESH:D001234' unchanged; bare 'D001234' -> 'MESH:D001234'."""
    s = str(x).strip()
    if "::" in s:
        s = s.split("::", 1)[1]
    s = re.sub(r"\s*\{.*\}\s*$", "", s)
    if re.fullmatch(r"[CD]\d{6,9}", s):
        return "MESH:" + s
    if re.fullmatch(r"\d+", s):
        return f"MONDO:{int(s):07d}"
    pfx, _, local = s.partition(":")
    return f"{pfx.upper()}:{local}" if local else s


# --------------------------------------------------------------------------- bridges
class Bridges:
    """MONDO <-> DOID <-> MeSH, from the same files notebook 07 uses."""

    def __init__(self, gs_dir: Path):
        do = pd.read_csv(gs_dir / "do_diseases.csv", dtype=str)
        self.doid_to_mondo = defaultdict(set)
        self.mondo_to_doid = defaultdict(set)
        for d, m in zip(do["doid"], do["mondo_id"]):
            if isinstance(d, str) and isinstance(m, str):
                self.doid_to_mondo[d].add(m)
                self.mondo_to_doid[m].add(d)
        self.mondo_label = dict(zip(do["mondo_id"], do["mondo_name"]))
        md = pd.read_csv(gs_dir / "mesh_to_doid.csv", dtype=str, on_bad_lines="skip")
        md["mesh"] = md["mesh_id"].astype(str).str.split(" ").str[0]
        self.mesh_to_doid = defaultdict(set)
        self.doid_to_mesh = defaultdict(set)
        for ms, d in zip(md["mesh"], md["doid"]):
            if isinstance(d, str) and ms.startswith("MESH:"):
                self.mesh_to_doid[ms].add(d)
                self.doid_to_mesh[d].add(ms)
        # SSSOM adds direct MONDO<->MeSH and MONDO<->DOID exact matches and MONDO labels.
        self.mondo_to_mesh = defaultdict(set)
        self.mesh_to_mondo = defaultdict(set)
        sssom = gs_dir / "mondo.sssom.tsv"
        self.sssom_used = sssom.exists()
        if self.sssom_used:
            with open(sssom, encoding="utf-8") as f:
                for line in f:
                    if line.startswith("#") or line.startswith("subject_id"):
                        continue
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 4:
                        continue
                    subj, label, pred, obj = p[0], p[1], p[2], p[3]
                    if not pred.endswith("exactMatch"):
                        continue
                    self.mondo_label.setdefault(subj, label)
                    if obj.startswith("mesh:") or obj.startswith("MESH:"):
                        ms = "MESH:" + obj.split(":", 1)[1]
                        self.mondo_to_mesh[subj].add(ms)
                        self.mesh_to_mondo[ms].add(subj)
                    elif obj.startswith("DOID:"):
                        self.mondo_to_doid[subj].add(obj)
                        self.doid_to_mondo[obj].add(subj)
        # is_a hierarchy for the distance rule
        hier = pd.read_csv(gs_dir / "mondo_hierarchy.csv", dtype=str)
        self.adj = defaultdict(set)
        for c, p in zip(hier["child_id"], hier["parent_id"]):
            if c.startswith("MONDO:") and p.startswith("MONDO:"):
                self.adj[c].add(p)
                self.adj[p].add(c)

    def doid_for(self, mondo: str) -> str:
        return sorted(self.mondo_to_doid.get(mondo, []))[0] if self.mondo_to_doid.get(mondo) else ""

    def all_mesh_for(self, mondo: str) -> set:
        via = set(self.mondo_to_mesh.get(mondo, set()))
        for d in self.mondo_to_doid.get(mondo, []):
            via |= self.doid_to_mesh.get(d, set())
        return via

    def mesh_for(self, mondo: str) -> str:
        direct = self.mondo_to_mesh.get(mondo)
        if direct:
            return sorted(direct)[0]
        via = self.all_mesh_for(mondo)
        return sorted(via)[0] if via else ""

    def to_mondo(self, curie: str) -> set:
        """Any disease CURIE a graph uses -> the MONDO terms it maps to."""
        if curie.startswith("MONDO:"):
            return {curie}
        if curie.startswith("DOID:"):
            return set(self.doid_to_mondo.get(curie, []))
        if curie.startswith("MESH:"):
            out = set(self.mesh_to_mondo.get(curie, []))
            for d in self.mesh_to_doid.get(curie, []):
                out |= self.doid_to_mondo.get(d, set())
            return out
        return set()

    def within_hops(self, a: str, b: str, max_hops: int) -> bool:
        """True when a and b are at most max_hops is_a steps apart (undirected):
        0 = same term, 1 = parent/child, 2 = grandparent, grandchild or sibling."""
        if a == b:
            return True
        seen, frontier = {a}, deque([(a, 0)])
        while frontier:
            n, d = frontier.popleft()
            if d >= max_hops:
                continue
            for nb in self.adj.get(n, ()):
                if nb == b:
                    return True
                if nb not in seen:
                    seen.add(nb)
                    frontier.append((nb, d + 1))
        return False


# --------------------------------------------------------------------------- global stages
def split_ids(s: str) -> list:
    return [x.strip() for x in str(s).split(";") if x.strip() and x.strip() != "nan"]


def build_pairs(src: pd.DataFrame):
    """S0 -> S3. Returns (pairs DataFrame, dropped DataFrame, funnel rows)."""
    funnel, dropped = [], []
    n0 = len(src)
    funnel.append(dict(scope="global", stage="S0", label="regulator rows", n=n0, dropped=0, reason=""))

    # S1: reviewed new-indication rows with both identifiers
    keep = src[src["kg_ready"].astype(str).str.lower() == "true"].copy()
    for r in src[~src.index.isin(keep.index)].itertuples():
        if str(r.in_scope).lower() != "true":
            why = "out of scope: " + (r.exclusion_reason or "")
        elif r.repurposing_class != "new_indication":
            why = f"class={r.repurposing_class}"
        elif not r.mondo_id_final:
            why = f"disease {r.disease_status}"
        else:
            why = f"drug {r.drug_status}"
        dropped.append(dict(candidate_id=r.candidate_id, stage="S1", reason=why))
    funnel.append(dict(scope="global", stage="S1", label="reviewed new-indication rows, both identifiers resolved",
                       n=len(keep), dropped=n0 - len(keep),
                       reason="out-of-scope class, no disease target, identifiers pending or unresolved"))

    # S2: KG-mappable rows
    n1 = len(keep)
    mask_v = keep["drug_class"].astype(str).str.contains("vaccine")
    for r in keep[mask_v].itertuples():
        dropped.append(dict(candidate_id=r.candidate_id, stage="S2", reason="vaccine (product-level identifier; no graph models it as a compound)"))
    keep = keep[~mask_v]
    n_ids = keep["drug_id_final"].map(lambda s: len(split_ids(s)))
    for r in keep[n_ids > 1].itertuples():
        dropped.append(dict(candidate_id=r.candidate_id, stage="S2", reason="combination product (composite drug identifier matches no graph node)"))
    keep = keep[n_ids == 1]
    is_db = keep["drug_id_final"].str.startswith("DRUGBANK:DB")
    for r in keep[~is_db].itertuples():
        dropped.append(dict(candidate_id=r.candidate_id, stage="S2", reason=f"no DrugBank accession ({r.drug_id_final})"))
    keep = keep[is_db].copy()
    # A diagnostic agent's "indication" is the disease it images, not one it treats; the review kept the
    # MONDO term (flagged "diagnostic" in the note) so the row is identifiable, and it leaves here.
    is_dx = keep["disease_note"].str.contains("diagnostic", case=False)
    for r in keep[is_dx].itertuples():
        dropped.append(dict(candidate_id=r.candidate_id, stage="S2", reason="diagnostic agent (imaging target, not a treated disease)"))
    keep = keep[~is_dx].copy()
    # Tumour-agnostic and other generic targets ("cancer", "neoplasm") are not a disease node a graph could
    # be asked to predict; every graph links most oncology drugs to the root cancer term already.
    is_gen = keep["mondo_id_final"].map(lambda s: all(m in GENERIC_DISEASES for m in split_ids(s)))
    for r in keep[is_gen].itertuples():
        dropped.append(dict(candidate_id=r.candidate_id, stage="S2", reason=f"generic disease target ({r.mondo_label_final}; tumour-agnostic)"))
    keep = keep[~is_gen].copy()
    funnel.append(dict(scope="global", stage="S2", label="KG-mappable rows", n=len(keep), dropped=n1 - len(keep),
                       reason=f"vaccines {int(mask_v.sum())}, combination products {int((n_ids > 1).sum())}, no DrugBank id {int((~is_db).sum())}, diagnostic agents {int(is_dx.sum())}, generic disease targets {int(is_gen.sum())}"))

    # S3: unique (DrugBank, MONDO) pairs; one approval may name several diseases (multi:)
    keep["drugbank_id"] = keep["drug_id_final"].str.replace("DRUGBANK:", "", regex=False)
    keep["event_date"] = pd.to_datetime(keep["event_date"], errors="coerce")
    rows = []
    for r in keep.itertuples():
        mondos = split_ids(r.mondo_id_final)
        labels = split_ids(r.mondo_label_final) if r.mondo_label_final else []
        setting = bool(SETTING_RE.search(str(r.disease_note)))
        for i, m in enumerate(mondos):
            rows.append(dict(drugbank_id=r.drugbank_id, mondo_id=m, setting_change_noted=setting,
                             mondo_label=labels[i] if i < len(labels) else "",
                             drug_name=clean_name(r.active_ingredient or r.drug_label_final or r.product_name),
                             drug_label=clean_name(r.drug_label_final or ""), product_name=r.product_name,
                             event_date=r.event_date, regulator=r.regulator, event_type=r.event_type,
                             candidate_id=r.candidate_id, document_url=r.document_url or r.source_url,
                             drug_class=r.drug_class, disease_status=r.disease_status, class_source=r.class_source))
    ev = pd.DataFrame(rows)
    n_rows_expanded = len(ev)
    pairs = []
    for (db, m), g in ev.sort_values("event_date").groupby(["drugbank_id", "mondo_id"], sort=True):
        first = g.iloc[0]
        agencies = sorted(set(g["regulator"]))
        # the shortest regulator-listed ingredient name is the plain INN (salts and co-formulants make it longer);
        # the Node Normalizer label is the fallback (it can be a code name such as BB-10901)
        names = [n for n in g["drug_name"] if n] + [n for n in g["drug_label"] if n and not re.search(r"\d", n)
                                                     and n not in {clean_name(p) for p in g["product_name"]}]
        pairs.append(dict(
            drugbank_id=db, drug_name=NAME_OVERRIDES.get(db, min(names, key=len)).lower(),
            disease_name=first["mondo_label"] or m, mondo_id=m,
            approval_date=first["event_date"].date().isoformat() if pd.notna(first["event_date"]) else "",
            agency=first["regulator"], agencies=";".join(agencies), n_agencies=len(agencies),
            date_type=("CHMP positive opinion" if first["regulator"] == "EMA" else "regulator decision"),
            evidence_ref=";".join(g["candidate_id"]), evidence_url=first["document_url"],
            setting_change_noted=bool(g["setting_change_noted"].any()),
            verification="regulator record; identifiers resolved via NCATS Translator; ai-first-pass review (status proposed)",
            combination_type="single agent", class_source=first["class_source"]))
    pairs = pd.DataFrame(pairs)
    funnel.append(dict(scope="global", stage="S3", label="unique drug-disease pairs", n=len(pairs),
                       dropped=n_rows_expanded - len(pairs),
                       reason=f"{n_rows_expanded} (row, disease) events collapse to one pair each; {int((pairs['n_agencies'] > 1).sum())} pairs seen by more than one regulator"))
    return pairs, pd.DataFrame(dropped), funnel


def add_disease_xrefs(pairs: pd.DataFrame, br: Bridges) -> pd.DataFrame:
    pairs = pairs.copy()
    pairs["doid"] = pairs["mondo_id"].map(br.doid_for)
    pairs["mesh_id"] = pairs["mondo_id"].map(br.mesh_for)
    pairs["disease_name"] = [n if n and n != m else br.mondo_label.get(m, m) for n, m in zip(pairs["disease_name"], pairs["mondo_id"])]
    return pairs


# --------------------------------------------------------------------------- KG views
def extract_view(kg_name: str, config: dict, out_dir: Path):
    """Load one graph with the biokgsuite loader and write the compact view the checks need."""
    from biokgsuite.io import load_kg
    from biokgsuite.grounding import KGDossiers, load_entrez_to_symbol
    cfg = config["knowledge_graphs"][kg_name]
    kg, nodes = load_kg(kg_name, config)
    dt = tuple(cfg["entity_types"]["Drug"]) if not isinstance(cfg["entity_types"]["Drug"], str) else (cfg["entity_types"]["Drug"],)
    dst = tuple(cfg["entity_types"]["Disease"]) if not isinstance(cfg["entity_types"]["Disease"], str) else (cfg["entity_types"]["Disease"],)
    rel = cfg["relations"]["drug_disease"]
    rels = list(rel["relations"]) if "relations" in rel else [rel["relation"]]
    deg = pd.concat([kg["x_index"], kg["y_index"]]).value_counts()
    sub = nodes[nodes["type"].isin(dt + dst)].copy()
    sub["degree"] = sub["idx"].map(deg).fillna(0).astype(int)
    sub["role"] = ["drug" if t in dt else "disease" for t in sub["type"]]
    # How many evidence lines a grounding dossier could render for the node,
    # by the dossier code's own rule (gene and pathway neighbours for a drug,
    # ADME relations excluded where labelled; gene, pathway and related
    # disease for a disease). The cap does not apply: this asks whether there
    # is anything to show, not how much is shown. Only the grounding arms
    # need it, and only they get the whole-graph adjacency built here.
    g = config["analysis_params"]["grounding"]
    if kg_name in g["kgs"]:
        dossiers = KGDossiers(kg_name, config, ROOT, cap=g["cap"], anonymize_genes=g["anonymize_genes"],
                              bridge_mode=g["bridge_mode"], entrez_to_symbol=load_entrez_to_symbol(ROOT), loaded=(kg, nodes))
        sub["n_mech"] = [dossiers.n_mechanistic(int(i), role) for i, role in zip(sub["idx"], sub["role"])]
        del dossiers
    else:
        sub["n_mech"] = pd.NA
    # Every edge between a drug node and a disease node, whatever its relation: any direct edge is a
    # leak for the prospective test, while only the indication relations count as prior indications.
    ntype = dict(zip(nodes["idx"], nodes["type"]))
    xt, yt = kg["x_index"].map(ntype), kg["y_index"].map(ntype)
    fwd = xt.isin(dt) & yt.isin(dst)
    rev = xt.isin(dst) & yt.isin(dt)
    edges = pd.DataFrame({"drug_idx": pd.concat([kg.loc[fwd, "x_index"], kg.loc[rev, "y_index"]]),
                          "disease_idx": pd.concat([kg.loc[fwd, "y_index"], kg.loc[rev, "x_index"]]),
                          "relation": pd.concat([kg.loc[fwd, "relation"], kg.loc[rev, "relation"]])}).drop_duplicates()
    edges["is_indication"] = edges["relation"].isin(rels)
    out_dir.mkdir(parents=True, exist_ok=True)
    sub[["idx", "id", "type", "name", "role", "degree", "n_mech"] + (["uid"] if "uid" in sub.columns else [])].to_csv(out_dir / f"{kg_name}_nodes.csv", index=False)
    edges.to_csv(out_dir / f"{kg_name}_edges.csv", index=False)
    meta = dict(kg=kg_name, extracted_at=datetime.now(timezone.utc).isoformat(timespec="seconds"),
                n_nodes_total=int(len(nodes)), n_edges_total=int(len(kg)), n_drug_nodes=int((sub["role"] == "drug").sum()),
                n_disease_nodes=int((sub["role"] == "disease").sum()), n_drug_disease_edges=int(len(edges)), n_indication_edges=int(edges["is_indication"].sum()),
                n_drug_nodes_grounded=(int(((sub["role"] == "drug") & (sub["n_mech"].fillna(0) > 0)).sum())
                                       if kg_name in g["kgs"] else None),
                indication_relations=rels, drug_types=list(dt), disease_types=list(dst))
    (out_dir / f"{kg_name}_meta.json").write_text(json.dumps(meta, indent=1))
    grounded = (f" ({meta['n_drug_nodes_grounded']:,} with a renderable gene or pathway neighbour)"
                if meta["n_drug_nodes_grounded"] is not None else "")
    print(f"{kg_name}: {meta['n_drug_nodes']:,} drug nodes{grounded}, {meta['n_disease_nodes']:,} disease nodes, "
          f"{meta['n_drug_disease_edges']:,} drug-disease edges -> {out_dir}")


def load_view(kg_name: str, view_dir: Path, config: dict):
    nodes = pd.read_csv(view_dir / f"{kg_name}_nodes.csv", dtype={"id": str, "name": str})
    if "n_mech" not in nodes.columns:
        nodes["n_mech"] = float("nan")      # a view from before the grounded check, or a graph that is not a grounding arm
    edges = pd.read_csv(view_dir / f"{kg_name}_edges.csv")
    meta = json.loads((view_dir / f"{kg_name}_meta.json").read_text())
    cfg = config["knowledge_graphs"][kg_name]
    # drug index: DrugBank accession -> idx
    drug_idx = {}
    drugs = nodes[nodes["role"] == "drug"]
    if kg_name == "openbilink":
        xw_path = cfg.get("pubchem_to_drugbank_path")
        xw = pd.read_csv(ROOT / xw_path, dtype=str) if xw_path and (ROOT / xw_path).exists() else pd.DataFrame()
        pc_col = next((c for c in xw.columns if "pubchem" in c.lower()), None)
        db_col = next((c for c in xw.columns if "drugbank" in c.lower()), None)
        pc_to_db = dict(zip(xw[pc_col].astype(str).str.replace("PUBCHEM.COMPOUND:", "", regex=False), xw[db_col])) if pc_col and db_col else {}
        for i, nid in zip(drugs["idx"], drugs["id"]):
            db = pc_to_db.get(str(nid).replace("PUBCHEM.COMPOUND:", ""))
            if db:
                drug_idx.setdefault(db, i)
    else:
        for i, nid in zip(drugs["idx"], drugs["id"]):
            m = DB_RE.search(str(nid))
            if m:
                drug_idx.setdefault(m.group(0), i)
    # disease index: normalised CURIE -> idx (a PrimeKG compound id "1234_5678" indexes both)
    dis_idx = {}
    sep = cfg.get("disease_id_separator")
    for i, nid in zip(nodes.loc[nodes["role"] == "disease", "idx"], nodes.loc[nodes["role"] == "disease", "id"]):
        parts = str(nid).split(sep) if sep and sep in str(nid) else [str(nid)]
        for p in parts:
            dis_idx.setdefault(norm_curie(p), i)
    degree = dict(zip(nodes["idx"], nodes["degree"]))
    n_mech = {i: (None if pd.isna(v) else int(v)) for i, v in zip(nodes["idx"], nodes["n_mech"])}
    node_curie = {i: norm_curie(nid) for i, nid in zip(nodes["idx"], nodes["id"])}
    edge_set = set(zip(edges["drug_idx"], edges["disease_idx"]))          # any relation: leakage
    prior = defaultdict(set)                                              # indication relations only
    for d, z in zip(edges.loc[edges["is_indication"], "drug_idx"], edges.loc[edges["is_indication"], "disease_idx"]):
        prior[d].add(z)
    return dict(meta=meta, drug_idx=drug_idx, dis_idx=dis_idx, degree=degree, n_mech=n_mech, node_curie=node_curie,
                edge_set=edge_set, prior=prior, n_drugs=len(drugs), n_diseases=int((nodes["role"] == "disease").sum()))


# --------------------------------------------------------------------------- per-KG checks
def check_pairs(pairs: pd.DataFrame, kg_name: str, view: dict, br: Bridges, snapshot: dict, max_hops: int):
    cd = snapshot.get("content_date")
    content_date = pd.Timestamp(cd) if cd else None
    rows = []
    for r in pairs.itertuples():
        d_idx = view["drug_idx"].get(r.drugbank_id)
        drug_present = d_idx is not None
        drug_ok = drug_present and view["degree"].get(d_idx, 0) >= MIN_DEGREE
        # every node the disease could be in this graph: the MONDO term, every DOID and every MeSH id it
        # maps to (a graph may carry the same disease under more than one identifier); the pair is covered
        # if any of them is present and leaked if the drug has an edge to any of them
        cands = [("mondo", r.mondo_id)]
        cands += [("doid", d) for d in sorted(br.mondo_to_doid.get(r.mondo_id, []))]
        cands += [("mesh", m) for m in sorted(br.all_mesh_for(r.mondo_id))]
        z_nodes, hows = [], []
        for kind, val in cands:
            zi = view["dis_idx"].get(val) if val else None
            if zi is not None and zi not in z_nodes:
                z_nodes.append(zi)
                hows.append(kind)
        z_idx, how = (z_nodes[0], "+".join(dict.fromkeys(hows))) if z_nodes else (None, "")
        covered = bool(drug_ok and z_idx is not None)
        leaked = bool(covered and any((d_idx, zi) in view["edge_set"] for zi in z_nodes))
        if content_date is None:
            date_ok = "unknown"
        else:
            date_ok = str(bool(pd.Timestamp(r.approval_date) > content_date)) if r.approval_date else "False"
        related, n_prior, nearest = False, 0, ""
        if drug_present:
            prior_mondo = set()
            for z in view["prior"].get(d_idx, ()):
                prior_mondo |= br.to_mondo(view["node_curie"].get(z, ""))
            n_prior = len(prior_mondo)
            for pm in sorted(prior_mondo):
                if br.within_hops(r.mondo_id, pm, max_hops):
                    related, nearest = True, pm
                    break
        evaluable = bool(covered and not leaked and date_ok == "True")
        # None when the view carries no neighbour count (not a grounding arm, or extracted before the check)
        drug_n_mech = view["n_mech"].get(d_idx) if drug_present else 0
        dz = [view["n_mech"].get(zi) for zi in z_nodes]
        disease_n_mech = max((v for v in dz if v is not None), default=(0 if z_nodes and None not in dz else None))
        grounded = None if (covered and drug_n_mech is None) else bool(covered and drug_n_mech > 0)
        rows.append(dict(kg=kg_name, drugbank_id=r.drugbank_id, mondo_id=r.mondo_id,
                         drug_present=drug_present, drug_degree=int(view["degree"].get(d_idx, 0)) if drug_present else 0,
                         disease_present=z_idx is not None, disease_resolved_by=how, covered=covered,
                         date_ok=date_ok, leaked=leaked, n_prior_indications=n_prior,
                         related_to_prior=related, nearest_prior_mondo=nearest,
                         evaluable=evaluable, evaluable_strict=bool(evaluable and not related),
                         drug_n_mech=drug_n_mech, disease_n_mech=disease_n_mech, grounded=grounded))
    st = pd.DataFrame(rows)
    n = len(st)
    known = st.grounded.notna()
    grd = st.grounded.map(lambda v: bool(v) if v is not None else False)
    funnel = [
        dict(scope=kg_name, stage="S3", label="unique pairs", n=n, dropped=0, reason=""),
        dict(scope=kg_name, stage="C2", label="covered (drug and disease nodes present)", n=int(st.covered.sum()),
             dropped=int(n - st.covered.sum()),
             reason=f"drug absent {int((~st.drug_present).sum())}, drug degree < {MIN_DEGREE}: {int((st.drug_present & ~(st.drug_degree >= MIN_DEGREE)).sum())}, disease absent {int((st.drug_present & (st.drug_degree >= MIN_DEGREE) & ~st.disease_present).sum())}"),
        dict(scope=kg_name, stage="C1", label="prospective (after content date, no direct edge)", n=int(st.evaluable.sum()),
             dropped=int(st.covered.sum() - st.evaluable.sum()),
             reason=f"leaked edge {int((st.covered & st.leaked).sum())}, approved before content date {int((st.covered & ~st.leaked & (st.date_ok == 'False')).sum())}, content date unknown {int((st.covered & ~st.leaked & (st.date_ok == 'unknown')).sum())}"),
        dict(scope=kg_name, stage="C3", label=f"repurposing (> {max_hops} is_a steps from every prior indication)", n=int(st.evaluable_strict.sum()),
             dropped=int(st.evaluable.sum() - st.evaluable_strict.sum()),
             reason=f"related to a disease the graph already links the drug to (line extension) {int((st.evaluable & st.related_to_prior).sum())}"),
        # Grounding experiment only (the notebook 07 figure reindexes to S3, C2, C1, C3 and ignores this
        # row): of the prospective pairs, those whose drug node has a neighbour a dossier would render.
        # On record so "the graph lacks the pair" and "the graph holds the pair but has nothing to say
        # about it" are counted apart.
        dict(scope=kg_name, stage="G", label="grounded (drug has a renderable gene or pathway neighbour)",
             n=int((st.evaluable & grd).sum()) if known.all() else None,
             dropped=int((st.evaluable & ~grd).sum()) if known.all() else None,
             reason=(f"covered but no gene or pathway neighbour {int((st.covered & ~grd).sum())}, "
                     f"of which prospective {int((st.evaluable & ~grd).sum())}") if known.all()
                    else "not computed (graph is not a grounding arm, or the view predates the check)"),
    ]
    return st, funnel


# --------------------------------------------------------------------------- figure
def draw_funnel(funnel: pd.DataFrame, out_stem: Path, kg_order: list, version: str):
    """Two panels: the global stages as a box funnel (a) and the per-graph checks as grouped bars (b).

    Drawn in the house style (FIGURE_STYLE_SPEC.md): 183 mm wide, 6 to 7 pt text, panel letters, graph
    colours from KG_PALETTE with the published graph names, and no title (the caption carries it).
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
    from biokgsuite.viz.style import (DOUBLE_COL_W, KG_NAMES, KG_ORDER, KG_PALETTE, MUTED, NEUTRAL, TEXT_COLOR,
                                      clean_ax, facet_title, panel_label, setup_style)
    setup_style()
    FILL, FINAL_FILL = "#F2F2F2", "#EAF3FA"   # one neutral fill for every box; the final set is tinted, not framed
    g = funnel[funnel.scope == "global"].reset_index(drop=True)
    present = set(funnel.scope)
    kgs = [k for k in KG_ORDER if k in present] + [k for k in kg_order if k in present and k not in KG_ORDER]
    fig_w = DOUBLE_COL_W if kgs else DOUBLE_COL_W / 2
    fig, axes = plt.subplots(1, 2 if kgs else 1, figsize=(fig_w, 3.3), squeeze=False,
                             gridspec_kw=dict(width_ratios=[1.15, 1]) if kgs else None, constrained_layout=True)
    axes = axes[0]
    ax = axes[0]
    ax.set_xlim(0, 10); ax.set_ylim(0, 10.6); ax.axis("off")
    top, h, gap = 10.5, 1.5, 1.15
    for i, r in g.iterrows():
        y = top - i * (h + gap)
        final = i == len(g) - 1
        ax.add_patch(FancyBboxPatch((0, y - h), 6.8, h, boxstyle="round,pad=0.02,rounding_size=0.06",
                                    linewidth=0.5, edgecolor=TEXT_COLOR, facecolor=FINAL_FILL if final else FILL, zorder=2))
        # Count and stage code stacked at the left; the label has the rest of the box and may take two lines.
        ax.text(0.3, y - h / 2 + 0.3, f"{r.n:,}", fontsize=7, color=TEXT_COLOR, va="center", zorder=3)
        ax.text(0.3, y - h / 2 - 0.35, r.stage, fontsize=6, color=MUTED, va="center", zorder=3)
        label = str(r.label)[:1].upper() + str(r.label)[1:]
        ax.text(1.7, y - h / 2, _wrap(label, 44), fontsize=7, color=TEXT_COLOR, va="center", zorder=3, linespacing=1.1)
        if i < len(g) - 1:
            nxt = g.iloc[i + 1]
            ax.add_patch(FancyArrowPatch((0.6, y - h - 0.05), (0.6, y - h - gap + 0.1), arrowstyle="-|>",
                                         mutation_scale=6, color=NEUTRAL, lw=0.6, zorder=3))
            ax.text(1.0, y - h - gap / 2 + 0.02, f"−{int(nxt.dropped):,}", fontsize=7, color=TEXT_COLOR, va="center")
            ax.text(2.1, y - h - gap / 2 + 0.02, _wrap(nxt.reason, 66), fontsize=6, color=MUTED, va="center", linespacing=1.15)
    facet_title(ax, "Global stages")
    footer = ""
    if kgs:
        ax2 = axes[1]
        stages = ["S3", "C2", "C1", "C3"]
        labels = ["Unique\npairs", "Covered\n(C2)", "Prospective\n(C1)", "Repurposing\n(C3)"]
        x = range(len(stages))
        width = 0.6 / max(len(kgs), 1)
        # A count on every bar only when the panel is sparse enough to read them (two graphs, eight bars).
        show_values = len(kgs) * len(stages) <= 8
        for j, k in enumerate(kgs):
            sub = funnel[funnel.scope == k].set_index("stage").reindex(stages)
            col = KG_PALETTE.get(k, NEUTRAL)
            bars = ax2.bar([xi + (j - (len(kgs) - 1) / 2) * width for xi in x], sub["n"].fillna(0), width=width,
                           color=col, label=KG_NAMES.get(k, k), zorder=2)
            if show_values:
                for b, v in zip(bars, sub["n"].fillna(0)):
                    ax2.text(b.get_x() + b.get_width() / 2, b.get_height() + 1, f"{int(v)}", ha="center", va="bottom",
                             fontsize=6, color=TEXT_COLOR)
        ax2.set_xticks(list(x)); ax2.set_xticklabels(labels)
        ax2.set_xlim(-0.5, len(stages) - 0.5)
        ax2.set_ylim(0, ax2.get_ylim()[1] * (1.06 if show_values else 1.0))
        clean_ax(ax2, ylabel="Pairs")
        facet_title(ax2, "Per-graph checks on the S3 pairs")
        # Every graph keeps every S3 pair, so the tallest group is at the left and the top right is free.
        ax2.legend(loc="upper right", ncol=1 if len(kgs) <= 3 else 2)
        cs = funnel[funnel.scope == "common_support"]
        if len(cs):
            counts = {r.label: int(r.n) for r in cs.itertuples()}
            ev = counts.get("evaluable in every graph"); st = counts.get("strict in every graph")
            footer = f"Pairs evaluable in all {len(kgs)} graphs: {ev}" + (f" (strict rule: {st})" if st is not None else "")
        for a, letter in zip(axes, "ab"):
            panel_label(a, letter, x=-0.02 if a is ax else -0.2)
    # Only the common-support count goes under the panels; the set version and
    # the script that derived the counts are in PROVENANCE.txt and the caption.
    if footer:
        fig.supxlabel(_wrap(footer, 150), fontsize=6, color=MUTED, x=0.0, ha="left", linespacing=1.2)
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(out_stem.with_suffix("." + ext), dpi=600, bbox_inches="tight")
    plt.close(fig)
    print("figure ->", out_stem.with_suffix(".pdf"))


def _wrap(s: str, width: int) -> str:
    words, lines, cur = str(s).split(), [], ""
    for w in words:
        if len(cur) + len(w) + 1 > width and cur:
            lines.append(cur)
            cur = w
        else:
            cur = (cur + " " + w).strip()
    if cur:
        lines.append(cur)
    return "\n".join(lines[:3])


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", type=Path, default=SOURCE_DEFAULT)
    ap.add_argument("--extract-views", nargs="*", metavar="KG", help="load these graphs and write compact views")
    ap.add_argument("--kgs", nargs="*", metavar="KG", help="run checks 1 to 3 on these graphs from their views ('all' = every view present)")
    ap.add_argument("--max-hops", type=int, default=None, help="is_a distance for check 3 (default: prospective_max_hops in config/params.yaml)")
    ap.add_argument("--figure", action="store_true")
    ap.add_argument("--out-dir", type=Path, default=GS_DIR)
    a = ap.parse_args()

    from biokgsuite.io import load_config
    config = load_config(ROOT, trial_subset=False)   # graphs are named on the command line, not by the trial file
    snaps = yaml.safe_load((ROOT / "config/kg_snapshots.yaml").read_text())["snapshots"]
    if a.max_hops is None:
        a.max_hops = int(config["analysis_params"]["prospective_max_hops"])   # no silent default

    if a.extract_views:
        for k in a.extract_views:
            extract_view(k, config, VIEW_DIR)

    src = pd.read_csv(a.source, dtype=str).fillna("")
    version_file = a.source.parent / "VERSION"
    version = version_file.read_text().strip() if version_file.exists() else "unversioned"
    br = Bridges(BRIDGE_DIR)
    pairs, dropped, funnel = build_pairs(src)
    pairs = add_disease_xrefs(pairs, br)
    n_doid = int((pairs["doid"] != "").sum())
    n_mesh = int((pairs["mesh_id"] != "").sum())
    print(f"S3: {len(pairs)} pairs, {pairs.drugbank_id.nunique()} drugs, {pairs.mondo_id.nunique()} diseases; "
          f"DOID for {n_doid}, MeSH for {n_mesh}")

    status_frames = []
    kgs = a.kgs or []
    if kgs == ["all"]:
        kgs = sorted(p.name[:-len("_meta.json")] for p in VIEW_DIR.glob("*_meta.json"))
    for k in kgs:
        if not (VIEW_DIR / f"{k}_meta.json").exists():
            print(f"{k}: no view in {VIEW_DIR}; run --extract-views {k} first", file=sys.stderr)
            continue
        snap = snaps.get(k, {})
        if not snap.get("verified"):
            print(f"note: {k} content date {snap.get('content_date')} is unverified (config/kg_snapshots.yaml)", file=sys.stderr)
        view = load_view(k, VIEW_DIR, config)
        st, f = check_pairs(pairs, k, view, br, snap, a.max_hops)
        status_frames.append(st)
        funnel.extend(f)
        print(f"{k}: covered {int(st.covered.sum())}, leaked {int(st.leaked.sum())}, date_ok {int((st.date_ok == 'True').sum())}, "
              f"evaluable {int(st.evaluable.sum())}, strict {int(st.evaluable_strict.sum())}")
    if status_frames:
        status = pd.concat(status_frames, ignore_index=True)
        done = sorted(status.kg.unique())
        wide = status.pivot_table(index=["drugbank_id", "mondo_id"], columns="kg", values=["evaluable", "evaluable_strict"], aggfunc="first")
        cs_eval = int(wide["evaluable"][done].all(axis=1).sum())
        cs_strict = int(wide["evaluable_strict"][done].all(axis=1).sum())
        funnel.append(dict(scope="common_support", stage="C1", label="evaluable in every graph", n=cs_eval, dropped=0, reason=",".join(done)))
        funnel.append(dict(scope="common_support", stage="C3", label="strict in every graph", n=cs_strict, dropped=0, reason=",".join(done)))
        status.to_csv(a.out_dir / "post_cutoff_pair_status.csv", index=False)
        print(f"common support over {done}: evaluable {cs_eval}, strict {cs_strict}")

    cols = ["drugbank_id", "drug_name", "disease_name", "doid", "approval_date", "agency", "mondo_id", "mesh_id",
            "evidence_ref", "evidence_url", "verification", "combination_type", "agencies", "n_agencies", "date_type", "class_source",
            "setting_change_noted"]
    a.out_dir.mkdir(parents=True, exist_ok=True)
    pairs[cols].to_csv(a.out_dir / "post_cutoff_approvals.csv", index=False)
    funnel_df = pd.DataFrame(funnel)
    funnel_df.to_csv(a.out_dir / "post_cutoff_funnel.csv", index=False)
    dropped.to_csv(a.out_dir / "post_cutoff_dropped.csv", index=False)
    prov = [
        "# Provenance for post_cutoff_approvals.csv, post_cutoff_pair_status.csv, post_cutoff_funnel.csv",
        f"# built:          {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
        "# built_by:       tools/build_post_cutoff_approvals.py",
        f"# source:         {a.source.relative_to(ROOT) if a.source.is_relative_to(ROOT) else a.source}",
        f"# source_version: {version}",
        f"# source_sha256:  {sha256(a.source)}",
        f"# bridges:        do_diseases.csv, mesh_to_doid.csv, mondo_hierarchy.csv, mondo.sssom.tsv ({'used' if br.sssom_used else 'absent'})",
        f"# min_degree:     {MIN_DEGREE}",
        f"# max_hops:       {a.max_hops}",
        f"# kgs_checked:    {','.join(kgs) if kgs else 'none'}",
        f"# n_pairs:        {len(pairs)}",
        "# approval_date is the earliest regulator decision in the source; for EMA-first pairs it is the CHMP positive-opinion date (date_type column).",
        "# Cite the source dataset version and the graph release dates in config/kg_snapshots.yaml, not the build date.",
    ]
    (a.out_dir / "PROVENANCE.txt").write_text("\n".join(prov) + "\n")
    if a.figure:
        draw_funnel(funnel_df, FIG_DIR / "07h_post_cutoff_funnel", list(config["knowledge_graphs"].keys()), version)
    print(f"-> {a.out_dir / 'post_cutoff_approvals.csv'} ({len(pairs)} pairs); funnel and dropped tables alongside")


if __name__ == "__main__":
    main()
