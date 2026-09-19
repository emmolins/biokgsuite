"""
Loader and validator for the post-cutoff approval set (approvals dated after every graph's content date).

WHY THIS EXISTS
---------------
The prospective arm used to be assembled inside the notebook from three
sources: a hand-curated FDA CSV, a Python module of literal EMA/PMDA/TGA
tuples, and a 269-line alias map translating free-text indication strings
("Ph+ CML-CP") into Disease Ontology names. Three consequences:

  * Disease resolution went through NAMES, so a term the alias map did not
    anticipate was silently dropped rather than failing loudly.
  * The FDA CSV was hand-curated and, when it went missing, the section
    degraded to a much weaker claim that nothing downstream surfaced.
  * Nothing checked the gold standard itself. A candidate file with composite
    DrugBank IDs ("DB08912,DB08911") matches no knowledge-graph node, so such
    a pair can never be flagged as leaked and never counts as covered; it
    passes every filter by failing to match anything.

This module replaces all three with one tidy table whose identifiers are
already resolved, and validates that table before the notebook uses it.

SCHEMA
------
Required: drugbank_id, drug_name, disease_name, doid, approval_date, agency
Optional but used when present:
    mondo_id   PrimeKG keys diseases by MONDO; without it that graph's disease
               lookup, and with it the per-graph leakage test, cannot resolve
    mesh_id    BioKG and DRKG key diseases by MeSH
    evidence_ref, evidence_url, verification, combination_type

VALIDATION
----------
Errors (refuse to load):
    composite drugbank_id, duplicate (drugbank_id, doid), unparseable date,
    missing required column
Warnings (load, but report):
    mondo_id absent and not derivable, approval_date outside the declared
    window, DOID label sharing no word with disease_name
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

REQUIRED = ["drugbank_id", "drug_name", "disease_name",
            "approval_date", "agency"]
# At least one of these must be present on every row. DOID alone was required
# originally, which was wrong: the Disease Ontology has no term for a good
# number of real indications (enthesitis-related juvenile idiopathic arthritis,
# isolated short stature, complement 3 glomerulopathy), and MONDO does. MONDO is
# also what PrimeKG keys diseases by, so requiring DOID rejected rows the
# pipeline could actually resolve.
DISEASE_ID_COLS = ["doid", "mondo_id", "mesh_id"]
OPTIONAL = ["doid", "mondo_id", "mesh_id", "evidence_ref", "evidence_url",
            "verification", "combination_type"]

_DB_RE = re.compile(r"^DB\d{5}$")
_DOID_RE = re.compile(r"^DOID:\d+$")


@dataclass
class LoadReport:
    n_rows: int = 0
    n_usable: int = 0
    n_with_doid: int = 0
    n_mondo_only: int = 0
    errors: list = field(default_factory=list)
    warnings: list = field(default_factory=list)
    mondo_backfilled: int = 0
    suspect_disease_ids: list = field(default_factory=list)

    def render(self) -> str:
        out = [f"post-cutoff approval set: {self.n_usable}/{self.n_rows} rows usable"]
        if self.mondo_backfilled:
            out.append(f"  MONDO IDs backfilled from DOID: {self.mondo_backfilled}")
        if self.n_mondo_only:
            out.append(f"  rows resolved by MONDO with no DOID available: "
                       f"{self.n_mondo_only}")
        for e in self.errors:
            out.append(f"  ERROR   {e}")
        for w in self.warnings:
            out.append(f"  WARNING {w}")
        if self.suspect_disease_ids:
            out.append(f"  {len(self.suspect_disease_ids)} row(s) where the DOID's "
                       f"ontology label shares no word with disease_name:")
            for d, lab, nm in self.suspect_disease_ids[:12]:
                out.append(f"      {d}  label='{lab}'  stated='{nm}'")
            out.append("  Check these before trusting the result: the DOID is what "
                       "the graph is queried with, so a wrong ID scores the model "
                       "against the wrong answer.")
        return "\n".join(out)


def _norm_tokens(s: str) -> set:
    return {w for w in re.sub(r"[^a-z0-9 ]+", " ", str(s).lower()).split()
            if len(w) > 3}


def load_post_cutoff_approvals(
    path, do_diseases_path=None, window=("2023-01-01", "2026-12-31"),
    strict: bool = True,
):
    """Read, validate and normalise the post-cutoff approval set.

    Returns (DataFrame, LoadReport). With strict=True a fatal validation error
    raises; with strict=False the offending rows are dropped and recorded.
    """
    path = Path(path)
    rep = LoadReport()
    if not path.exists():
        rep.errors.append(f"file not found: {path}")
        if strict:
            raise FileNotFoundError(path)
        return pd.DataFrame(columns=REQUIRED + OPTIONAL), rep

    df = pd.read_csv(path, dtype=str)
    df.columns = [c.strip().lower() for c in df.columns]
    rep.n_rows = len(df)

    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        rep.errors.append(f"missing required column(s): {missing}")
        if strict:
            raise ValueError(f"post-cutoff approval set is missing {missing}; "
                             f"found {list(df.columns)}")
        return df, rep
    for c in OPTIONAL:
        if c not in df.columns:
            df[c] = pd.NA

    df = df.copy()
    for c in REQUIRED + OPTIONAL:
        df[c] = df[c].astype("string").str.strip()

    # --- fatal: composite drug identifiers --------------------------------
    composite = df[df.drugbank_id.fillna("").str.contains(r"[,;+/]", regex=True)]
    if len(composite):
        detail = ", ".join(f"{r.drug_name}={r.drugbank_id}"
                           for r in composite.head(5).itertuples())
        msg = (f"{len(composite)} row(s) carry a composite drugbank_id ({detail}). "
               f"A composite ID matches no graph node, so the pair can never be "
               f"flagged as leaked and never counts as covered. Resolve each to "
               f"the single component the approval is held by.")
        rep.errors.append(msg)
        if strict:
            raise ValueError(msg)
        df = df.drop(composite.index)

    bad_db = df[~df.drugbank_id.fillna("").str.match(_DB_RE)]
    if len(bad_db):
        rep.warnings.append(f"{len(bad_db)} row(s) have a drugbank_id that is not "
                            f"DBnnnnn, e.g. {list(bad_db.drugbank_id.head(3))}")

    # Warn only on a value that is present and wrong, not on an absent one.
    _has_doid = df.doid.notna() & df.doid.fillna("").str.len().gt(0)
    bad_doid = df[_has_doid & ~df.doid.fillna("").str.match(_DOID_RE)]
    if len(bad_doid):
        rep.warnings.append(f"{len(bad_doid)} row(s) have a malformed doid, "
                            f"e.g. {list(bad_doid.doid.head(3))}")

    # Every row needs at least one disease identifier to reach a graph node.
    _no_id = df[df[DISEASE_ID_COLS].isna().all(axis=1)]
    if len(_no_id):
        msg = (f"{len(_no_id)} row(s) carry no disease identifier in any of "
               f"{DISEASE_ID_COLS}, so they cannot be resolved to a graph node "
               f"by anything but the disease name: "
               f"{list(_no_id.disease_name.head(4))}")
        rep.errors.append(msg)
        if strict:
            raise ValueError(msg)
        df = df.drop(_no_id.index)
    rep.n_with_doid = int(df.doid.notna().sum()) if "doid" in df else 0
    rep.n_mondo_only = int((df.doid.isna() & df.mondo_id.notna()).sum())

    # --- fatal: duplicate pairs -------------------------------------------
    # Key on the best identifier the row actually carries, falling back to the
    # disease name. Keying on DOID alone made every DOID-less row collide with
    # every other DOID-less row for the same drug, so two genuinely different
    # indications (psoriasis-related vs enthesitis-related juvenile idiopathic
    # arthritis) looked like one duplicated pair and one of them was dropped.
    def _disease_key(row):
        for col in DISEASE_ID_COLS:
            v = row.get(col)
            if pd.notna(v) and str(v).strip():
                return str(v).strip()
        return "name:" + re.sub(r"\s+", " ", str(row.get("disease_name", "")).lower()).strip()

    df["_pair"] = (df.drugbank_id.fillna("") + "|"
                   + df.apply(_disease_key, axis=1))
    dup = df[df._pair.duplicated(keep=False)]
    if len(dup):
        msg = (f"{dup._pair.nunique()} (drug, disease) pair(s) appear more than "
               f"once ({len(dup)} rows). The evaluation asks one ranking question "
               f"per pair, so duplicates inflate n and are perfectly correlated. "
               f"e.g. {list(dup.drug_name.head(4))}")
        rep.errors.append(msg)
        if strict:
            raise ValueError(msg)
        df = df.drop_duplicates("_pair", keep="first")

    # --- dates -------------------------------------------------------------
    parsed = pd.to_datetime(df.approval_date, errors="coerce")
    if parsed.isna().any():
        n = int(parsed.isna().sum())
        msg = f"{n} row(s) have an unparseable approval_date"
        rep.errors.append(msg)
        if strict:
            raise ValueError(msg)
        df = df[parsed.notna()]
        parsed = parsed[parsed.notna()]
    df["approval_date_parsed"] = parsed
    lo, hi = pd.Timestamp(window[0]), pd.Timestamp(window[1])
    out_of_window = df[(df.approval_date_parsed < lo) | (df.approval_date_parsed > hi)]
    if len(out_of_window):
        rep.warnings.append(
            f"{len(out_of_window)} row(s) fall outside the declared window "
            f"{window[0]}..{window[1]}, e.g. "
            f"{[(r.drug_name, str(r.approval_date)[:10]) for r in out_of_window.head(3).itertuples()]}. "
            f"A pair approved before the graph snapshot is not a prospective test.")

    # --- MONDO backfill and disease-label sanity ---------------------------
    if do_diseases_path and Path(do_diseases_path).exists():
        do = pd.read_csv(do_diseases_path, dtype=str)
        doid_to_mondo = dict(zip(do.doid, do.mondo_id))
        doid_to_label = dict(zip(do.doid, do.mondo_name))

        fill = df.mondo_id.isna() & df.doid.map(lambda d: d in doid_to_mondo)
        rep.mondo_backfilled = int(fill.sum())
        df.loc[fill, "mondo_id"] = df.loc[fill, "doid"].map(doid_to_mondo)

        still = df[df.mondo_id.isna()]
        if len(still):
            rep.warnings.append(
                f"{len(still)} row(s) still have no MONDO ID, so PrimeKG cannot "
                f"resolve their disease and they are invisible to that graph's "
                f"leakage and coverage tests: "
                f"{[r.disease_name for r in still.head(4).itertuples()]}")

        for r in df.itertuples():
            lab = doid_to_label.get(r.doid) or doid_to_label.get(r.mondo_id)
            if lab and not (_norm_tokens(lab) & _norm_tokens(r.disease_name)):
                rep.suspect_disease_ids.append((r.doid, lab, r.disease_name))

    df = df.drop(columns=["_pair"])
    rep.n_usable = len(df)
    return df.reset_index(drop=True), rep


def to_new_indication_approvals(df: pd.DataFrame) -> dict:
    """Shape the table into the structure the notebook's downstream cells expect:
    {drugbank_id: (drug_name, [(disease_term, approval_date, note), ...])}.

    The disease TERM carried here is the ontology-resolved name where one is
    available, not the free-text indication, so name-based lookup no longer
    depends on the alias map guessing correctly.
    """
    out: dict = {}
    for r in df.itertuples():
        dbid = str(r.drugbank_id)
        term = str(r.disease_name)
        note = f"{r.agency} — {term}"
        if dbid not in out:
            out[dbid] = (str(r.drug_name), [])
        if term.lower() not in {t[0].lower() for t in out[dbid][1]}:
            out[dbid][1].append((term, str(r.approval_date), note))
    return out


def disease_id_lookup(df: pd.DataFrame) -> dict:
    """(drugbank_id, disease_name) -> {'doid':.., 'mondo':.., 'mesh':..}

    Lets the per-KG resolver try identifiers first and fall back to the name,
    instead of relying on name matching alone.
    """
    out = {}
    for r in df.itertuples():
        out[(str(r.drugbank_id), str(r.disease_name).lower())] = {
            "doid": None if pd.isna(r.doid) else str(r.doid),
            "mondo": None if pd.isna(r.mondo_id) else str(r.mondo_id),
            "mesh": None if pd.isna(r.mesh_id) else str(r.mesh_id),
        }
    return out


# --------------------------------------------------------------------------
# Identifier spelling
# --------------------------------------------------------------------------
# Kept here rather than in the notebook so it can be tested. Graphs disagree
# about whether the vocabulary prefix is part of a disease node's id.

def norm_curie(x) -> str:
    """'Disease::DOID:8778' -> 'DOID8778';  'MESH:D005317' -> 'MESHD005317'."""
    s = str(x)
    if "::" in s:
        s = s.split("::", 1)[1]
    return re.sub(r"[^A-Z0-9]", "", s.upper())


def curie_variants(v) -> list:
    """Every spelling a graph might use for one identifier, most specific first.

    Hetionet writes `Disease::DOID:8778`, BioKG writes a bare MeSH accession
    `D005317`, PrimeKG writes bare MONDO digits. Trying only the prefixed form
    matched 0 of the 92 MeSH-carrying gold rows against BioKG, which holds
    11,347 disease nodes; adding the prefix-stripped form matches 90.

    The bare form is offered only when it still carries a letter (`D005317`),
    never when it reduces to digits alone: `DOID:684` and a MONDO node numbered
    684 are different diseases, and a numeric match across vocabularies would
    bind an indication to the wrong node while looking like a clean hit. Bare
    numeric node ids are reached from the index side instead, where they are
    registered under their own vocabulary prefix.
    """
    k = norm_curie(v)
    out = [k, re.sub(r"^(MONDO|DOID|MESH)0*", r"\1", k)]
    bare = re.sub(r"^(MONDO|DOID|MESH|UMLS|OMIM)", "", k)
    if bare and not bare.isdigit():
        out.append(bare)
    seen, ordered = set(), []
    for x in out:
        if x and x not in seen:
            seen.add(x)
            ordered.append(x)
    return ordered
