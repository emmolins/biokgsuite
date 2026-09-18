"""Guards on the post-cutoff approval set loader.

Each test corresponds to a defect found in the 2026-09 audit of the previous
gold standard, so a regression reintroduces a known-bad state loudly.
"""
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from biokgsuite.gold.post_cutoff_approvals import (  # noqa: E402
    load_post_cutoff_approvals, to_new_indication_approvals, disease_id_lookup,
)

BASE = Path(__file__).resolve().parents[1]
DO = BASE / "data/gold_standards/bridges/do_diseases.csv"

GOOD = [
    dict(drugbank_id="DB11652", drug_name="tucatinib",
         disease_name="colorectal carcinoma", doid="DOID:0080199",
         approval_date="2023-01-19", agency="FDA", mondo_id="MONDO:0024331"),
    dict(drugbank_id="DB15035", drug_name="zanubrutinib",
         disease_name="chronic myeloid leukemia", doid="DOID:0081088",
         approval_date="2023-03-07", agency="FDA", mondo_id=""),
]


def _write(tmp_path, rows, name="gs.csv"):
    p = tmp_path / name
    pd.DataFrame(rows).to_csv(p, index=False)
    return p


def test_loads_a_conforming_file(tmp_path):
    df, rep = load_post_cutoff_approvals(_write(tmp_path, GOOD), DO)
    assert len(df) == 2
    assert rep.n_usable == 2
    assert not rep.errors


def test_composite_drugbank_id_is_fatal(tmp_path):
    """A composite ID matches no graph node, so it can never be flagged as
    leaked and never counts as covered: it passes filters by matching nothing."""
    rows = GOOD + [dict(GOOD[0], drugbank_id="DB08912,DB08911",
                        drug_name="dabrafenib+trametinib", doid="DOID:3908")]
    with pytest.raises(ValueError, match="composite"):
        load_post_cutoff_approvals(_write(tmp_path, rows), DO)
    df, rep = load_post_cutoff_approvals(_write(tmp_path, rows), DO, strict=False)
    assert len(df) == 2 and any("composite" in e for e in rep.errors)


def test_duplicate_pair_is_fatal(tmp_path):
    """One ranking question per pair; duplicates inflate n and correlate."""
    rows = GOOD + [dict(GOOD[0], drug_name="tucatinib (EMA)", agency="EMA")]
    with pytest.raises(ValueError, match="more than"):
        load_post_cutoff_approvals(_write(tmp_path, rows), DO)


def test_unparseable_date_is_fatal(tmp_path):
    rows = GOOD + [dict(GOOD[1], drugbank_id="DB00001", approval_date="soon")]
    with pytest.raises(ValueError, match="approval_date"):
        load_post_cutoff_approvals(_write(tmp_path, rows), DO)


def test_out_of_window_warns_not_fails(tmp_path):
    """A 2018 approval cannot postdate the graph snapshot."""
    rows = GOOD + [dict(GOOD[1], drugbank_id="DB00002",
                        approval_date="2018-11-21", doid="DOID:4398")]
    df, rep = load_post_cutoff_approvals(_write(tmp_path, rows), DO)
    assert len(df) == 3
    assert any("outside the declared window" in w for w in rep.warnings)


def test_missing_required_column_is_fatal(tmp_path):
    """`doid` is deliberately NOT in this list: MONDO alone is sufficient."""
    rows = [{k: v for k, v in GOOD[0].items() if k != "approval_date"}]
    with pytest.raises(ValueError, match="missing"):
        load_post_cutoff_approvals(_write(tmp_path, rows), DO)


def test_a_missing_doid_column_is_not_fatal(tmp_path):
    rows = [{k: v for k, v in r.items() if k != "doid"} for r in MONDO_ONLY]
    df, rep = load_post_cutoff_approvals(_write(tmp_path, rows), DO,
                                             strict=False)
    assert len(df) == 2 and not rep.errors


@pytest.mark.skipif(not DO.exists(), reason="ontology table absent")
def test_mondo_is_backfilled_from_doid(tmp_path):
    """Without a MONDO ID, PrimeKG cannot resolve the disease, which silently
    disables that graph's leakage and coverage tests for the row."""
    df, rep = load_post_cutoff_approvals(_write(tmp_path, GOOD), DO)
    assert rep.mondo_backfilled >= 1
    assert df.mondo_id.notna().all()


@pytest.mark.skipif(not DO.exists(), reason="ontology table absent")
def test_mismatched_disease_id_is_surfaced(tmp_path):
    """DOID:684 is hepatocellular carcinoma; labelling it hepatitis C scores the
    model against the wrong answer."""
    rows = GOOD + [dict(GOOD[0], drugbank_id="DB13879",
                        drug_name="glecaprevir", disease_name="acute hepatitis C",
                        doid="DOID:684")]
    df, rep = load_post_cutoff_approvals(_write(tmp_path, rows), DO)
    assert any(d == "DOID:684" for d, _, _ in rep.suspect_disease_ids)


def test_downstream_shapes(tmp_path):
    df, _ = load_post_cutoff_approvals(_write(tmp_path, GOOD), DO)
    nia = to_new_indication_approvals(df)
    assert set(nia) == {"DB11652", "DB15035"}
    assert nia["DB11652"][0] == "tucatinib"
    assert len(nia["DB11652"][1]) == 1
    lut = disease_id_lookup(df)
    assert lut[("DB11652", "colorectal carcinoma")]["doid"] == "DOID:0080199"


# --- regressions from the 2026-09-18 dataset ---------------------------------

MONDO_ONLY = [
    dict(drugbank_id="DB11569", drug_name="ixekizumab",
         disease_name="psoriasis-related juvenile idiopathic arthritis",
         doid="", mondo_id="MONDO:0019436",
         approval_date="2025-07-25", agency="EMA"),
    dict(drugbank_id="DB11569", drug_name="ixekizumab",
         disease_name="enthesitis-related juvenile idiopathic arthritis",
         doid="", mondo_id="MONDO:0019437",
         approval_date="2025-07-25", agency="EMA"),
]


def test_two_diseases_without_doid_are_not_a_duplicate(tmp_path):
    """Keying the pair on DOID alone made every DOID-less row for one drug
    collide with every other, so two genuinely different indications looked
    like one duplicated pair and one was silently dropped."""
    df, rep = load_post_cutoff_approvals(_write(tmp_path, MONDO_ONLY), DO)
    assert len(df) == 2
    assert not rep.errors
    assert set(df.disease_name) == {r["disease_name"] for r in MONDO_ONLY}


def test_a_mondo_id_is_sufficient_when_the_disease_ontology_has_no_term(tmp_path):
    """The Disease Ontology has no term for a number of real indications that
    MONDO does carry, and MONDO is what PrimeKG keys diseases by, so requiring
    DOID rejected rows the pipeline could resolve."""
    df, rep = load_post_cutoff_approvals(_write(tmp_path, MONDO_ONLY), DO)
    assert len(df) == 2
    assert rep.n_mondo_only == 2
    assert not any("malformed doid" in w for w in rep.warnings)


def test_a_row_with_no_disease_identifier_at_all_is_fatal(tmp_path):
    rows = [dict(MONDO_ONLY[0], doid="", mondo_id="", mesh_id="")]
    with pytest.raises(ValueError, match="no disease identifier"):
        load_post_cutoff_approvals(_write(tmp_path, rows), DO)


def test_a_genuine_duplicate_is_still_caught_without_doid(tmp_path):
    """Widening the key must not stop the duplicate check working."""
    rows = MONDO_ONLY + [dict(MONDO_ONLY[0])]
    with pytest.raises(ValueError, match="more than"):
        load_post_cutoff_approvals(_write(tmp_path, rows), DO)


# --- identifier spelling across graphs ---------------------------------------

from biokgsuite.gold.post_cutoff_approvals import norm_curie, curie_variants  # noqa: E402


def test_norm_curie_strips_the_node_type_prefix():
    assert norm_curie("Disease::DOID:8778") == "DOID8778"
    assert norm_curie("MONDO:0019436") == "MONDO0019436"


def test_bare_mesh_accession_is_offered():
    """BioKG writes diseases as a bare MeSH accession (D005317). Offering only
    the prefixed form matched 0 of 92 MeSH-carrying rows against its 11,347
    disease nodes."""
    assert "D005317" in curie_variants("MESH:D005317")
    assert curie_variants("MESH:D005317")[0] == "MESHD005317"   # prefixed first


def test_a_bare_numeric_form_is_never_offered():
    """DOID:684 and a MONDO node numbered 684 are different diseases. A numeric
    match across vocabularies would bind an indication to the wrong node while
    looking like a clean hit."""
    assert "684" not in curie_variants("DOID:684")
    assert "0019436" not in curie_variants("MONDO:0019436")


def test_variants_are_unique_and_ordered_most_specific_first():
    v = curie_variants("DOID:0080199")
    assert len(v) == len(set(v))
    assert v[0] == "DOID0080199"
