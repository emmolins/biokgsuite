"""Guards on tools/build_post_cutoff_approvals.py: the funnel accounting, the distance rule and the
per-graph checks, on small synthetic inputs (no graph files needed)."""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tools"))
import build_post_cutoff_approvals as B  # noqa: E402


class FakeBridges:
    """MONDO:1 is_a MONDO:0; MONDO:2 is_a MONDO:0; MONDO:3 is_a MONDO:2 (so 1 and 3 are 3 hops apart)."""
    def __init__(self):
        self.adj = {"MONDO:0": {"MONDO:1", "MONDO:2"}, "MONDO:1": {"MONDO:0"},
                    "MONDO:2": {"MONDO:0", "MONDO:3"}, "MONDO:3": {"MONDO:2"}}
        self.doid_to_mondo = {"DOID:10": {"MONDO:1"}, "DOID:20": {"MONDO:2"}}
        self.mesh_to_mondo, self.mesh_to_doid = {}, {}
        self.mondo_to_doid = {"MONDO:1": {"DOID:10"}, "MONDO:2": {"DOID:20"}}
        self.mondo_to_mesh, self.doid_to_mesh = {}, {}
    within_hops = B.Bridges.within_hops
    to_mondo = B.Bridges.to_mondo
    all_mesh_for = B.Bridges.all_mesh_for


def _source_rows():
    base = dict(in_scope="True", exclusion_reason="", repurposing_class="new_indication", kg_ready="True",
                mondo_id_final="MONDO:1", mondo_label_final="disease one", drug_id_final="DRUGBANK:DB00001",
                drug_label_final="drug a", active_ingredient="", product_name="", drug_class="", disease_status="reviewed_id",
                drug_status="reviewed_id", class_source="pipeline", event_type="", document_url="u", source_url="",
                disease_note="")
    rows = [
        dict(base, candidate_id="FDA:1", regulator="FDA", event_date="2024-03-01"),
        dict(base, candidate_id="EMA:1", regulator="EMA", event_date="2023-11-01"),           # same pair, earlier
        dict(base, candidate_id="TGA:1", regulator="TGA", event_date="2024-01-01", mondo_id_final="MONDO:1; MONDO:2",
             mondo_label_final="disease one; disease two"),                                  # multi-disease row
        dict(base, candidate_id="HC:1", regulator="HC", event_date="2024-01-01", drug_class="vaccine"),
        dict(base, candidate_id="HC:2", regulator="HC", event_date="2024-01-01", drug_id_final="DRUGBANK:DB00001; DRUGBANK:DB00002"),
        dict(base, candidate_id="HC:3", regulator="HC", event_date="2024-01-01", kg_ready="False", repurposing_class="expanded_population"),
        dict(base, candidate_id="PMDA:1", regulator="PMDA", event_date="2024-01-01", disease_note="flag as diagnostic"),
    ]
    return pd.DataFrame(rows)


def test_funnel_accounting_and_collapse():
    pairs, dropped, funnel = B.build_pairs(_source_rows())
    f = {r["stage"]: r for r in funnel}
    assert f["S0"]["n"] == 7 and f["S1"]["n"] == 6 and f["S2"]["n"] == 3
    assert len(pairs) == 2                                  # (DB00001, MONDO:1) and (DB00001, MONDO:2)
    p1 = pairs[pairs.mondo_id == "MONDO:1"].iloc[0]
    assert p1.approval_date == "2023-11-01" and p1.agency == "EMA" and p1.n_agencies == 3
    assert p1.date_type == "CHMP positive opinion"
    assert set(dropped.stage) == {"S1", "S2"} and len(dropped) == 4
    # every drop is accounted for: S0 - S2 == rows dropped at S1 and S2
    assert f["S0"]["n"] - f["S2"]["n"] == len(dropped)


def test_distance_rule():
    br = FakeBridges()
    assert br.within_hops("MONDO:1", "MONDO:1", 2)
    assert br.within_hops("MONDO:1", "MONDO:2", 2)         # siblings
    assert not br.within_hops("MONDO:1", "MONDO:3", 2)     # three steps apart
    assert br.within_hops("MONDO:1", "MONDO:3", 3)


def test_per_graph_checks():
    br = FakeBridges()
    pairs = pd.DataFrame([
        dict(drugbank_id="DB00001", mondo_id="MONDO:1", doid="DOID:10", mesh_id="", approval_date="2024-01-01"),
        dict(drugbank_id="DB00001", mondo_id="MONDO:3", doid="", mesh_id="", approval_date="2024-01-01"),
        dict(drugbank_id="DB00009", mondo_id="MONDO:1", doid="DOID:10", mesh_id="", approval_date="2024-01-01"),
    ])
    view = dict(drug_idx={"DB00001": 1}, dis_idx={"DOID:10": 10, "DOID:20": 20}, degree={1: 9, 10: 3, 20: 3},
                n_mech={1: 4, 10: 2, 20: 0},
                node_curie={1: "DB00001", 10: "DOID:10", 20: "DOID:20"},
                edge_set={(1, 20)}, prior={1: {20}})
    st, funnel = B.check_pairs(pairs, "toy", view, br, {"content_date": "2020-01-01"}, 2)
    r1 = st.iloc[0]
    assert r1.covered and not r1.leaked and r1.date_ok == "True"
    assert r1.related_to_prior and r1.nearest_prior_mondo == "MONDO:2"    # sibling of the prior indication
    assert r1.evaluable and not r1.evaluable_strict
    assert r1.grounded and r1.drug_n_mech == 4 and r1.disease_n_mech == 2
    assert [f["stage"] for f in funnel] == ["S3", "C2", "C1", "C3", "G"]
    assert funnel[-1]["n"] == 1 and funnel[-1]["dropped"] == 0
    # a drug that clears the degree floor on edges a dossier never renders is covered, not grounded;
    # evaluable (notebook 07) is unchanged by it
    blank = dict(view, n_mech={1: 0, 10: 2, 20: 0})
    st_b, funnel_b = B.check_pairs(pairs, "toy", blank, br, {"content_date": "2020-01-01"}, 2)
    assert st_b.iloc[0].covered and st_b.iloc[0].evaluable and not st_b.iloc[0].grounded
    assert funnel_b[-1]["n"] == 0 and funnel_b[-1]["dropped"] == 1
    # a view with no neighbour counts (not a grounding arm) leaves grounded empty and the G row unfilled
    old = dict(view, n_mech={})
    st_o, funnel_o = B.check_pairs(pairs, "toy", old, br, {"content_date": "2020-01-01"}, 2)
    assert st_o.iloc[0].grounded is None and st_o.iloc[0].evaluable
    assert funnel_o[-1]["n"] is None and "not computed" in funnel_o[-1]["reason"]
    assert not st.iloc[1].covered                                         # MONDO:3 has no DOID, graph keys by DOID
    assert not st.iloc[2].drug_present
    # a snapshot with no content date leaves date_ok unknown and nothing evaluable
    st2, _ = B.check_pairs(pairs, "toy", view, br, {}, 2)
    assert st2.iloc[0].date_ok == "unknown" and not st2.evaluable.any()
    # an approval before the content date is not prospective
    st3, _ = B.check_pairs(pairs, "toy", view, br, {"content_date": "2025-01-01"}, 2)
    assert st3.iloc[0].date_ok == "False" and not st3.evaluable.any()


def test_norm_curie():
    assert B.norm_curie("Disease::DOID:14227") == "DOID:14227"
    assert B.norm_curie("D001234") == "MESH:D001234"
    assert B.norm_curie("4975") == "MONDO:0004975"
    assert B.norm_curie("Disease::MESH:C000600608") == "MESH:C000600608"
