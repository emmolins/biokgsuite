"""The benchmark frame, query construction and the mock ranking loop, on a
synthetic post-cutoff set written to a temporary repo layout."""
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))

from biokgsuite.grounding import benchmark as B  # noqa: E402
from biokgsuite.grounding import runner as R  # noqa: E402

KGS = ['kga', 'kgb']


@pytest.fixture
def base(tmp_path):
    root = tmp_path
    (root / B.SET_DIR).mkdir(parents=True)
    (root / B.BRIDGE_DIR).mkdir(parents=True)
    gold = pd.DataFrame({
        'drugbank_id': ['DB1', 'DB2', 'DB3', 'DB4'], 'mondo_id': ['MONDO:1', 'MONDO:2', 'MONDO:3', 'MONDO:2'],
        'drug_name': ['alpha', 'beta', 'gamma', 'delta'], 'disease_name': ['dis one', 'dis two', 'dis three', 'dis two'],
        'doid': [None] * 4, 'mesh_id': [None] * 4,
        'approval_date': ['2023-02-01', '2024-03-01', '2025-04-01', '2026-05-01'],
        'agencies': ['FDA', 'EMA;FDA', 'PMDA', 'TGA'], 'n_agencies': [1, 2, 1, 1],
        'combination_type': ['single agent'] * 4})
    gold.to_csv(root / B.SET_DIR / 'post_cutoff_approvals.csv', index=False)
    rows = []
    for kg in KGS:
        for db, m in zip(gold.drugbank_id, gold.mondo_id):
            cov = not (kg == 'kgb' and db == 'DB3')          # kgb misses pair 3
            leak = kg == 'kga' and db == 'DB4'               # kga already holds pair 4
            rows.append(dict(kg=kg, drugbank_id=db, mondo_id=m, covered=cov, leaked=leak,
                             evaluable=cov and not leak, evaluable_strict=cov and not leak and db != 'DB2'))
    pd.DataFrame(rows).to_csv(root / B.SET_DIR / 'post_cutoff_pair_status.csv', index=False)
    (root / B.SET_DIR / 'PROVENANCE.txt').write_text('# source_sha256:  abc\n')
    # bridges the classifier and the pool need
    pd.DataFrame({'doid': ['DOID:1'], 'mondo_id': ['MONDO:1'], 'mondo_name': ['dis one']}).to_csv(root / B.BRIDGE_DIR / 'do_diseases.csv', index=False)
    pd.DataFrame({'mesh_id': ['MESH:D1'], 'doid': ['DOID:1'], 'source': ['DO']}).to_csv(root / B.BRIDGE_DIR / 'mesh_to_doid.csv', index=False)
    pd.DataFrame({'child_id': ['MONDO:20', 'MONDO:21'], 'parent_id': ['MONDO:2', 'MONDO:20']}).to_csv(root / B.BRIDGE_DIR / 'mondo_hierarchy.csv', index=False)
    pool = pd.DataFrame({'drugbank_id': [f'DB{i}' for i in range(1, 31)], 'drug_name': [f'drug{i}' for i in range(1, 31)],
                         'chembl_id': [f'CHEMBL{i}' for i in range(1, 31)],
                         'drug_type': ['Small molecule'] * 20 + ['Antibody'] * 10, 'id_source': ['x'] * 30})
    pool.to_csv(root / B.BRIDGE_DIR / 'approved_drug_pool.csv', index=False)
    # DB10 is indicated for dis two; DB11 for a descendant of dis two (MONDO:21)
    pd.DataFrame({'drugbank_id': ['DB10', 'DB11'], 'mondo_id': ['MONDO:2', 'MONDO:21'], 'source_disease_id': ['MONDO_2', 'MONDO_21']}
                 ).to_csv(root / B.BRIDGE_DIR / 'approved_indications.csv', index=False)
    return root


def test_benchmark_flags_and_common(base):
    b = B.load_benchmark(base, kgs=KGS)
    assert list(b.pair_id) == [1, 2, 3, 4]
    assert b.evaluable_kgb.tolist() == [True, True, False, True]
    assert b.evaluable_kga.tolist() == [True, True, True, False]
    assert b.common.tolist() == [True, True, False, False]
    assert b.n_kg_covered.tolist() == [2, 2, 1, 2]
    assert b.approval_year.tolist() == [2023, 2024, 2025, 2026]
    assert b.modality.tolist() == ['small molecule'] * 4
    strict = B.load_benchmark(base, kgs=KGS, strict=True)
    assert strict.common.tolist() == [True, False, False, False]


def test_queries_exclude_indicated_and_are_deterministic(base):
    b = B.load_benchmark(base, kgs=KGS)
    pool = B.load_distractor_pool(base)
    children = B.mondo_children(base)
    q = B.build_queries(b, pool, family='main', pool_size=8, seed=3, children=children)
    assert len(q) == 4 and all(len(x['candidates']) == 8 for x in q)
    assert all(x['candidates'][0]['is_pos'] for x in q)
    for x in q:
        distractors = {c['drug_id'] for c in x['candidates'][1:]}
        assert x['drugbank_id'] not in distractors
        if x['mondo_id'] == 'MONDO:2':
            # indicated for the disease, indicated for its descendant, and the other positive
            assert distractors.isdisjoint({'DB10', 'DB11', 'DB2', 'DB4'})
    again = B.build_queries(b, pool, family='main', pool_size=8, seed=3, children=children)
    assert [x['candidates'] for x in q] == [x['candidates'] for x in again]
    other = B.build_queries(b, pool, family='main', pool_size=8, seed=4, children=children)
    assert [x['candidates'] for x in q] != [x['candidates'] for x in other]
    common = B.build_queries(b, pool, family='common', pool_size=8, seed=3, children=children)
    assert sorted(x['pair_id'] for x in common) == [1, 2]
    with pytest.raises(ValueError):
        B.build_queries(b, pool, family='other', pool_size=8, seed=0)


def test_mock_run_rows_and_idempotence(base, tmp_path):
    b = B.load_benchmark(base, kgs=KGS)
    q = B.build_queries(b, B.load_distractor_pool(base), family='main', pool_size=8, seed=0)
    out = tmp_path / 'runs'
    df, aborted = R.run_family(family='main', queries=q, models=['m1'], kgs=KGS, dossiers={}, mock=True, seed=0,
                               shuffles=2, pool_size=8, anonymize=True, temperature=0.7, api_delay=0, out_dir=out,
                               example_dir=tmp_path / 'ex')
    assert not aborted
    assert list(df.columns) == R.RUN_COLUMNS
    assert len(df) == 4 * 2 * (1 + len(KGS))
    assert set(df.arm) == {'none', *KGS} and set(df.pair_id) == {1, 2, 3, 4}
    assert (df.parsed == 1).all() and df['rank'].between(1, 8).all()
    assert (df[df.arm == 'none'].n_with_evidence == 0).all()
    assert (tmp_path / 'ex' / 'main__kga.txt').exists()
    path = R.run_path(out, 'main', 'm1', 0)
    assert path.exists()
    mtime = path.stat().st_mtime
    df2, _ = R.run_family(family='main', queries=q, models=['m1'], kgs=KGS, dossiers={}, mock=True, seed=0,
                          shuffles=2, pool_size=8, anonymize=True, temperature=0.7, api_delay=0, out_dir=out)
    assert path.stat().st_mtime == mtime and len(df2) == len(df)


def test_shuffle_streams_are_independent_of_order():
    a = R._shuffle_rng(0, 7, 1, 'kga').permutation(8).tolist()
    b = R._shuffle_rng(0, 7, 1, 'kga').permutation(8).tolist()
    c = R._shuffle_rng(0, 7, 1, 'kgb').permutation(8).tolist()
    assert a == b and a != c
