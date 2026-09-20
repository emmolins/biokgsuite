"""Unit tests for biokgsuite.grounding.analysis and biokgsuite.grounding.bias.

    python -m pytest tests/test_grounding_analysis.py -q

The run frame is synthetic: three models, four arms, two seeds, two shuffles
and twenty pairs, with the first KG arm drawn to rank the positive first far
more often than the others so its lift is known to be positive.
"""
import numpy as np
import pandas as pd
import pytest

from biokgsuite.grounding.analysis import (absence_penalty_table, fabrication_rate, variance_decomposition, RUN_KEY, arm_table, bh_fdr, boot_ci,
                                           capacity_table, chance_mrr, contrast_table,
                                           cross_model_corr, diagnostics,
                                           difficulty_table, load_runs, per_query,
                                           scaling_table, significance_tests,
                                           validate_runs)
from biokgsuite.grounding.bias import odds_ratio, run_bias_tests

MODELS = ['alpha-1', 'vendor:beta-2', 'gamma:3b']
KGS = ['primekg', 'drkg', 'biokg']
ARMS = ['none'] + KGS
LABELS = {'alpha-1': 'Alpha', 'vendor:beta-2': 'Beta', 'gamma:3b': 'Gamma',
          'primekg': 'PrimeKG', 'drkg': 'DRKG', 'biokg': 'BioKG'}
PARAMS_B = {'alpha-1': 1.0, 'vendor:beta-2': 8.0, 'gamma:3b': 70.0}
N_PAIRS, POOL, SEEDS, SHUFFLES = 20, 8, (0, 1), (0, 1)
AREAS = ['oncology', 'neurology', 'cardiology', 'immunology']


ABSENT_IN = {'drkg': {3, 6, 9}}   # pairs whose positive carries no dossier in that arm


def make_runs(family='main', seed=0, lifted_kg='primekg', p_top=0.7):
    """One row per (model, arm, seed, shuffle, pair). ``lifted_kg`` ranks the
    positive first with probability ``p_top``; every other arm is uniform,
    except that a positive with no dossier (``ABSENT_IN``) is ranked last."""
    rng = np.random.default_rng(seed)
    rows = []
    for model in MODELS:
        for arm in ARMS:
            for s in SEEDS:
                for sh in SHUFFLES:
                    for pid in range(1, N_PAIRS + 1):
                        absent = pid in ABSENT_IN.get(arm, set())
                        if absent:
                            rank = POOL
                        elif arm == lifted_kg and rng.random() < p_top:
                            rank = 1
                        else:
                            rank = int(rng.integers(1, POOL + 1))
                        rows.append(dict(
                            family=family, model=model, arm=arm,
                            condition='no_kg' if arm == 'none' else 'kg',
                            seed=s, shuffle=sh, pair_id=pid,
                            drugbank_id=f'DB{pid:05d}', mondo_id=f'MONDO:{pid:07d}',
                            disease=f'disease {pid}', drug=f'drug {pid}',
                            pool_size=POOL, positive_letter='A', rank=rank,
                            reciprocal_rank=1.0 / rank, hit1=int(rank <= 1),
                            hit3=int(rank <= 3), hit5=int(rank <= 5), parsed=1, error=None,
                            pos_has_evidence=arm != 'none' and not absent,
                            n_with_evidence=0 if arm == 'none' else 3,
                            pos_basis='prior' if arm == 'none' else 'kg',
                            pos_agreement='insufficient' if arm == 'none'
                            else rng.choice(['consistent', 'conflicting', 'insufficient']),
                            pos_confidence=int(rng.integers(1, 6)), pos_rationale='r',
                            frac_basis_kg=0.0 if arm == 'none' else 0.5,
                            mean_confidence=3.0))
    return pd.DataFrame(rows)


def make_benchmark(seed=0):
    rng = np.random.default_rng(seed)
    pids = np.arange(1, N_PAIRS + 1)
    bench = pd.DataFrame(dict(
        pair_id=pids, drugbank_id=[f'DB{p:05d}' for p in pids],
        mondo_id=[f'MONDO:{p:07d}' for p in pids], drug_name=[f'drug {p}' for p in pids],
        disease_name=[f'disease {p}' for p in pids],
        area=[AREAS[p % len(AREAS)] for p in pids],
        modality=rng.choice(['Small molecule', 'Antibody'], size=N_PAIRS),
        approval_year=rng.integers(1980, 2024, size=N_PAIRS)))
    for kg in KGS:
        bench[f'covered_{kg}'] = rng.random(N_PAIRS) < 0.7
        bench[f'grounded_{kg}'] = bench[f'covered_{kg}'] & (rng.random(N_PAIRS) < 0.8)
    bench['n_kg_covered'] = bench[[f'covered_{kg}' for kg in KGS]].sum(axis=1)
    bench['n_kg_grounded'] = bench[[f'grounded_{kg}' for kg in KGS]].sum(axis=1)
    return bench


def grounded_by_kg(bench):
    return {kg: set(bench.loc[bench[f'grounded_{kg}'], 'pair_id']) for kg in KGS}


def test_validate_runs_passes_on_clean_frame():
    df = make_runs()
    assert validate_runs(df) is df


def test_validate_runs_raises_on_duplicate_row():
    df = make_runs()
    dup = pd.concat([df, df.iloc[[5]]], ignore_index=True)
    with pytest.raises(ValueError, match='duplicated'):
        validate_runs(dup)


def test_validate_runs_raises_when_an_arm_misses_a_pair():
    df = make_runs()
    drop = df[(df['model'] == MODELS[0]) & (df['arm'] == 'primekg') & (df['seed'] == 0)
              & (df['shuffle'] == 0) & (df['pair_id'] == 7)].index
    with pytest.raises(ValueError, match='same pairs'):
        validate_runs(df.drop(drop))


def test_load_runs_concatenates_and_filters(tmp_path):
    df = make_runs()
    for (family, model, seed), part in df.groupby(['family', 'model', 'seed']):
        part.to_csv(tmp_path / f'{family}__{model}__seed{seed}.csv', index=False)
    loaded = load_runs(tmp_path)
    assert len(loaded) == len(df)
    assert set(loaded.columns) >= set(RUN_KEY)
    only = load_runs(tmp_path, models=[MODELS[0]])
    assert set(only['model']) == {MODELS[0]}
    with pytest.raises(FileNotFoundError):
        load_runs(tmp_path / 'empty')


def test_chance_mrr():
    assert chance_mrr(8) == pytest.approx(sum(1 / k for k in range(1, 9)) / 8)


def test_per_query_shape():
    pq = per_query(make_runs())
    assert len(pq) == len(MODELS) * N_PAIRS
    assert set(pq.columns) == set(ARMS)
    assert pq.index.names == ['model', 'pair_id']


def test_boot_ci_brackets_the_mean():
    rng = np.random.default_rng(1)
    v = rng.normal(size=50)
    lo, hi = boot_ci(v, 500, rng)
    assert lo <= v.mean() <= hi
    assert all(np.isnan(boot_ci([1.0], 100, rng)))


def test_arm_table_lift_sign_and_interval():
    df = make_runs()
    rng = np.random.default_rng(0)
    arms = arm_table(df, 'main', ARMS, 300, rng, LABELS)
    assert set(arms['family']) == {'main'}
    assert len(arms) == len(MODELS) * len(ARMS)
    prime = arms[arms['arm'] == 'primekg']
    assert (prime['lift'] > 0).all()
    assert (prime['lift_lo'] > 0).all()
    kg = arms[arms['arm'] != 'none']
    assert ((kg['lift_lo'] <= kg['lift']) & (kg['lift'] <= kg['lift_hi'])).all()
    assert arms['lift'].isna().equals(arms['arm'] == 'none')
    assert set(arms['arm_label']) == {'No KG', 'PrimeKG', 'DRKG', 'BioKG'}
    assert set(arms['model_label']) == {'Alpha', 'Beta', 'Gamma'}
    with pytest.raises(ValueError):
        arm_table(df, 'absent', ARMS, 10, rng, LABELS)


def test_contrast_table_separates_flag():
    df = make_runs()
    rng = np.random.default_rng(0)
    pairs = [('primekg', 'drkg'), ('primekg', 'biokg'), ('biokg', 'drkg')]
    con = contrast_table(df, 'main', pairs, 300, rng, LABELS)
    assert len(con) == len(MODELS) * len(pairs)
    expected = (con['lo'] > 0) | (con['hi'] < 0)
    assert con['separates'].equals(expected)
    assert (con.loc[con['kg_a'] == 'primekg', 'separates']).all()
    assert 'PrimeKG - DRKG' in set(con['contrast'])


def test_bh_fdr_known_vector():
    q = bh_fdr([0.01, 0.04, 0.03, 0.2])
    assert q == pytest.approx([0.04, 0.16 / 3, 0.16 / 3, 0.2])
    assert bh_fdr([]).size == 0


def test_significance_tests_rows_and_fdr_monotone():
    df = make_runs()
    bench = make_benchmark()
    common = set.intersection(*grounded_by_kg(bench).values())
    subsets = {'all pairs': df, 'evaluable in every graph': df[df['pair_id'].isin(common)]}
    pairs = [('primekg', 'drkg'), ('primekg', 'biokg'), ('biokg', 'drkg')]
    tests = significance_tests(subsets, pairs, LABELS)
    assert len(tests) == len(subsets) * len(MODELS) * len(pairs)
    assert not tests.duplicated(['subset', 'model', 'contrast']).any()
    assert (tests['q_bh'] >= tests['p_wilcoxon'] - 1e-12).all()
    assert (tests['q_bh'] <= 1).all()
    assert tests['sig_fdr'].equals(tests['q_bh'] < 0.05)


def test_diagnostics_no_correction_when_everything_parses():
    df = make_runs()
    diag = diagnostics(df, POOL)
    assert len(diag) == len(MODELS) * len(ARMS)
    assert (diag['parse_rate'] == 1).all()
    assert (diag['n_failed'] == 0).all()
    assert (diag['n_error'] == 0).all()
    assert diag['shift'].abs().max() == 0
    broken = df.copy()
    idx = broken.index[:10]
    broken.loc[idx, ['parsed', 'rank', 'reciprocal_rank']] = [0, POOL + 1, 1 / (POOL + 1)]
    broken.loc[idx[:3], 'error'] = 'timeout'
    d2 = diagnostics(broken, POOL)
    assert d2['n_failed'].sum() == 10
    assert d2['n_error'].sum() == 3
    assert (d2.loc[d2['n_failed'] > 0, 'shift'] > 0).all()


def test_capacity_table_sorted_by_params():
    df = make_runs(family='capacity')
    rng = np.random.default_rng(0)
    cap = capacity_table(df, PARAMS_B, ARMS, 200, rng, LABELS)
    assert cap['params_b'].is_monotonic_increasing
    assert len(cap) == len(MODELS) * len(ARMS)
    assert list(cap['arm'][:len(ARMS)]) == ARMS
    subset = capacity_table(df, {MODELS[2]: 70.0}, ARMS, 50, rng, LABELS)
    assert set(subset['model']) == {MODELS[2]}


def test_scaling_table_has_pooled_and_grounded_lift():
    df = make_runs()
    bench = make_benchmark()
    rng = np.random.default_rng(0)
    sc = scaling_table(df, PARAMS_B, grounded_by_kg(bench), 200, rng)
    assert list(sc['model']) == sorted(PARAMS_B, key=PARAMS_B.get)
    for col in ('lift_pooled', 'lift_pooled_lo', 'lift_pooled_hi', 'lift_grounded', 'lift_grounded_lo',
                'lift_grounded_hi', 'mrr_none', 'mrr_kg', 'mrr_kg_grounded'):
        assert col in sc.columns
    assert (sc['n_pairs'] == N_PAIRS).all()
    assert (sc['n_pairs_grounded'] <= N_PAIRS).all()
    assert sc[['lift_pooled', 'lift_grounded']].notna().all().all()


def test_absence_penalty_reads_what_the_model_saw():
    df = make_runs()
    rng = np.random.default_rng(0)
    ap = absence_penalty_table(df, 200, rng, LABELS)
    assert set(ap['positive']) == {'shown', 'absent'}
    assert set(ap['arm']) == set(ARMS) - {'none'}
    drkg = ap[ap['arm'] == 'drkg'].set_index(['model', 'positive'])
    other = ap[ap['arm'] != 'drkg']
    # Only DRKG has blank positives, and it ranked them last: a clear negative lift.
    assert (other.loc[other['positive'] == 'absent', 'n_pairs'] == 0).all()
    assert (drkg.xs('absent', level='positive')['n_pairs'] == len(ABSENT_IN['drkg'])).all()
    assert (drkg.xs('absent', level='positive')['lift_hi'] < 0).all()
    assert (drkg.xs('shown', level='positive')['n_pairs'] == N_PAIRS - len(ABSENT_IN['drkg'])).all()


def test_difficulty_table_one_row_per_pair():
    df = make_runs()
    bench = make_benchmark()
    per = difficulty_table(df, bench)
    assert len(per) == N_PAIRS
    assert per['pair_id'].is_unique
    assert per['mrr'].is_monotonic_increasing
    for col in ('mrr', 'mrr_none', 'mrr_grounded', 'hit1', 'n_kg_covered', 'n_kg_grounded',
                'pos_confidence', 'frac_consistent', 'area', 'modality', 'approval_year',
                'drug_name', 'disease_name', 'model_spread'):
        assert col in per.columns
    expected = bench.set_index('pair_id')[[f'grounded_{k}' for k in KGS]].sum(axis=1)
    assert per.set_index('pair_id')['n_kg_grounded'].astype(int).equals(expected.loc[per['pair_id']])
    assert per.loc[per['n_kg_grounded'] == 0, 'mrr_grounded'].isna().all()
    assert per.loc[per['n_kg_grounded'] > 0, 'mrr_grounded'].notna().all()
    model_cols = [c for c in per.columns if c.startswith('mrr:')]
    assert sorted(c[4:] for c in model_cols) == sorted(MODELS)
    assert (per['n_runs'] == len(KGS) * len(MODELS) * len(SEEDS) * len(SHUFFLES)).all()


def test_cross_model_corr_is_square_over_models():
    per = difficulty_table(make_runs(), make_benchmark())
    corr = cross_model_corr(per)
    assert corr.shape == (len(MODELS), len(MODELS))
    assert sorted(corr.index) == sorted(MODELS)
    assert np.allclose(np.diag(corr.to_numpy()), 1.0)


def test_odds_ratio_and_bias_tests():
    rng = np.random.default_rng(0)
    n = 200
    df = pd.DataFrame(dict(
        area=rng.choice(['a', 'b', 'c'], size=n),
        year=rng.integers(1990, 2020, size=n)))
    # Rows in area 'a' are kept far more often, so its odds ratio is above one.
    df['kept'] = np.where(df['area'] == 'a', rng.random(n) < 0.9, rng.random(n) < 0.3)
    orr, lo, hi, p = odds_ratio(df['area'] == 'a', df['kept'])
    assert lo < orr < hi
    assert orr > 1 and p < 0.05
    res = run_bias_tests(df, 'kept', ['area'], ['year'], min_n=8)
    assert set(res['kind']) == {'categorical', 'continuous'}
    assert (res['covariate'] == 'area').sum() == 3
    cont = res[res['kind'] == 'continuous'].iloc[0]
    assert cont['effect'] == pytest.approx(cont['median_kept'] - cont['median_dropped'])
    assert res.loc[res['trait'] == 'area = a', 'effect'].iloc[0] == pytest.approx(orr)


def test_variance_decomposition_shares_sum_to_one_and_find_the_lifted_arm():
    df = make_runs(lifted_kg='primekg', p_top=0.9)
    rng = np.random.default_rng(0)
    vd = variance_decomposition(df, 'main', ARMS, 100, rng).set_index('quantity')
    shares = vd.loc[['share_arm', 'share_model', 'share_interaction'], 'value']
    assert abs(shares.sum() - 1) < 1e-9 and (shares >= 0).all()
    # Every model is uniform except on the lifted arm, so the arm lever dominates.
    assert vd.loc['share_arm', 'value'] > vd.loc['share_model', 'value']
    assert vd.loc['ratio', 'value'] > 1
    assert (vd['lo'] <= vd['value']).all() and (vd['value'] <= vd['hi']).all()
    with pytest.raises(ValueError):
        variance_decomposition(df[df.model == MODELS[0]], 'main', ARMS, 10, rng)


def test_fabrication_rate_reads_the_none_arm_only():
    df = make_runs()
    df.loc[(df.arm == 'none') & (df.pair_id <= 2), 'pos_basis'] = 'both'
    fr = fabrication_rate(df)
    assert set(fr.model) == set(MODELS) and (fr.n == N_PAIRS * len(SEEDS) * len(SHUFFLES)).all()
    assert np.allclose(fr.claimed_kg, 2 / N_PAIRS)
