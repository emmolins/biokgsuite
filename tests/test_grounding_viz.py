"""Smoke tests for biokgsuite.viz.grounding: every figure function returns a
Figure on the synthetic frames and closes cleanly.

    python -m pytest tests/test_grounding_viz.py -q
"""
import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from biokgsuite.grounding.analysis import (arm_table, capacity_table, chance_mrr,  # noqa: E402
                                           contrast_table, difficulty_table,
                                           scaling_table)
from biokgsuite.viz import grounding as viz  # noqa: E402
from test_grounding_analysis import (ARMS, KGS, LABELS, MODELS, PARAMS_B, POOL,  # noqa: E402
                                     grounded_by_kg, make_benchmark, make_runs)

KG_PAIRS = [('primekg', 'drkg'), ('primekg', 'biokg'), ('biokg', 'drkg')]


@pytest.fixture(scope='module')
def runs():
    return make_runs()


@pytest.fixture(scope='module')
def bench():
    return make_benchmark()


@pytest.fixture(scope='module')
def tables(runs, bench):
    rng = np.random.default_rng(0)
    common = set.intersection(*grounded_by_kg(bench).values())
    cov = runs[runs['pair_id'].isin(common)].assign(family='covered')
    return dict(
        arms=arm_table(runs, 'main', ARMS, 100, rng, LABELS),
        con_main=contrast_table(runs, 'main', KG_PAIRS, 100, rng, LABELS),
        con_cov=contrast_table(cov, 'covered', KG_PAIRS, 100, rng, LABELS),
        cap=capacity_table(runs, PARAMS_B, ARMS, 100, rng, LABELS),
        scaling=scaling_table(runs, PARAMS_B, grounded_by_kg(bench), 100, rng),
        per=difficulty_table(runs, bench),
        n_cov=len(common))


@pytest.fixture(scope='module')
def selection():
    rng = np.random.default_rng(2)
    n = 120
    df = pd.DataFrame(dict(
        area=rng.choice(['oncology', 'neurology', 'cardiology'], size=n),
        modality=rng.choice(['Small molecule', 'Antibody'], size=n),
        approval_year=rng.integers(1980, 2024, size=n),
        kept=rng.random(n) < 0.5))
    return df


def _check(fig):
    assert isinstance(fig, Figure)
    plt.close(fig)


def test_fig_grounding_gain(tables):
    _check(viz.fig_grounding_gain(tables['arms'], tables['con_main'], tables['con_cov'],
                                  MODELS, KGS, 20, tables['n_cov']))


def test_fig_capacity(tables):
    _check(viz.fig_capacity(tables['cap'], chance_mrr(POOL), KGS))


def test_fig_scaling(tables, runs):
    _check(viz.fig_lift_vs_size(tables['scaling']))
    _check(viz.fig_mrr_vs_size(tables['scaling']))
    _check(viz.fig_kg_by_model(runs, PARAMS_B, KGS, LABELS))


def test_fig_rank_distribution_and_evidence(runs):
    _check(viz.fig_rank_distribution(runs, MODELS, ARMS, LABELS))
    _check(viz.fig_rank_distribution(runs, MODELS[:1], ARMS))
    _check(viz.fig_evidence_quality(runs))


def test_fig_difficulty(tables):
    per = tables['per']
    _check(viz.fig_hardest_pairs(per, n=10, labels=LABELS))
    _check(viz.fig_difficulty_drivers(per))
    _check(viz.fig_difficulty_by_area(per, min_n=2))
    _check(viz.fig_intrinsic_difficulty(per, labels=LABELS))
    _check(viz.fig_coverage_fair(per))


def test_fig_bias(selection):
    traits = {'oncology': selection['area'] == 'oncology',
              'antibody': selection['modality'] == 'Antibody',
              'post-2010': selection['approval_year'] > 2010}
    _check(viz.fig_bias_forest(selection, 'kept', traits))
    _check(viz.fig_bias_distribution(selection, 'kept', 'approval_year'))
    _check(viz.fig_bias_composition(selection, 'kept', ['area', 'modality']))


def test_fig_benchmark_overview(bench):
    # The real frame carries approval_date, leaked_<kg> and evaluable_<kg>;
    # the mock has the covered and grounded flags, so derive the rest as the loader does.
    rng = np.random.default_rng(3)
    b = bench.copy()
    b['approval_date'] = pd.to_datetime('2023-01-01') + pd.to_timedelta(rng.integers(0, 900, len(b)), unit='D')
    for kg in KGS:
        b[f'leaked_{kg}'] = b[f'covered_{kg}'] & (rng.random(len(b)) < 0.1)
        b[f'evaluable_{kg}'] = b[f'grounded_{kg}'] & ~b[f'leaked_{kg}']
    dates = {'primekg': '2022-05-01', 'drkg': '2020-04-19', 'biokg': '2020-06-01', 'hetionet': None}
    _check(viz.fig_benchmark_overview(b, dates, KGS))
