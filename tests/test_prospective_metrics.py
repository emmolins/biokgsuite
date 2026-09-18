"""Guards on the three prospective numbers and the pairwise comparison.

The failure these exist to prevent is a small, easy evaluable subset producing a
flattering conditional score that then gets compared against a graph evaluated on
four times as many, harder pairs.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from biokgsuite.stats.prospective import (  # noqa: E402
    coverage, hits_at_k_all, hits_at_k_evaluable, pairwise_agreement,
    prospective_summary,
)

N_GOLD = 128


def test_coverage_is_over_the_whole_gold_set():
    assert coverage(31, N_GOLD) == pytest.approx(31 / 128)
    assert coverage(128, N_GOLD) == 1.0
    assert np.isnan(coverage(5, 0))


def test_hits_at_k_all_counts_the_inexpressible_as_a_miss():
    # A graph that evaluates 31 pairs and ranks every one of them first still
    # only retrieved 31 of 128 approvals.
    assert hits_at_k_all([1] * 31, N_GOLD, 10) == pytest.approx(31 / 128)
    # The same graph looks perfect on the conditional metric.
    assert hits_at_k_evaluable([1] * 31, 10) == 1.0


def test_the_two_denominators_diverge_exactly_by_coverage():
    ranks = [1, 3, 9, 40, 77]
    cond = hits_at_k_evaluable(ranks, 10)
    allv = hits_at_k_all(ranks, N_GOLD, 10)
    assert allv == pytest.approx(cond * len(ranks) / N_GOLD)


def test_nothing_is_imputed_for_pairs_that_were_never_evaluated():
    # Only the evaluated ranks are passed in; the metric never sees a zero, a
    # NaN or a placeholder for the rest.
    assert hits_at_k_all([1, 2], N_GOLD, 10) == pytest.approx(2 / 128)


def test_hits_at_k_respects_k():
    ranks = [1, 11, 50]
    assert hits_at_k_evaluable(ranks, 10) == pytest.approx(1 / 3)
    assert hits_at_k_evaluable(ranks, 50) == 1.0


def test_no_evaluable_pairs_is_nan_not_zero_conditionally():
    assert np.isnan(hits_at_k_evaluable([], 10))
    # But over the whole gold set, retrieving nothing IS zero.
    assert hits_at_k_all([], N_GOLD, 10) == 0.0


def test_summary_reports_all_three_numbers_together():
    df = prospective_summary(
        {'small_easy': {('A', '1'): 1, ('A', '2'): 2},
         'big_mixed': {('A', '1'): 5, ('A', '2'): 30, ('B', '3'): 1, ('B', '4'): 2}},
        n_gold=N_GOLD, k=10)
    assert set(df.columns) == {'KG', 'n_evaluable', 'n_gold', 'coverage',
                               'hits@10_all', 'hits@10_evaluable'}
    small = df.set_index('KG').loc['small_easy']
    big = df.set_index('KG').loc['big_mixed']
    # The small graph wins on the conditional metric and loses on the honest one.
    assert small['hits@10_evaluable'] > big['hits@10_evaluable']
    assert small['hits@10_all'] < big['hits@10_all']


def test_pairwise_compares_only_the_shared_pairs():
    ranks = {
        'a': {('D1', 'X'): 1, ('D2', 'Y'): 2, ('D3', 'Z'): 3},
        'b': {('D1', 'X'): 5, ('D2', 'Y'): 1},          # shares 2, wins 1 of them
    }
    row = pairwise_agreement(ranks, min_shared=1).iloc[0]
    assert row['n_shared'] == 2
    assert row['a_better'] == 1 and row['b_better'] == 1
    assert row['a_win_rate'] == pytest.approx(0.5)
    assert bool(row['comparable']) is True


def test_pairwise_flags_a_comparison_with_too_few_shared_pairs():
    ranks = {'a': {('D1', 'X'): 1}, 'b': {('D1', 'X'): 2, ('D9', 'Q'): 1}}
    row = pairwise_agreement(ranks, min_shared=10).iloc[0]
    assert row['n_shared'] == 1
    assert bool(row['comparable']) is False


def test_pairwise_excludes_ties_from_the_win_rate():
    ranks = {'a': {('D1', 'X'): 4, ('D2', 'Y'): 2},
             'b': {('D1', 'X'): 4, ('D2', 'Y'): 9}}
    row = pairwise_agreement(ranks, min_shared=1).iloc[0]
    assert row['ties'] == 1
    assert row['a_win_rate'] == 1.0          # 1 decided comparison, a won it


def test_pairwise_handles_no_overlap_without_dividing_by_zero():
    ranks = {'a': {('D1', 'X'): 1}, 'b': {('D2', 'Y'): 1}}
    row = pairwise_agreement(ranks, min_shared=1).iloc[0]
    assert row['n_shared'] == 0
    assert np.isnan(row['a_win_rate'])
