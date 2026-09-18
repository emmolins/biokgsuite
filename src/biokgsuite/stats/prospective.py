"""Prospective-generalisation metrics that differ in what they count.

A prospective pair is *evaluable* for a graph only when the graph holds both the
drug and the disease. Graphs differ enormously in how many of the gold pairs
they can express at all, and the three numbers below answer three different
questions about that. Keeping them apart is the point of this module.

    coverage                how often can the graph answer?
    hits_at_k_all           can it deliver a future approval, over the whole
                            gold set, counting what it cannot express as a miss
    conditional AUROC       when it can answer, is the ranking any good
                            (computed elsewhere, on the evaluable set only)

The distinction that matters: a rank-based metric such as AUROC is defined over
a ranked candidate list, so a pair the graph cannot express has no list and no
defensible value to insert. Imputing zero there asserts that the true positive
ranked last, which is a strong claim rather than a neutral one, and it puts the
graph below AUROC's 0.5 chance floor. A count-based metric such as Hits@K has
the gold pairs as its denominator, so a pair the graph cannot rank is a genuine
miss and needs no imputation at all. Use each for what it can carry.

Comparing conditional AUROCs across graphs is comparing estimates on different
samples, and the samples are not missing at random: an older graph tends to hold
the better-studied drugs, which are also the easier ones. `pairwise_agreement`
is the complete-case answer to that, comparing two graphs only on the pairs both
can evaluate.
"""
from __future__ import annotations

import itertools

import numpy as np
import pandas as pd


def coverage(n_evaluable: int, n_gold: int) -> float:
    """Fraction of the gold set the graph can express at all."""
    if n_gold <= 0:
        return float('nan')
    return n_evaluable / n_gold


def hits_at_k_all(ranks, n_gold: int, k: int) -> float:
    """Hits@K over the WHOLE gold set, not over the evaluable subset.

    ``ranks`` are the ranks of the positives the graph could evaluate, 1-based,
    against whatever candidate set they were scored in. Every gold pair the
    graph could not evaluate is counted as a miss, which is what it is: the
    graph failed to retrieve an approval that happened. Nothing is imputed, so
    this is a description of what occurred rather than a guess at what would
    have.
    """
    if n_gold <= 0:
        return float('nan')
    r = np.asarray(list(ranks), dtype=float)
    r = r[~np.isnan(r)]
    return float((r <= k).sum()) / n_gold


def hits_at_k_evaluable(ranks, k: int) -> float:
    """Hits@K over the evaluable subset, the conditional counterpart above."""
    r = np.asarray(list(ranks), dtype=float)
    r = r[~np.isnan(r)]
    if r.size == 0:
        return float('nan')
    return float((r <= k).mean())


def prospective_summary(ranks_by_kg: dict, n_gold: int, k: int = 10) -> pd.DataFrame:
    """One row per graph: coverage, full-denominator and conditional Hits@K.

    ``ranks_by_kg`` maps a KG name to {pair_key: rank}. Reporting all three
    together is the point; any one of them alone invites the reader to mistake a
    small easy sample for a good graph.
    """
    rows = []
    for kg, ranks in ranks_by_kg.items():
        vals = list(ranks.values())
        rows.append({
            'KG': kg,
            'n_evaluable': len(vals),
            'n_gold': n_gold,
            'coverage': round(coverage(len(vals), n_gold), 4),
            f'hits@{k}_all': round(hits_at_k_all(vals, n_gold, k), 4),
            f'hits@{k}_evaluable': round(hits_at_k_evaluable(vals, k), 4),
        })
    return pd.DataFrame(rows).sort_values(f'hits@{k}_all', ascending=False)


def pairwise_agreement(ranks_by_kg: dict, min_shared: int = 10) -> pd.DataFrame:
    """Compare every pair of graphs on the pairs BOTH can evaluate.

    The global intersection across a whole panel collapses to almost nothing as
    soon as one graph has poor coverage, so it answers nothing. Pairwise
    intersections stay large enough to be informative and are the standard
    complete-case comparison.

    For each unordered pair of graphs this counts, over the shared pairs, how
    often each ranked the true approval better. Ties are excluded from the win
    rate and reported. It assumes nothing about the distribution of the scores,
    which is what makes it readable at this sample size; a win rate near 0.5
    means the two graphs order the shared pairs about equally well.
    """
    rows = []
    for a, b in itertools.combinations(sorted(ranks_by_kg), 2):
        ra, rb = ranks_by_kg[a], ranks_by_kg[b]
        shared = sorted(set(ra) & set(rb))
        a_better = sum(1 for key in shared if ra[key] < rb[key])
        b_better = sum(1 for key in shared if rb[key] < ra[key])
        ties = len(shared) - a_better - b_better
        decided = a_better + b_better
        rows.append({
            'kg_a': a, 'kg_b': b,
            'n_shared': len(shared),
            'a_better': a_better, 'b_better': b_better, 'ties': ties,
            'a_win_rate': round(a_better / decided, 4) if decided else float('nan'),
            'comparable': len(shared) >= min_shared,
        })
    return pd.DataFrame(rows).sort_values('n_shared', ascending=False)
