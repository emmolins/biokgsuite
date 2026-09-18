"""Evaluation metrics for link prediction and ranking tasks.

References:
    Hanley, J.A. & McNeil, B.J. "The meaning and use of the area under
    a receiver operating characteristic (ROC) curve."
    Radiology 143(1), 29–36 (1982).

    Efron, B. & Tibshirani, R.J. "An Introduction to the Bootstrap."
    Chapman & Hall/CRC (1993).
"""

import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve


def compute_metrics(scores, labels):
    """Compute AUROC, AUPRC, MRR, and Hits@10.

    Parameters
    ----------
    scores : array-like of float
        Predicted scores (higher = more likely positive).
    labels : array-like of int
        Binary labels (1 = positive, 0 = negative).

    Returns
    -------
    dict with keys: auroc, auprc, mrr, hits@10
    """
    scores = np.array(scores, dtype=float)
    labels = np.array(labels, dtype=int)

    auroc = roc_auc_score(labels, scores)
    auprc = average_precision_score(labels, scores)

    # MRR and Hits@10: for each positive, find its rank among all scores
    pos_scores = scores[labels == 1]
    reciprocal_ranks, hits = [], []
    for ps in pos_scores:
        rank = int(np.sum(scores >= ps))      # 1-indexed; ties broken conservatively
        reciprocal_ranks.append(1.0 / rank if rank > 0 else 0.0)
        hits.append(1 if rank <= 10 else 0)

    return {
        'auroc':   float(auroc),
        'auprc':   float(auprc),
        'mrr':     float(np.mean(reciprocal_ranks)),
        'hits@10': float(np.mean(hits)),
    }


# ── Bootstrap confidence intervals ────────────────────────────────────────────
# Unified bootstrap CI functions used by notebooks 06, 07, and 08.
# Replaces per-notebook local definitions with a single shared implementation.

def _bootstrap_resample(scores_arr, labels_arr, metric_fn, n_boot, rng,
                        stratified=True, min_valid=10):
    """Core bootstrap resampling loop.

    Parameters
    ----------
    scores_arr, labels_arr : np.ndarray
        Prediction scores and binary labels.
    metric_fn : callable
        Function(scores, labels) → float.
    n_boot : int
        Number of bootstrap resamples.
    rng : np.random.RandomState
        RNG for reproducibility.
    stratified : bool
        If True, resample positives and negatives separately to preserve
        class balance (recommended for imbalanced datasets).
    min_valid : int
        Minimum number of valid bootstrap samples required.

    Returns
    -------
    list of float
        Bootstrap metric estimates, or empty list if insufficient data.
    """
    n = len(scores_arr)
    if labels_arr.sum() < 2 or (1 - labels_arr).sum() < 2:
        return []

    pos_idx = np.where(labels_arr == 1)[0]
    neg_idx = np.where(labels_arr == 0)[0]

    boots = []
    for _ in range(n_boot):
        if stratified and len(pos_idx) > 0 and len(neg_idx) > 0:
            b_pos = pos_idx[rng.randint(0, len(pos_idx), len(pos_idx))]
            b_neg = neg_idx[rng.randint(0, len(neg_idx), len(neg_idx))]
            idx = np.concatenate([b_pos, b_neg])
        else:
            idx = rng.randint(0, n, n)

        b_labels = labels_arr[idx]
        if b_labels.sum() == 0 or b_labels.sum() == len(b_labels):
            continue
        try:
            boots.append(metric_fn(labels_arr[idx], scores_arr[idx]))
        except Exception:
            continue

    return boots


def bootstrap_metric_ci(scores_arr, labels_arr, metric_fn, n_boot=1000,
                        rng=None, stratified=True, alpha=0.05):
    """Percentile bootstrap 95% CI for any metric.

    Parameters
    ----------
    scores_arr : array-like of float
        Prediction scores.
    labels_arr : array-like of int
        Binary labels (1 = positive, 0 = negative).
    metric_fn : callable
        Function(labels, scores) → float (sklearn convention).
    n_boot : int
        Number of bootstrap resamples.
    rng : np.random.RandomState, optional
        RNG for reproducibility. If None, creates one with seed 42.
    stratified : bool
        Stratified resampling to preserve class balance.
    alpha : float
        Significance level (default 0.05 → 95% CI).

    Returns
    -------
    tuple of (float, float, float)
        (point_estimate, ci_lo, ci_hi) or (nan, nan, nan) if insufficient data.
    """
    scores_arr = np.asarray(scores_arr, dtype=float)
    labels_arr = np.asarray(labels_arr, dtype=int)

    if rng is None:
        rng = np.random.RandomState(42)

    boots = _bootstrap_resample(scores_arr, labels_arr, metric_fn,
                                n_boot, rng, stratified)
    if len(boots) < 10:
        return np.nan, np.nan, np.nan

    lo, hi = np.percentile(boots, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(np.mean(boots)), float(lo), float(hi)


def bootstrap_auroc_ci(scores_arr, labels_arr, n_boot=1000, rng_=None,
                       stratified=True):
    """Percentile bootstrap 95% CI for AUROC.

    Returns (point_estimate, ci_lo, ci_hi). Returns (nan, nan, nan) if there
    are fewer than 2 positives or 2 negatives, or if fewer than 10 valid
    bootstrap samples are generated.

    Parameters
    ----------
    scores_arr : array-like of float
        Prediction scores.
    labels_arr : array-like of int
        Binary labels (1 = positive, 0 = negative).
    n_boot : int
        Number of bootstrap samples.
    rng_ : np.random.RandomState, optional
        RNG for reproducibility. If None, creates one with seed 42.
    stratified : bool
        Stratified resampling to preserve class balance.

    Returns
    -------
    tuple of (float, float, float)
        (point_estimate, ci_lo, ci_hi) or (nan, nan, nan) if insufficient data.
    """
    return bootstrap_metric_ci(scores_arr, labels_arr, roc_auc_score,
                               n_boot=n_boot, rng=rng_, stratified=stratified)


def bootstrap_auprc_ci(scores_arr, labels_arr, n_boot=1000, rng_=None,
                       stratified=True):
    """Percentile bootstrap 95% CI for AUPRC (Average Precision).

    Returns (point_estimate, ci_lo, ci_hi). Returns (nan, nan, nan) if there
    are fewer than 2 positives or 2 negatives, or if fewer than 10 valid
    bootstrap samples are generated.

    Parameters
    ----------
    scores_arr : array-like of float
        Prediction scores.
    labels_arr : array-like of int
        Binary labels (1 = positive, 0 = negative).
    n_boot : int
        Number of bootstrap samples.
    rng_ : np.random.RandomState, optional
        RNG for reproducibility. If None, creates one with seed 42.
    stratified : bool
        Stratified resampling to preserve class balance.

    Returns
    -------
    tuple of (float, float, float)
        (point_estimate, ci_lo, ci_hi) or (nan, nan, nan) if insufficient data.
    """
    return bootstrap_metric_ci(scores_arr, labels_arr, average_precision_score,
                               n_boot=n_boot, rng=rng_, stratified=stratified)


def bootstrap_metric_ci_clustered(scores_arr, labels_arr, cluster_ids,
                                  metric_fn, n_boot=1000, rng=None, alpha=0.05):
    """Percentile bootstrap CI that resamples CLUSTERS, not individual rows.

    Why this exists
    ---------------
    The post-cutoff set contains several drug-disease pairs that share a
    disease, and each positive carries its own generated negatives. Rows within
    a disease are therefore not independent: they share a candidate pool, a
    neighbourhood in the graph, and a difficulty. Resampling rows treats them as
    independent and produces an interval that is too narrow, which is the
    failure mode that makes two graphs look separated when they are not.

    Resampling whole diseases with replacement keeps the within-cluster
    correlation intact, so the interval reflects how much the estimate would
    move if a different set of diseases had been approved.

    Parameters
    ----------
    scores_arr, labels_arr : array-like
        Scores and binary labels, one entry per scored candidate.
    cluster_ids : array-like
        Cluster label per entry, normally the disease node. Entries sharing a
        label are resampled together.
    metric_fn : callable
        metric_fn(labels, scores) -> float, e.g. roc_auc_score.

    Note on the point estimate
    --------------------------
    This returns the OBSERVED metric, which is what compute_metrics reports and
    therefore what the notebook prints beside the interval. The row-resampling
    bootstrap_metric_ci above returns the MEAN OF THE REPLICATES instead, which
    is a different quantity and is biased for a rank metric. The notebook
    discards that value, so the two have never disagreed in output, but they are
    not interchangeable if a caller starts using the first return value.

    Returns
    -------
    (point_estimate, ci_lo, ci_hi), or (nan, nan, nan) when the data cannot
    support an estimate (fewer than 2 clusters, or a degenerate label set).
    """
    import numpy as _np

    s = _np.asarray(scores_arr, dtype=float)
    y = _np.asarray(labels_arr, dtype=int)
    c = _np.asarray(cluster_ids)
    if s.size != y.size or s.size != c.size:
        raise ValueError("scores, labels and cluster_ids must be the same length")

    finite = _np.isfinite(s)
    s, y, c = s[finite], y[finite], c[finite]
    if y.sum() < 2 or (1 - y).sum() < 2:
        return (float("nan"),) * 3

    uniq = _np.unique(c)
    if uniq.size < 2:
        # One cluster: a cluster bootstrap has nothing to resample over. Say so
        # by returning nan rather than silently falling back to the row
        # bootstrap, which would report an interval the design cannot support.
        return (float(metric_fn(y, s)), float("nan"), float("nan"))

    idx_by_cluster = {u: _np.flatnonzero(c == u) for u in uniq}
    point = float(metric_fn(y, s))

    rng = rng if rng is not None else _np.random.RandomState(42)
    stats = []
    for _ in range(n_boot):
        pick = rng.choice(uniq.size, uniq.size, replace=True)
        idx = _np.concatenate([idx_by_cluster[uniq[p]] for p in pick])
        yy, ss = y[idx], s[idx]
        if yy.sum() < 1 or (1 - yy).sum() < 1:
            continue
        try:
            stats.append(float(metric_fn(yy, ss)))
        except ValueError:
            continue
    if len(stats) < 10:
        return (point, float("nan"), float("nan"))
    lo = float(_np.percentile(stats, 100 * alpha / 2))
    hi = float(_np.percentile(stats, 100 * (1 - alpha / 2)))
    return (point, lo, hi)


def bootstrap_auroc_ci_clustered(scores_arr, labels_arr, cluster_ids,
                                 n_boot=1000, rng_=None):
    """Cluster bootstrap CI for AUROC. See bootstrap_metric_ci_clustered."""
    from sklearn.metrics import roc_auc_score as _roc
    return bootstrap_metric_ci_clustered(scores_arr, labels_arr, cluster_ids,
                                         _roc, n_boot=n_boot, rng=rng_)


def bootstrap_auprc_ci_clustered(scores_arr, labels_arr, cluster_ids,
                                 n_boot=1000, rng_=None):
    """Cluster bootstrap CI for AUPRC. See bootstrap_metric_ci_clustered."""
    from sklearn.metrics import average_precision_score as _ap
    return bootstrap_metric_ci_clustered(scores_arr, labels_arr, cluster_ids,
                                         _ap, n_boot=n_boot, rng=rng_)
