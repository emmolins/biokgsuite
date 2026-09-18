"""The disease-clustered bootstrap used by the prospective arm of notebook 07.

Prospective pairs are not independent: several share a disease, and each
positive carries its own generated negatives, so rows within a disease share a
candidate pool, a neighbourhood and a difficulty. Resampling rows treats them
as independent and produces an interval that is too narrow.
"""
import importlib.util
from pathlib import Path

import numpy as np
import pytest
from sklearn.metrics import roc_auc_score

_spec = importlib.util.spec_from_file_location(
    "ev", Path(__file__).resolve().parents[1] / "src/biokgsuite/stats/evaluation.py")
ev = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ev)


def _clustered_data(seed=0, n_disease=14, per=10):
    """Diseases differ in how cleanly their positives separate. Which diseases
    happen to be approved is the uncertainty a cluster bootstrap measures."""
    rng = np.random.RandomState(seed)
    c, s, y = [], [], []
    for d in range(n_disease):
        delta = rng.normal(1.0, 1.4)
        for j in range(per):
            lab = 1 if j < 2 else 0
            c.append(d); y.append(lab)
            s.append(rng.normal(delta if lab else 0.0, 0.4))
    return np.array(s), np.array(y), np.array(c)


def test_point_estimate_is_the_observed_metric():
    s, y, c = _clustered_data()
    p, _, _ = ev.bootstrap_auroc_ci_clustered(s, y, c, n_boot=500,
                                              rng_=np.random.RandomState(1))
    assert p == pytest.approx(roc_auc_score(y, s))


def test_clustering_widens_the_interval():
    s, y, c = _clustered_data()
    _, l1, h1 = ev.bootstrap_auroc_ci(s, y, n_boot=3000,
                                      rng_=np.random.RandomState(1))
    _, l2, h2 = ev.bootstrap_auroc_ci_clustered(s, y, c, n_boot=3000,
                                                rng_=np.random.RandomState(1))
    assert (h2 - l2) > (h1 - l1), (
        "with disease-level difficulty the row bootstrap is overconfident")


def test_single_cluster_gives_no_interval():
    """One cluster cannot support a cluster bootstrap. Return no interval rather
    than silently falling back to the row bootstrap."""
    s, y, c = _clustered_data()
    p, lo, hi = ev.bootstrap_auroc_ci_clustered(s, y, np.zeros_like(c), n_boot=200)
    assert not np.isnan(p) and np.isnan(lo) and np.isnan(hi)


def test_degenerate_labels_return_nan():
    s, y, c = _clustered_data()
    assert np.isnan(ev.bootstrap_auroc_ci_clustered(s, np.ones_like(y), c,
                                                    n_boot=200)[0])


def test_length_mismatch_raises():
    s, y, c = _clustered_data()
    with pytest.raises(ValueError):
        ev.bootstrap_auroc_ci_clustered(s, y, c[:-1])


def test_auprc_variant_runs():
    s, y, c = _clustered_data()
    p, lo, hi = ev.bootstrap_auprc_ci_clustered(s, y, c, n_boot=400,
                                                rng_=np.random.RandomState(2))
    assert 0.0 <= p <= 1.0 and lo <= p <= hi
