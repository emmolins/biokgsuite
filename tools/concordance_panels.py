"""Kendall's W concordance for the six-graph and MATRIX-excluded panels.

Reproduces the concordance cell of biokgsuite_dimensions/00_summary/notebook.ipynb
from results/tables/00_summary/benchmark_summary.csv and writes
results/tables/00_summary/concordance_specifications.csv with a `panel`
column, so the five-graph statistics quoted in Results point 2 and
Supplementary Table S4d sit in the same file as the six-graph ones.

Usage (from the repo root):
    python tools/concordance_panels.py
"""
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'results' / 'tables' / '00_summary' / 'benchmark_summary.csv'
OUT = ROOT / 'results' / 'tables' / '00_summary' / 'concordance_specifications.csv'

N_PERM = 20_000
DISC_TOL = 0.05
DIMS = ['Coverage', 'Annotation accuracy', 'Trustworthiness', 'Topology',
        'Stability', 'Task performance', 'Generalisation']
# metric -> dimension, in the order of Table 1
METRIC_DIM = {
    'Entity coverage': 'Coverage', 'Relation coverage': 'Coverage',
    'Entity validity': 'Annotation accuracy', 'Relational consistency': 'Annotation accuracy',
    'Edge traceability': 'Trustworthiness', 'Source diversity': 'Trustworthiness',
    'Uncertainty quantification': 'Trustworthiness',
    'Graph cohesion': 'Topology', 'Clustering coefficient': 'Topology',
    'Graph reachability': 'Topology', 'Community purity': 'Topology',
    'Random-dropout stability': 'Stability', 'Periphery-dropout stability': 'Stability',
    'Link prediction': 'Task performance', 'Neighbourhood retrieval': 'Task performance',
    'Multi-hop reasoning': 'Task performance',
    'Data-sparse generalisation': 'Generalisation', 'Cross-domain generalisation': 'Generalisation',
    'Prospective generalisation': 'Generalisation',
}


def rank_matrix(df, cols):
    return df[cols].rank(ascending=False, method='average')


def kendall_w(df, cols):
    R = rank_matrix(df, cols).dropna(axis=1)
    Rt = R.T.to_numpy()
    k, n = Rt.shape
    s = Rt.sum(axis=0)
    W = 12 * ((s - s.mean()) ** 2).sum() / (k ** 2 * (n ** 3 - n))
    rho = (k * W - 1) / (k - 1)
    return W, rho, Rt, k


def perm_p(Rt, W, rng):
    k, n = Rt.shape
    null = np.empty(N_PERM)
    for i in range(N_PERM):
        P = np.apply_along_axis(rng.permutation, 1, Rt)
        s = P.sum(axis=0)
        null[i] = 12 * ((s - s.mean()) ** 2).sum() / (k ** 2 * (n ** 3 - n))
    return float((null >= W).mean()), float(null.mean())


def specifications(df, panel):
    rng = np.random.default_rng(0)
    metrics = [m for m in METRIC_DIM if m in df.columns]
    disc = [m for m in metrics if df[m].notna().sum() >= 2 and (df[m].max() - df[m].min()) > DISC_TOL]
    dimcols = [f'{d} [dim]' for d in DIMS]
    rows = []
    for label, cols in [('All metrics', metrics), ('Discriminating metrics', disc), ('Dimension scores', dimcols)]:
        W, rho, Rt, k = kendall_w(df, cols)
        p, nullW = perm_p(Rt, W, rng)
        rows.append(dict(panel=panel, specification=label, n_judges=k, W=round(W, 3),
                         mean_rho=round(rho, 3), null_W=round(nullW, 3), perm_p=p, note=''))
    by_dim = {}
    for m in disc:
        by_dim.setdefault(METRIC_DIM[m], []).append(m)
    combos = list(product(*by_dim.values()))
    Wc = np.array([kendall_w(df, list(c))[0] for c in combos])
    nd = len(by_dim)
    rows.append(dict(panel=panel, specification='One metric per dimension', n_judges=nd,
                     W=round(float(np.median(Wc)), 3),
                     mean_rho=round(float((nd * np.median(Wc) - 1) / (nd - 1)), 3),
                     null_W=round(1 / nd, 3), perm_p=np.nan,
                     note=f'median of {len(combos)} combinations, range {Wc.min():.3f}-{Wc.max():.3f}'))
    Wj = np.array([kendall_w(df, [x for x in disc if x != m])[0] for m in disc])
    nj = len(disc) - 1
    rows.append(dict(panel=panel, specification='Leave-one-metric-out', n_judges=nj,
                     W=round(float(np.median(Wj)), 3),
                     mean_rho=round(float((nj * np.median(Wj) - 1) / (nj - 1)), 3),
                     null_W=round(1 / nj, 3), perm_p=np.nan,
                     note=f'median over {len(disc)} refits, range {Wj.min():.3f}-{Wj.max():.3f}'))
    return rows


def main():
    df = pd.read_csv(SRC).set_index('KG').dropna(how='all')
    rows = specifications(df, 'all six graphs')
    if 'matrix' in df.index:
        rows += specifications(df.drop('matrix'), 'MATRIX excluded')
    out = pd.DataFrame(rows)
    out.to_csv(OUT, index=False)
    print(out.to_string(index=False))


if __name__ == '__main__':
    main()
