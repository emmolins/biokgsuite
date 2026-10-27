"""Supplementary Table S8d: pairwise contrasts between graphs on link prediction
and neighbourhood retrieval, and which reversals between the two tasks survive
multiple-testing correction.

Reads the scored specification of each task from
results/tables/06_task_performance/ (Adamic-Adar with type-constrained
negatives; Recall@K with K = 0.5% of the candidate pool) and writes
results/tables/06_task_performance/task_contrasts.csv. Each graph is scored on
its own held-out set, so a contrast is a two-sample z-test with each standard
error taken from the 1,000-replicate bootstrap interval ((hi - lo) / 3.92).
q values are Benjamini-Hochberg over the 15 graph pairs within each task.
"""
import csv, itertools, math
from pathlib import Path
from statistics import NormalDist

ROOT = Path(__file__).resolve().parents[1]
T = ROOT / 'results/tables/06_task_performance'
OUT = T / 'task_contrasts.csv'
N = NormalDist()


def bh(ps):
    m = len(ps); order = sorted(range(m), key=lambda i: ps[i]); q = [0.0] * m; prev = 1.0
    for rank, i in reversed(list(enumerate(order, 1))):
        prev = min(prev, ps[i] * m / rank); q[i] = prev
    return q


def main():
    lp = {r['kg']: (float(r['auroc']), (float(r['auroc_ci_hi']) - float(r['auroc_ci_lo'])) / 3.92)
          for r in csv.DictReader(open(T / 'link_prediction.csv'))
          if r['heuristic'] == 'Adamic-Adar' and r['strategy'] == 'type-constrained'}
    rc = {r['KG']: (float(r['recall@frac']), (float(r['recall@frac_ci_hi']) - float(r['recall@frac_ci_lo'])) / 3.92)
          for r in csv.DictReader(open(T / 'neighbourhood_retrieval.csv'))}
    kgs = [k for k in lp if k in rc]
    pairs = list(itertools.combinations(kgs, 2))

    def contrasts(d):
        ps = []
        for a, b in pairs:
            se = math.hypot(d[a][1], d[b][1])
            z = (d[a][0] - d[b][0]) / se if se > 0 else math.inf
            ps.append(2 * (1 - N.cdf(abs(z))))
        return ps, bh(ps)

    lp_p, lp_q = contrasts(lp); rc_p, rc_q = contrasts(rc)
    with open(OUT, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['kg_a', 'kg_b', 'lp_diff', 'lp_p', 'lp_q', 'retrieval_diff', 'retrieval_p', 'retrieval_q',
                    'opposite_order', 'significant_reversal'])
        for (a, b), p1, q1, p2, q2 in zip(pairs, lp_p, lp_q, rc_p, rc_q):
            dl = lp[a][0] - lp[b][0]; dr = rc[a][0] - rc[b][0]
            opp = dl * dr < 0
            w.writerow([a, b, round(dl, 4), f'{p1:.3g}', f'{q1:.3g}', round(dr, 4), f'{p2:.3g}', f'{q2:.3g}',
                        int(opp), int(opp and q1 < 0.05 and q2 < 0.05)])
    print('wrote', OUT.relative_to(ROOT))
    for line in open(OUT): print(' ', line.rstrip())


if __name__ == '__main__':
    main()
