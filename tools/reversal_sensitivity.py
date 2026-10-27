"""Supplementary Table S4d: rank-reversal rates between dimensions under near-tie
restrictions, with and without MATRIX.

Reads results/tables/00_summary/dimension_scores.csv and writes
results/tables/00_summary/reversal_sensitivity.csv. A reversal is a pair of
graphs ordered oppositely on two dimensions; tied comparisons are excluded, and
the gap-restricted rows keep a comparison only when both dimension differences
are at least the stated gap. Annotation accuracy is constant across the five
earlier graphs, so its pairs contribute no untied comparisons once MATRIX is
excluded; the 'discriminating' columns drop it for the six-graph panel too.
"""
import csv, itertools
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'results/tables/00_summary/dimension_scores.csv'
OUT = ROOT / 'results/tables/00_summary/reversal_sensitivity.csv'
GAPS = (0.0, 0.02, 0.05, 0.10)


def rates(rows, dims, gap):
    rev = cmp = 0
    for d1, d2 in itertools.combinations(dims, 2):
        for a, b in itertools.combinations(rows, 2):
            x = a[d1] - b[d1]; y = a[d2] - b[d2]
            if x == 0 or y == 0 or abs(x) < gap or abs(y) < gap:
                continue
            cmp += 1; rev += x * y < 0
    return rev, cmp


def main():
    raw = list(csv.DictReader(open(SRC)))
    dims = [k for k in raw[0] if k != 'KG']
    rows = [{d: float(r[d]) for d in dims} | {'KG': r['KG']} for r in raw]
    panels = {
        'all six graphs': rows,
        'MATRIX excluded': [r for r in rows if r['KG'] != 'matrix'],
    }
    dimsets = {'all dimensions': dims,
               'discriminating dimensions': [d for d in dims if d != 'Annotation accuracy']}
    with open(OUT, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['panel', 'dimension_set', 'min_gap', 'reversals', 'comparisons', 'rate'])
        for pname, prow in panels.items():
            for dname, ds in dimsets.items():
                for g in GAPS:
                    rev, cmp = rates(prow, ds, g)
                    w.writerow([pname, dname, g, rev, cmp, round(rev / cmp, 3) if cmp else ''])
    print('wrote', OUT.relative_to(ROOT))
    for line in open(OUT): print(' ', line.rstrip())


if __name__ == '__main__':
    main()
