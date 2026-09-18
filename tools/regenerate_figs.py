#!/usr/bin/env python3
"""Redraw every figure from cached checkpoints, without re-running the analysis.

WHY THIS EXISTS. Rendering is kept in `biokgsuite.viz.dimension*` rather than
inline in the notebooks for one concrete reason: a journal revision asks for a
font change, a colour change for print, or a different panel order, and every
figure has to be regenerated. Re-running the notebooks to do that means hours of
graph building, a Louvain pass expensive enough to be cached, and two graphs that
exhaust memory. Reading the checkpoints instead takes seconds.

Each notebook writes the distilled figure inputs to
``results/checkpoints/<nb>.pkl``; this reads them and re-renders. If a
checkpoint is absent, that notebook's figures are skipped with a clear message
rather than the whole run failing.

    python tools/regenerate_figs.py                # all
    python tools/regenerate_figs.py --only 03,04   # some
    python tools/regenerate_figs.py --check        # report only, write nothing
"""
from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from biokgsuite.io import Outputs                      # noqa: E402
from biokgsuite.viz.style import setup_style           # noqa: E402

CKPT = ROOT / 'results' / 'checkpoints'


def _load(stem: str):
    p = CKPT / f'{stem}.pkl'
    if not p.exists():
        return None
    with p.open('rb') as fh:
        return pickle.load(fh)


# Each entry: notebook -> (checkpoint stem, renderer). A renderer takes the
# checkpoint and an Outputs, and writes through `out.figure` so regenerated
# figures land exactly where the notebook puts them, under the same names.
def _render_03(ck, out) -> int:
    from biokgsuite.viz.dimension03 import fig03_source_diversity
    fig = _load('03_figdata')
    if not fig:
        print('    03_figdata.pkl missing; run the notebook once to cache it')
        return 0
    out.figure(fig03_source_diversity(fig['kg_counts'], fig['panel_config']),
               'source_diversity')
    return 1


def _render_04(ck, out) -> int:
    from biokgsuite.viz.dimension04 import (fig04_clustering, fig04_cohesion,
                                            fig04_recovery, fig04_resilience)
    import pandas as pd
    n = 0
    if ck.get('conn_records'):
        out.figure(fig04_cohesion(pd.DataFrame(ck['conn_records']).set_index('KG')),
                   'graph_cohesion'); n += 1
    if ck.get('clust_records'):
        out.figure(fig04_clustering(pd.DataFrame(ck['clust_records']).set_index('KG')),
                   'clustering_coefficient'); n += 1
    if ck.get('recovery_results'):
        out.figure(fig04_recovery(ck['recovery_results']), 'hop_recovery'); n += 1
    if ck.get('diversity_scalars'):
        out.figure(fig04_resilience(ck['diversity_scalars']),
                   'differential_resilience'); n += 1
    return n


def _render_01(ck, out) -> int:
    """Coverage heatmaps, rebuilt from the checkpoint's stored fractions."""
    import pandas as pd
    from biokgsuite.viz.primitives import coverage_heatmap
    not_scored = set(ck.get('reported_not_scored', ()))
    kgs = list(ck['entity_cov'])

    ent = pd.DataFrame({kg: {t: (None if v is None else 100 * v)
                             for t, v in types.items()}
                        for kg, types in ck['entity_cov'].items()})[kgs]
    out.figure(coverage_heatmap(
        ent, 'Entity coverage against task-derived reference sets',
        subtitle='share of the reference entities each graph contains',
        not_scored=not_scored), 'entity_coverage')

    rel = pd.DataFrame({kg: {r: d['pct'] for r, d in rels.items()}
                        for kg, rels in ck['relation_coverage'].items()})[kgs]
    out.figure(coverage_heatmap(
        rel, 'Relation coverage against task-derived reference sets',
        subtitle='share of the reference pairs each graph links',
        not_scored=not_scored), 'relation_coverage')
    return 2


def _render_02(ck, out) -> int:
    """Entity validity heatmap, pooled to one row per canonical type."""
    import pandas as pd
    from biokgsuite.stats.annotation import pool_by_canonical
    from biokgsuite.viz.primitives import coverage_heatmap
    pooled = pool_by_canonical(pd.DataFrame(ck['entity_validity']))
    pivot = pooled.pivot(index='Canonical', columns='KG', values='Validity (%)')
    out.figure(coverage_heatmap(
        pivot, 'Entity validity',
        subtitle='share of identifiers whose shape fits a scheme the declared type accepts',
        cbar_label='Entity validity (%)'), 'entity_validity')
    return 1


def _render_05(ck, out) -> int:
    """Stability curves against the degree-matched null, and scorer agreement."""
    import pandas as pd
    from biokgsuite.viz.dimension05 import fig05_scorers, fig05_stability
    if not ck.get('raw_r'):
        return 0
    out.figure(fig05_stability(ck['raw_r'], ck['dropout_rates'],
                               ck.get('raw_r_null')), 'stability')
    n = 1
    if ck.get('stab_records') and ck.get('scorers'):
        out.figure(fig05_scorers(pd.DataFrame(ck['stab_records']),
                                 ck['primary_dropout'], ck['scorers'],
                                 ck['scored_by']), 'scorer_agreement')
        n += 1
    return n


def _render_06(ck, out) -> int:
    """Task-performance figures, from the 06 figure-data cache."""
    from biokgsuite.viz.dimension06 import (fig06_link_heatmap, fig06_link_pr,
                                            fig06_multihop, fig06_retrieval)
    fig = _load('06_figdata')
    if not fig:
        print('    06_figdata.pkl missing; run the notebook once to cache it')
        return 0
    n = 0
    if fig.get('link_records'):
        out.figure(fig06_link_heatmap(fig['link_records'], fig['primary_heuristic']),
                   'link_prediction_heatmap'); n += 1
    if fig.get('pr_curves'):
        out.figure(fig06_link_pr(fig['pr_curves'], fig['anchor_strategy']),
                   'link_prediction_pr'); n += 1
    if fig.get('target_results'):
        out.figure(fig06_retrieval(fig['target_results']), 'neighbourhood_retrieval'); n += 1
    if fig.get('multihop_results'):
        out.figure(fig06_multihop(fig['multihop_results']), 'multihop_reasoning'); n += 1
    return n


RENDERERS = {
    '01_coverage':        ('01_coverage',        _render_01),
    '02_annotation_accuracy': ('02_annotation_accuracy', _render_02),
    '03_trustworthiness': ('03_trustworthiness', _render_03),
    '04_topology':        ('04_topology',        _render_04),
    '05_stability':       ('05_stability',       _render_05),
    '06_task_performance': ('06_task_performance', _render_06),
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--only', help='comma-separated notebook names, e.g. 03,04')
    ap.add_argument('--check', action='store_true',
                    help='report which checkpoints are available; write nothing')
    a = ap.parse_args()

    wanted = list(RENDERERS)
    if a.only:
        keys = [k.strip() for k in a.only.split(',') if k.strip()]
        wanted = [n for n in RENDERERS if any(n.startswith(k) or k in n for k in keys)]
        unknown = [k for k in keys
                   if not any(n.startswith(k) or k in n for n in RENDERERS)]
        if unknown:
            raise SystemExit(f'unknown notebook(s): {unknown}. '
                             f'Known: {sorted(RENDERERS)}')

    setup_style()
    total = 0
    for nb in wanted:
        stem, render = RENDERERS[nb]
        ck = _load(stem)
        if ck is None:
            print(f'  {nb:22s} SKIP  ({stem}.pkl not cached — run the notebook once)')
            continue
        if a.check:
            print(f'  {nb:22s} ready')
            continue
        n = render(ck, Outputs(nb, ROOT))
        print(f'  {nb:22s} {n} figure(s) redrawn')
        total += n
    if not a.check:
        print(f'\n{total} figure(s) regenerated from cached checkpoints.')


if __name__ == '__main__':
    main()
