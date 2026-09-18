#!/usr/bin/env python3
"""Does notebook 05's stability result survive a learned readout?

Notebook 05 measures stability with neighbourhood heuristics, which are a
deterministic function of the graph: the only thing that varies across the sweep
is the graph, which is what makes the standard deviations interpretable.
Practitioners use knowledge-graph embeddings, so the open question is whether
the heuristic conclusion transfers. This answers it at the scored dropout rate.

The design point that matters. Embedding training is stochastic in
initialisation, negative sampling and SGD, so a drop in rank agreement after
retraining mixes real signal loss with run-to-run training variance. This script
measures the variance FLOOR first, by training several seeds on the UNPERTURBED
graph and correlating their rankings with each other. A perturbed result only
means something if it falls below that floor. Without it the numbers are
uninterpretable, which is why this check does not belong inside notebook 05.

What is compared, and why it is comparable: the same therapeutic pairs, drawn by
the same seeded rule from `biokgsuite.stats.stability`, perturbed by the same
triple-level dropout, and summarised by the same Spearman statistic as the
heuristic sweep. Only the scorer changes.

    python supplementary/embeddings/run_stability_embeddings.py --kg hetionet
    python supplementary/embeddings/run_stability_embeddings.py --kg hetionet --smoke

Full runs belong on the GPU; see hpc/README.md. Results land in
results/tables/supp_stability_embeddings.csv, appended per (KG, model) so the
script can be run one graph at a time and resumed.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = next(p for p in Path(__file__).resolve().parents
            if (p / 'config' / 'kgs.yaml').exists())
sys.path.insert(0, str(ROOT / 'src'))

from biokgsuite.io import find_config, load_config, load_kg, require   # noqa: E402
from biokgsuite.stats import stability as stab                         # noqa: E402

OUT = ROOT / 'results' / 'tables' / 'supp_stability_embeddings.csv'


def _indication_relation_index(rel_to_idx, rels):
    """Relation index used to score a (drug, disease) pair, and its inverse."""
    for r in rels:
        if r in rel_to_idx:
            return rel_to_idx[r], rel_to_idx.get(f'{r}_inv')
    return None, None


def rank_pairs(model, pairs, rel_idx, rel_idx_inv):
    """Model score per (drug, disease) pair, max over the relation and its inverse.

    Mirrors `compute_embedding_metrics`: indication edges are stored in either
    direction across the panel, so scoring one direction only would penalise
    whichever graphs chose the other.
    """
    fwd = np.asarray(model.score_pairs(pairs, rel_idx), dtype=float)
    if rel_idx_inv is None:
        return fwd
    rev = np.asarray(model.score_pairs([(t, h) for h, t in pairs], rel_idx_inv),
                     dtype=float)
    return np.maximum(fwd, rev)


def train_and_rank(kg_df, pairs, model_name, dim, epochs, seed, n_entities):
    """Train one model on one edge table and return its ranking of `pairs`.

    The therapeutic pairs are held out of training. Leaving them in would let
    the model read the answer off the edge it is being asked to score, and the
    ranking would measure memorisation rather than structure.
    """
    from biokgsuite.models.embedding import MODELS, build_train_triples

    pair_set = set(pairs)
    rels = REL_CACHE['rels']
    triples, rel_to_idx, _ = build_train_triples(kg_df, pair_set, list(rels))
    rel_idx, rel_idx_inv = _indication_relation_index(rel_to_idx, rels)
    if rel_idx is None:
        raise RuntimeError(
            f'none of the configured drug_disease relations {rels} survive in the '
            f'training triples; nothing to score against')

    cls = MODELS[model_name]
    kw = dict(dim=dim, seed=seed)
    if model_name == 'RotatE':
        kw.update(dim=max(dim // 2, 8), margin=6.0)
    else:
        kw.update(margin=1.0)
    model = cls(n_entities=n_entities, n_relations=len(rel_to_idx), **kw)
    model.fit(triples, n_epochs=epochs, verbose=False)
    return rank_pairs(model, pairs, rel_idx, rel_idx_inv)


REL_CACHE: dict = {}


def run_kg(kg_name, models, dim, epochs, n_trials, rate, seed, sample_size, smoke):
    config = load_config(find_config(ROOT))
    kg_df, nodes_df = load_kg(kg_name, config)
    n_entities = int(max(kg_df['x_index'].max(), kg_df['y_index'].max())) + 1

    import networkx as nx
    from biokgsuite.graph.graph_utils import build_graph
    lcc = max(nx.connected_components(build_graph(kg_df, nodes_df)), key=len)

    REL_CACHE['rels'] = stab.relation_labels(config, kg_name, 'drug_disease')
    all_pairs = stab.therapeutic_pairs(kg_df, config, kg_name, lcc)
    pairs = stab.sample_pairs(all_pairs, sample_size, seed)
    if len(pairs) < 11:
        raise RuntimeError(f'{kg_name}: only {len(pairs)} therapeutic pairs; '
                           f'too few for a rank correlation')
    print(f'{kg_name}: {len(pairs)} of {len(all_pairs):,} therapeutic pairs, '
          f'{len(kg_df):,} triples, {n_entities:,} entities')

    rows = []
    for model_name in models:
        t0 = time.time()

        # 1. Variance floor: same graph, different seeds. Anything the
        #    perturbation does has to beat this to count as signal loss.
        base = [train_and_rank(kg_df, pairs, model_name, dim, epochs, seed + t, n_entities)
                for t in range(n_trials)]
        # drop_zero=False: TransE and RotatE score by negative distance, so the
        # heuristic zero-mask would discard every pair.
        floor = [stab.spearman(base[i], base[j], drop_zero=False)[0]
                 for i in range(len(base)) for j in range(i + 1, len(base))]
        floor_mean = float(np.mean(floor)) if floor else float('nan')
        print(f'  {model_name}: retraining floor r={floor_mean:.4f} '
              f'over {len(floor)} seed pairs ({time.time() - t0:.0f}s)')
        rows.append({'KG': kg_name, 'Model': model_name, 'Comparison': 'retraining floor',
                     'Strategy': '-', 'Dropout': 0.0,
                     'Spearman r': round(floor_mean, 4), 'Spearman SD': round(
                         float(np.std(floor, ddof=1)), 4) if len(floor) > 1 else 0.0,
                     'Trials': n_trials, 'Pairs': len(pairs),
                     'dim': dim, 'epochs': epochs, 'seed': seed})

        # 2. Perturbed: each trial drops the same fraction of triples, then
        #    retrains. Compared against the matching baseline seed, so training
        #    variance enters both sides the same way.
        for strategy in ('random', 'periphery'):
            vals = []
            for t in range(n_trials):
                rng = np.random.RandomState(seed + t)
                sub, n_dropped = stab.drop_rows(kg_df, rate, strategy, rng)
                pert = train_and_rank(sub, pairs, model_name, dim, epochs,
                                      seed + t, n_entities)
                vals.append(stab.spearman(base[t], pert, drop_zero=False)[0])
            m, sd = float(np.mean(vals)), (float(np.std(vals, ddof=1))
                                           if len(vals) > 1 else 0.0)
            print(f'  {model_name} {strategy:>9s} @ {rate:.0%}: r={m:.4f} (sd {sd:.4f})'
                  f'   {"BELOW" if m < floor_mean else "AT OR ABOVE"} the floor')
            rows.append({'KG': kg_name, 'Model': model_name, 'Comparison': 'perturbed',
                         'Strategy': strategy, 'Dropout': rate,
                         'Spearman r': round(m, 4), 'Spearman SD': round(sd, 4),
                         'Trials': n_trials, 'Pairs': len(pairs),
                         'dim': dim, 'epochs': epochs, 'seed': seed})
        print(f'  {model_name}: {time.time() - t0:.0f}s total')
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--kg', required=True, help='knowledge graph name from config/kgs.yaml')
    ap.add_argument('--models', default='TransE,RotatE')
    ap.add_argument('--dim', type=int, default=32)
    ap.add_argument('--epochs', type=int, default=5)
    ap.add_argument('--smoke', action='store_true',
                    help='tiny config to prove the path runs; NOT a result')
    a = ap.parse_args()

    params = load_config(find_config(ROOT))['analysis_params']
    seed = require(params, 'random_seed')
    n_trials = require(params, 'n_trials', 'robustness')
    rate = require(params, 'primary_dropout', 'robustness')
    sample_size = require(params, 'sample_size', 'robustness')

    dim, epochs = a.dim, a.epochs
    if a.smoke:
        dim, epochs, n_trials, sample_size = 8, 1, 2, 60
        print('SMOKE RUN: dim=8 epochs=1 trials=2 pairs=60. Proves the path runs; '
              'the numbers are not a result.')

    df = run_kg(a.kg, a.models.split(','), dim, epochs, n_trials, rate, seed,
                sample_size, a.smoke)
    if a.smoke:
        print('\n' + df.to_string(index=False))
        return
    OUT.parent.mkdir(parents=True, exist_ok=True)
    if OUT.exists():
        prev = pd.read_csv(OUT)
        prev = prev[~((prev.KG == a.kg) & (prev.Model.isin(a.models.split(','))))]
        df = pd.concat([prev, df], ignore_index=True)
    df.to_csv(OUT, index=False)
    print(f'\n-> {OUT.relative_to(ROOT)}')


if __name__ == '__main__':
    main()
