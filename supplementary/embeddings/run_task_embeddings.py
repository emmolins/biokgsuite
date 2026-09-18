#!/usr/bin/env python3
"""Embedding cross-check for notebook 06's OTHER two tasks.

`run_link_prediction_embeddings.py` established that the link-prediction
ranking survives a scorer from outside the shared-neighbour family. Retrieval
and multi-hop were still Adamic-Adar alone, so the same objection applied to
two thirds of the dimension. This script closes that.

Tasks
-----
retrieval : given a disease, rank every gene. Adamic-Adar scores by shared
            neighbours; the embedding scores the (disease, r, gene) triple
            directly, taking the max over whichever relations actually carry
            disease-gene edges in that graph (they are selected by entity type
            in 06, so the relation vocabulary differs per KG and comes from the
            checkpoint).
multihop  : given (drug, ?, disease), rank the intermediate gene. Adamic-Adar
            uses AA(drug,g) * AA(g,disease). The embedding analogue composes
            the two hops additively, score(drug,r1,g) + score(g,r2,disease),
            because these models score by NEGATIVE distance, so a sum in score
            space is a product in likelihood space.

Leakage
-------
06 hides the answer edge at scoring time and puts it back afterwards. A trained
model cannot do that: by scoring time it has already seen the edge. So for each
task the held-out pairs are removed from the training triples BEFORE training,
which is the stricter of the two and never the looser.

Epochs are calibrated per task on a validation split carved from pairs that are
NOT in the evaluation set. The evaluation set never influences selection.

Usage
-----
    python supplementary/embeddings/run_task_embeddings.py \
        --task retrieval --models transe,rotate --epochs 1,2,5,10,20,40,80
"""
import argparse
import json
import pickle
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))

from biokgsuite.io.config import load_config                     # noqa: E402
from biokgsuite.io.loading import load_kg                        # noqa: E402
from biokgsuite.models.embedding import MODELS, kg_to_triples    # noqa: E402

CKPT = ROOT / 'results' / 'checkpoints' / '06_task_performance.pkl'
OUT_DIR = ROOT / 'results' / 'tables' / '06_task_performance'
ALIAS = {'transe': 'TransE', 'rotate': 'RotatE', 'gemma': 'Gemma'}


def drop_pairs(kg_df, hidden):
    """Remove every row whose unordered index pair is hidden.

    Relation-agnostic on purpose: 06 selects disease-gene edges by entity type,
    so the same pair can be asserted under several relation names and hiding
    only one of them would leave the answer visible under another.
    """
    x = kg_df['x_index'].to_numpy()
    y = kg_df['y_index'].to_numpy()
    lo = np.minimum(x, y)
    hi = np.maximum(x, y)
    keep = ~pd.MultiIndex.from_arrays([lo, hi]).isin(hidden)
    return kg_df[keep]


def rel_indices(rel_to_idx, names):
    """Forward and inverse relation indices for a set of relation names."""
    out = []
    for n in names:
        if n in rel_to_idx:
            out.append(rel_to_idx[n])
        if f'{n}_inv' in rel_to_idx:
            out.append(rel_to_idx[f'{n}_inv'])
    return sorted(set(out))


def score_matrix(model, head, rels, tails):
    """Best score from `head` to every tail, maximised over `rels`."""
    pairs = np.column_stack([np.full(len(tails), head, dtype=np.int64), tails])
    best = None
    for r in rels:
        s = model.score_pairs(pairs, r)
        best = s if best is None else np.maximum(best, s)
    return best


def eval_retrieval(model, queries, genes, rels, k_fixed, k_norm):
    gene_arr = np.asarray(genes, dtype=np.int64)
    pos = {int(g): j for j, g in enumerate(gene_arr)}
    r100, rfrac, rr = [], [], []
    for disease, held in queries:
        if int(held) not in pos:
            continue
        s = score_matrix(model, int(disease), rels, gene_arr)
        t = s[pos[int(held)]]
        rank = int(np.sum(s >= t))          # same conservative tie rule as 06
        r100.append(1 if rank <= k_fixed else 0)
        rfrac.append(1 if rank <= k_norm else 0)
        rr.append(1.0 / rank)
    if not r100:
        return {}
    return {'recall@100': float(np.mean(r100)),
            'recall@frac': float(np.mean(rfrac)),
            'mrr': float(np.mean(rr)), 'n': len(r100)}


def eval_multihop(model, triples, genes, dg_rels, gd_rels, k_fixed, k_norm):
    gene_arr = np.asarray(genes, dtype=np.int64)
    pos = {int(g): j for j, g in enumerate(gene_arr)}
    h100, hfrac, rr = [], [], []
    for drug, gene, disease in triples:
        if int(gene) not in pos:
            h100.append(0); hfrac.append(0); rr.append(0.0)
            continue
        a = score_matrix(model, int(drug), dg_rels, gene_arr)
        b = score_matrix(model, int(disease), gd_rels, gene_arr)
        s = a + b            # negative distances: sum in score space
        t = s[pos[int(gene)]]
        rank = int(np.sum(s >= t))
        h100.append(1 if rank <= k_fixed else 0)
        hfrac.append(1 if rank <= k_norm else 0)
        rr.append(1.0 / rank)
    return {'hits@100': float(np.mean(h100)),
            'hits@frac': float(np.mean(hfrac)),
            'mrr': float(np.mean(rr)), 'n': len(h100)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--task', required=True, choices=['retrieval', 'multihop'])
    ap.add_argument('--models', default='transe,rotate')
    ap.add_argument('--epochs', default='1,2,5,10,20,40,80')
    ap.add_argument('--dim', type=int, default=64)
    ap.add_argument('--val-frac', type=float, default=0.15)
    ap.add_argument('--seed', type=int, default=42)
    args = ap.parse_args()

    ckpt = pickle.load(open(CKPT, 'rb'))
    need = ['retrieval_queries', 'gene_lists', 'disease_gene_rels',
            'multihop_triples']
    missing = [k for k in need if k not in ckpt]
    if missing:
        sys.exit(f'Checkpoint missing {missing}; re-run notebook 06.')

    cfg = load_config()
    KGC = cfg['knowledge_graphs']
    grid = [int(e) for e in args.epochs.split(',')]
    models = [ALIAS.get(m.strip().lower(), m.strip())
              for m in args.models.split(',')]
    k_frac = cfg['analysis_params']['k_fraction']

    aa_ref = {}
    for n, d in ckpt.get('nbhd_scalars', {}).items():
        aa_ref.setdefault(n, {})['retrieval'] = (d.get('recall@100'),
                                                 d.get('recall@frac'))
    for n, d in ckpt.get('multihop_scalars', {}).items():
        aa_ref.setdefault(n, {})['multihop'] = (d.get('hits@100'),
                                                d.get('hits@frac'))

    rows = []
    for name in ckpt['gene_lists']:
        genes = ckpt['gene_lists'][name]
        k_norm = max(1, int(np.ceil(k_frac * len(genes))))
        kg_df, _nodes = load_kg(name, cfg)

        if args.task == 'retrieval':
            items = [tuple(q) for q in ckpt['retrieval_queries'][name]]
            hidden = pd.MultiIndex.from_tuples(
                [(min(a, b), max(a, b)) for a, b in items]) if items else None
            rels_names = ckpt['disease_gene_rels'][name]
        else:
            items = [tuple(t) for t in ckpt['multihop_triples'][name]]
            hidden = pd.MultiIndex.from_tuples(
                [(min(d, s), max(d, s)) for d, _g, s in items]) if items else None
            rels_names = ckpt['disease_gene_rels'][name]

        if not items:
            print(f'{name}: no {args.task} items, skipped')
            continue

        train_df = drop_pairs(kg_df, hidden) if hidden is not None else kg_df
        triples, rel_to_idx, _ = kg_to_triples(train_df)
        dg_cfg = KGC[name]['relations'].get('drug_target', {})
        dg_names = (list(dg_cfg['relations']) if 'relations' in dg_cfg
                    else [dg_cfg['relation']] if 'relation' in dg_cfg else [])
        r_dis_gene = rel_indices(rel_to_idx, rels_names)
        r_drug_gene = rel_indices(rel_to_idx, dg_names)
        if not r_dis_gene or (args.task == 'multihop' and not r_drug_gene):
            print(f'{name}: required relations absent after hiding, skipped')
            continue

        # ---- Validation items, disjoint from the evaluation items -------
        # Epochs must not be chosen on the set being reported. Validation items
        # have the SAME structure as evaluation items and are hidden from
        # training in the same way, so the model never sees either, and only
        # validation drives selection.
        rng = np.random.RandomState(args.seed)
        gene_set = set(int(g) for g in genes)
        typed = kg_df[['x_index', 'y_index', 'relation']].to_numpy()

        if args.task == 'retrieval':
            dis_gene = defaultdict(set)
            keep_rel = set(rels_names)
            for a, b, r in typed:
                if str(r) not in keep_rel:
                    continue
                a, b = int(a), int(b)
                if b in gene_set:
                    dis_gene[a].add(b)
                elif a in gene_set:
                    dis_gene[b].add(a)
            eval_heads = {int(q[0]) for q in items}
            eval_pairs = {(int(a), int(b)) for a, b in items}
            cands = sorted(d for d, gs in dis_gene.items()
                           if d not in eval_heads and len(gs) >= 5)
            n_val = min(len(cands), max(100, int(args.val_frac * len(items))))
            val_items = []
            for d in (rng.choice(cands, size=n_val, replace=False)
                      if cands else []):
                gs = sorted(dis_gene[int(d)])
                val_items.append((int(d), int(gs[rng.randint(len(gs))])))

            if not val_items:
                # Hetionet has 121 eligible diseases and evaluates all 121, so
                # there is no spare disease to validate on. Fall back to a
                # SECOND, different gene from each evaluated disease: still not
                # an evaluation answer, and hidden from training just the same.
                # Knowing one partner of a disease does not reveal another, and
                # hiding two partners instead of one leaves the model with LESS
                # signal than Adamic-Adar had, not more, so the fallback is
                # conservative. Stated as a limitation: here validation and
                # evaluation queries share their head entities.
                for d, held in items:
                    gs = [g for g in sorted(dis_gene[int(d)])
                          if (int(d), g) not in eval_pairs]
                    if gs:
                        val_items.append((int(d), int(gs[rng.randint(len(gs))])))
                if val_items:
                    print(f'{name}: no spare diseases; validating on a second '
                          f'held-out gene from each evaluated disease '
                          f'({len(val_items):,} queries)')
        else:
            # Validation triples come from notebook 06, drawn from the wider CTD
            # pool with the same resolution, leakage filter and entity caps as
            # the benchmark, and disjoint from it.
            #
            # NOT constructed here. A locally-built proxy (pick any gene that
            # already bridges drug and disease in the graph) measures an easier
            # task than the benchmark: in testing it rose from 0.18 to 0.70
            # while the real metric FELL from 0.057 to 0.049, so epochs chosen
            # on it would have been the worst available for the actual task.
            # Validation has to match the evaluation's difficulty, not just its
            # shape.
            val_items = [tuple(t) for t in
                         ckpt.get('multihop_validation', {}).get(name, [])]
            eval_pairs = {(int(t[0]), int(t[2])) for t in items}
            val_items = [v for v in val_items
                         if (int(v[0]), int(v[2])) not in eval_pairs]

        if not val_items:
            print(f'{name}: could not build validation items, skipped')
            continue
        val_hidden = pd.MultiIndex.from_tuples(
            [(min(a, b), max(a, b)) for a, b in
             ([(v[0], v[1]) for v in val_items] if args.task == 'retrieval'
              else [(v[0], v[2]) for v in val_items])])
        train_df = drop_pairs(train_df, val_hidden)
        triples, rel_to_idx, _ = kg_to_triples(train_df)
        r_dis_gene = rel_indices(rel_to_idx, rels_names)
        r_drug_gene = rel_indices(rel_to_idx, dg_names)
        n_ent = int(max(triples[:, 0].max(), triples[:, 2].max())) + 1
        print(f'{name}: {len(items):,} eval items, {len(val_items):,} validation '
              f'items, both hidden from {len(triples):,} training triples')

        for mname in models:
            cls = MODELS[mname]
            m = cls(n_ent, len(rel_to_idx), dim=args.dim, seed=args.seed)
            trace, done = [], 0
            for n_ep in grid:
                if n_ep <= done:
                    continue
                t0 = time.time()
                m.fit(triples, n_epochs=n_ep - done, verbose=False)
                done = n_ep
                key = 'recall@frac' if args.task == 'retrieval' else 'hits@frac'
                if args.task == 'retrieval':
                    vres = eval_retrieval(m, val_items, genes, r_dis_gene, 100, k_norm)
                    eres = eval_retrieval(m, items, genes, r_dis_gene, 100, k_norm)
                else:
                    vres = eval_multihop(m, val_items, genes, r_drug_gene,
                                         r_dis_gene, 100, k_norm)
                    eres = eval_multihop(m, items, genes, r_drug_gene,
                                         r_dis_gene, 100, k_norm)
                trace.append({'epochs': n_ep, 'val': vres.get(key),
                              **{f'eval_{k}': v for k, v in eres.items()}})
                print(f'  {name} {mname} {args.task} epochs={n_ep:>3}  '
                      f'val={vres.get(key, float("nan")):.4f}  '
                      f'eval={eres.get(key, float("nan")):.4f}  '
                      f'[{time.time()-t0:.0f}s]', flush=True)
            # Selection on validation only.
            best = max(trace, key=lambda r: (r.get('val') if r.get('val') is not None else -1))
            ref = aa_ref.get(name, {}).get(args.task, (None, None))
            rows.append({'kg': name, 'model': mname, 'task': args.task,
                         'dim': args.dim, 'epochs_reported': best['epochs'],
                         'val_metric': best.get('val'),
                         'metric_fixedK': best.get(
                             'eval_recall@100' if args.task == 'retrieval'
                             else 'eval_hits@100'),
                         'metric_frac': best.get(f'eval_{key}'),
                         'adamic_adar_fixedK': ref[0],
                         'adamic_adar_frac': ref[1],
                         'k_norm': k_norm, 'n': best.get('eval_n'),
                         'trace': json.dumps(trace)})

    if not rows:
        sys.exit('No results.')
    df = pd.DataFrame(rows)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = OUT_DIR / f'{args.task}_embeddings.csv'
    df.to_csv(out, index=False)
    print(f'\nwrote {out}')
    for mname in models:
        sub = df[df.model == mname]
        if len(sub) < 2:
            continue
        aa_order = list(sub.sort_values('adamic_adar_frac', ascending=False).kg)
        em_order = list(sub.sort_values('metric_frac', ascending=False).kg)
        print(f'  {mname}: AA {aa_order} / {mname} {em_order}  '
              f'{"PRESERVED" if aa_order == em_order else "FLIPPED"}')


if __name__ == '__main__':
    main()
