#!/usr/bin/env python3
"""Embedding cross-check for notebook 06's link-prediction ranking.

Why this exists
---------------
All three task metrics in 06 are scored with Adamic-Adar, and the two
"alternative" scorers reported beside it, Common Neighbors and Jaccard, are
members of the same family: every one of them counts shared neighbours. So the
whole dimension is conditional on local neighbourhood overlap being the right
notion of proximity, and a ranking that only holds under that assumption is a
ranking about one family of heuristics, not about the graphs.

This script scores the SAME test pairs and the SAME negatives with two learned
embedding models whose inductive bias is entirely different: TransE (translation
in real space) and RotatE (rotation in complex space). Neither counts shared
neighbours. If the ranking of the graphs survives, the Adamic-Adar result is
substantially better supported. If it flips, that is a finding and belongs in
the paper.

Epoch calibration
-----------------
An untrained or undertrained embedding sits at chance, and reporting a
chance-level AUROC as "the embedding disagrees" would be a statement about the
training budget, not about the graph. Epochs are therefore calibrated on a
VALIDATION split carved out of the training positives, never on the test pairs.
The test set is scored exactly once, at the epoch count validation chose.

Inputs come from notebook 06's checkpoint so the split, the negatives and the
test pairs are identical to the ones Adamic-Adar was scored on. The script
refuses to run if the checkpoint does not carry them.

Usage
-----
    python supplementary/embeddings/run_link_prediction_embeddings.py \
        [--models transe,rotate] [--epochs 20,50,100,200] [--dim 128]
"""
import argparse
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'src'))

from biokgsuite.io.config import load_config                     # noqa: E402
from biokgsuite.io.loading import load_kg                        # noqa: E402
from biokgsuite.models.embedding import (                        # noqa: E402
    MODELS, build_train_triples)

CKPT = ROOT / 'results' / 'checkpoints' / '06_task_performance.pkl'
OUT_DIR = ROOT / 'results' / 'tables' / '06_task_performance'
ANCHOR_STRAT = 'type-constrained'


def _auroc(model, pos, neg, rel_idx, rel_idx_inv):
    pairs = np.array(list(pos) + list(neg))
    labels = np.array([1] * len(pos) + [0] * len(neg))
    s = model.score_pairs(pairs, rel_idx)
    if rel_idx_inv is not None:
        s = np.maximum(s, model.score_pairs(pairs, rel_idx_inv))
    return float(roc_auc_score(labels, s)), s, labels


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--models', default='transe,rotate')
    ap.add_argument('--epochs', default='1,2,5,10,20,40,80',
                    help='candidate epoch counts, chosen on validation AUROC')
    ap.add_argument('--dim', type=int, default=128)
    ap.add_argument('--val-frac', type=float, default=0.1)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--kgs', default='',
                    help='comma-separated graphs to run (default: every graph in the checkpoint)')
    ap.add_argument('--force', action='store_true',
                    help='recompute graphs already present in the output table')
    args = ap.parse_args()

    if not CKPT.exists():
        sys.exit(f'No checkpoint at {CKPT}. Run notebook 06 first.')
    ckpt = pickle.load(open(CKPT, 'rb'))

    needed = ('splits', 'negatives')
    missing = [k for k in needed if k not in ckpt]
    if missing:
        sys.exit(
            f'Checkpoint is missing {missing}. Notebook 06 must export the split '
            f'and the negatives so this cross-check scores exactly the pairs '
            f'Adamic-Adar was scored on. Re-run 06 after the export is added; '
            f'reconstructing the split here would not be verifiably the same.')

    cfg = load_config()
    KG_CONFIG = cfg['knowledge_graphs']

    def indication_rels(name):
        c = KG_CONFIG[name]['relations']['drug_disease']
        return list(c['relations']) if 'relations' in c else [c['relation']]

    epoch_grid = [int(e) for e in args.epochs.split(',')]
    model_names = [m.strip() for m in args.models.split(',')]
    _alias = {'transe': 'TransE', 'rotate': 'RotatE', 'gemma': 'Gemma'}
    model_names = [_alias.get(m.lower(), m) for m in model_names]

    aa = {r['kg']: r['auroc'] for r in ckpt['link_records']
          if r['heuristic'] == 'Adamic-Adar' and r['strategy'] == ANCHOR_STRAT}

    # The table is written after every graph and read back at the start, so a
    # job that stops at the wall clock keeps what it finished and the next
    # launch skips it. MATRIX alone can take most of a day.
    out_path = OUT_DIR / 'link_prediction_embeddings.csv'
    prev = pd.read_csv(out_path) if out_path.exists() and not args.force else pd.DataFrame()
    rows = prev.to_dict('records')
    have = {(r['kg'], r['model']) for r in rows}
    wanted = [k.strip() for k in args.kgs.split(',') if k.strip()] or list(ckpt['splits'])
    unknown = [k for k in wanted if k not in ckpt['splits']]
    if unknown:
        sys.exit(f'not in the 06 checkpoint: {unknown}; known: {list(ckpt["splits"])}')

    def write():
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(out_path, index=False)

    for name in wanted:
        if all((name, m) in have for m in model_names):
            print(f'{name}: already in {out_path.name}, skipped (use --force to redo)')
            continue
        sp = ckpt['splits'][name]
        neg = ckpt['negatives'][name][ANCHOR_STRAT]
        kg_df, _nodes = load_kg(name, cfg)
        ind_rels = indication_rels(name)
        train_triples, rel_to_idx, _ = build_train_triples(
            kg_df, set(map(tuple, sp['test_pos'])), ind_rels)

        # Validation carved from TRAIN positives only. The test pairs are not
        # touched until the epoch count is fixed.
        rng = np.random.RandomState(args.seed)
        tr = list(sp['train_pos'])
        perm = rng.permutation(len(tr))
        n_val = max(10, int(args.val_frac * len(tr)))
        val_pos = [tr[i] for i in perm[:n_val]]
        # Validation negatives drawn from the same pool shape as the test
        # negatives, at the same ratio.
        ratio = max(1, len(neg) // max(1, len(sp['test_pos'])))
        val_neg_idx = rng.choice(len(neg), size=min(len(neg), n_val * ratio),
                                 replace=False)
        val_neg = [neg[i] for i in sorted(val_neg_idx)]

        rel_key = next((r for r in rel_to_idx if r in set(ind_rels)), None)
        rel_idx = rel_to_idx.get(rel_key)
        rel_idx_inv = rel_to_idx.get(f'{rel_key}_inv') if rel_key else None
        if rel_idx is None:
            print(f'{name}: no indication relation in the training triples, skipped')
            continue

        n_ent = int(max(train_triples[:, 0].max(), train_triples[:, 2].max())) + 1
        n_rel = len(rel_to_idx)

        for mname in model_names:
            if (name, mname) in have:
                continue
            cls = MODELS[mname]
            m = cls(n_ent, n_rel, dim=args.dim, seed=args.seed)
            # Incremental training: fit() continues from the current weights, so
            # the sweep costs the MAX epoch count rather than their sum. Test
            # AUROC is recorded at every checkpoint but never used to choose one;
            # selection is on validation alone.
            trace, done = [], 0
            for n_ep in epoch_grid:
                if n_ep <= done:
                    continue
                t0 = time.time()
                m.fit(train_triples, n_epochs=n_ep - done, verbose=False)
                done = n_ep
                v, _, _ = _auroc(m, val_pos, val_neg, rel_idx, rel_idx_inv)
                t, _, _ = _auroc(m, sp['test_pos'], neg, rel_idx, rel_idx_inv)
                trace.append({'epochs': n_ep, 'val': v, 'test': t})
                print(f'  {name} {mname} epochs={n_ep:>4}  val={v:.4f}  '
                      f'test={t:.4f}  [{time.time()-t0:.0f}s]', flush=True)

            best = max(trace, key=lambda r: r['val'])
            calibrated = best['val'] >= 0.55
            if not calibrated:
                print(f'  {name} {mname}: validation AUROC peaked at '
                      f'{best["val"]:.4f} across {epoch_grid}, i.e. at chance. '
                      f'Reported as UNCALIBRATED, not as disagreement with '
                      f'Adamic-Adar: an undertrained model tells you about the '
                      f'training budget, not about the graph.')
            rows.append({
                'kg': name, 'model': mname, 'dim': args.dim,
                'epochs_selected': best['epochs'],
                'val_auroc': round(best['val'], 4),
                'test_auroc': round(best['test'], 4),
                'adamic_adar_auroc': aa.get(name),
                'calibrated': bool(calibrated),
                'epoch_trace': json.dumps(trace),
                'n_test_pos': len(sp['test_pos']), 'n_neg': len(neg),
            })
            print(f'{name:12s} {mname:8s}: test AUROC={best["test"]:.4f} '
                  f'at {best["epochs"]} epochs  (AA={aa.get(name)})', flush=True)
            have.add((name, mname))
            write()

    if not rows:
        sys.exit('No results produced.')
    df = pd.DataFrame(rows)
    write()
    print(f'\nwrote {out_path}')

    print('\nDoes the Adamic-Adar ranking survive?')
    for mname in model_names:
        sub = df[(df.model == mname) & df.calibrated]
        if len(sub) < 2:
            print(f'  {mname}: fewer than two calibrated graphs, no ranking to compare')
            continue
        aa_order = list(sub.sort_values('adamic_adar_auroc', ascending=False).kg)
        em_order = list(sub.sort_values('test_auroc', ascending=False).kg)
        print(f'  {mname}: AA order {aa_order} / {mname} order {em_order}  '
              f'{"PRESERVED" if aa_order == em_order else "FLIPPED"}')


if __name__ == '__main__':
    main()
