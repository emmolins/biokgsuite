#!/usr/bin/env python3
"""Run EmbeddingGemma-300m on MATRIX (Every Cure KG).

MATRIX is the long pole. After config-based filtering to the canonical
{Drug, Disease, Gene/Protein, Pathway, Phenotype} subset, it still has
~500K-3M entities depending on what's kept. Even encoding 500K short
strings with a 308M-param model on CPU takes several hours; on GPU,
~10-30 min.

This script offers two modes:

  --mode full        Encode every entity. Use on a workstation with a GPU.
  --mode sampled     Subsample to N entities (default 100K, all drugs +
                     diseases retained, random gene/pathway/phenotype
                     subsample). Gives a defensible-but-noisy number in
                     ~30 min on CPU.

Prereqs: same as scripts/run_gemma_benchmark.sh

Output: results/cache/embedding_matrix.json (Gemma key added; preserves
existing TransE/RotatE entries if present).
"""
import argparse, json, sys, time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from biokgsuite.models.embedding import GemmaNameEmbedder, compute_embedding_metrics
from biokgsuite.io import find_config, load_config, load_kg
from biokgsuite.graph.negative_sampling import generate_negatives
from biokgsuite.models.embedding import build_train_triples


def prepare_matrix(seed=42, neg_ratio=5, subsample_to=None):
    """Like run_emb_model.prepare_kg but with MATRIX-specific knobs."""
    config = load_config(ROOT)
    kg_df, nodes_df = load_kg('matrix', config)
    kg_cfg = config['knowledge_graphs']['matrix']
    etypes = kg_cfg['entity_types']

    type_map = dict(zip(nodes_df['idx'], nodes_df['type']))
    drug_idx = {i for i, t in type_map.items() if t == etypes.get('Drug', 'Drug')}
    disease_idx = {i for i, t in type_map.items() if t == etypes.get('Disease', 'Disease')}
    gene_idx = {i for i, t in type_map.items() if t == etypes.get('Gene/Protein', 'Gene')}

    dd = kg_cfg.get('relations', {}).get('drug_disease', {})
    ind_rels = [dd['relation']] if 'relation' in dd else dd.get('relations', [])
    dt = kg_cfg.get('relations', {}).get('drug_target', {})
    dt_rels = [dt['relation']] if 'relation' in dt else dt.get('relations', [])

    # Drug-disease indication pairs
    mask = kg_df['relation'].isin(ind_rels)
    pairs = set()
    for _, row in kg_df[mask].iterrows():
        h, t = int(row['x_index']), int(row['y_index'])
        if h in drug_idx and t in disease_idx:
            pairs.add((h, t))
        elif t in drug_idx and h in disease_idx:
            pairs.add((t, h))
    pairs = list(pairs)
    print(f'MATRIX: {len(pairs):,} drug-disease indication pairs')

    rng = np.random.RandomState(seed)
    perm = rng.permutation(len(pairs))
    split = int(0.9 * len(pairs))
    test_pos = [pairs[i] for i in perm[split:]]
    all_pos = set(pairs)

    train_triples, rel_to_idx, _ = build_train_triples(
        kg_df, set(test_pos), ind_rels)
    n_ent_full = int(nodes_df['idx'].max()) + 1
    node_name_map = dict(zip(nodes_df['idx'], nodes_df['name']))

    drug_targets = {}
    dt_mask = kg_df['relation'].isin(dt_rels)
    for _, row in kg_df[dt_mask].iterrows():
        h, t = int(row['x_index']), int(row['y_index'])
        if h in drug_idx and t in gene_idx:
            drug_targets.setdefault(h, set()).add(t)
        elif t in drug_idx and h in gene_idx:
            drug_targets.setdefault(t, set()).add(h)

    neg_pairs = generate_negatives(
        test_pos, len(test_pos) * neg_ratio, 'type-constrained',
        drug_idx, disease_idx, drug_targets, node_name_map, all_pos, rng)

    rel_idx = rel_to_idx[ind_rels[0]]

    # Subsampling: keep all drugs+diseases (needed for the test pairs) +
    # a random subsample of other entities.
    if subsample_to is not None and subsample_to < n_ent_full:
        keep = set(drug_idx) | set(disease_idx)
        # Also keep all idx that appear in test_pos / neg_pairs
        for h, t in test_pos + neg_pairs:
            keep.add(int(h)); keep.add(int(t))
        remaining_budget = max(0, subsample_to - len(keep))
        other = list(set(range(n_ent_full)) - keep)
        if remaining_budget > 0 and other:
            sample = rng.choice(other, size=min(remaining_budget, len(other)),
                                replace=False)
            keep.update(sample.tolist())
        # Renumber: idx -> compact
        old_to_new = {old: new for new, old in enumerate(sorted(keep))}
        new_to_name = {old_to_new[old]: node_name_map.get(old, '')
                       for old in keep}
        # Remap test / neg pairs
        test_pos = [(old_to_new[h], old_to_new[t]) for h, t in test_pos
                    if h in old_to_new and t in old_to_new]
        neg_pairs = [(old_to_new[h], old_to_new[t]) for h, t in neg_pairs
                     if h in old_to_new and t in old_to_new]
        n_ent = len(keep)
        names = [new_to_name.get(i, '') for i in range(n_ent)]
        print(f'MATRIX subsampled: {n_ent_full:,} → {n_ent:,} entities; '
              f'kept {len(test_pos):,} test_pos, {len(neg_pairs):,} neg_pairs')
    else:
        n_ent = n_ent_full
        names = [node_name_map.get(i, '') for i in range(n_ent)]

    return {
        'names': names,
        'n_ent': n_ent,
        'n_rels': len(rel_to_idx),
        'test_pos': test_pos,
        'neg_pairs': neg_pairs,
        'rel_idx': rel_idx,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', choices=['full', 'sampled'], default='sampled')
    ap.add_argument('--subsample', type=int, default=100_000,
                    help='entity budget for --mode sampled (default 100K)')
    ap.add_argument('--dim', type=int, default=768,
                    help='Matryoshka dim (128/256/512/768)')
    ap.add_argument('--batch', type=int, default=128)
    ap.add_argument('--seed', type=int, default=42)
    args = ap.parse_args()

    sub = None if args.mode == 'full' else args.subsample
    p = prepare_matrix(seed=args.seed, subsample_to=sub)

    model = GemmaNameEmbedder(n_entities=p['n_ent'], n_relations=p['n_rels'],
                              dim=args.dim, batch_size=args.batch, seed=args.seed)
    cache_path = ROOT / 'results' / 'cache' / (
        f'gemma_emb_matrix_d{args.dim}'
        f'{"_sampled"+str(args.subsample) if sub else ""}.npz')

    t0 = time.time()
    if cache_path.exists():
        try:
            model.load_embeddings(cache_path)
            print(f'Loaded cached embeddings from {cache_path.name}')
        except Exception as e:
            print(f'Cache load failed ({e}) — re-encoding')
            model.encode_entities(p['names'])
            model.save_embeddings(cache_path)
    else:
        model.encode_entities(p['names'])
        model.save_embeddings(cache_path)
    t_enc = time.time() - t0
    print(f'Encoding: {t_enc:.1f}s')

    m = compute_embedding_metrics(model, p['test_pos'], p['neg_pairs'],
                                  p['rel_idx'], rel_idx_inv=None)
    m['train_time_s'] = t_enc
    m['n_epochs'] = 0
    m['dim'] = args.dim
    m['matrix_mode'] = args.mode
    m['matrix_subsample'] = sub
    print(f"AUROC={m['auroc']:.4f}  AUPRC={m['auprc']:.4f}  "
          f"MRR={m['mrr']:.4f}  H@10={m['hits@10']:.4f}")

    # Save to embedding_matrix.json under 'Gemma'
    out = ROOT / 'results' / 'cache' / 'embedding_matrix.json'
    data = {}
    if out.exists():
        with open(out) as f:
            data = json.load(f)
    data.setdefault('kg', 'matrix')
    data.setdefault('models', {})
    data['models']['Gemma'] = m
    data['n_test_gemma'] = len(p['test_pos'])
    data['n_neg_gemma'] = len(p['neg_pairs'])
    data['n_entities_gemma'] = p['n_ent']
    with open(out, 'w') as f:
        json.dump(data, f, indent=2)
    print(f'Saved to {out}')


if __name__ == '__main__':
    main()
