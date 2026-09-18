#!/usr/bin/env python3
"""Run a single embedding model on a single KG. Saves to results/cache/."""
import json, sys, time, numpy as np, pickle
from pathlib import Path

ROOT = next(p for p in Path(__file__).resolve().parents
            if (p / 'config' / 'kgs.yaml').exists())
sys.path.insert(0, str(ROOT / 'src'))

from biokgsuite.models.embedding import (TransE, RotatE, GemmaNameEmbedder,
                           build_train_triples, kg_to_triples,
                           compute_embedding_metrics)
from biokgsuite.graph.negative_sampling import generate_negatives


def prepare_kg(kg_name, seed=42, neg_ratio=5):
    """Load and prepare KG data, cache to pkl."""
    cache = ROOT / 'results' / 'cache' / f'{kg_name}_prep.pkl'
    if cache.exists():
        with open(cache, 'rb') as f:
            data = pickle.load(f)
        if data is not None and isinstance(data, dict) and 'train_triples' in data:
            return data

    from biokgsuite.io import find_config, load_config, load_kg
    config = load_config(find_config(ROOT))
    kg_df, nodes_df = load_kg(kg_name, config)
    kg_cfg = config['knowledge_graphs'][kg_name]
    etypes = kg_cfg['entity_types']
    type_map = dict(zip(nodes_df['idx'], nodes_df['type']))
    drug_idx = {i for i, t in type_map.items() if t == etypes.get('Drug', 'Drug')}
    disease_idx = {i for i, t in type_map.items() if t == etypes.get('Disease', 'Disease')}
    gene_idx = {i for i, t in type_map.items() if t == etypes.get('Gene/Protein', 'Gene')}

    dd = kg_cfg.get('relations', {}).get('drug_disease', {})
    ind_rels = [dd['relation']] if 'relation' in dd else dd.get('relations', [])
    dt = kg_cfg.get('relations', {}).get('drug_target', {})
    dt_rels = [dt['relation']] if 'relation' in dt else dt.get('relations', [])

    mask = kg_df['relation'].isin(ind_rels)
    pairs = set()
    for _, row in kg_df[mask].iterrows():
        h, t = int(row['x_index']), int(row['y_index'])
        if h in drug_idx and t in disease_idx: pairs.add((h, t))
        elif t in drug_idx and h in disease_idx: pairs.add((t, h))
    pairs = list(pairs)

    rng = np.random.RandomState(seed)
    perm = rng.permutation(len(pairs))
    split = int(0.9 * len(pairs))
    test_pos = [pairs[i] for i in perm[split:]]
    all_pos = set(pairs)

    train_triples, rel_to_idx, idx_to_rel = build_train_triples(
        kg_df, set(test_pos), ind_rels)
    n_ent = int(nodes_df['idx'].max()) + 1
    node_name_map = dict(zip(nodes_df['idx'], nodes_df['name']))

    drug_targets = {}
    dt_mask = kg_df['relation'].isin(dt_rels)
    for _, row in kg_df[dt_mask].iterrows():
        h, t = int(row['x_index']), int(row['y_index'])
        if h in drug_idx and t in gene_idx: drug_targets.setdefault(h, set()).add(t)
        elif t in drug_idx and h in gene_idx: drug_targets.setdefault(t, set()).add(h)

    neg_pairs = generate_negatives(test_pos, len(test_pos) * neg_ratio,
        'type-constrained', drug_idx, disease_idx, drug_targets,
        node_name_map, all_pos, rng)

    rel_idx = rel_to_idx[ind_rels[0]]
    inv_name = f'{ind_rels[0]}_inv'
    rel_idx_inv = rel_to_idx.get(inv_name)

    prep = {'train_triples': train_triples, 'rel_to_idx': rel_to_idx,
            'n_ent': n_ent, 'n_rels': len(rel_to_idx),
            'test_pos': test_pos, 'neg_pairs': neg_pairs,
            'rel_idx': rel_idx, 'rel_idx_inv': rel_idx_inv,
            'n_train': split, 'n_test': len(test_pos)}

    with open(cache, 'wb') as f:
        pickle.dump(prep, f)
    return prep


def _run_gemma(kg_name, p, dim=768, batch_size=64, seed=42):
    """Gemma branch: encode entity names with EmbeddingGemma, score pairs by
    cosine similarity. No training. Reuses prepare_kg's test/neg splits so
    metrics are directly comparable to TransE/RotatE.

    Caches the (n_entities × dim) embedding matrix to disk so re-runs across
    negative-sampling strategies don't re-encode.
    """
    from biokgsuite.io import find_config, load_config, load_kg
    config = load_config(find_config(ROOT))
    _, nodes_df = load_kg(kg_name, config)

    n_ent = p['n_ent']
    # Build name list in idx order; missing idx → '' (rare, but possible)
    name_by_idx = dict(zip(nodes_df['idx'].astype(int), nodes_df['name']))
    names = [name_by_idx.get(i, '') for i in range(n_ent)]
    # Audit: how many names look like opaque IDs?
    import re
    id_like = sum(1 for n in names if not n or re.fullmatch(r'[A-Z0-9:_.\-]+', str(n)))
    pct_id = 100.0 * id_like / max(n_ent, 1)
    print(f'{kg_name}/Gemma: {n_ent} entities | ~{pct_id:.1f}% look ID-only '
          f'(opaque alphanumerics)', flush=True)
    if pct_id > 50:
        print(f'  WARNING: majority of names look like IDs, not human-readable.\n'
              f'  Word-priors result for this KG will be near-noise without '
              f'external name resolution.', flush=True)

    # Cache the embedding matrix per (kg, model, dim) so repeated runs don't
    # re-encode (encoding is the slow part).
    emb_cache = ROOT / 'results' / 'cache' / f'gemma_emb_{kg_name}_d{dim}.npz'
    model = GemmaNameEmbedder(n_entities=n_ent, n_relations=p['n_rels'],
                              dim=dim, batch_size=batch_size, seed=seed)

    t0 = time.time()
    if emb_cache.exists():
        try:
            model.load_embeddings(emb_cache)
            print(f'  Loaded cached embeddings from {emb_cache.name}', flush=True)
        except Exception as e:
            print(f'  Cache load failed ({e}) — re-encoding', flush=True)
            model.encode_entities(names)
            model.save_embeddings(emb_cache)
    else:
        model.encode_entities(names)
        model.save_embeddings(emb_cache)
    ts = time.time() - t0
    print(f'Gemma encoded in {ts:.1f}s', flush=True)
    return model, ts


def run_model(kg_name, model_name, n_epochs=5, dim=32, seed=42):
    p = prepare_kg(kg_name, seed=seed)
    print(f'{kg_name}/{model_name}: {len(p["train_triples"])} triples, '
          f'{p["n_ent"]} entities, {p["n_rels"]} rels', flush=True)

    if model_name == 'Gemma':
        # Use 768-d by default for Gemma (full Matryoshka). Override via CLI.
        gemma_dim = 768 if dim in (32, 0) else dim
        model, ts = _run_gemma(kg_name, p, dim=gemma_dim, seed=seed)
        m = compute_embedding_metrics(model, p['test_pos'], p['neg_pairs'],
                                      p['rel_idx'], rel_idx_inv=None)
        m['train_time_s'] = ts
        m['n_epochs'] = 0
        m['dim'] = gemma_dim
        print(f'AUROC={m["auroc"]:.4f} AUPRC={m["auprc"]:.4f} MRR={m["mrr"]:.4f} '
              f'H@10={m["hits@10"]:.4f} H@50={m["hits@50"]:.4f} '
              f'H@100={m["hits@100"]:.4f}', flush=True)
        # Save into the same per-KG JSON the KGE models use, under 'Gemma' key
        out = ROOT / 'results' / 'cache' / f'embedding_{kg_name}.json'
        data = {}
        if out.exists():
            with open(out) as f:
                data = json.load(f)
        data.setdefault('kg', kg_name)
        data.setdefault('models', {})
        data['models']['Gemma'] = m
        data['n_test'] = len(p['test_pos'])
        data['n_neg'] = len(p['neg_pairs'])
        data['n_entities'] = p['n_ent']
        data['n_relations'] = p['n_rels']
        data['n_train_triples'] = len(p['train_triples'])
        with open(out, 'w') as f:
            json.dump(data, f, indent=2)
        print(f'Saved to {out}', flush=True)
        return m

    Cls = TransE if model_name == 'TransE' else RotatE
    kw = dict(n_entities=p['n_ent'], n_relations=p['n_rels'],
              dim=dim, seed=seed, lr=0.01)
    if model_name == 'RotatE':
        kw['dim'] = max(dim // 2, 8)
        kw['margin'] = 6.0
    else:
        kw['margin'] = 1.0

    t0 = time.time()
    model = Cls(**kw)
    model.fit(p['train_triples'], n_epochs=n_epochs, batch_size=32768, verbose=True)
    ts = time.time() - t0
    print(f'{model_name} trained in {ts:.1f}s', flush=True)

    m = compute_embedding_metrics(model, p['test_pos'], p['neg_pairs'],
                                  p['rel_idx'], rel_idx_inv=p['rel_idx_inv'])
    m['train_time_s'] = ts
    m['n_epochs'] = n_epochs
    m['dim'] = dim
    print(f'AUROC={m["auroc"]:.4f} AUPRC={m["auprc"]:.4f} MRR={m["mrr"]:.4f} '
          f'H@10={m["hits@10"]:.4f} H@50={m["hits@50"]:.4f} H@100={m["hits@100"]:.4f}',
          flush=True)

    # Save/update JSON
    out = ROOT / 'results' / 'cache' / f'embedding_{kg_name}.json'
    data = {}
    if out.exists():
        with open(out) as f:
            data = json.load(f)
    data.setdefault('kg', kg_name)
    data.setdefault('models', {})
    data['models'][model_name] = m
    data['n_test'] = len(p['test_pos'])
    data['n_neg'] = len(p['neg_pairs'])
    data['n_entities'] = p['n_ent']
    data['n_relations'] = p['n_rels']
    data['n_train_triples'] = len(p['train_triples'])
    with open(out, 'w') as f:
        json.dump(data, f, indent=2)
    print(f'Saved to {out}', flush=True)
    return m


if __name__ == '__main__':
    kg = sys.argv[1]
    model = sys.argv[2]
    epochs = int(sys.argv[3]) if len(sys.argv) > 3 else 5
    dim = int(sys.argv[4]) if len(sys.argv) > 4 else 32
    run_model(kg, model, n_epochs=epochs, dim=dim)
