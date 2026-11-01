"""Topology: does a graph's structure reflect biological organisation?

    graph cohesion          share of nodes in the largest connected component
    clustering coefficient  local clustering relative to a density-matched
                            random graph, sampled identically for every graph
    graph reachability      whether validated drug-disease pairs stay connected
                            under edge removal, relative to type-matched random
                            pairs
    community purity        agreement between Louvain communities and entity
                            types, by AMI

The reachability gold standard is the reference set notebook 01 scores against,
resolved through the shared :class:`~biokgsuite.io.resolve.IdResolver`, so
the result does not depend on how much mapping work has been done for each
graph.

Used by ``biokgsuite_dimensions/04_topology``.
"""
from __future__ import annotations

import math
import random
from collections import Counter

import networkx as nx
import numpy as np
import pandas as pd

NA = float('nan')

# Clustering is estimated from a fixed sample of this many nodes in EVERY graph.
# Sampling k nodes costs k/n of the full computation, so it grows cheaper as
# graphs grow (~60 s for a 47k-node graph, ~3 s for a 7M-node one) while full
# computation does not scale at all. One sample size for every graph also keeps
# the comparison like for like rather than exact for some and estimated for
# others. Clustering lies in [0, 1], so its standard deviation cannot exceed
# 0.5, which bounds the standard error on the score below 0.03 for any graph.
CLUSTER_SAMPLE_N = 10_000

# Ceiling for the small-world ratio. The log reflects that 1x -> 10x is a more
# meaningful step than 100x -> 1000x; a lower ceiling saturates several graphs
# at 1.0 and loses all discrimination.
SW_CEILING = 1000


# ------------------------------------------------------------------- cohesion
def connectivity(graphs: dict) -> tuple:
    """Per graph: component structure. Returns (frame, largest component per KG)."""
    rows, lcc_sets = [], {}
    for name, G in graphs.items():
        components = sorted(nx.connected_components(G), key=len, reverse=True)
        lcc = components[0] if components else set()
        lcc_sets[name] = lcc
        n = G.number_of_nodes()
        if len(components) > 1:
            outside = Counter(G.nodes[v].get('node_type')
                              for c in components[1:] for v in c)
            top = ', '.join(f'{t} ({c:,})' for t, c in outside.most_common(3))
        else:
            top = '-'
        rows.append({
            'KG': name, 'Nodes': n, 'LCC nodes': len(lcc),
            # Both, deliberately. 'LCC (%)' is the display column and is
            # rounded; 'LCC fraction' is what the sub-metric is scored from.
            # Scoring off the rounded percent wrote a twice-rounded value into
            # raw_statistics.csv, which Outputs.table promises is full precision.
            'LCC fraction': len(lcc) / n if n else NA,
            'LCC (%)': round(100 * len(lcc) / n, 2) if n else NA,
            'Components': len(components),
            'Isolated': sum(1 for v in G if G.degree(v) == 0),
            'Top types outside LCC': top,
        })
    return pd.DataFrame(rows).set_index('KG'), lcc_sets


# ----------------------------------------------------------------- clustering
def small_world_score(ratio: float, ceiling: int = SW_CEILING) -> float:
    """C/C_ER mapped to [0, 1] on a log scale. Defined ONCE and used everywhere."""
    if not ratio or ratio <= 1 or ratio != ratio:
        return 0.0
    return round(min(1.0, math.log(ratio) / math.log(ceiling)), 4)


def clustering(graphs: dict, seed: int, k: int = CLUSTER_SAMPLE_N) -> pd.DataFrame:
    """Per graph: sampled clustering, its standard error, and the scored ratio.

    Per-node values are kept rather than discarded so the sampling error is
    measured from the sample itself and can be reported, instead of asserted.
    """
    rows = []
    for name, G in graphs.items():
        n, m = G.number_of_nodes(), G.number_of_edges()
        mean_degree = 2 * m / n if n else 0.0
        # Erdos-Renyi expectation: <k>/(n-1). This matches the graph's DENSITY,
        # not its degree sequence, and assumes a Poisson degree distribution.
        # Biomedical KGs are heavy-tailed, where the configuration-model
        # expectation is far larger, so C/C_ER overstates enrichment relative
        # to a degree-preserving null. Reported as an ER ratio for that reason.
        c_er = mean_degree / (n - 1) if n > 1 else 0.0

        sample = (random.Random(seed).sample(list(G.nodes()), k)
                  if n > k else list(G.nodes()))
        per_node = np.fromiter(nx.clustering(G, nodes=sample).values(), dtype=float)
        c = float(per_node.mean()) if per_node.size else 0.0
        se = float(per_node.std(ddof=1) / math.sqrt(per_node.size)) if per_node.size > 1 else 0.0
        ratio = c / c_er if c_er else NA

        score = small_world_score(ratio)
        hi = small_world_score((c + se) / c_er) if c_er else NA
        lo = small_world_score(max(c - se, 1e-9) / c_er) if c_er else NA
        rows.append({
            'KG': name, 'Nodes sampled': len(sample),
            'Sampled (%)': round(100 * len(sample) / n, 1) if n else NA,
            'C (empirical)': round(c, 4), 'C SE': round(se, 4),
            'C_ER (random)': round(c_er, 6),
            'C / C_ER': round(ratio, 1) if ratio == ratio else NA,
            'SW score': score,
            'SW score SE': round(abs(hi - lo) / 2, 4) if hi == hi else NA,
        })
    return pd.DataFrame(rows).set_index('KG')


# --------------------------------------------------------------- reachability
def canonical_index_maps(nodes_df, type_map: dict, resolver) -> dict:
    """canonical type -> {canonical id: node index}, via the shared resolver.

    Replaces per-graph identifier parsing: every graph's drugs resolve to ChEMBL
    and its diseases to MONDO, the namespaces the reference set is keyed in.
    """
    maps = {}
    for canon in ('Drug', 'Disease'):
        sub = nodes_df[nodes_df['type'].isin(type_map.get(canon, ()))]
        m = {}
        for raw, idx in zip(sub['id'].astype(str), sub['idx'].astype(int)):
            for cid in resolver.resolve(canon, [raw])[0]:
                m.setdefault(cid, idx)
        maps[canon] = m
    return maps


def map_gold_pairs(gold_pairs, index_maps: dict, lcc: set) -> list:
    """Gold (drug, disease) pairs as node-index pairs, both endpoints in the LCC.

    Restricting to the largest component matters: a pair with an endpoint
    outside it is unreachable by construction, which would be recorded as a
    topology failure rather than as the absence it is.
    """
    dmap, smap = index_maps['Drug'], index_maps['Disease']
    return sorted({(dmap[d], smap[s]) for d, s in gold_pairs
                   if d in dmap and s in smap
                   and dmap[d] in lcc and smap[s] in lcc})


def negative_pairs(G, lcc: set, type_map: dict, positives: set,
                   n: int, rng) -> list:
    """Type-matched random drug-disease pairs that are not gold-standard pairs.

    Type-matched, so the comparison isolates whether GOLD pairs are specifically
    resilient rather than whether drug-disease pairs in general are.
    """
    drug_labels = set(type_map.get('Drug', ()))
    disease_labels = set(type_map.get('Disease', ()))
    drugs = [v for v in lcc if G.nodes[v].get('node_type') in drug_labels]
    diseases = [v for v in lcc if G.nodes[v].get('node_type') in disease_labels]
    if not drugs or not diseases:
        return []
    out, seen, attempts, cap = [], set(), 0, n * 50
    while len(out) < n and attempts < cap:
        attempts += 1
        d, s = drugs[rng.randint(len(drugs))], diseases[rng.randint(len(diseases))]
        if (d, s) in positives or (s, d) in positives or (d, s) in seen:
            continue
        seen.add((d, s))
        out.append((d, s))
    return out


def indexed_edges(G) -> tuple:
    """(edge array, node -> row index) for fast repeated edge dropout."""
    node_to_idx = {v: i for i, v in enumerate(G.nodes())}
    edges = np.fromiter((i for u, v in G.edges()
                         for i in (node_to_idx[u], node_to_idx[v])),
                        dtype=np.int32).reshape(-1, 2)
    return edges, node_to_idx


def _csr(edges, n_nodes: int) -> tuple:
    """Symmetric adjacency as CSR arrays (indptr, indices).

    The previous implementation rebuilt a Python dict of sets on every dropout
    trial: 36 rebuilds per graph, each boxing ~4M integers into set objects.
    That dominated both time and memory and is what made this cell unrunnable on
    a small machine. CSR holds the same adjacency in two int32 arrays, roughly
    8 bytes per directed edge instead of ~100, and is built with numpy rather
    than a Python loop.
    """
    src = np.concatenate([edges[:, 0], edges[:, 1]])
    dst = np.concatenate([edges[:, 1], edges[:, 0]])
    order = np.argsort(src, kind='stable')
    indices = dst[order].astype(np.int32, copy=False)
    indptr = np.zeros(n_nodes + 1, dtype=np.int64)
    np.cumsum(np.bincount(src, minlength=n_nodes), out=indptr[1:])
    return indptr, indices


def _reachable(indptr, indices, src: int, tgt: int, max_hops: int,
               seen: np.ndarray, stamp: int) -> bool:
    """Breadth-first reachability within max_hops over CSR adjacency.

    `seen` is a reusable int array stamped per query, so no allocation happens
    inside the loop.
    """
    if src == tgt:
        return True
    seen[src] = stamp
    frontier = [src]
    for _ in range(max_hops):
        nxt = []
        for node in frontier:
            for j in range(indptr[node], indptr[node + 1]):
                nbr = indices[j]
                if nbr == tgt:
                    return True
                if seen[nbr] != stamp:
                    seen[nbr] = stamp
                    nxt.append(nbr)
        if not nxt:
            return False
        frontier = nxt
    return False


def dropout_curve(edges, pairs: list, drop_rates, rng, max_hops: int,
                  n_trials: int, n_nodes: int) -> tuple:
    """(mean recovery per rate, per-trial values per rate) under random dropout."""
    means, trials_by_rate = [], {}
    n_edges = len(edges)
    seen = np.zeros(n_nodes, dtype=np.int32)
    stamp = 0
    for rate in drop_rates:
        n_remove = int(round(rate * n_edges))
        trials = []
        for _ in range(n_trials):
            if n_remove:
                keep = np.ones(n_edges, dtype=bool)
                keep[rng.choice(n_edges, size=n_remove, replace=False)] = False
                kept = edges[keep]
            else:
                kept = edges
            indptr, indices = _csr(kept, n_nodes)
            hits = 0
            for s, t in pairs:
                stamp += 1
                if _reachable(indptr, indices, s, t, max_hops, seen, stamp):
                    hits += 1
            trials.append(hits / len(pairs) if pairs else 0.0)
            del indptr, indices
        trials_by_rate[rate] = trials
        means.append(float(np.mean(trials)))
    return means, trials_by_rate


def resilience_index(drop_rates, recoveries) -> float:
    """Normalised area under the recovery-vs-dropout curve, in [0, 1]."""
    x = np.asarray(drop_rates, dtype=float)
    y = np.asarray(recoveries, dtype=float)
    span = float(x[-1] - x[0])
    return round(float(np.trapezoid(y, x)) / span, 4) if span > 0 else 0.0


def differential_resilience(ri_gold: float, ri_rand: float) -> float:
    """Share of the headroom above random that gold-standard pairs occupy.

        DR = (RI_gold - RI_rand) / max(eps, 1 - RI_rand), clipped at 0

    Density-invariant: a larger or denser graph that preserves ALL pairs equally
    gains nothing. Defined once, so the checkpoint and the score carry the same
    value.
    """
    return round(max(0.0, (ri_gold - ri_rand) / max(1e-6, 1.0 - ri_rand)), 4)


# ------------------------------------------------------------------ community
# Louvain on graphs above this size is restricted to the largest component.
# Nodes outside it form their own trivial communities and shift the agreement
# scores while telling you nothing about how the graph is organised.
LOUVAIN_LCC_THRESHOLD = 1_000_000


def canonical_labels(raw_types, type_map: dict) -> tuple:
    """Raw KG-internal type strings mapped to the benchmark's canonical types.

    A raw label the config does not map keeps its own name and is returned in
    the second element, rather than being dropped or folded into one bucket:
    an unmapped type is a schema-coverage fact the notebook should report.
    """
    lookup = {}
    for canon, labels in (type_map or {}).items():
        for raw in (labels if isinstance(labels, (list, tuple, set)) else (labels,)):
            lookup[raw] = canon
    mapped = [lookup.get(t, t) for t in raw_types]
    unmapped = sorted({t for t in raw_types if t not in lookup})
    return mapped, unmapped


def communities(G, seed: int, resolution: float = 1.0,
                lcc_threshold: int = 0,
                type_map: dict | None = None) -> dict:
    """Louvain communities and their agreement with entity types.

    AMI is the scored statistic. NMI rises with the number of communities even
    when the partition carries no information, so it would reward a graph for
    being fragmented; AMI corrects for chance (Vinh et al. 2010). NMI and ARI
    are returned beside it for interpretation.

    Agreement is measured against the CANONICAL types from `type_map`
    (`config/kgs.yaml: entity_types`), not the raw strings a graph happens to
    use. BioKG splits diseases into `Disease` (MeSH) and `GeneticDisorder`
    (OMIM); the config adjudicates both as `Disease` and notebook 01 counts
    them as one, so grading AMI on the raw split would penalise the graph for a
    labelling artefact this benchmark has already resolved -- and would let a
    graph's own schema decide what counts as a pure community. Raw labels the
    config does not map are kept as their own class and reported in
    `unmapped_types`. Passing `type_map=None` falls back to raw labels.
    """
    from sklearn.metrics import (adjusted_mutual_info_score, adjusted_rand_score,
                                 normalized_mutual_info_score)

    # Always partition the largest connected component, not the whole graph.
    #
    # Louvain assigns every isolated node its own community, so on a graph with
    # a long tail of singletons the partition being graded is mostly singletons:
    # Hetionet carries 1,873 isolated nodes and the full-graph run returned 1,882
    # communities, and AMI, NMI, ARI and Q would be computed against that.
    # Restricting every graph to its largest component (the default,
    # `lcc_threshold` = 0) keeps the scored AMI comparable across the panel.
    # `nodes_clustered` reports what was actually partitioned.
    target = G
    if G.number_of_nodes() > lcc_threshold and G.number_of_nodes() > 0:
        components = list(nx.connected_components(G))
        if components:
            target = G.subgraph(max(components, key=len)).copy()

    found = nx.community.louvain_communities(target, resolution=resolution, seed=seed)
    node_to_community = {v: i for i, comm in enumerate(found) for v in comm}

    nodes = list(target.nodes())
    raw = [target.nodes[v].get('node_type', 'unknown') for v in nodes]
    entity, unmapped = canonical_labels(raw, type_map) if type_map else (raw, [])
    label = [node_to_community.get(v, -1) for v in nodes]

    return {
        'n_communities': len(found),
        'n_nodes': G.number_of_nodes(),
        'n_edges': G.number_of_edges(),
        # Stamped so a cached partition can be invalidated when the library
        # that produced it changes. environment.yml pins networkx only to
        # >=3.2,<4, and the Louvain implementation has changed inside that
        # range, so node and edge counts alone do not identify a partition.
        'networkx_version': nx.__version__,
        'nodes_clustered': target.number_of_nodes(),
        'n_type_classes': len(set(entity)),
        'unmapped_types': unmapped,
        # What the agreement scores were graded against. Part of the cache key:
        # a partition is still valid when the labelling basis changes, but its
        # AMI/NMI/ARI are not.
        'ami_basis': 'canonical' if type_map else 'raw',
        'ami': float(adjusted_mutual_info_score(entity, label)),
        'nmi': float(normalized_mutual_info_score(entity, label)),
        'ari': float(adjusted_rand_score(entity, label)),
        'modularity': float(nx.community.modularity(target, found, resolution=resolution)),
        'node_to_community': node_to_community,
        'communities': found,
    }
