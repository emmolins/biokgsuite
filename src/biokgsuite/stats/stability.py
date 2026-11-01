"""Stability: how well does a graph's predictive signal survive edge loss?

    random-dropout stability      Spearman r of Adamic-Adar rankings of known
                                  drug-disease pairs under uniform edge removal
    periphery-dropout stability   the same, removing lowest-degree edges first

Both are reported against a DEGREE-MATCHED NULL: the same measurement on random
drug-disease pairs whose endpoint degrees match the therapeutic pairs. Adamic-Adar
is a degree-driven statistic, so a high correlation on its own does not show that
a graph encodes robust *therapeutic* structure; it may only show that rankings
driven by hub nodes survive random edge loss. The null is what separates the two,
and on the graphs measured so far most of the stability is degree-driven. It is
reported, not scored, in the same way pathway evidence is reported in notebook 01
and NMI/ARI in notebook 04.

Used by ``biokgsuite_dimensions/05_stability``.
"""
from __future__ import annotations

import bisect

import networkx as nx
import numpy as np
from scipy import stats

NA = float('nan')

# Nodes either side of the matched degree that a null endpoint may be drawn from.
DEGREE_MATCH_WINDOW = 10


# ------------------------------------------------------------------- scoring
def cles(r: float) -> float:
    """Linear rescaling of a correlation to [0, 1]: (r + 1) / 2.

    The name is historical and is kept because callers and stored checkpoints
    use it. It is NOT the Common Language Effect Size, and the value is not a
    probability. (r + 1) / 2 is the probability of concordance for Kendall's
    tau; for Spearman's rho it carries no probabilistic reading, and CLES is a
    different quantity again (the probability that a draw from one group
    exceeds a draw from another), so the value must not be read as
    P(rank preserved).

    What the transform does provide: a common [0, 1] scale with a chance floor
    of 0.5, so a graph whose rankings are destroyed scores 0.5 and only perfect
    rank inversion reaches 0. Monotonic in r, so it does not reorder graphs.
    """
    # float(), not the numpy scalar scipy hands back: the checkpoint is a
    # pickle other notebooks and tools read, and half its scores arriving as
    # np.float64 and half as float is a portability wart in a published artefact.
    return round(float((r + 1) / 2), 4) if r == r else NA


def spearman(a: np.ndarray, b: np.ndarray, min_pairs: int = 10,
             drop_zero: bool = True) -> tuple:
    """(Spearman r, n pairs actually correlated) between two rankings of pairs.

    With ``drop_zero`` (the default, for neighbourhood heuristics) pairs scoring
    0 both before and after are dropped. They were tied at the bottom before the
    perturbation and are tied at the bottom after it, so they carry no rank
    information and would pad the correlation with identical values on both
    sides. Adamic-Adar, Common Neighbors and Jaccard are all non-negative and
    cannot rise when edges are removed, so that set is exactly the pairs with a
    zero baseline: the same set at every dropout rate, which keeps the
    correlation comparable across the sweep.

    **Embedding scores need ``drop_zero=False``.** TransE and RotatE score a
    triple by a negative distance, so their scores are routinely below zero and
    zero carries no meaning at all. Applying the heuristic mask to them masks
    out every pair and returns nan, which is how this argument came to exist.

    The count is returned because it is not the sample size. Of 500 sampled
    therapeutic pairs, 455 carry Adamic-Adar signal in Hetionet and 419 in
    BioKG; reporting r beside "500 pairs" overstates what was measured.
    """
    if drop_zero:
        mask = (a > 0) | (b > 0)
        a, b = a[mask], b[mask]
    n = len(a)
    r = stats.spearmanr(a, b).statistic if n > min_pairs else NA
    return float(r), n


def top_k_overlap(a: np.ndarray, b: np.ndarray, k: int) -> float:
    """Overlap of the top-k ranked pairs; nan when k covers the whole sample.

    With k >= len(a) both top-k sets are 'every pair' and the overlap is 1.0 by
    construction, which reads as perfect stability rather than as the
    non-measurement it is.
    """
    if k >= len(a):
        return NA
    return len(set(np.argsort(-a)[:k]) & set(np.argsort(-b)[:k])) / k


# --------------------------------------------------------------- perturbation
def periphery_order(G, rng=None) -> list:
    """Edges ordered lowest-degree-first. Deterministic, or replicated.

    With ``rng=None`` the order is a pure function of the graph: the key is
    (minimum endpoint degree, sum of endpoint degrees, sorted endpoint indices).
    Sorting on minimum degree alone leaves tie groups that Python's stable sort
    then resolves by graph insertion order, so at the primary rate on BioKG the
    boundary degree is 4 and its tie group holds 31,301 edges, and most of what
    "periphery-first" removed was decided by the order the loader happened to
    build the graph in.

    With an ``rng``, ties are broken by a seeded permutation INSIDE each
    minimum-degree band instead. Two edges with the same minimum endpoint degree
    are equally peripheral by this metric's own definition, so which of them is
    removed first is arbitrary, and the deterministic key settles it by a rule
    that is reproducible without being meaningful. Repeating the sweep across
    seeds turns that arbitrariness into a standard deviation: the deterministic
    ordering is one draw from this space, and the spread says how much the
    periphery score depends on the draw rather than on the graph.

    The endpoint tuples are canonicalised to (min, max) either way, so nothing
    downstream inherits the insertion order this function exists to remove.
    """
    deg = dict(G.degree())
    edges = [(min(u, v), max(u, v)) for u, v in G.edges()]
    if rng is None:
        return sorted(edges, key=lambda e: (min(deg[e[0]], deg[e[1]]),
                                            deg[e[0]] + deg[e[1]],
                                            e[0], e[1]))
    band = np.fromiter((min(deg[u], deg[v]) for u, v in edges),
                       dtype=np.int64, count=len(edges))
    # lexsort takes the PRIMARY key last: band first, then the permutation.
    order = np.lexsort((rng.permutation(len(edges)), band))
    return [edges[i] for i in order]


def without(G, drop_edges, view_above: int = 1_000_000):
    """G with `drop_edges` removed; a filtered view for very large graphs."""
    if G.number_of_nodes() > view_above:
        dropped = frozenset((min(u, v), max(u, v)) for u, v in drop_edges)
        return nx.subgraph_view(
            G, filter_edge=lambda u, v: (min(u, v), max(u, v)) not in dropped)
    H = G.copy()
    H.remove_edges_from(drop_edges)
    return H


# ------------------------------------------------------- therapeutic pair set
def relation_labels(config, name: str, key: str) -> tuple:
    """Relation labels config/kgs.yaml declares for `key` on this KG.

    The config writes a single relation as `relation:` and several as
    `relations:`; both shapes are accepted.
    """
    spec = (config['knowledge_graphs'][name].get('relations') or {}).get(key)
    if not spec:
        return ()
    return tuple(spec['relations']) if 'relations' in spec else (spec['relation'],)


def entity_labels(config, name: str, canon: str) -> set:
    """KG-internal labels for a canonical entity type (always a tuple in config)."""
    return set(config['knowledge_graphs'][name].get('entity_types', {}).get(canon, ()))


def therapeutic_pairs(kg_df, config, name: str, lcc: set) -> list:
    """(drug, disease) node-index pairs inside the LCC, oriented drug-first.

    Defined here rather than in a notebook cell because the heuristic sweep in
    notebook 05 and the embedding arm in supplementary/ must measure the SAME
    pairs for their results to be comparable at all. Two copies of this rule
    would drift, which is exactly how the hardcoded KG_THERAPEUTIC table drifted
    from the config before it was removed.

    Restricted to the largest component: a pair with an endpoint outside it has
    no path to lose, so its ranking cannot degrade and it would make the graph
    look more stable than it is.
    """
    rels = relation_labels(config, name, 'drug_disease')
    if not rels:
        return []
    drug = entity_labels(config, name, 'Drug')
    disease = entity_labels(config, name, 'Disease')
    sub = (kg_df[kg_df['relation'].isin(rels)][['x_index', 'y_index', 'x_type', 'y_type']]
           .astype({'x_index': int, 'y_index': int}))
    sub = sub[sub['x_index'].isin(lcc) & sub['y_index'].isin(lcc)]
    fwd = sub[sub['x_type'].isin(drug) & sub['y_type'].isin(disease)]
    rev = sub[sub['y_type'].isin(drug) & sub['x_type'].isin(disease)]
    # sorted, not list(set(...)): the pair order decides which pairs the sample
    # draws, so it has to be deterministic rather than set-iteration order.
    return sorted(set(list(zip(fwd['x_index'], fwd['y_index'])) +
                      list(zip(rev['y_index'], rev['x_index']))))


def sample_pairs(pairs: list, n: int, seed: int) -> list:
    """The scored subset. One draw, its own generator, so changing the number of
    dropout trials or of scorers cannot change WHICH pairs are measured."""
    if not pairs:
        return []
    rng = np.random.RandomState(seed)
    return [pairs[i] for i in rng.choice(len(pairs), min(n, len(pairs)), replace=False)]


# ------------------------------------------------- triple loss on an edge table
def row_multiplicity(kg_df):
    """Per-row triple weight for an edge table: every row IS one triple.

    The graph path counts triples per collapsed node pair; on the edge table the
    two coincide, so dropping a fraction of rows is the same perturbation the
    graph path applies. Keeping them the same operation is what makes the
    heuristic and embedding arms comparable.
    """
    return np.ones(len(kg_df), dtype=np.int64)


def edge_table_degrees(kg_df):
    """Node degree counted from the edge table, without building a graph."""
    x = kg_df['x_index'].values.astype(np.int64)
    y = kg_df['y_index'].values.astype(np.int64)
    n = int(max(x.max(initial=0), y.max(initial=0))) + 1
    return np.bincount(x, minlength=n) + np.bincount(y, minlength=n)


def drop_rows(kg_df, rate: float, strategy: str, rng=None):
    """Edge table with `rate` of its triples removed, random or periphery-first.

    Operates on rows, which are triples, so this is the same perturbation as
    `triple_dropout` applies through the graph. Periphery order is
    (minimum endpoint degree, sum of endpoint degrees); with an `rng`, ties
    inside a minimum-degree band are broken by a seeded permutation rather than
    by row order, matching `periphery_order`.
    """
    n = len(kg_df)
    n_drop = int(round(rate * n))
    if n_drop <= 0:
        return kg_df, 0
    if strategy == 'random':
        if rng is None:
            raise ValueError('random dropout needs an rng')
        gone = rng.choice(n, n_drop, replace=False)
    elif strategy == 'periphery':
        deg = edge_table_degrees(kg_df)
        x = kg_df['x_index'].values.astype(np.int64)
        y = kg_df['y_index'].values.astype(np.int64)
        band = np.minimum(deg[x], deg[y])
        tie = (rng.permutation(n) if rng is not None
               else np.lexsort((np.arange(n), deg[x] + deg[y])))
        gone = np.lexsort((tie, band))[:n_drop]
    else:
        raise ValueError(f'unknown strategy {strategy!r}')
    keep = np.ones(n, dtype=bool)
    keep[gone] = False
    return kg_df[keep], int(n_drop)


# ----------------------------------------------------------- triple-level loss
def indexed_edges(G) -> tuple:
    """(edge list, per-edge triple count) in one pass over G.edges().

    `n_triples` is the number of rows of the edge table that asserted the pair.
    Edges without the attribute fall back to 1, which gives pair-level dropout.
    """
    edges, mult = [], []
    for u, v, d in G.edges(data=True):
        edges.append((u, v))
        mult.append(int(d.get('n_triples', 1)))
    return edges, np.asarray(mult, dtype=np.int64)


def triple_dropout(edges, mult, rate: float, rng) -> tuple:
    """Edges to remove when `rate` of the underlying TRIPLES is dropped.

    Incompleteness in a real KG is missing ASSERTIONS, not missing node pairs,
    and a pair supported by several predicates should survive losing one of
    them. Dropping whole pairs made a triply-asserted pair exactly as easy to
    sever as a singly-asserted one: severed with probability `rate` instead of
    roughly `rate ** 3`.

    Triples are drawn uniformly without replacement across the whole graph, and
    an edge disappears only when every triple supporting it is gone. Returns
    (edges to remove, triples dropped).
    """
    n_triples = int(mult.sum())
    n_drop = int(round(rate * n_triples))
    if n_drop <= 0:
        return [], 0
    owner = np.repeat(np.arange(len(edges), dtype=np.int64), mult)
    hit = np.bincount(owner[rng.choice(n_triples, n_drop, replace=False)],
                      minlength=len(edges))
    return [edges[i] for i in np.flatnonzero(hit == mult)], n_drop


def periphery_triples(order, mult_by_edge: dict, rate: float,
                      n_triples: int) -> tuple:
    """Lowest-degree edges first, spending a budget counted in triples.

    The budget is in triples so the two strategies remove comparable amounts of
    evidence at the same nominal rate. A pair asserted twice costs two of the
    budget, exactly as it would cost two draws under random dropout.
    """
    budget = int(round(rate * n_triples))
    if budget <= 0:
        return [], 0
    cost = np.fromiter((mult_by_edge.get(e, 1) for e in order),
                       dtype=np.int64, count=len(order))
    cum = np.cumsum(cost)
    k = int(np.searchsorted(cum, budget, side='left')) + 1
    k = min(k, len(order))
    return list(order[:k]), int(cum[k - 1])


# ------------------------------------------------------------------ null pairs
def degree_matched_pairs(G, gold, drug_pool, disease_pool, rng,
                         window: int = DEGREE_MATCH_WINDOW,
                         exclude=None) -> list:
    """Random drug-disease pairs whose endpoint degrees match the gold pairs.

    For each gold pair, each endpoint is replaced by one of the `window` nodes
    nearest it in degree within the same entity pool, drawn uniformly. Matching
    on degree is what makes the null informative: Adamic-Adar is a function of
    shared-neighbour degrees, so an unmatched random pair would differ from a
    therapeutic pair in degree as well as in biology, and could not separate the
    two.

    `exclude` is the set of pairs that must not be drawn, and it has to be every
    known therapeutic pair in the graph, not the sample being matched. Banning
    only `gold` leaves the rest eligible: with 500 pairs sampled from a graph
    holding tens of thousands, almost every real indication can be drawn into
    the null. The null then contains genuine therapeutic pairs, which pulls it
    toward gold and shrinks the gold-minus-null gap that carries the biological
    claim. Callers pass the full pair set; the default falls back to `gold` so
    an older caller is no worse off than before, but it is not the right
    argument.
    """
    deg = dict(G.degree())
    _ban_src = gold if exclude is None else exclude
    banned = set(_ban_src) | {(v, u) for u, v in _ban_src}

    def sorted_pool(pool):
        pairs = sorted((deg[v], v) for v in pool)
        return [d for d, _ in pairs], [v for _, v in pairs]

    ddeg, dnode = sorted_pool(drug_pool)
    sdeg, snode = sorted_pool(disease_pool)
    if not dnode or not snode:
        return []

    def near(degrees, nodes, target):
        i = bisect.bisect_left(degrees, target)
        lo = max(0, min(i - window // 2, len(nodes) - window))
        hi = min(len(nodes), lo + window)
        return nodes[lo:hi]

    out = []
    for u, v in gold:
        du, dv = near(ddeg, dnode, deg[u]), near(sdeg, snode, deg[v])
        for _ in range(20):
            cand = (du[rng.randint(len(du))], dv[rng.randint(len(dv))])
            if cand[0] != cand[1] and cand not in banned:
                out.append(cand)
                break
    return out


def degree_match_quality(G, gold, null) -> dict:
    """How closely the null's endpoint degrees track the gold pairs'.

    Reported so the null can be judged rather than trusted: a null that failed to
    match degrees would make the comparison meaningless in the flattering
    direction.
    """
    deg = dict(G.degree())
    if not null:
        return {'n': 0, 'median_abs_log2_ratio': NA}
    ratios = [abs(np.log2(max(deg[b], 1) / max(deg[a], 1)))
              for (a, _), (b, _) in zip(gold, null)]
    ratios += [abs(np.log2(max(deg[b], 1) / max(deg[a], 1)))
               for (_, a), (_, b) in zip(gold, null)]
    return {'n': len(null), 'median_abs_log2_ratio': round(float(np.median(ratios)), 4)}
