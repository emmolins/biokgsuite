"""Unit tests for biokgsuite.stats.stability.

These cover the properties the dimension's claims rest on, on graphs small
enough to reason about by hand: that the chance floor is where the notebook says
it is, that the number of pairs actually correlated is reported rather than
assumed, that top-K refuses to report a degenerate 1.0, that periphery order is
a function of the graph and not of its construction, and that the null matches
degrees.

    python -m pytest tests/test_stability.py -q
"""
import networkx as nx
import numpy as np
import pytest

from biokgsuite.stats import stability as st


def test_cles_chance_floor_is_half():
    # The scale's floor is 0.5, not 0. A run that scored a failure as 0.0 would
    # place it below every possible real result.
    assert st.cles(0.0) == 0.5
    assert st.cles(1.0) == 1.0
    assert st.cles(-1.0) == 0.0
    assert st.cles(float('nan')) != st.cles(float('nan'))     # nan propagates


def test_spearman_reports_pairs_correlated_not_pairs_passed():
    base = np.array([3.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    pert = np.array([3.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    r, n = st.spearman(base, pert, min_pairs=2)
    assert n == 3                       # nine both-zero pairs carry no rank information
    assert r == pytest.approx(1.0)


def test_spearman_undefined_below_threshold_rather_than_zero():
    a = np.array([1.0, 2.0, 0.0, 0.0])
    r, n = st.spearman(a, a, min_pairs=10)
    assert r != r                       # nan, so the notebook can raise
    assert n == 2


def test_top_k_overlap_is_undefined_when_k_covers_the_sample():
    a = np.array([3.0, 2.0, 1.0])
    assert st.top_k_overlap(a, a, k=3) != st.top_k_overlap(a, a, k=3)   # nan
    assert st.top_k_overlap(a, a, k=2) == 1.0


def test_top_k_overlap_detects_reordering():
    a = np.array([4.0, 3.0, 2.0, 1.0])
    b = np.array([1.0, 2.0, 3.0, 4.0])
    assert st.top_k_overlap(a, b, k=2) == 0.0


def _diamond():
    G = nx.Graph()
    G.add_edges_from([(1, 2), (1, 3), (2, 4), (3, 4), (4, 5), (5, 6)])
    return G


def test_periphery_order_is_independent_of_insertion_order():
    G = _diamond()
    H = nx.Graph()
    H.add_edges_from(reversed(list(G.edges())))
    assert st.periphery_order(G) == st.periphery_order(H)


def test_periphery_order_takes_lowest_degree_first():
    G = _diamond()
    deg = dict(G.degree())
    mins = [min(deg[u], deg[v]) for u, v in st.periphery_order(G)]
    assert mins == sorted(mins)


def test_without_removes_only_the_named_edges():
    G = _diamond()
    H = st.without(G, [(1, 2)])
    assert not H.has_edge(1, 2)
    assert H.number_of_edges() == G.number_of_edges() - 1
    assert G.has_edge(1, 2)             # the original is not mutated


def test_degree_matched_null_matches_degrees_and_avoids_gold():
    G = nx.Graph()
    for d in range(1, 9):                       # drugs d1..d8 of varying degree
        for k in range(d):
            G.add_edge(f'd{d}', f'g{k}')
    for s in range(1, 9):                       # diseases s1..s8 likewise
        for k in range(s):
            G.add_edge(f's{s}', f'g{k}')
    gold = [('d5', 's5')]
    null = st.degree_matched_pairs(G, gold, [f'd{i}' for i in range(1, 9)],
                                   [f's{i}' for i in range(1, 9)],
                                   np.random.RandomState(0), window=3)
    assert len(null) == 1
    assert null[0] not in set(gold)
    q = st.degree_match_quality(G, gold, null)
    assert q['median_abs_log2_ratio'] < 1.0     # within a factor of two on degree


def test_degree_matched_null_is_deterministic_given_a_seed():
    G = nx.complete_bipartite_graph(['d1', 'd2', 'd3'], ['s1', 's2', 's3'])
    args = (G, [('d1', 's1')], ['d1', 'd2', 'd3'], ['s1', 's2', 's3'])
    a = st.degree_matched_pairs(*args, np.random.RandomState(7))
    b = st.degree_matched_pairs(*args, np.random.RandomState(7))
    assert a == b


# ---------------------------------------------------------- triple-level loss
def _multi_graph():
    """Pair (1,2) asserted three times, every other pair once."""
    import pandas as pd
    from biokgsuite.graph.graph_utils import build_graph
    kg = pd.DataFrame({
        'x_index': [1, 1, 1, 2, 3, 4],
        'y_index': [2, 2, 2, 3, 4, 5],
        'relation': ['a', 'b', 'c', 'd', 'e', 'f'],
    })
    nodes = pd.DataFrame({'idx': [1, 2, 3, 4, 5], 'type': ['T'] * 5,
                          'name': list('vwxyz')})
    return build_graph(kg, nodes)


def test_build_graph_keeps_every_predicate_and_counts_triples():
    G = _multi_graph()
    assert G[1][2]['relations'] == frozenset({'a', 'b', 'c'})
    assert G[1][2]['n_triples'] == 3
    assert sum(d['n_triples'] for *_, d in G.edges(data=True)) == 6


def test_build_graph_counts_assertions_not_rows():
    """A predicate stored in both directions (PrimeKG) or repeated is one
    assertion; two different predicates on one pair are two."""
    import pandas as pd
    from biokgsuite.graph.graph_utils import build_graph
    kg = pd.DataFrame({
        'x_index': [1, 2, 1, 3, 3, 4],
        'y_index': [2, 1, 2, 4, 4, 3],
        'relation': ['a', 'a', 'a', 'd', 'e', 'd'],
    })
    nodes = pd.DataFrame({'idx': [1, 2, 3, 4], 'type': ['T'] * 4, 'name': list('vwxy')})
    G = build_graph(kg, nodes)
    assert G[1][2]['n_triples'] == 1 and G[1][2]['relations'] == frozenset({'a'})
    assert G[3][4]['n_triples'] == 2 and G[3][4]['relations'] == frozenset({'d', 'e'})


def test_indexed_edges_defaults_to_one_triple_when_absent():
    G = nx.Graph()
    G.add_edge(1, 2)
    _, mult = st.indexed_edges(G)
    assert list(mult) == [1]


def test_triple_dropout_spares_a_multiply_asserted_pair():
    G = _multi_graph()
    edges, mult = st.indexed_edges(G)
    i12 = edges.index((1, 2))
    # Drop half the triples, many times: the triple pair should survive far
    # more often than a singly-asserted one.
    survived_multi = survived_single = 0
    i45 = edges.index((4, 5))
    for seed in range(60):
        gone = set(st.triple_dropout(edges, mult, 0.5, np.random.RandomState(seed))[0])
        survived_multi += (edges[i12] not in gone)
        survived_single += (edges[i45] not in gone)
    assert survived_multi > survived_single


def test_triple_dropout_budget_is_counted_in_triples():
    G = _multi_graph()
    edges, mult = st.indexed_edges(G)
    _, dropped = st.triple_dropout(edges, mult, 0.5, np.random.RandomState(0))
    assert dropped == 3                      # 50% of six triples, not of four edges


def test_triple_dropout_removes_an_edge_only_when_all_its_triples_go():
    G = _multi_graph()
    edges, mult = st.indexed_edges(G)
    removed, dropped = st.triple_dropout(edges, mult, 1.0, np.random.RandomState(0))
    assert dropped == 6 and set(removed) == set(edges)


def test_periphery_budget_is_also_counted_in_triples():
    G = _multi_graph()
    edges, mult = st.indexed_edges(G)
    by_edge = dict(zip(edges, mult))
    order = st.periphery_order(G)
    chosen, spent = st.periphery_triples(order, by_edge, 0.5, int(mult.sum()))
    assert spent >= 3                        # at least the budget, in triples
    assert all(e in set(edges) for e in chosen)


# ------------------------------------------------- replicated periphery order
def test_periphery_order_with_rng_still_takes_lowest_degree_first():
    G = _diamond()
    deg = dict(G.degree())
    order = st.periphery_order(G, rng=np.random.RandomState(3))
    mins = [min(deg[u], deg[v]) for u, v in order]
    assert mins == sorted(mins)          # bands are respected, only ties shuffle


def test_periphery_order_with_rng_is_reproducible_per_seed():
    G = _diamond()
    a = st.periphery_order(G, rng=np.random.RandomState(11))
    b = st.periphery_order(G, rng=np.random.RandomState(11))
    assert a == b


def test_periphery_order_seeds_disagree_inside_a_tie_band():
    # A star plus a clique gives one large band of degree-1 edges to shuffle.
    G = nx.Graph()
    for i in range(2, 12):
        G.add_edge(1, i)                 # ten leaves, all min-degree 1
    orders = {tuple(st.periphery_order(G, rng=np.random.RandomState(s)))
              for s in range(8)}
    assert len(orders) > 1               # the tie-break is genuinely arbitrary


def test_periphery_order_covers_the_same_edges_however_ties_break():
    G = _diamond()
    det = set(st.periphery_order(G))
    rnd = set(st.periphery_order(G, rng=np.random.RandomState(5)))
    assert det == rnd


def test_spearman_drop_zero_is_wrong_for_signed_scores():
    # Embedding scores are negative distances; the heuristic zero-mask would
    # discard every pair and return nan.
    a = np.linspace(-5, -1, 40)
    b = a + np.random.RandomState(0).normal(0, 0.01, 40)
    r_masked, n_masked = st.spearman(a, b)
    assert n_masked == 0 and r_masked != r_masked        # nan, all masked out
    r_all, n_all = st.spearman(a, b, drop_zero=False)
    assert n_all == 40 and r_all > 0.9


# ── graph_utils.find_node: entity types arrive as tuples ─────────────────────
# load_config normalises every configured entity type to a tuple, and find_node
# compared with `==`. A tuple caller therefore got None for every node, with no
# exception: notebook 07's prospective section resolved zero drugs in Hetionet
# for this reason and reported it as the KG lacking them.

def test_find_node_accepts_plain_string_type():
    from biokgsuite.graph.graph_utils import find_node
    nm, tm = {1: 'Aspirin'}, {1: 'Compound'}
    assert find_node('aspirin', nm, tm, 'Compound') == 1


def test_find_node_accepts_tuple_type():
    from biokgsuite.graph.graph_utils import find_node
    nm, tm = {1: 'Aspirin'}, {1: 'Compound'}
    assert find_node('aspirin', nm, tm, ('Compound',)) == 1


def test_find_node_matches_every_label_in_a_multi_label_type():
    """BioKG's Disease is ('Disease', 'GeneticDisorder'); both carry nodes, so
    taking the first label would silently drop the second."""
    from biokgsuite.graph.graph_utils import find_node
    nm = {1: 'Asthma', 2: 'Rare thing'}
    tm = {1: 'Disease', 2: 'GeneticDisorder'}
    types = ('Disease', 'GeneticDisorder')
    assert find_node('asthma', nm, tm, types) == 1
    assert find_node('rare thing', nm, tm, types) == 2


def test_find_node_still_rejects_the_wrong_type():
    from biokgsuite.graph.graph_utils import find_node
    nm, tm = {1: 'Aspirin'}, {1: 'Compound'}
    assert find_node('aspirin', nm, tm, ('Disease',)) is None
    assert find_node('aspirin', nm, tm, None) is None
