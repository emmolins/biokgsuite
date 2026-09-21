"""Graph construction utilities for BioKGSuite."""

import networkx as nx
import numpy as np
import pandas as pd


def build_graph(kg: pd.DataFrame, nodes: pd.DataFrame) -> nx.Graph:
    """Build an undirected NetworkX graph from edge and node tables.

    The graph is SIMPLE, not a multigraph, because every metric computed on it
    is defined that way: Adamic-Adar weights each shared neighbour by
    1/log(degree), and a neighbour joined by two predicates is not twice as
    promiscuous; clustering coefficient, Louvain modularity and largest-connected
    -component size are likewise defined on simple graphs, and networkx's
    ``clustering`` refuses a MultiGraph outright.

    Collapsing pairs does not mean discarding what supported them. Each edge
    carries:

    ``relations``   frozenset of every predicate asserted between the pair
    ``n_triples``   how many DISTINCT predicates assert it (len(relations))

    The previous version stored a single ``relation`` and let the LAST row win,
    silently dropping the others: 142,473 Hetionet pairs carry two predicates
    and 15 carry three, so one of each was lost. Keeping the full set also lets
    a caller perturb the graph at the level of ASSERTIONS rather than of pairs,
    which is what incompleteness in a real KG looks like — a pair supported by
    three predicates should survive losing one of them.

    ``n_triples`` counts distinct predicates, not rows. Rows overcount: PrimeKG
    stores every undirected edge in both directions, so counting rows gave every
    PrimeKG pair a multiplicity of exactly 2.000 and a 30% triple loss severed
    only 9% of its pairs, against 27 to 30% for graphs that store each edge
    once. A predicate repeated in the reverse direction, or repeated outright,
    is one assertion; a pair carried by two different predicates (DRKG's
    DrugBank and GNBR relations, say) is two.

    The frozensets are interned: distinct predicate sets number in the hundreds
    while edges number in the millions, so one shared object per distinct set
    keeps this to a pointer per edge instead of a new set per edge.

    Parameters
    ----------
    kg : pd.DataFrame
        Edge table with columns x_index, y_index, and (optionally) relation.
    nodes : pd.DataFrame
        Node table with columns idx, type, name.

    Returns
    -------
    nx.Graph
        Undirected. Node attributes: 'node_type', 'name'. Edge attributes (when
        the kg DataFrame has a 'relation' column): 'relations', 'n_triples'.
    """
    G = nx.Graph()

    # Add nodes with attributes (itertuples is ~50x faster than iterrows)
    G.add_nodes_from(
        (row.idx, {'node_type': row.type, 'name': row.name})
        for row in nodes.itertuples(index=False)
    )

    x = kg['x_index'].values.astype(np.int64)
    y = kg['y_index'].values.astype(np.int64)
    lo = np.minimum(x, y)
    hi = np.maximum(x, y)

    if 'relation' not in kg.columns:
        G.add_edges_from(zip(lo, hi))
        return G

    # One code per unordered pair, so multiplicity and predicate sets can be
    # aggregated with numpy rather than a Python loop over millions of rows.
    span = int(max(hi.max(initial=0), lo.max(initial=0))) + 1
    codes, uniq = pd.factorize(lo * span + hi, sort=False)
    n_pairs = len(uniq)
    mult = np.bincount(codes, minlength=n_pairs)

    rel = kg['relation'].values
    first = np.full(n_pairs, -1, dtype=np.int64)
    first[codes[::-1]] = np.arange(len(codes) - 1, -1, -1)   # first row per pair

    intern: dict = {}

    def share(labels) -> frozenset:
        fs = frozenset(labels)
        return intern.setdefault(fs, fs)

    # Singletons are the overwhelming majority and need no grouping.
    sets = np.empty(n_pairs, dtype=object)
    singles = np.flatnonzero(mult == 1)
    for i in singles:
        sets[i] = share((rel[first[i]],))

    multi = np.flatnonzero(mult > 1)
    if multi.size:
        keep = np.isin(codes, multi)
        grouped = pd.Series(rel[keep]).groupby(codes[keep]).agg(lambda s: share(s))
        for code, fs in grouped.items():
            sets[code] = fs

    pair_lo = uniq // span
    pair_hi = uniq % span
    G.add_edges_from(
        (int(u), int(v), {'relations': sets[i], 'n_triples': len(sets[i])})
        for i, (u, v) in enumerate(zip(pair_lo, pair_hi))
    )
    return G


def build_lookup_maps(nodes: pd.DataFrame) -> dict:
    """Build quick-lookup dictionaries from the node table.

    Parameters
    ----------
    nodes : pd.DataFrame
        Node table with columns idx, type, name (and optionally uid).

    Returns
    -------
    dict with keys:
        'node_type_map' : {idx -> type}
        'node_name_map' : {idx -> name}
        'uid_to_idx'    : {uid -> idx}  (only if ``uid`` column exists)
    """
    maps = {
        'node_type_map': dict(zip(nodes['idx'], nodes['type'])),
        'node_name_map': dict(zip(nodes['idx'], nodes['name'])),
    }
    if 'uid' in nodes.columns:
        maps['uid_to_idx'] = dict(zip(nodes['uid'], nodes['idx']))
    return maps


def node_id_lookup(nodes: pd.DataFrame, node_type: str, node_id: str):
    """Look up a node index by type and external ID — collision-safe.

    This is the recommended way to map external IDs (e.g. DrugBank, MONDO)
    to node indices in KGs where the same external ID string can refer to
    different entity types (notably PrimeKG).

    Parameters
    ----------
    nodes : pd.DataFrame
        Node table with columns idx, id, type.
    node_type : str
        Entity type (e.g. 'drug', 'disease', 'gene/protein').
    node_id : str
        External identifier (e.g. 'DB15119', '9891').

    Returns
    -------
    int or None — node index if found, else None.
    """
    mask = (nodes['type'] == node_type) & (nodes['id'].astype(str) == str(node_id))
    match = nodes.loc[mask, 'idx']
    if len(match) == 1:
        return int(match.iloc[0])
    elif len(match) > 1:
        # Multiple matches (shouldn't happen but be safe) — return first.
        return int(match.iloc[0])
    return None


def find_node(name: str, node_name_map: dict, node_type_map: dict,
              node_type):
    """Find a node index by name and type (case-insensitive).

    Parameters
    ----------
    name : str
        Node name to search for.
    node_name_map : dict {idx -> name}
    node_type_map : dict {idx -> type}
    node_type : str or tuple/list/set of str
        Expected node type(s), e.g. 'drug' or ('Disease', 'GeneticDisorder').

        `load_config` normalises every configured entity type to a TUPLE of
        labels, so callers reading types from config pass a tuple here. This
        compared with `==`, which is never true against a tuple, so such a call
        returned None for every node and did so silently: no exception, just an
        empty lookup that reads as "the KG does not contain this entity".
        Membership is required rather than indexing [0], because a type like
        BioKG's ('Disease', 'GeneticDisorder') carries real nodes under both
        labels and taking the first would drop the second.

    Returns
    -------
    int or None — node index if found, else None.
    """
    types = ({node_type} if isinstance(node_type, str)
             else set(node_type) if node_type is not None else set())
    name_lower = name.lower()
    for idx, n in node_name_map.items():
        if str(n).lower() == name_lower and node_type_map.get(idx) in types:
            return idx
    return None
