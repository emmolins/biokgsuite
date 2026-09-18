"""Graph heuristic scorers for link prediction.

Each scorer follows the signature:
    scorer(G, pairs) -> np.ndarray of float scores

where pairs is a list of (u, v) node index tuples and higher scores
indicate a stronger predicted link.

Reference:
    Liben-Nowell, D. & Kleinberg, J. "The link-prediction problem for
    social networks." JASIST 58(7), 1019–1031 (2007).
    Adamic, L.A. & Adar, E. "Friends and neighbors on the Web."
    Social Networks 25(3), 211–230 (2003).
"""

import numpy as np
import networkx as nx


def score_common_neighbors(G, pairs):
    """Number of shared neighbours between u and v.

    CN(u, v) = |N(u) ∩ N(v)|
    """
    scores = np.zeros(len(pairs), dtype=float)
    for i, (u, v) in enumerate(pairs):
        if G.has_node(u) and G.has_node(v):
            scores[i] = sum(1 for _ in nx.common_neighbors(G, u, v))
    return scores


def score_jaccard(G, pairs):
    """Jaccard coefficient: common neighbours normalised by union size.

    J(u, v) = |N(u) ∩ N(v)| / |N(u) ∪ N(v)|
    """
    scores = np.zeros(len(pairs), dtype=float)
    for i, (u, v) in enumerate(pairs):
        if G.has_node(u) and G.has_node(v):
            nu, nv = set(G.neighbors(u)), set(G.neighbors(v))
            denom  = len(nu | nv)
            scores[i] = len(nu & nv) / denom if denom > 0 else 0.0
    return scores


def score_adamic_adar(G, pairs):
    """Adamic-Adar score: down-weights high-degree shared neighbours.

    AA(u, v) = Σ_{w ∈ N(u) ∩ N(v)} 1 / log|N(w)|
    """
    # Precompute inverse log-degree once for all nodes
    inv_log_deg = {}
    for n in G.nodes():
        d = G.degree(n)
        inv_log_deg[n] = 1.0 / np.log(d) if d > 1 else 0.0

    scores = np.zeros(len(pairs), dtype=float)
    for i, (u, v) in enumerate(pairs):
        if G.has_node(u) and G.has_node(v):
            shared    = set(G.neighbors(u)) & set(G.neighbors(v))
            scores[i] = sum(inv_log_deg.get(w, 0.0) for w in shared)
    return scores


# ── Registry ──────────────────────────────────────────────────────────────────

SCORERS = {
    'Common Neighbors': score_common_neighbors,
    'Jaccard':          score_jaccard,
    'Adamic-Adar':      score_adamic_adar,
}
