"""Negative sampling strategies for link prediction evaluation.

Three strategies of increasing difficulty are implemented, following:
    Kotnis, B. & Nastase, V. "Analysis of the impact of negative sampling
    on link prediction in knowledge graphs." AKBC Workshop, 2017.
"""

import numpy as np


def generate_negatives(pos_pairs, n_neg, strategy,
                       drug_idx, disease_idx, drug_targets,
                       node_name_map, all_pos_set, rng):
    """Generate negative drug-disease pairs.

    Parameters
    ----------
    pos_pairs : list of (drug_idx, disease_idx)
        Known positive pairs (used as seeds for shared-target strategy).
    n_neg : int
        Number of negative pairs to generate.
    strategy : str
        One of 'random', 'type-constrained', 'shared-target'.
    drug_idx : set of int
        All drug node indices in the KG.
    disease_idx : set of int
        All disease node indices in the KG.
    drug_targets : dict {drug_idx: set of gene_idx}
        Maps each drug to its known protein targets.
    node_name_map : dict {idx: name}
        Node index to name lookup (unused here; kept for API compatibility).
    all_pos_set : set of (drug_idx, disease_idx)
        Full positive set — negatives are drawn strictly outside this set.
    rng : np.random.RandomState

    Returns
    -------
    list of (drug_idx, disease_idx) — length min(n_neg, feasible pairs found)
    """
    drug_list    = np.array(sorted(drug_idx),    dtype=int)
    disease_list = np.array(sorted(disease_idx), dtype=int)
    negatives    = set()
    max_attempts = n_neg * 30

    if strategy in ('random', 'type-constrained'):
        # Random: uniform sample over drug × disease; type-constrained is
        # equivalent here since we already restrict to drug/disease types.
        attempts = 0
        while len(negatives) < n_neg and attempts < max_attempts:
            d   = int(drug_list[rng.randint(len(drug_list))])
            dis = int(disease_list[rng.randint(len(disease_list))])
            if (d, dis) not in all_pos_set:
                negatives.add((d, dis))
            attempts += 1

    elif strategy == 'shared-target':
        # Hard negatives: drug shares ≥1 protein target with a positive drug
        # but is NOT known to treat the sampled disease — tests whether the KG
        # can distinguish therapeutically similar but non-indicated pairs.
        drug_targets_inv: dict[int, set] = {}
        for d, targets in drug_targets.items():
            for t in targets:
                drug_targets_inv.setdefault(t, set()).add(d)

        attempts = 0
        while len(negatives) < n_neg and attempts < max_attempts:
            src_d, _ = pos_pairs[rng.randint(len(pos_pairs))]
            src_tgts = list(drug_targets.get(src_d, set()))
            if not src_tgts:
                attempts += 1
                continue
            t     = src_tgts[rng.randint(len(src_tgts))]
            peers = list(drug_targets_inv.get(t, set()) - {src_d})
            if not peers:
                attempts += 1
                continue
            d   = peers[rng.randint(len(peers))]
            dis = int(disease_list[rng.randint(len(disease_list))])
            if (d, dis) not in all_pos_set:
                negatives.add((d, dis))
            attempts += 1

    else:
        raise ValueError(f"Unknown strategy '{strategy}'. "
                         "Choose from: 'random', 'type-constrained', 'shared-target'.")

    return list(negatives)[:n_neg]
