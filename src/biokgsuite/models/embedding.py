"""TransE and RotatE implementations for KG link prediction.

PyTorch implementations that mirror the standard formulations and are
GPU-native (auto-detects CUDA > MPS > CPU). The class interfaces, checkpoint
formats, and all downstream code are unchanged from the original NumPy version.

TransE:
    Bordes, A. et al. "Translating embeddings for modeling multi-relational
    data." NeurIPS 2013.
    Score: -||h + r - t||

RotatE:
    Sun, Z. et al. "RotatE: Knowledge Graph Embedding by Relational Rotation
    in Complex Space." ICLR 2019.
    Score: -||h o r - t||  (complex Hadamard product)

Both use margin-based ranking loss (hinge) with uniform negative sampling
and vanilla SGD with sparse gradient updates (only touched entity rows are
updated each step), identical in spirit to the original per-entity SGD.
"""

import math
import os

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from pathlib import Path


# ── Device selection ──────────────────────────────────────────────────────────

def _get_device() -> torch.device:
    """Auto-detect the best available device: CUDA > MPS > CPU."""
    if torch.cuda.is_available():
        return torch.device("cuda")
    if getattr(torch.backends, "mps", None) and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


# ── Helpers ───────────────────────────────────────────────────────────────────

def _l2_normalize_rows(x, eps=1e-12):
    """L2-normalise each row of a numpy array in-place (kept for callers outside
    this module that may still import it)."""
    norms = np.linalg.norm(x, axis=1, keepdims=True)
    np.maximum(norms, eps, out=norms)
    x /= norms
    return x


def _save_ckpt_atomic(path, **arrays):
    """Atomic-rename save so a crash mid-write doesn't corrupt the checkpoint.

    Accepts numpy arrays as kwargs and writes them with np.savez (uncompressed).
    Uses a temp file + os.replace so the checkpoint is never half-written.
    """
    path = Path(path)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "wb") as f:
        np.savez(f, **arrays)
    os.replace(tmp, path)


def _normalize_emb_inplace(emb: nn.Embedding):
    """L2-normalise entity embedding rows in-place (no gradient tracking)."""
    with torch.no_grad():
        F.normalize(emb.weight.data, p=2, dim=1, out=emb.weight.data)


# ── TransE ────────────────────────────────────────────────────────────────────

class TransE:
    """TransE with sparse SGD on the best available device (CUDA / MPS / CPU).

    The public interface is identical to the original NumPy version:
      - ``score(h_idx, r_idx, t_idx)``  → numpy array
      - ``score_pairs(pairs, rel_idx)``  → numpy array
      - ``fit(triples, ...)``
      - ``ent_emb`` / ``rel_emb`` properties return float32 numpy arrays
      - checkpoint files use the same .npz keys (ent_emb, rel_emb,
        epoch_completed) so existing checkpoints are loaded transparently.

    Parameters
    ----------
    n_entities, n_relations : int
    dim : int
        Embedding dimension.
    margin : float
    lr : float
        Learning rate for vanilla SGD.
    seed : int
    """

    name = "TransE"

    def __init__(self, n_entities, n_relations, dim=128, margin=1.0,
                 lr=0.01, seed=42):
        self.dim = dim
        self.margin = float(margin)
        self.lr = lr
        self.n_entities = n_entities
        self.n_relations = n_relations
        self._device = _get_device()
        self._rng = np.random.RandomState(seed)

        # Initialise embeddings with the same Xavier-uniform distribution as
        # the original NumPy version.
        bound = 6.0 / math.sqrt(dim)
        ent_w = torch.from_numpy(
            self._rng.uniform(-bound, bound, (n_entities, dim)).astype(np.float32))
        ent_w = F.normalize(ent_w, p=2, dim=1)
        rel_w = torch.from_numpy(
            self._rng.uniform(-bound, bound, (n_relations, dim)).astype(np.float32))

        # sparse=True → SGD only touches the rows that appear in each batch,
        # which is crucial for large entity tables (e.g. Matrix at 4.8M ents).
        self._ent = nn.Embedding(n_entities, dim, sparse=True, _weight=ent_w.clone())
        self._rel = nn.Embedding(n_relations, dim, sparse=True, _weight=rel_w.clone())
        self._ent.to(self._device)
        self._rel.to(self._device)

        self._opt = torch.optim.SGD(
            list(self._ent.parameters()) + list(self._rel.parameters()), lr=lr)

    # ── Numpy-compatible properties ──────────────────────────────────────────

    @property
    def ent_emb(self) -> np.ndarray:
        """Entity embedding matrix as a (n_entities, dim) float32 numpy array."""
        return self._ent.weight.detach().cpu().numpy()

    @property
    def rel_emb(self) -> np.ndarray:
        """Relation embedding matrix as a (n_relations, dim) float32 numpy array."""
        return self._rel.weight.detach().cpu().numpy()

    # ── Scoring ──────────────────────────────────────────────────────────────

    def score(self, h_idx, r_idx, t_idx) -> np.ndarray:
        """Return -||h + r - t||₂  (higher = more plausible)."""
        with torch.no_grad():
            dev = self._device
            h = self._ent.weight[torch.as_tensor(np.asarray(h_idx, dtype=np.int64), device=dev)]
            r = self._rel.weight[torch.as_tensor(np.asarray(r_idx, dtype=np.int64), device=dev)]
            t = self._ent.weight[torch.as_tensor(np.asarray(t_idx, dtype=np.int64), device=dev)]
            return -(h + r - t).norm(p=2, dim=-1).cpu().numpy()

    def score_pairs(self, pairs, rel_idx) -> np.ndarray:
        pairs = np.asarray(pairs)
        h_idx, t_idx = pairs[:, 0], pairs[:, 1]
        if np.isscalar(rel_idx):
            rel_idx = np.full(len(pairs), rel_idx, dtype=np.int64)
        return self.score(h_idx, rel_idx, t_idx)

    # ── Training ─────────────────────────────────────────────────────────────

    def fit(self, triples, n_epochs=50, batch_size=16384, verbose=True,
            checkpoint_path=None, checkpoint_every=10):
        """Train with hinge (margin) ranking loss and sparse SGD.

        Semantics are identical to the original NumPy ``fit()``:
          - Uniform negative sampling (corrupt head or tail 50/50).
          - Only violated triplets contribute to the loss (F.relu clamps
            non-violated samples to zero, so their gradient is 0).
          - Entity embeddings are L2-normalised after every epoch.
          - Checkpoint resume: if a .npz checkpoint with matching shapes
            exists at ``checkpoint_path``, training resumes from the saved
            epoch.
        """
        if verbose:
            print(f"  [{self.name}] device={self._device}  "
                  f"entities={self.n_entities:,}  dim={self.dim}  "
                  f"epochs={n_epochs}  batch={batch_size}", flush=True)

        triples = np.asarray(triples, dtype=np.int32)
        n = len(triples)
        dev = self._device

        # ── Resume from checkpoint ──────────────────────────────────────────
        start_epoch = 0
        if checkpoint_path and Path(checkpoint_path).exists():
            try:
                ckpt = np.load(checkpoint_path)
                if (ckpt["ent_emb"].shape == (self.n_entities, self.dim) and
                        ckpt["rel_emb"].shape == (self.n_relations, self.dim)):
                    self._ent.weight.data.copy_(
                        torch.from_numpy(ckpt["ent_emb"].astype(np.float32)))
                    self._rel.weight.data.copy_(
                        torch.from_numpy(ckpt["rel_emb"].astype(np.float32)))
                    start_epoch = int(ckpt["epoch_completed"])
                    if verbose:
                        print(f"  Resumed from checkpoint at epoch "
                              f"{start_epoch}/{n_epochs}", flush=True)
                else:
                    if verbose:
                        print("  Checkpoint shape mismatch — starting fresh",
                              flush=True)
            except Exception as e:
                if verbose:
                    print(f"  Checkpoint load failed ({e}) — starting fresh",
                          flush=True)

        # ── Training loop ───────────────────────────────────────────────────
        for epoch in range(start_epoch, n_epochs):
            perm = self._rng.permutation(n)
            epoch_loss = 0.0
            n_batches = 0

            for start in range(0, n, batch_size):
                batch = triples[perm[start:start + batch_size]]
                bs = len(batch)
                heads = batch[:, 0].astype(np.int64)
                rels  = batch[:, 1].astype(np.int64)
                tails = batch[:, 2].astype(np.int64)

                # Corrupt head or tail (50/50)
                neg_heads = heads.copy()
                neg_tails = tails.copy()
                flip = self._rng.random(bs) < 0.5
                n_h = flip.sum()
                if n_h > 0:
                    neg_heads[flip] = self._rng.randint(0, self.n_entities, n_h)
                if bs - n_h > 0:
                    neg_tails[~flip] = self._rng.randint(
                        0, self.n_entities, bs - n_h)

                # Transfer index tensors to device
                heads_t     = torch.as_tensor(heads,     device=dev)
                rels_t      = torch.as_tensor(rels,      device=dev)
                tails_t     = torch.as_tensor(tails,     device=dev)
                neg_heads_t = torch.as_tensor(neg_heads, device=dev)
                neg_tails_t = torch.as_tensor(neg_tails, device=dev)

                # Forward — embeddings (with gradient tracking via nn.Embedding)
                h  = self._ent(heads_t)
                r  = self._rel(rels_t)
                t  = self._ent(tails_t)
                nh = self._ent(neg_heads_t)
                nt = self._ent(neg_tails_t)

                d_pos = (h + r - t).norm(p=2, dim=1)
                d_neg = (nh + r - nt).norm(p=2, dim=1)
                # F.relu gives zero gradient for non-violated pairs,
                # matching the original ``mask = violation > 0`` logic.
                # NOTE: .sum() — NOT .mean(). The validated NumPy reference
                # applies a per-sample update of scale ~lr to each touched
                # entity. With .mean(), autograd divides every entity's
                # gradient by batch_size, shrinking the effective step
                # ~batch_size× and badly undertraining the model (AUROC
                # collapses toward / below chance on the harder KGs). Do not
                # "simplify" this back to .mean() without rescaling lr by bs.
                loss = F.relu(self.margin + d_pos - d_neg).sum()

                self._opt.zero_grad()
                loss.backward()
                self._opt.step()

                epoch_loss += loss.item()
                n_batches += 1

            # Re-normalise entity embeddings once per epoch
            _normalize_emb_inplace(self._ent)

            if verbose and (epoch + 1) % max(1, n_epochs // 10) == 0:
                avg = epoch_loss / max(n_batches, 1)
                print(f"  Epoch {epoch + 1:>4d}/{n_epochs}  loss={avg:.4f}",
                      flush=True)

            # Periodic checkpoint (same .npz format as original)
            if checkpoint_path and (epoch + 1) % checkpoint_every == 0:
                _save_ckpt_atomic(
                    checkpoint_path,
                    ent_emb=self.ent_emb,
                    rel_emb=self.rel_emb,
                    epoch_completed=np.int32(epoch + 1),
                )


# ── RotatE ────────────────────────────────────────────────────────────────────

class RotatE:
    """RotatE with sparse SGD on the best available device (CUDA / MPS / CPU).

    Complex-space rotational embeddings: each entity is a complex vector of
    dimension ``dim`` (stored as real and imaginary parts separately) and each
    relation is a phase vector that rotates entity embeddings via Hadamard
    product in C^dim.

    The public interface and checkpoint format (.npz with keys ent_re, ent_im,
    rel_phase, epoch_completed) are identical to the original NumPy version.

    Parameters
    ----------
    n_entities, n_relations : int
    dim : int
        Complex embedding dimension. Each entity has 2*dim float parameters.
    margin : float
    lr : float
    seed : int
    """

    name = "RotatE"

    def __init__(self, n_entities, n_relations, dim=64, margin=6.0,
                 lr=0.01, seed=42):
        self.dim = dim
        self.margin = float(margin)
        self.lr = lr
        self.n_entities = n_entities
        self.n_relations = n_relations
        self._device = _get_device()
        self._rng = np.random.RandomState(seed)

        bound = 6.0 / math.sqrt(dim)
        ent_re_w = torch.from_numpy(
            self._rng.uniform(-bound, bound, (n_entities, dim)).astype(np.float32))
        ent_im_w = torch.from_numpy(
            self._rng.uniform(-bound, bound, (n_entities, dim)).astype(np.float32))
        rel_ph_w = torch.from_numpy(
            self._rng.uniform(-math.pi, math.pi, (n_relations, dim)).astype(np.float32))

        self._ent_re  = nn.Embedding(n_entities, dim, sparse=True,
                                     _weight=ent_re_w.clone())
        self._ent_im  = nn.Embedding(n_entities, dim, sparse=True,
                                     _weight=ent_im_w.clone())
        self._rel_ph  = nn.Embedding(n_relations, dim, sparse=True,
                                     _weight=rel_ph_w.clone())
        self._ent_re.to(self._device)
        self._ent_im.to(self._device)
        self._rel_ph.to(self._device)

        self._opt = torch.optim.SGD(
            list(self._ent_re.parameters()) +
            list(self._ent_im.parameters()) +
            list(self._rel_ph.parameters()),
            lr=lr)

    # ── Numpy-compatible properties ──────────────────────────────────────────

    @property
    def ent_re(self) -> np.ndarray:
        return self._ent_re.weight.detach().cpu().numpy()

    @property
    def ent_im(self) -> np.ndarray:
        return self._ent_im.weight.detach().cpu().numpy()

    @property
    def rel_phase(self) -> np.ndarray:
        return self._rel_ph.weight.detach().cpu().numpy()

    # ── Internal distance helper ─────────────────────────────────────────────

    @staticmethod
    def _dist_tensors(h_re, h_im, r_phase, t_re, t_im) -> torch.Tensor:
        """||h ∘ r - t||₂ in complex space (all inputs are tensors)."""
        r_re = torch.cos(r_phase)
        r_im = torch.sin(r_phase)
        hr_re = h_re * r_re - h_im * r_im
        hr_im = h_re * r_im + h_im * r_re
        d_re = hr_re - t_re
        d_im = hr_im - t_im
        return torch.sqrt((d_re ** 2 + d_im ** 2).sum(dim=-1) + 1e-12)

    # ── Scoring ──────────────────────────────────────────────────────────────

    def score(self, h_idx, r_idx, t_idx) -> np.ndarray:
        """Return -||h ∘ r - t||₂  (higher = more plausible)."""
        with torch.no_grad():
            dev = self._device
            hi = torch.as_tensor(np.asarray(h_idx, dtype=np.int64), device=dev)
            ri = torch.as_tensor(np.asarray(r_idx, dtype=np.int64), device=dev)
            ti = torch.as_tensor(np.asarray(t_idx, dtype=np.int64), device=dev)
            dist = self._dist_tensors(
                self._ent_re.weight[hi], self._ent_im.weight[hi],
                self._rel_ph.weight[ri],
                self._ent_re.weight[ti], self._ent_im.weight[ti])
            return -dist.cpu().numpy()

    def score_pairs(self, pairs, rel_idx) -> np.ndarray:
        pairs = np.asarray(pairs)
        h_idx, t_idx = pairs[:, 0], pairs[:, 1]
        if np.isscalar(rel_idx):
            rel_idx = np.full(len(pairs), rel_idx, dtype=np.int64)
        return self.score(h_idx, rel_idx, t_idx)

    # ── Training ─────────────────────────────────────────────────────────────

    def fit(self, triples, n_epochs=50, batch_size=16384, verbose=True,
            checkpoint_path=None, checkpoint_every=10):
        """Train with hinge (margin) ranking loss and sparse SGD.

        Checkpoint format: .npz with keys ent_re, ent_im, rel_phase,
        epoch_completed — identical to the original NumPy version.
        """
        if verbose:
            print(f"  [{self.name}] device={self._device}  "
                  f"entities={self.n_entities:,}  dim={self.dim}  "
                  f"epochs={n_epochs}  batch={batch_size}", flush=True)

        triples = np.asarray(triples, dtype=np.int32)
        n = len(triples)
        dev = self._device

        # ── Resume from checkpoint ──────────────────────────────────────────
        start_epoch = 0
        if checkpoint_path and Path(checkpoint_path).exists():
            try:
                ckpt = np.load(checkpoint_path)
                if (ckpt["ent_re"].shape    == (self.n_entities,  self.dim) and
                        ckpt["ent_im"].shape    == (self.n_entities,  self.dim) and
                        ckpt["rel_phase"].shape == (self.n_relations, self.dim)):
                    self._ent_re.weight.data.copy_(
                        torch.from_numpy(ckpt["ent_re"].astype(np.float32)))
                    self._ent_im.weight.data.copy_(
                        torch.from_numpy(ckpt["ent_im"].astype(np.float32)))
                    self._rel_ph.weight.data.copy_(
                        torch.from_numpy(ckpt["rel_phase"].astype(np.float32)))
                    start_epoch = int(ckpt["epoch_completed"])
                    if verbose:
                        print(f"  Resumed from checkpoint at epoch "
                              f"{start_epoch}/{n_epochs}", flush=True)
                else:
                    if verbose:
                        print("  Checkpoint shape mismatch — starting fresh",
                              flush=True)
            except Exception as e:
                if verbose:
                    print(f"  Checkpoint load failed ({e}) — starting fresh",
                          flush=True)

        # ── Training loop ───────────────────────────────────────────────────
        for epoch in range(start_epoch, n_epochs):
            perm = self._rng.permutation(n)
            epoch_loss = 0.0
            n_batches = 0

            for start in range(0, n, batch_size):
                batch = triples[perm[start:start + batch_size]]
                bs = len(batch)
                heads = batch[:, 0].astype(np.int64)
                rels  = batch[:, 1].astype(np.int64)
                tails = batch[:, 2].astype(np.int64)

                neg_heads = heads.copy()
                neg_tails = tails.copy()
                flip = self._rng.random(bs) < 0.5
                n_h = flip.sum()
                if n_h > 0:
                    neg_heads[flip] = self._rng.randint(0, self.n_entities, n_h)
                if bs - n_h > 0:
                    neg_tails[~flip] = self._rng.randint(
                        0, self.n_entities, bs - n_h)

                heads_t     = torch.as_tensor(heads,     device=dev)
                rels_t      = torch.as_tensor(rels,      device=dev)
                tails_t     = torch.as_tensor(tails,     device=dev)
                neg_heads_t = torch.as_tensor(neg_heads, device=dev)
                neg_tails_t = torch.as_tensor(neg_tails, device=dev)

                h_re  = self._ent_re(heads_t);  h_im  = self._ent_im(heads_t)
                t_re  = self._ent_re(tails_t);  t_im  = self._ent_im(tails_t)
                nh_re = self._ent_re(neg_heads_t); nh_im = self._ent_im(neg_heads_t)
                nt_re = self._ent_re(neg_tails_t); nt_im = self._ent_im(neg_tails_t)
                rp    = self._rel_ph(rels_t)

                d_pos = self._dist_tensors(h_re,  h_im,  rp, t_re,  t_im)
                d_neg = self._dist_tensors(nh_re, nh_im, rp, nt_re, nt_im)
                # .sum() not .mean() — see TransE.fit: .mean() divides each
                # entity's gradient by batch_size and undertrains the model.
                loss  = F.relu(self.margin + d_pos - d_neg).sum()

                self._opt.zero_grad()
                loss.backward()
                self._opt.step()

                epoch_loss += loss.item()
                n_batches += 1

            if verbose and (epoch + 1) % max(1, n_epochs // 10) == 0:
                avg = epoch_loss / max(n_batches, 1)
                print(f"  Epoch {epoch + 1:>4d}/{n_epochs}  loss={avg:.4f}",
                      flush=True)

            if checkpoint_path and (epoch + 1) % checkpoint_every == 0:
                _save_ckpt_atomic(
                    checkpoint_path,
                    ent_re=self.ent_re,
                    ent_im=self.ent_im,
                    rel_phase=self.rel_phase,
                    epoch_completed=np.int32(epoch + 1),
                )


# ── Triple preparation ────────────────────────────────────────────────────────

def kg_to_triples(kg_df, add_inverse=True):
    """Convert a BioKGSuite edge DataFrame to integer-indexed triples.

    Parameters
    ----------
    kg_df : pd.DataFrame
        Must have columns: x_index, y_index, relation.
    add_inverse : bool
        If True, add inverse triples (t, r_inv, h) for each (h, r, t).
        This mirrors the undirected-graph assumption of heuristic scorers
        and makes TransE/RotatE direction-agnostic.

    Returns
    -------
    triples : np.ndarray of shape (N, 3) — (head, relation, tail)
    rel_to_idx : dict
    idx_to_rel : dict
    """
    relations = list(kg_df['relation'].unique())
    if add_inverse:
        inv_relations = [f'{r}_inv' for r in relations]
        all_relations = relations + inv_relations
    else:
        all_relations = relations
    rel_to_idx = {r: i for i, r in enumerate(all_relations)}
    idx_to_rel = {i: r for r, i in rel_to_idx.items()}

    h_vals = kg_df['x_index'].values
    r_vals = kg_df['relation'].map(rel_to_idx).values
    t_vals = kg_df['y_index'].values
    fwd = np.column_stack([h_vals, r_vals, t_vals]).astype(np.int32)

    if add_inverse:
        r_inv_vals = kg_df['relation'].map(
            lambda r: rel_to_idx[f'{r}_inv']).values
        rev = np.column_stack([t_vals, r_inv_vals, h_vals]).astype(np.int32)
        triples = np.vstack([fwd, rev])
    else:
        triples = fwd

    return triples, rel_to_idx, idx_to_rel


def build_train_triples(kg_df, test_pairs_set, indication_rels):
    """Build training triples by removing held-out drug-disease test edges.

    Parameters
    ----------
    kg_df : pd.DataFrame
        Full edge table.
    test_pairs_set : set of (int, int)
        Drug-disease pairs held out for testing.
    indication_rels : list of str
        Relation names corresponding to drug-disease indications.

    Returns
    -------
    train_triples : np.ndarray of shape (M, 3)
    rel_to_idx : dict
    idx_to_rel : dict
    """
    is_indication = kg_df['relation'].isin(set(indication_rels))
    x_vals = kg_df['x_index'].values
    y_vals = kg_df['y_index'].values

    test_both = test_pairs_set | {(b, a) for a, b in test_pairs_set}

    ind_idx = np.where(is_indication)[0]
    drop_mask = np.zeros(len(kg_df), dtype=bool)
    if len(ind_idx) > 0:
        pairs_check = list(zip(x_vals[ind_idx].astype(int),
                               y_vals[ind_idx].astype(int)))
        is_test = np.array([p in test_both for p in pairs_check])
        drop_mask[ind_idx[is_test]] = True

    train_df = kg_df[~drop_mask]
    return kg_to_triples(train_df)


def compute_embedding_metrics(model, test_pairs, neg_pairs, rel_idx,
                              rel_idx_inv=None, ks=(10, 50, 100)):
    """Evaluate embedding model on test pairs with negative sampling.

    Scores each pair using both the forward relation and its inverse
    (if provided), taking the maximum score. This handles KGs where
    indication edges may be stored in either direction.

    Parameters
    ----------
    model : TransE or RotatE instance
    test_pairs : list of (head, tail) tuples — positive pairs
    neg_pairs : list of (head, tail) tuples — negative pairs
    rel_idx : int
        Forward relation index for scoring (drug, rel, disease).
    rel_idx_inv : int, optional
        Inverse relation index. If given, also scores (disease, rel_inv, drug)
        and takes element-wise max.
    ks : tuple of int
        Values of K for Hits@K.

    Returns
    -------
    dict with: auroc, auprc, mrr, hits@K for each K
    """
    from sklearn.metrics import roc_auc_score, average_precision_score

    all_pairs = np.array(test_pairs + neg_pairs)
    labels = np.array([1] * len(test_pairs) + [0] * len(neg_pairs))

    scores_fwd = model.score_pairs(all_pairs, rel_idx)

    if rel_idx_inv is not None:
        scores_inv = model.score_pairs(all_pairs, rel_idx_inv)
        scores = np.maximum(scores_fwd, scores_inv)
    else:
        scores = scores_fwd

    auroc = roc_auc_score(labels, scores)
    auprc = average_precision_score(labels, scores)

    pos_scores = scores[:len(test_pairs)]
    results = {'auroc': float(auroc), 'auprc': float(auprc)}

    reciprocal_ranks = []
    hits_at_k = {k: [] for k in ks}
    for ps in pos_scores:
        rank = int(np.sum(scores >= ps))
        reciprocal_ranks.append(1.0 / rank if rank > 0 else 0.0)
        for k in ks:
            hits_at_k[k].append(1 if rank <= k else 0)

    results['mrr'] = float(np.mean(reciprocal_ranks))
    for k in ks:
        results[f'hits@{k}'] = float(np.mean(hits_at_k[k]))

    results['scores'] = scores.astype(float).tolist()
    results['labels'] = labels.astype(int).tolist()

    return results


# ── Text-embedding baseline: EmbeddingGemma (word-priors only) ───────────────

class GemmaNameEmbedder:
    """Score drug-disease pairs using ONLY a text encoder's prior over entity names.

    This is intentionally NOT a KG embedding method. It exists to answer a
    specific question: "how much of the drug-disease indication signal is
    already latent in a pretrained language model's word priors, with zero
    knowledge of the graph?"

    The model embeds each entity's bare name (no type prefix, no task prefix,
    no neighborhood context) and scores a (drug, disease) pair by cosine
    similarity of their two embeddings. There are no relation parameters.

    Designed to plug into the same eval harness as TransE/RotatE — exposes
    `score_pairs(pairs, rel_idx)` (rel_idx is ignored) and a `name` attr,
    so ``compute_embedding_metrics`` works unchanged.

    Parameters
    ----------
    n_entities : int
        Size of the entity space (so we can allocate the embedding matrix).
    n_relations : int
        Accepted but ignored — keeps the constructor signature symmetric
        with TransE/RotatE.
    dim : int
        Matryoshka-truncated embedding dimension (EmbeddingGemma supports
        128, 256, 512, 768). Lower = less RAM, slightly lower fidelity.
    model_name : str
        HuggingFace model id. Default 'google/embeddinggemma-300m'.
    batch_size : int
        Encoding batch size on CPU/GPU.
    seed : int
        Accepted for interface parity; the encoder itself is deterministic.
    """

    name = 'EmbeddingGemma-300m'

    def __init__(self, n_entities, n_relations=0, dim=768,
                 model_name='google/embeddinggemma-300m',
                 batch_size=64, seed=42, token=None):
        """
        token : str or None
            HuggingFace token for accessing gated models. If None, the
            sentence-transformers library falls back to env vars
            ``HF_TOKEN`` / ``HUGGING_FACE_HUB_TOKEN`` or the cached
            credential at ``~/.cache/huggingface/token``.
        """
        self.n_entities = n_entities
        self.n_relations = n_relations  # ignored
        self.dim = dim
        self.model_name = model_name
        self.batch_size = batch_size
        self.seed = seed
        self.token = token
        self._encoder = None
        # Lazy-allocate; encode_entities fills this in
        self.ent_emb = None  # (n_entities, dim) float32, L2-normalized

    def _load_encoder(self):
        if self._encoder is not None:
            return
        from sentence_transformers import SentenceTransformer
        import torch
        device = ('cuda' if torch.cuda.is_available()
                  else 'mps' if getattr(torch.backends, 'mps', None)
                                  and torch.backends.mps.is_available()
                  else 'cpu')
        kwargs = {}
        if self.dim < 768:
            kwargs['truncate_dim'] = self.dim
        token = (self.token
                 or os.environ.get('HF_TOKEN')
                 or os.environ.get('HUGGING_FACE_HUB_TOKEN'))
        if token is not None:
            kwargs['token'] = token
        self._encoder = SentenceTransformer(self.model_name, device=device, **kwargs)
        try:
            torch.manual_seed(self.seed)
        except Exception:
            pass

    def encode_entities(self, names, verbose=True):
        """Encode an ordered list of names into self.ent_emb.

        Parameters
        ----------
        names : list of str
            Length must equal n_entities. Order corresponds to entity idx.
            Empty / NaN names are replaced with the literal "" to keep
            indices aligned (those entities will get the model's default
            "empty string" embedding — distinguishable from missing).

        Returns
        -------
        self.ent_emb : np.ndarray of shape (n_entities, dim), L2-normalized.
        """
        assert len(names) == self.n_entities, (
            f"got {len(names)} names but n_entities={self.n_entities}")
        self._load_encoder()
        clean = [('' if (n is None or (isinstance(n, float) and np.isnan(n)))
                  else str(n)) for n in names]
        if verbose:
            print(f'  Encoding {len(clean)} names with {self.model_name} '
                  f'(dim={self.dim}, batch={self.batch_size})...', flush=True)
        emb = self._encoder.encode(
            clean,
            batch_size=self.batch_size,
            show_progress_bar=verbose,
            convert_to_numpy=True,
            normalize_embeddings=True,
        ).astype(np.float32)
        self.ent_emb = emb
        return self.ent_emb

    def score(self, h_idx, r_idx, t_idx):
        """Higher = more plausible. Ignores r_idx (no relation concept)."""
        if self.ent_emb is None:
            raise RuntimeError("Call encode_entities(names) before scoring.")
        h = self.ent_emb[h_idx]
        t = self.ent_emb[t_idx]
        return (h * t).sum(axis=-1)

    def score_pairs(self, pairs, rel_idx):
        """Score a list of (head, tail) pairs. rel_idx is accepted but ignored."""
        pairs = np.asarray(pairs)
        return self.score(pairs[:, 0], None, pairs[:, 1])

    def fit(self, *args, **kwargs):
        """No-op: there's nothing to train. encode_entities() does the work."""
        return self

    # ── Cache helpers ──────────────────────────────────────────────────────
    def save_embeddings(self, path):
        """Save self.ent_emb to disk (atomic)."""
        if self.ent_emb is None:
            raise RuntimeError("Nothing to save — encode_entities first.")
        _save_ckpt_atomic(path, ent_emb=self.ent_emb,
                          dim=np.int32(self.dim),
                          n_entities=np.int32(self.n_entities))

    def load_embeddings(self, path):
        """Load self.ent_emb from a prior save_embeddings() call."""
        ckpt = np.load(path)
        if int(ckpt['n_entities']) != self.n_entities:
            raise ValueError(
                f"cache n_entities={int(ckpt['n_entities'])} != "
                f"current n_entities={self.n_entities}")
        if int(ckpt['dim']) != self.dim:
            raise ValueError(
                f"cache dim={int(ckpt['dim'])} != current dim={self.dim}")
        self.ent_emb = ckpt['ent_emb'].astype(np.float32, copy=True)
        return self.ent_emb


# ── Model registry ────────────────────────────────────────────────────────────

MODELS = {
    'TransE': TransE,
    'RotatE': RotatE,
    'Gemma':  GemmaNameEmbedder,
}
