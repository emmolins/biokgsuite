# Supplementary: embedding cross-checks

| Path | Contents |
|---|---|
| `embeddings/` | TransE and RotatE link prediction on the same held-out pairs as notebook 06, and an EmbeddingGemma name prior; Supplementary Table S7. `EMBEDDING_CROSSCHECK.md` describes the design |
| `embeddings/run_stability_embeddings.py` | whether notebook 05's stability result survives a learned readout; not run at full scale for the manuscript |

The embedding code is `src/biokgsuite/models/embedding.py` and its figures
`src/biokgsuite/viz/embeddings.py`, so the supplementary uses the same
implementation as the rest of the repository. `torch` is imported only by
`biokgsuite.models.embedding`; the seven dimension notebooks do not need it.

## Stability under a learned readout

Notebook 05 measures stability with neighbourhood heuristics, which are a
deterministic function of the graph: the only thing varying across its sweep is
the graph, which is what makes its standard deviations mean anything. Embedding
training is stochastic, so a drop in rank agreement after retraining mixes real
signal loss with run-to-run training variance.

`run_stability_embeddings.py` measures the variance FLOOR first, by training
`n_trials` seeds on the unperturbed graph and correlating their rankings with
each other. A perturbed result only counts as signal loss if it falls below that
floor. It uses the same therapeutic pairs (`stats.stability.therapeutic_pairs`,
shared with notebook 05 so the two cannot drift), the same triple-level dropout,
and the same Spearman statistic. Only the scorer changes.

```
python supplementary/embeddings/run_stability_embeddings.py --kg hetionet --smoke
python supplementary/embeddings/run_stability_embeddings.py --kg hetionet
```

The smoke config (dim 8, 1 epoch) exists to prove the path runs and is not a
result. It is also instructive: at that size TransE's retraining floor is about
r = -0.02, meaning two runs on the identical graph agree no better than chance,
and every perturbed number sits above the floor and means nothing. A real run
needs enough epochs for the floor to be high before any perturbation reading is
interpretable.
