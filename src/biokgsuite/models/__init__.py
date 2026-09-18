"""Knowledge-graph embedding models.

Kept in the package rather than beside the supplementary notebook so the
dimension notebooks, the supplementary notebooks and the HPC drivers all import
the same implementation from one place. `torch` is imported lazily by
`biokgsuite.models.embedding`, so importing `biokgsuite` costs nothing when the
embedding arms are not being run.
"""
