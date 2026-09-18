"""Evaluation sets and the resolution they need.

These three modules build and interpret what notebook 07 scores against: the
post-cutoff approval set and its validator, the regulator approval records it
is derived from, and the ontology classifier that assigns the therapeutic
domains defining the cross-domain strata.

They lived in `biokgsuite_dimensions/07_generalization/` and were importable
only because Jupyter puts the notebook's own directory on `sys.path`. That made
them invisible to `pip install -e .`, unavailable to any other notebook or
tool, and reachable from the tests only through a hand-written path insert. They
are ordinary library code and belong with the rest of it.
"""
from .ontology_classifier import ONTOLOGY_DOMAINS, OntologyClassifier  # noqa: F401
from .post_cutoff_approvals import (  # noqa: F401
    curie_variants, disease_id_lookup, load_post_cutoff_approvals, norm_curie,
    to_new_indication_approvals,
)
from .regulatory_approvals import get_non_fda_pairs, load_expanded_approvals  # noqa: F401
