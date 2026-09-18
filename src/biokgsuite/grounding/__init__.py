"""KG-grounded drug repurposing ranking: does knowledge-graph evidence help an
LLM place the true post-cutoff repurposed drug above distractors?

    benchmark.py  positives from the post-cutoff approval set, per-graph flags,
                  the Open Targets distractor pool, query construction
    prompt.py     prompt text, output schema, response parsing
    clients.py    mock / hosted API / Ollama dispatch
    dossiers.py   KG evidence blocks under the entity-type policy
    runner.py     the ranking loop and the tidy run-row files
    analysis.py   per-query frame, paired bootstrap, contrasts, FDR, capacity
    bias.py       kept-versus-dropped selection audit

The stage notebooks under grounding_experiment/ are thin callers of these.
"""
from .analysis import (arm_table, bh_fdr, capacity_table, chance_mrr, contrast_table,  # noqa: F401
                       cross_model_corr, diagnostics, difficulty_table, load_runs,
                       mean_pairwise_corr, per_query, scaling_table, significance_tests,
                       validate_runs)
from .benchmark import (build_queries, load_benchmark, load_distractor_pool,  # noqa: F401
                        mondo_children, selection_frame)
from .bias import run_bias_tests  # noqa: F401
from .clients import api_key_for, is_api_model, load_dotenv, ollama_models_available, rank_call  # noqa: F401
from .dossiers import KGDossiers, load_entrez_to_symbol  # noqa: F401
from .prompt import build_prompt, parse_response, rank_of  # noqa: F401
from .runner import run_family, run_path  # noqa: F401
