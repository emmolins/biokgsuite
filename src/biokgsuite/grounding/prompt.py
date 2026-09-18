"""Prompt text, output schema and response parsing for the grounding experiment.

The prompt asks the model to weigh knowledge-graph evidence against its own
pharmacological knowledge, flag conflicts, abstain when the evidence is
insufficient, and then emit one JSON object. Every reliance field is required
so that the analysis can measure how much the model leaned on the graph, not
only where it placed the positive.

The wording is kept verbatim from the run that produced the reported results:
changing a word here changes the experiment.
"""
import json
import re

PROMPT_HEADER = (
    "You are assisting a drug repurposing scientist. For the target disease below,\n"
    "you are given candidate drugs. Some include evidence retrieved from a\n"
    "biomedical knowledge graph, reported as mechanistic links\n"
    "(drug -> target -> pathway -> disease) and related annotations.\n"
)
PROMPT_INSTRUCTIONS = (
    "Assess each candidate's mechanistic plausibility for this disease using BOTH\n"
    "the knowledge-graph evidence and your own pharmacological knowledge. The\n"
    "knowledge graph reports associations; treat them as one source, not as\n"
    "established truth. Where its evidence is sparse, or conflicts with what you\n"
    "know, say so and weight it accordingly -- do not construct a rationale to fit\n"
    "weak evidence. If the evidence is insufficient to judge a candidate, mark it\n"
    "\"uncertain\" rather than forcing a confident rank.\n\n"
    "Respond with ONLY a single JSON object (no text before or after) of the form:\n"
    "{ \"reasoning\": \"<a few sentences weighing the candidates and flagging any "
    "KG-vs-knowledge conflicts>\",\n"
    "  \"ranking\": [ { \"rank\": 1, \"drug\": \"<letter>\",\n"
    "    \"basis\": \"kg\" | \"prior\" | \"both\",\n"
    "    \"evidence_agreement\": \"consistent\" | \"conflicting\" | \"insufficient\",\n"
    "    \"confidence\": 1-5, \"rationale\": \"<one sentence>\" }, ... ] }\n"
    "Rank every candidate. ALL fields are required for EVERY candidate, even when "
    "no evidence is shown: set basis=\"prior\" when you used only your own "
    "knowledge (\"kg\"/\"both\" when you used the provided evidence), "
    "evidence_agreement=\"insufficient\" when no evidence was given; confidence "
    "(1-5) and a one-sentence rationale are ALWAYS required."
)

# Schema for providers that enforce structured output (Ollama `format=`, OpenAI
# json_schema). The parser below does not rely on it: a model that ignores the
# schema is still scored from whatever it returned.
RANKING_SCHEMA = {
    "type": "object",
    "properties": {
        "reasoning": {"type": "string"},
        "ranking": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "rank": {"type": "integer"},
                    "drug": {"type": "string"},
                    "basis": {"type": "string", "enum": ["kg", "prior", "both"]},
                    "evidence_agreement": {"type": "string",
                                           "enum": ["consistent", "conflicting", "insufficient"]},
                    "confidence": {"type": "integer", "minimum": 1, "maximum": 5},
                    "rationale": {"type": "string"},
                },
                "required": ["rank", "drug", "basis", "evidence_agreement",
                             "confidence", "rationale"],
            },
        },
    },
    "required": ["reasoning", "ranking"],
}


def build_prompt(disease_name, disease_id, candidates, disease_profile=None):
    """Render the prompt for one query.

    ``candidates`` is the candidate view: dicts with ``letter``, ``drug_name``
    and ``evidence`` (a string, or None). ``disease_profile`` is the KG disease
    dossier, shown once in the header on KG arms only.
    """
    lines = [PROMPT_HEADER, f"\nTarget disease: {disease_name} ({disease_id})"]
    if disease_profile:
        lines.append(f"Disease profile (from the knowledge graph): {disease_profile}")
    lines.append("\nCandidates (listed order is arbitrary):")
    for c in candidates:
        lines.append(f"{c['letter']}. {c['drug_name']}")
        # A candidate without evidence gets no profile line at all, so the no-KG
        # arm and an unresolved candidate on a KG arm look identical to the model.
        if c['evidence']:
            lines.append(f"   Drug profile (from the knowledge graph): {c['evidence']}")
    lines.append("")
    lines.append(PROMPT_INSTRUCTIONS)
    return "\n".join(lines)


def _candidate_letter(item, valid_letters):
    """The candidate letter an entry names, or None when it names none."""
    d = str(item.get('drug', '')).strip().upper()
    return d[0] if d and d[0] in valid_letters else None


def _parse_json_ranking(text, valid_letters):
    """Ordered letters and per-letter fields from the first JSON object in
    ``text``. Returns ([], {}) when there is none or it carries no ranking."""
    m = re.search(r'\{.*\}', text, re.DOTALL)
    if not m:
        return [], {}
    fields, ranked = {}, []
    try:
        obj = json.loads(m.group(0))
        for it in obj.get('ranking', []):
            d = _candidate_letter(it, valid_letters)
            if not d:
                continue
            if d not in fields:
                fields[d] = {
                    'basis': (str(it.get('basis', '')).lower() or None),
                    'evidence_agreement': (str(it.get('evidence_agreement', '')).lower() or None),
                    'confidence': _as_int(it.get('confidence')),
                    'rationale': (str(it.get('rationale', '')).strip() or None),
                }
            # A letter listed twice keeps its first position; a missing or
            # non-numeric rank falls back to list order.
            if d not in [r[0] for r in ranked]:
                rk = it.get('rank')
                ranked.append((d, rk if isinstance(rk, (int, float)) else len(ranked) + 1))
        ranked.sort(key=lambda t: t[1])
    except (json.JSONDecodeError, TypeError, AttributeError):
        # Fields read before a malformed entry still describe the model's
        # reliance on the evidence, so they survive even when the order is lost.
        return [], fields
    return [d for d, _ in ranked], fields


def parse_ranking(text, valid_letters):
    """Ordered list of candidate letters, best first. [] if unparseable.

    Tries the JSON object first; when that yields nothing, falls back to the
    order in which capital letters first appear in the text, so a model that
    answered in prose ("B, then A, then C") is still scored.
    """
    if not text:
        return []
    ordered, _ = _parse_json_ranking(text, valid_letters)
    if ordered:
        return ordered
    seen = []
    for tok in re.findall(r'\b([A-Z])\b', text.upper()):
        if tok in valid_letters and tok not in seen:
            seen.append(tok)
    return seen


def _as_int(x):
    """Confidence clipped to the 1-5 scale, or None when it is not a number."""
    try:
        return max(1, min(5, int(x)))
    except (ValueError, TypeError):
        return None


def parse_response(text, valid_letters):
    """Like parse_ranking, plus per-letter reliance fields
    {letter: {basis, evidence_agreement, confidence, rationale}}.

    The fields are kept even when the JSON ranking is empty and the order
    comes from the free-text fallback, because they say how the model used
    the evidence regardless of where it placed each candidate.
    """
    ordered, fields = ([], {}) if not text else _parse_json_ranking(text, valid_letters)
    if not ordered:
        ordered = parse_ranking(text, valid_letters)
    return ordered, fields


def rank_of(positive_letter, ordered_letters, pool_size):
    """1-based rank of the positive; pool_size + 1 (a miss) if absent or unparsed."""
    if positive_letter in ordered_letters:
        return ordered_letters.index(positive_letter) + 1
    return pool_size + 1
