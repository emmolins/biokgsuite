"""Model dispatch for the grounding experiment: mock, hosted APIs and Ollama.

A model name selects the route. ``provider:model`` names a hosted provider from
API_PROVIDERS explicitly; bare ``claude*`` and ``gpt-*`` / ``o1`` / ``o3`` /
``o4`` names go to Anthropic and OpenAI; anything else is a local Ollama tag.
Every route returns the raw response text, or a string starting with
``__ERROR__`` so that a failed call is recorded in the run rows rather than
aborting a sweep that may have taken hours.
"""
import os
import random
import time

import numpy as np

# The Ollama server this process talks to. OLLAMA_HOST is the variable the
# ollama binary itself reads, so one setting names the server for both sides;
# a Slurm job sets it to a port of its own so two jobs on one node never share
# (and stall) a single server.
def _ollama_base():
    host = os.environ.get('OLLAMA_HOST', '127.0.0.1:11434').strip()
    if not host.startswith('http'):
        host = 'http://' + host
    return host.rstrip('/')


OLLAMA_URL = _ollama_base() + '/api/generate'
TAGS_URL = _ollama_base() + '/api/tags'

API_PROVIDERS = {
    'anthropic':  dict(url='https://api.anthropic.com/v1/messages',
                       key='ANTHROPIC_API_KEY', style='anthropic'),
    'openai':     dict(url='https://api.openai.com/v1/chat/completions',
                       key='OPENAI_API_KEY', style='openai'),
    'groq':       dict(url='https://api.groq.com/openai/v1/chat/completions',
                       key='GROQ_API_KEY', style='openai'),
    'openrouter': dict(url='https://openrouter.ai/api/v1/chat/completions',
                       key='OPENROUTER_API_KEY', style='openai'),
    'gemini':     dict(url='https://generativelanguage.googleapis.com/v1beta/openai/chat/completions',
                       key='GEMINI_API_KEY', style='openai'),
    'deepseek':   dict(url='https://api.deepseek.com/v1/chat/completions',
                       key='DEEPSEEK_API_KEY', style='openai'),
}

# Reasoning models spend tokens thinking before the JSON; a 1500-token budget
# truncates them mid-answer, which parses as a miss and would be scored as one.
_REASONING_HINTS = ('gpt-oss', 'o1', 'o3', 'o4', '-r1', 'deepseek-r', 'qwq', 'reason', 'think')

# Whether a provider accepted response_format=json_object is learnt on the
# first call and remembered, so a provider that rejects it costs one failed
# request per process rather than one per query.
_API_FORMAT_CACHE = {}


def provider_model(name):
    """(provider, model) for a hosted-API name, (None, None) for an Ollama tag."""
    prefix = name.split(':', 1)[0]
    if prefix in API_PROVIDERS and ':' in name:
        return prefix, name.split(':', 1)[1]
    low = name.lower()
    if low.startswith('claude'):
        return 'anthropic', name
    if low.startswith(('gpt-', 'o1', 'o3', 'o4')):
        return 'openai', name
    return None, None


def is_api_model(name):
    return provider_model(name)[0] is not None


def _max_tokens(name):
    n = name.lower()
    if any(h in n for h in _REASONING_HINTS):
        return 8000
    # Gemini 2.5 and 3.x Flash think by default, so the budget covers the
    # thinking as well as the JSON answer.
    if 'gemini' in n:
        return 8000
    return 1500


def api_key_for(name):
    prov, _ = provider_model(name)
    return os.environ.get(API_PROVIDERS[prov]['key']) if prov else None


def load_dotenv(base):
    """Read ``<base>/.env`` into os.environ without overriding variables that
    are already set, so a key exported in the shell wins over the file."""
    envf = base / '.env'
    if not envf.exists():
        return
    for line in envf.read_text().splitlines():
        line = line.strip()
        if line and not line.startswith('#') and '=' in line:
            k, v = line.split('=', 1)
            os.environ.setdefault(k.strip(), v.strip().strip('"').strip("'"))


def ollama_models_available():
    """Names of the models the local Ollama server has pulled."""
    import requests
    try:
        data = requests.get(TAGS_URL, timeout=3).json()
        return {m['name'] for m in data.get('models', [])}
    except Exception as e:
        raise SystemExit(f"Ollama not reachable at {TAGS_URL}: {e}")


def ollama_rank(model, prompt, seed, temperature):
    import requests
    # Plain JSON mode rather than the strict schema: the schema made some
    # local models loop or truncate, and the parser tolerates loose JSON.
    # num_predict is raised to at least 4000 because 70B models are verbose
    # enough to hit the default budget mid-ranking.
    payload = {'model': model, 'prompt': prompt, 'stream': False,
               'format': 'json',
               'options': {'temperature': temperature,
                           'num_predict': max(_max_tokens(model), 4000),
                           'num_ctx': 8192, 'seed': int(seed)}}
    try:
        r = requests.post(OLLAMA_URL, json=payload, timeout=600)
        r.raise_for_status()
        return r.json().get('response', '')
    except Exception as e:
        return f'__ERROR__ {e}'


def _post_with_backoff(url, headers, body, timeout=180, max_retries=8):
    """POST with exponential backoff on transient errors and 429s.

    Returns (response, None) or (None, error string). Retry-After is honoured
    up to 60 s with jitter so concurrent runs do not retry in lockstep.
    """
    import requests
    delay = 3.0
    detail = ''
    for _ in range(max_retries):
        try:
            r = requests.post(url, headers=headers, json=body, timeout=timeout)
        except Exception:
            time.sleep(min(delay, 30))
            delay = min(delay * 2, 60)
            continue
        if r.status_code == 429:
            detail = ' '.join(r.text.split())[:200]
            if is_daily_quota(detail):
                # A per-day quota does not reopen within any retry ladder;
                # report it at once so the arm aborts instead of idling.
                return None, f'__ERROR__ 429 daily quota exhausted: {detail}'
            ra = r.headers.get('retry-after', '')
            try:
                wait = float(ra)
            except ValueError:
                wait = delay
            time.sleep(min(wait, 60) + random.uniform(0, 1.5))
            delay = min(delay * 2, 60)
            continue
        return r, None
    return None, f'__ERROR__ 429 rate-limited (gave up after retries): {detail}'


def is_daily_quota(text):
    """True when a 429 body names a per-day quota (Gemini's
    GenerateRequestsPerDay..., Groq's tokens per day), which a retry ladder
    cannot wait out."""
    t = text.lower()
    return 'perday' in t or 'per day' in t or 'daily' in t


def api_rank(name, prompt, temperature):
    """Hosted-API completion for ``name``; ``__ERROR__ ...`` on any failure."""
    prov, model = provider_model(name)
    cfg = API_PROVIDERS.get(prov)
    if not cfg:
        return '__ERROR__ not an API model'
    key = os.environ.get(cfg['key'])
    if not key:
        return f'__ERROR__ no key ({cfg["key"]})'
    mx = _max_tokens(name)
    if cfg['style'] == 'anthropic':
        r, err = _post_with_backoff(
            cfg['url'],
            {'x-api-key': key, 'anthropic-version': '2023-06-01', 'content-type': 'application/json'},
            {'model': model, 'max_tokens': mx, 'temperature': temperature,
             'messages': [{'role': 'user', 'content': prompt}]})
        if err:
            return err
        if r.status_code >= 400:
            return f'__ERROR__ {r.status_code} {r.text[:160]}'
        try:
            return ''.join(b.get('text', '') for b in r.json().get('content', []))
        except Exception as e:
            return f'__ERROR__ parse {e}'
    headers = {'Authorization': f'Bearer {key}', 'content-type': 'application/json'}
    base = {'model': model, 'max_tokens': mx, 'temperature': temperature,
            'messages': [{'role': 'user', 'content': prompt}]}
    ck = (prov, model)
    formats = ([_API_FORMAT_CACHE[ck]] if ck in _API_FORMAT_CACHE
               else [{'type': 'json_object'}, None])
    last = '__ERROR__ no response'
    for fmt in formats:
        body = {**base, **({'response_format': fmt} if fmt else {})}
        r, err = _post_with_backoff(cfg['url'], headers, body)
        if err:
            return err
        if r.status_code >= 400:
            last = f'__ERROR__ {r.status_code} {r.text[:160]}'
            continue
        _API_FORMAT_CACHE[ck] = fmt
        try:
            return r.json()['choices'][0]['message']['content']
        except Exception as e:
            return f'__ERROR__ parse {e}'
    return last


def mock_rank(candidates, seed):
    """Simulate a discerning model: candidates with evidence float up, with
    noise, so the mock pipeline shows a KG lift without being deterministic."""
    rng = np.random.default_rng(seed)
    scored = [(1.0 * bool(c['evidence']) + rng.normal(0, 0.5), c['letter']) for c in candidates]
    scored.sort(key=lambda t: -t[0])
    ordered = [letter for _, letter in scored]
    items = ", ".join(
        f'{{"rank": {i + 1}, "drug": "{letter}", "basis": "both", '
        f'"evidence_agreement": "consistent", "confidence": 3, "rationale": "mock"}}'
        for i, letter in enumerate(ordered))
    return f'Reasoning: mock.\n{{"ranking": [{items}]}}'


def rank_call(model, prompt, candidates, seed, *, mock, temperature, api_delay):
    """Single dispatch: mock, hosted API, or local Ollama tag."""
    if mock:
        return mock_rank(candidates, seed)
    if is_api_model(model):
        if api_delay:
            time.sleep(api_delay)
        return api_rank(model, prompt, temperature)
    return ollama_rank(model, prompt, seed, temperature)
