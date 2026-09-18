"""Unit tests for biokgsuite.grounding.clients.

    python -m pytest tests/test_grounding_clients.py -q
"""
import pytest

from biokgsuite.grounding import clients
from biokgsuite.grounding.prompt import parse_response

LETTERS = ['A', 'B', 'C', 'D', 'E']


def _candidates(evidence_for=('B', 'D')):
    return [{'letter': letter, 'drug_name': f'drug{letter}',
             'evidence': 'some evidence' if letter in evidence_for else None}
            for letter in LETTERS]


def test_mock_rank_is_parseable_and_ranks_every_letter():
    text = clients.mock_rank(_candidates(), seed=3)
    ordered, fields = parse_response(text, LETTERS)
    assert sorted(ordered) == LETTERS
    assert set(fields) == set(LETTERS)
    assert all(f['basis'] == 'both' and f['confidence'] == 3 for f in fields.values())


def test_mock_rank_is_seeded():
    assert clients.mock_rank(_candidates(), seed=7) == clients.mock_rank(_candidates(), seed=7)


def test_rank_call_mock_never_touches_the_network(monkeypatch):
    def boom(*args, **kwargs):
        raise AssertionError('network route called in mock mode')
    monkeypatch.setattr(clients, 'api_rank', boom)
    monkeypatch.setattr(clients, 'ollama_rank', boom)
    monkeypatch.setattr(clients, '_post_with_backoff', boom)
    monkeypatch.setattr(clients.time, 'sleep', boom)
    for model in ('gpt-4.1-mini', 'llama3.3:70b'):
        text = clients.rank_call(model, 'prompt', _candidates(), 0,
                                 mock=True, temperature=0.7, api_delay=1.0)
        ordered, _ = parse_response(text, LETTERS)
        assert sorted(ordered) == LETTERS


def test_rank_call_routes_api_and_ollama(monkeypatch):
    calls = []
    monkeypatch.setattr(clients, 'api_rank',
                        lambda name, prompt, temperature: calls.append(('api', name, temperature)) or 'x')
    monkeypatch.setattr(clients, 'ollama_rank',
                        lambda model, prompt, seed, temperature: calls.append(('ollama', model, seed)) or 'y')
    monkeypatch.setattr(clients.time, 'sleep', lambda s: calls.append(('sleep', s)))
    assert clients.rank_call('gpt-4.1-mini', 'p', _candidates(), 5,
                             mock=False, temperature=0.2, api_delay=1.5) == 'x'
    assert clients.rank_call('llama3.3:70b', 'p', _candidates(), 5,
                             mock=False, temperature=0.2, api_delay=1.5) == 'y'
    assert calls == [('sleep', 1.5), ('api', 'gpt-4.1-mini', 0.2), ('ollama', 'llama3.3:70b', 5)]


@pytest.mark.parametrize('name, expected', [
    ('gpt-4.1-mini', ('openai', 'gpt-4.1-mini')),
    ('groq:llama-3.3-70b-versatile', ('groq', 'llama-3.3-70b-versatile')),
    ('gemini:gemini-3.5-flash-lite', ('gemini', 'gemini-3.5-flash-lite')),
    ('llama3.3:70b', (None, None)),
    ('claude-sonnet-4-5', ('anthropic', 'claude-sonnet-4-5')),
])
def test_provider_model(name, expected):
    assert clients.provider_model(name) == expected
    assert clients.is_api_model(name) == (expected[0] is not None)


def test_api_key_for_reads_provider_variable(monkeypatch):
    monkeypatch.setenv('GROQ_API_KEY', 'secret')
    monkeypatch.delenv('OPENAI_API_KEY', raising=False)
    assert clients.api_key_for('groq:llama-3.3-70b-versatile') == 'secret'
    assert clients.api_key_for('gpt-4.1-mini') is None
    assert clients.api_key_for('llama3.3:70b') is None


def test_api_rank_without_key_reports_error_not_exception(monkeypatch):
    monkeypatch.delenv('OPENAI_API_KEY', raising=False)
    assert clients.api_rank('gpt-4.1-mini', 'p', 0.7).startswith('__ERROR__ no key')
    assert clients.api_rank('llama3.3:70b', 'p', 0.7) == '__ERROR__ not an API model'


def test_load_dotenv_does_not_override_existing(tmp_path, monkeypatch):
    (tmp_path / '.env').write_text('# comment\nNEW_KEY="abc"\nOLD_KEY=from_file\n')
    monkeypatch.delenv('NEW_KEY', raising=False)
    monkeypatch.setenv('OLD_KEY', 'from_shell')
    clients.load_dotenv(tmp_path)
    assert clients.os.environ['NEW_KEY'] == 'abc'
    assert clients.os.environ['OLD_KEY'] == 'from_shell'
