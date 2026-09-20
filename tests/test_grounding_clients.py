"""Rate-limit handling in the hosted-API client."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))

from biokgsuite.grounding import clients as C  # noqa: E402


def test_daily_quota_is_recognised_from_the_429_body():
    gemini = ('{"error": {"code": 429, "status": "RESOURCE_EXHAUSTED", "details": [{"violations": '
              '[{"quotaMetric": "generativelanguage.googleapis.com/generate_content_free_tier_requests", '
              '"quotaId": "GenerateRequestsPerDayPerProjectPerModel-FreeTier"}]}]}}')
    assert C.is_daily_quota(gemini)
    assert C.is_daily_quota('Rate limit reached for tokens per day (TPD)')
    assert not C.is_daily_quota('{"error": {"message": "Rate limit reached for requests per minute"}}')
    assert not C.is_daily_quota('')


def test_backoff_returns_at_once_on_a_daily_quota(monkeypatch):
    calls = []

    class R:
        status_code = 429
        headers = {}
        text = '{"quotaId": "GenerateRequestsPerDayPerProjectPerModel-FreeTier"}'

    class Req:
        @staticmethod
        def post(url, headers=None, json=None, timeout=None):
            calls.append(1)
            return R()

    monkeypatch.setitem(sys.modules, 'requests', Req)
    monkeypatch.setattr(C.time, 'sleep', lambda s: None)
    r, err = C._post_with_backoff('http://x', {}, {})
    assert r is None and err.startswith('__ERROR__ 429 daily quota exhausted')
    assert len(calls) == 1


def test_backoff_keeps_the_429_body_after_retries(monkeypatch):
    class R:
        status_code = 429
        headers = {'retry-after': '0'}
        text = 'Rate limit reached for requests per minute'

    class Req:
        @staticmethod
        def post(url, headers=None, json=None, timeout=None):
            return R()

    monkeypatch.setitem(sys.modules, 'requests', Req)
    monkeypatch.setattr(C.time, 'sleep', lambda s: None)
    r, err = C._post_with_backoff('http://x', {}, {}, max_retries=3)
    assert r is None and '429' in err and 'per minute' in err
