"""Unit tests for biokgsuite.grounding.prompt.

    python -m pytest tests/test_grounding_prompt.py -q
"""
import json

from biokgsuite.grounding.prompt import (PROMPT_HEADER, PROMPT_INSTRUCTIONS, build_prompt,
                                         parse_response, rank_of)

LETTERS = ['A', 'B', 'C', 'D']


def _candidates(evidence_for=()):
    return [{'letter': letter, 'drug_name': f'drug{letter}',
             'evidence': f'drug{letter} targets GENE{letter}' if letter in evidence_for else None}
            for letter in LETTERS]


def test_build_prompt_lists_every_letter_and_only_given_evidence():
    text = build_prompt('asthma', 'MONDO:0004979', _candidates(evidence_for=('B',)),
                        disease_profile='asthma is associated with gene IL13')
    for letter in LETTERS:
        assert f'\n{letter}. drug{letter}' in text
    assert text.count('Drug profile (from the knowledge graph)') == 1
    assert 'drugB targets GENEB' in text
    assert 'Disease profile (from the knowledge graph): asthma is associated with gene IL13' in text
    assert text.startswith(PROMPT_HEADER)
    assert text.endswith(PROMPT_INSTRUCTIONS)


def test_build_prompt_omits_profile_lines_when_evidence_is_none():
    text = build_prompt('asthma', 'MONDO:0004979', _candidates())
    assert 'Drug profile' not in text
    assert 'Disease profile' not in text


def _ranking_json(order, **overrides):
    items = []
    for i, letter in enumerate(order):
        item = {'rank': i + 1, 'drug': letter, 'basis': 'prior',
                'evidence_agreement': 'insufficient', 'confidence': 3, 'rationale': 'r'}
        item.update(overrides.get(letter, {}))
        items.append(item)
    return json.dumps({'reasoning': 'because', 'ranking': items})


def test_parse_response_well_formed_json():
    text = 'Some preamble.\n' + _ranking_json(['C', 'A', 'D', 'B'], A={'basis': 'KG', 'confidence': 9})
    ordered, fields = parse_response(text, LETTERS)
    assert ordered == ['C', 'A', 'D', 'B']
    assert set(fields) == set(LETTERS)
    assert fields['A']['basis'] == 'kg'            # normalised to lower case
    assert fields['A']['confidence'] == 5          # clipped to the 1-5 scale
    assert fields['B']['evidence_agreement'] == 'insufficient'
    assert fields['B']['rationale'] == 'r'


def test_parse_response_orders_by_rank_field_and_tolerates_bad_ranks():
    # Ranks are given out of list order, one is a string and one is missing:
    # numeric ranks decide the order, and the others fall back to list position.
    text = json.dumps({'reasoning': '', 'ranking': [
        {'rank': 3, 'drug': 'A', 'basis': 'prior', 'evidence_agreement': 'insufficient',
         'confidence': 'high', 'rationale': ''},
        {'rank': 'first', 'drug': 'B'},
        {'drug': 'C'},
        {'rank': 1, 'drug': 'd'},
        {'rank': 2, 'drug': 'Z'},
    ]})
    ordered, fields = parse_response(text, LETTERS)
    assert ordered == ['D', 'B', 'A', 'C']
    assert fields['A']['confidence'] is None
    assert fields['A']['rationale'] is None
    assert fields['B'] == {'basis': None, 'evidence_agreement': None,
                           'confidence': None, 'rationale': None}
    assert 'Z' not in fields


def test_parse_response_free_text_fallback():
    text = 'I would rank B first, then D, then A, and finally C. B is best.'
    ordered, fields = parse_response(text, LETTERS)
    assert ordered == ['B', 'D', 'A', 'C']
    assert fields == {}


def test_parse_response_empty_and_malformed():
    assert parse_response('', LETTERS) == ([], {})
    assert parse_response(None, LETTERS) == ([], {})
    ordered, fields = parse_response('{"ranking": [not json', LETTERS)
    assert ordered == []
    assert fields == {}


def test_rank_of_position_and_miss():
    assert rank_of('B', ['C', 'B', 'A'], pool_size=8) == 2
    assert rank_of('A', ['C', 'B', 'A'], pool_size=8) == 3
    assert rank_of('D', ['C', 'B', 'A'], pool_size=8) == 9
    assert rank_of('A', [], pool_size=8) == 9
