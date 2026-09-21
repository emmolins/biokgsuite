"""The shape check behind annotation accuracy, and the coverage arithmetic,
on the identifier shapes the six graphs actually write.

    python -m pytest tests/test_annotation_coverage.py -q
"""
import pandas as pd
import pytest

from biokgsuite.io.schemes import IdSchemes
from biokgsuite.stats.annotation import (AGNOSTIC, CORRECT, UNKNOWN, WRONG_SCHEME, bare_number_types,
                                         verdict)
from biokgsuite.stats.coverage import coverage, schema_coverage


@pytest.fixture(scope='module')
def schemes():
    return IdSchemes()


@pytest.mark.parametrize('raw, declared, expected', [
    ('DB00590', 'Drug', CORRECT),
    ('Compound::DB00590', 'Drug', CORRECT),
    ('PUBCHEM.COMPOUND:2244', 'Drug', CORRECT),          # OpenBioLink's CURIE form
    ('Compound::MESH:D001241', 'Drug', CORRECT),          # DRKG literature-mined substances
    ('7254', 'Disease', CORRECT),                          # PrimeKG bare MONDO
    ('1200_1134_15512_5080_100078', 'Disease', CORRECT),   # PrimeKG grouped MONDO
    ('Disease::DOID:14227', 'Disease', CORRECT),
    ('Disease::MESH:D001943', 'Disease', CORRECT),         # DRKG
    ('UMLS:C0006142', 'Disease', CORRECT),                 # MATRIX, unbridged
    ('ICD9:804.61', 'Disease', CORRECT),                   # MATRIX
    ('D001943', 'Disease', CORRECT),
    ('MIM:139320', 'Disease', CORRECT),
    ('Gene::1956', 'Gene/Protein', CORRECT),
    ('NCBIGENE:1956', 'Gene/Protein', CORRECT),            # OpenBioLink
    ('ENSEMBL:ENST00000489698', 'Gene/Protein', UNKNOWN),  # a transcript is not a gene
    ('REACTOME:R-HSA-109582', 'Pathway', CORRECT),         # OpenBioLink
    ('CL:0000057', 'Cell type', CORRECT),
    ('P00533', 'Gene/Protein', CORRECT),
    ('8150', 'Biological Process', CORRECT),               # PrimeKG bare GO
    ('GO:0008150', 'Biological Process', CORRECT),
    ('2', 'Anatomy', CORRECT),                             # PrimeKG bare UBERON
    ('118', 'Phenotype', CORRECT),                         # PrimeKG bare HPO
    ('R-HSA-1006173', 'Complex', CORRECT),
    ('Tax::9606', 'Taxon', CORRECT),
    ('Atc::L01XE', 'ATC class', CORRECT),
    ('DB00590', 'Disease', WRONG_SCHEME),                  # a drug accession filed as a disease
    ('R-HSA-109582', 'Drug', WRONG_SCHEME),
    ('not-an-id', 'Drug', UNKNOWN),
    # A UMLS CUI names any kind of concept: under Disease it is a recognised
    # shape, under Anatomy it is not a mislabelling but an identifier the
    # check cannot judge (MATRIX files 110,327 of them under AnatomicalEntity).
    ('UMLS:C0006142', 'Disease', CORRECT),
    ('UMLS:C0006142', 'Anatomy', AGNOSTIC),
    ('C0000727', 'Gene/Protein', AGNOSTIC),
    ('MESH:D001943', 'Anatomy', AGNOSTIC),
    ('EFO:0009491', 'Phenotype', AGNOSTIC),
])
def test_verdict_on_real_shapes(schemes, raw, declared, expected):
    assert verdict(raw, declared, schemes) == expected


def test_schemes_examples_match_their_patterns(schemes):
    assert schemes.self_check() == []


def test_bare_number_types_reports_every_bare_type():
    nodes = pd.DataFrame({
        'type': ['disease', 'disease', 'gene/protein', 'drug', 'anatomy', 'Gene'],
        'id':   ['7254', '21223_385', '1956', 'DB00001', '2', 'Gene::1956']})
    tm = {'x': {'Disease': ('disease',), 'Gene/Protein': ('gene/protein', 'Gene'),
                'Drug': ('drug',), 'Anatomy': ('anatomy',)}}
    out = bare_number_types({'x': {'nodes': nodes}}, tm).set_index('Declared type')
    assert set(out.index) == {'disease', 'gene/protein', 'Gene', 'anatomy'}   # drug has none
    assert out.loc['disease', 'Bare numbers'] == 1 and out.loc['disease', 'Grouped numbers'] == 1
    assert out.loc['disease', 'Read as'] == 'MONDO'
    assert out.loc['anatomy', 'Read as'] == 'UBERON'
    assert out.loc['Gene', 'Bare numbers'] == 1                                 # Gene::1956 counts


def test_coverage_arithmetic():
    res = coverage({'a', 'b', 'c', 'z'}, {'a', 'b', 'd', 'e'})
    assert res['overlap_n'] == 2 and res['gold_n'] == 4 and res['kg_n'] == 4
    assert res['coverage_pct'] == 50.0
    assert coverage(set(), {'a'})['coverage_pct'] == 0.0


def test_schema_coverage_flags_undeclared_types():
    nodes = pd.DataFrame({'type': ['drug'] * 3 + ['exposure'] * 2, 'id': list('abcde')})
    kgs = {'x': {'nodes': nodes}}
    row = schema_coverage(kgs, {'x': {'Drug': ('drug',)}}).iloc[0]
    assert row['Mapped'] == 1 and row['Unmapped'] == 1 and row['Nodes mapped (%)'] == 60.0
    assert row['Unmapped types'] == 'exposure x2'
    row = schema_coverage(kgs, {'x': {'Drug': ('drug',), 'Exposure': ('exposure',)}}).iloc[0]
    assert row['Unmapped'] == 0 and row['Unmapped types'] == '-'
