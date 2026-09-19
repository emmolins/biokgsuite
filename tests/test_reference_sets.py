"""The rules that define the coverage reference sets, on hand-made inputs.

    python -m pytest tests/test_reference_sets.py -q
"""
import importlib.util
import sys
from pathlib import Path

import pandas as pd

spec = importlib.util.spec_from_file_location(
    'build_reference_sets', Path(__file__).resolve().parents[1] / 'tools' / 'build_reference_sets.py')
brs = importlib.util.module_from_spec(spec)
sys.modules['build_reference_sets'] = brs
spec.loader.exec_module(brs)


def test_descendants_walks_the_is_a_graph():
    isa = {('2', '1'), ('3', '2'), ('4', '1'), ('9', '8')}
    assert brs.descendants('1', isa) == {'1', '2', '3', '4'}
    assert brs.descendants('8', isa) == {'8', '9'}


def test_approved_therapeutics_rules():
    dm = pd.DataFrame([
        # id, parentId, drugType, stage
        ('CHEMBL707',     None,        'Small molecule',    'APPROVAL'),   # doxazosin: indicated
        ('CHEMBL1200561', 'CHEMBL707', 'Small molecule',    'APPROVAL'),   # its mesylate salt
        ('CHEMBL1200415', None,        'Small molecule',    'APPROVAL'),   # imaging agent, no indication
        ('CHEMBL9001',    None,        'Vaccine component', 'APPROVAL'),   # indicated vaccine component
        ('CHEMBL9002',    None,        'Small molecule',    'PHASE_III'),  # not approved
        ('CHEMBL9003',    'CHEMBL9004', 'Small molecule',   'APPROVAL'),   # salt indicated, parent is not
        ('CHEMBL9004',    None,        'Small molecule',    'APPROVAL'),
    ], columns=['id', 'parentId', 'drugType', 'maximumClinicalStage'])
    ci = pd.DataFrame([
        ('CHEMBL707',  'MONDO_0000001', 'APPROVAL'),
        ('CHEMBL9001', 'MONDO_0000002', 'APPROVAL'),
        ('CHEMBL9002', 'MONDO_0000003', 'PHASE_III'),
        ('CHEMBL9003', 'MONDO_0000004', 'APPROVAL'),
        ('CHEMBL1200415', 'MONDO_0000005', 'PHASE_II'),
    ], columns=['drugId', 'diseaseId', 'maxClinicalStage'])
    drugs, parent_of = brs.approved_therapeutics(dm, ci)
    # one row per parent; the salt collapses onto doxazosin
    assert drugs == {'CHEMBL707', 'CHEMBL9004'}
    # the salt's indication carries to the parent
    assert parent_of['CHEMBL9003'] == 'CHEMBL9004' and parent_of['CHEMBL707'] == 'CHEMBL707'
    # imaging agent without an approved indication, vaccine component, and a
    # phase III molecule are all out
    assert not {'CHEMBL1200415', 'CHEMBL9001', 'CHEMBL9002', 'CHEMBL1200561'} & drugs
