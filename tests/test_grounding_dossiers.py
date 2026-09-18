"""Unit tests for biokgsuite.grounding.dossiers on a synthetic graph.

The graph is small enough to read by hand and is shaped like DRKG (``Type::id``
node IDs, DrugBank relation names) so that the entity-type policy, the ADME
exclusion and the DOID/MeSH disease scheme are all exercised without data.

    python -m pytest tests/test_grounding_dossiers.py -q
"""
from collections import defaultdict

import pandas as pd
import pytest

from biokgsuite.grounding import dossiers
from biokgsuite.grounding.dossiers import KGDossiers, build_edge_index, gene_verb, make_resolvers

KG = 'drkg'
CONFIG = {'knowledge_graphs': {KG: {
    'entity_types': {'Drug': ('Compound',), 'Disease': ('Disease',),
                     'Gene/Protein': ('Gene',), 'Pathway': ('Pathway',)},
    'disease_id_scheme': 'doid_mesh',
}}}

NODES = pd.DataFrame([
    (0, 'Compound::DB001', 'Compound', 'aspirin'),
    (1, 'Disease::DOID:2841', 'Disease', 'asthma'),
    (2, 'Gene::3716', 'Gene', '3716'),
    (3, 'Gene::1576', 'Gene', '1576'),
    (4, 'Pathway::R-HSA-1', 'Pathway', 'Cytokine signalling'),
    (5, 'Gene::7124', 'Gene', '7124'),
    (6, 'Gene::3717', 'Gene', '3717'),
    (7, 'Gene::99999', 'Gene', '99999'),
    (8, 'Compound::DB002', 'Compound', 'ibuprofen'),
    (9, 'Disease::MESH:D001249', 'Disease', 'asthma (MeSH)'),
    (10, 'Gene::5', 'Gene', '5'),
], columns=['idx', 'id', 'type', 'name'])

EDGES = pd.DataFrame([
    (0, 2, 'DRUGBANK::target::Compound:Gene'),
    (0, 3, 'DRUGBANK::enzyme::Compound:Gene'),        # ADME protein, excluded
    (0, 4, 'PATHWAY::member'),
    (0, 1, 'DRUGBANK::treats::Compound:Disease'),     # the leak edge
    (0, 6, 'GNBR::N::Compound:Gene'),                 # inhibition modality
    (0, 7, 'DRUGBANK::target::Compound:Gene'),        # unreadable gene, dropped
    (0, 10, 'DRUGBANK::target::Compound:Gene'),
    (1, 2, 'DGN::associated'),
    (1, 5, 'GNBR::U::Gene:Disease'),
    (1, 9, 'Disease::related'),
    (2, 6, 'PPI'),                                    # raises JAK1's degree above SYM5's
], columns=['x_index', 'y_index', 'relation'])

ENTREZ = {'3716': 'JAK1', '1576': 'CYP3A4', '7124': 'TNF', '3717': 'JAK2', '5': 'SYM5'}


def _crosswalk(**extra):
    xw = {'mondo_to_doid_all': defaultdict(set, {'4979': {'2841'}}),
          'doid_to_mondo_all': defaultdict(set, {'2841': {'4979'}}),
          'mesh_to_doid': {'MESH:D001249': '2841'},
          'doid_to_mesh': defaultdict(set, {'2841': {'MESH:D001249'}}),
          'drug_equiv': defaultdict(set),
          'disease_equiv': defaultdict(set)}
    xw.update(extra)
    return xw


@pytest.fixture
def make(monkeypatch):
    monkeypatch.setattr(dossiers, 'load_kg', lambda kg_name, config: (EDGES.copy(), NODES.copy()))
    monkeypatch.setattr(dossiers, 'load_crosswalks', lambda base: _crosswalk())

    def _make(**kw):
        kw.setdefault('cap', 8)
        kw.setdefault('entrez_to_symbol', ENTREZ)
        return KGDossiers(KG, CONFIG, base=None, **kw)
    return _make


def test_drug_dossier_lists_genes_and_pathways_never_the_disease(make):
    d = make()
    text = d.drug_dossier('DB001', 'aspirin', 'MONDO:0004979')
    assert text == ('aspirin inhibits JAK2; '
                    'aspirin modulates pathway Cytokine signalling; '
                    'aspirin targets SYM5; aspirin targets JAK1')
    assert 'asthma' not in text
    assert 'CYP3A4' not in text                       # PK_EXCLUDE
    assert '99999' not in text                        # unreadable gene ID


def test_drug_dossier_honours_cap(make):
    text = make(cap=1).drug_dossier('DB001', 'aspirin', 'MONDO:0004979')
    # One neighbour per phrase group, the lowest-degree one first.
    assert text == ('aspirin inhibits JAK2; '
                    'aspirin modulates pathway Cytokine signalling; '
                    'aspirin targets SYM5')


def test_drug_dossier_none_when_unresolved_or_empty(make):
    d = make()
    assert d.drug_dossier('DB999', 'nothing', 'MONDO:0004979') is None
    assert d.drug_dossier('DB002', 'ibuprofen', 'MONDO:0004979') is None
    assert d.resolves('DB001') and d.resolves('DB002') and not d.resolves('DB999')


def test_disease_profile_lists_associated_genes(make):
    d = make()
    assert d.resolves_disease('MONDO:0004979')
    assert not d.resolves_disease('MONDO:9999999')
    text = d.disease_profile('MONDO:0004979', 'asthma')
    assert text == ('asthma is associated with gene TNF; asthma is associated with gene JAK1; '
                    'asthma is related to asthma (MeSH)')
    assert 'aspirin' not in text
    assert d.disease_profile('MONDO:0004979') .startswith('asthma is associated')
    assert d.disease_profile('MONDO:9999999') == ''


def test_disease_profile_empty_in_mechanism_only_mode(make):
    assert make(bridge_mode='mechanism_only').disease_profile('MONDO:0004979') == ''


def test_anonymised_genes_are_coded_consistently_per_disease(make):
    d = make(anonymize_genes=True)
    drug = d.drug_dossier('DB001', 'aspirin', 'MONDO:0004979')
    disease = d.disease_profile('MONDO:0004979', 'asthma')
    assert 'JAK1' not in drug and 'JAK1' not in disease
    assert 'Cytokine signalling' in drug           # only genes are coded
    jak1_code = [t.split()[-1] for t in drug.split('; ') if t.startswith('aspirin targets')][-1]
    assert jak1_code.startswith('Gene-')
    assert f'asthma is associated with gene {jak1_code}' in disease


def test_numeric_gene_ids_render_raw_without_symbol_map(make):
    text = make(entrez_to_symbol=None).drug_dossier('DB001', 'aspirin', 'MONDO:0004979')
    # Without a symbol map every gene name is digits, so all gene neighbours
    # are dropped as unreadable and only the pathway remains.
    assert text == 'aspirin modulates pathway Cytokine signalling'


def test_gene_verb_modalities():
    assert gene_verb('drug', {'DRUGBANK::target::Compound:Gene'}) == 'targets'
    assert gene_verb('drug', {'GNBR::N::Compound:Gene'}) == 'inhibits'
    assert gene_verb('drug', {'GNBR::A+::Compound:Gene'}) == 'activates'
    assert gene_verb('drug', {'Hetionet::CuG::Compound:Gene'}) == 'upregulates'
    assert gene_verb('drug', {'Hetionet::CbG::Compound:Gene'}) == 'binds'
    assert gene_verb('disease', {'Hetionet::DuG::Disease:Gene'}) == 'upregulates'
    assert gene_verb('disease', {'DGN::associated'}) == 'is associated with gene'


def test_build_edge_index_is_undirected():
    ix = build_edge_index(EDGES)
    assert (2, 'DRUGBANK::target::Compound:Gene') in ix[0]
    assert (0, 'DRUGBANK::target::Compound:Gene') in ix[2]
    assert len(ix[0]) == 7


def _resolver_config(scheme):
    return {'knowledge_graphs': {KG: {'disease_id_scheme': scheme}}}


def test_resolver_mondo_scheme():
    _, to_dis = make_resolvers(_resolver_config('mondo'), KG, _crosswalk())
    ids = to_dis('MONDO:0004979')
    assert {'4979', '0004979', 'MONDO:4979', 'MONDO_4979'} <= ids
    assert 'DOID:2841' not in ids


def test_resolver_doid_scheme_from_any_input_form():
    _, to_dis = make_resolvers(_resolver_config('doid'), KG, _crosswalk())
    for given in ('MONDO:0004979', 'DOID:2841', 'D001249', '0004979'):
        ids = to_dis(given)
        assert {'2841', 'DOID:2841', 'DOID_2841', 'doid:2841'} <= ids, given
        assert 'MESH:D001249' not in ids


def test_resolver_doid_mesh_scheme():
    _, to_dis = make_resolvers(_resolver_config('doid_mesh'), KG, _crosswalk())
    ids = to_dis('MONDO:0004979')
    assert {'DOID:2841', '2841', 'MESH:D001249', 'D001249'} <= ids


def test_resolver_mesh_scheme():
    _, to_dis = make_resolvers(_resolver_config('mesh'), KG, _crosswalk())
    ids = to_dis('MONDO:0004979')
    assert {'MESH:D001249', 'D001249', 'MESH_D001249'} <= ids
    assert 'DOID:2841' not in ids and '2841' not in ids


def test_resolver_widens_by_sssom_equivalence_and_drug_xref():
    xw = _crosswalk(disease_equiv=defaultdict(set, {'MONDO:4979': {'UMLS:C0004096'}}),
                    drug_equiv=defaultdict(set, {'DB001': {'PUBCHEM.COMPOUND:2244'}}))
    to_drug, to_dis = make_resolvers(_resolver_config('mondo'), KG, xw)
    assert {'UMLS:C0004096', 'C0004096'} <= to_dis('MONDO:0004979')
    assert {'DB001', 'PUBCHEM.COMPOUND:2244', 'pubchem.compound:2244', '2244'} <= to_drug('DB001')
    assert to_drug('DB777') == {'DB777'}


def test_load_crosswalks_from_files_without_drugbank_xref(tmp_path):
    bridges = tmp_path / 'data' / 'gold_standards' / 'bridges'
    bridges.mkdir(parents=True)
    (bridges / 'do_diseases.csv').write_text(
        'doid,mondo_id,mondo_name\nDOID:2841,MONDO:0004979,asthma\nDOID:0050156,,unmapped\n')
    (bridges / 'mesh_to_doid.csv').write_text(
        'mesh_id,doid,source\nMESH:D001249,DOID:0050156,MONDO\nMESH:D001249 {x},DOID:2841,DO\n')
    (bridges / 'mondo.sssom.tsv').write_text(
        '# curie_map: {}\nsubject_id\tpredicate_id\tobject_id\n'
        'MONDO:4979\tskos:exactMatch\tUMLS:C0004096\n'
        'UMLS:C0004096\tskos:exactMatch\tNCIT:C2985\n'
        'MONDO:4979\tskos:closeMatch\tEFO:0000270\n')
    xw = dossiers.load_crosswalks(tmp_path)
    assert xw['drug_equiv'] == {}
    assert xw['mondo_to_doid_all']['4979'] == {'2841'}
    assert xw['doid_to_mondo_all']['2841'] == {'4979'}
    assert xw['mesh_to_doid'] == {'MESH:D001249': '2841'}        # DO source wins over MONDO
    assert xw['doid_to_mesh']['2841'] == {'MESH:D001249'}
    assert xw['disease_equiv']['MONDO:4979'] == {'UMLS:C0004096', 'NCIT:C2985'}
    assert 'EFO:0000270' not in xw['disease_equiv']['MONDO:4979']


def test_load_crosswalks_requires_sssom(tmp_path):
    bridges = tmp_path / 'data' / 'gold_standards' / 'bridges'
    bridges.mkdir(parents=True)
    (bridges / 'do_diseases.csv').write_text('doid,mondo_id,mondo_name\n')
    (bridges / 'mesh_to_doid.csv').write_text('mesh_id,doid\n')
    with pytest.raises(FileNotFoundError):
        dossiers.load_crosswalks(tmp_path)


def test_load_entrez_to_symbol_required(tmp_path):
    with pytest.raises(FileNotFoundError):
        dossiers.load_entrez_to_symbol(tmp_path)
    maps = tmp_path / 'data' / 'gold_standards' / 'bridges'
    maps.mkdir(parents=True)
    (maps / 'entrez_to_symbol.tsv').write_text('3716\tJAK1\nbad line\n7124\tTNF\n')
    assert dossiers.load_entrez_to_symbol(tmp_path) == {'3716': 'JAK1', '7124': 'TNF'}
