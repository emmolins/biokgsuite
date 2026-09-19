"""Identifier resolution on every id shape the six graphs emit.

Guards the failure that zeroed PrimeKG's disease coverage: a bare MONDO number
('7254'), or a grouped PrimeKG node ('21223_385'), carries no namespace prefix
and used to resolve to nothing. Runs on a small synthetic crosswalk directory.

    python -m pytest tests/test_resolve.py -q
"""
import pytest

from biokgsuite.io.resolve import IdResolver, split_curie

MONDO_OBO = """format-version: 1.2

[Term]
id: MONDO:0007254
name: breast cancer

[Term]
id: MONDO:0000385
name: some disease

[Term]
id: MONDO:0021223
name: another disease

[Term]
id: MONDO:0000001
name: retired term
is_obsolete: true

[Term]
id: MONDO:0000002
name: merged term
is_obsolete: true
replaced_by: MONDO:0007254
"""

SSSOM = """# curie_map: {}
subject_id\tpredicate_id\tobject_id
MONDO:0007254\tskos:exactMatch\tDOID:1612
MONDO:0007254\tskos:exactMatch\tmesh:D001943
MONDO:0007254\tskos:exactMatch\tUMLS:C0006142
MONDO:0007254\tskos:exactMatch\tOMIM:114480
MONDO:0000385\tskos:exactMatch\tOrphanet:1234
MONDO:0021223\tskos:closeMatch\tDOID:9999
"""

GO_OBO = """format-version: 1.2

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process

[Term]
id: GO:0005575
name: cellular_component
namespace: cellular_component
"""

HP_OBO = """format-version: 1.2

[Term]
id: HP:0000118
name: Phenotypic abnormality

[Term]
id: HP:0000001
name: All
is_obsolete: true

[Term]
id: HP:0100237
name: obsolete Proximal foot symphalangism
is_obsolete: true
replaced_by: HP:0000118
"""

UBERON_OBO = """format-version: 1.2

[Term]
id: UBERON:0000002
name: uterine cervix
"""

HGNC = "\t".join(['hgnc_id', 'symbol', 'entrez_id', 'ensembl_gene_id', 'uniprot_ids']) + "\n" + \
       "\t".join(['HGNC:3236', 'EGFR', '1956', 'ENSG00000146648', 'P00533']) + "\n" + \
       "\t".join(['HGNC:11998', 'TP53', '7157', 'ENSG00000141510', 'P04637|Q9Y573']) + "\n"


@pytest.fixture(scope='module')
def resolver(tmp_path_factory):
    d = tmp_path_factory.mktemp('crosswalks')
    (d / 'mondo.obo').write_text(MONDO_OBO)
    (d / 'mondo.sssom.tsv').write_text(SSSOM)
    (d / 'hgnc_complete_set.txt').write_text(HGNC)
    # Open Targets drug_molecule: doxazosin (with a DrugBank xref) and its
    # mesylate salt (no xref); varenicline with no DrugBank xref at all.
    pd = pytest.importorskip('pandas')
    pytest.importorskip('pyarrow')
    dm = pd.DataFrame([
        {'id': 'CHEMBL707', 'parentId': None, 'childChemblIds': ['CHEMBL1200561'],
         'crossReferences': [{'source': 'drugbank', 'ids': ['DB00590']}]},
        {'id': 'CHEMBL1200561', 'parentId': 'CHEMBL707', 'childChemblIds': None, 'crossReferences': None},
        {'id': 'CHEMBL1076903', 'parentId': None, 'childChemblIds': None, 'crossReferences': None},
        {'id': 'CHEMBL1201497', 'parentId': None, 'childChemblIds': None, 'crossReferences': None,
         'name': 'INSULIN GLARGINE', 'synonyms': ['Lantus'], 'tradeNames': None},
    ])
    (d / 'ot_drug_molecule').mkdir()
    dm.to_parquet(d / 'ot_drug_molecule' / 'part-0.parquet', index=False)
    (d / 'unichem_drugbank_chembl.tsv').write_text('CHEMBL1076903\tDB01273\nCHEMBL707\tDB00590\n')
    # Name bridge (tools/build_drug_name_xref.py output): insulin glargine has
    # no structure mapping anywhere, so only its name can reach ChEMBL.
    (d / 'name_bridge.csv').write_text(
        'name_normalised,name,drugbank_id,sources,ambiguous\n'
        'insulin glargine,Insulin glargine,DB00047,hetionet,0\n'
        'doxazosin,Doxazosin,DB00590,hetionet,0\n'
        'lantus,Lantus,DB99001,biokg,1\n')
    (d / 'go.obo').write_text(GO_OBO)
    (d / 'uberon.obo').write_text(UBERON_OBO)
    (d / 'hp.obo').write_text(HP_OBO)
    return IdResolver(d, verbose=False, aux_registries=True, name_bridge=d / 'name_bridge.csv')


def _one(res, raw):
    out, st = res(raw if isinstance(raw, list) else [raw])
    return out, st


def test_split_curie_shapes():
    assert split_curie('DOID:1612') == ('DOID', '1612')
    assert split_curie('Disease::MESH:D001943') == ('MESH', 'D001943')
    assert split_curie('7254') == (None, '7254')


def test_bare_mondo_number_is_a_disease(resolver):
    # PrimeKG's native disease id and MATRIX's bridged id are both bare numbers.
    out, st = _one(resolver.diseases_to_mondo, ['7254', '0007254'])
    assert out == {'7254'} and st.n_resolved == 2


def test_bare_number_outside_the_registry_does_not_resolve(resolver):
    out, st = _one(resolver.diseases_to_mondo, ['9999999', '1'])   # unknown; obsolete
    assert out == set() and st.n_resolved == 0


def test_retired_mondo_follows_replaced_by(resolver):
    # PrimeKG carries 2022 MONDO ids; a merged one resolves to its successor.
    out, st = _one(resolver.diseases_to_mondo, ['2', 'MONDO:0000002'])
    assert out == {'7254'} and st.n_resolved == 2


def test_grouped_primekg_node_resolves_every_member(resolver):
    out, st = _one(resolver.diseases_to_mondo, ['21223_385'])
    assert out == {'21223', '385'} and st.n_resolved == 1
    assert dict(st.by_namespace) == {'grouped': 1}


def test_prefixed_disease_ids(resolver):
    out, _ = _one(resolver.diseases_to_mondo,
                  ['MONDO:0007254', 'DOID:1612', 'Disease::MESH:D001943', 'UMLS:C0006142',
                   'MIM:114480', 'ORPHANET:1234', 'orphanet:1234'])
    assert out == {'7254', '385'}
    assert resolver.diseases_to_mondo(['ICD9:804.61'])[1].n_resolved == 0


def test_gene_curie_shapes(resolver):
    out, st = _one(resolver.genes_to_ensembl,
                   ['ENSG00000146648', 'ENSEMBL:ENSG00000141510', 'NCBIGene:1956', '7157',
                    'UniProtKB:Q9Y573-2', 'HGNC:3236', 'EGFR'])
    assert out == {'ENSG00000146648', 'ENSG00000141510'} and st.n_resolved == 7


def test_non_gene_ensembl_ids_stay_unresolved(resolver):
    _, st = _one(resolver.genes_to_ensembl, ['ENSEMBL:ENST00000489698', 'ENSEMBL:ENSMUSG00000000001'])
    assert st.n_resolved == 0
    assert dict(st.failed_by_namespace) == {'ENSEMBL': 2}


def test_unprefixed_ids_are_reported_by_shape(resolver):
    # BioKG's non-human proteins must read as UniProt accessions, not 'bare'.
    _, st = _one(resolver.genes_to_ensembl, ['Q9DB72', 'A0A0B4J2D5', '999999999', 'P00533'])
    assert dict(st.failed_by_namespace) == {'UniProt': 2, 'integer': 1}
    assert dict(st.by_namespace) == {'UniProt': 1}
    _, st = _one(resolver.drugs_to_chembl, ['DB99999'])
    assert dict(st.failed_by_namespace) == {'DrugBank': 1}


def test_drugbank_resolves_to_the_chembl_family(resolver):
    # The parent carries the DrugBank xref; the salt is scored as held too.
    out, st = _one(resolver.drugs_to_chembl, ['DB00590', 'Compound::DB00590'])
    assert out == {'CHEMBL707', 'CHEMBL1200561'} and st.n_resolved == 2
    # The salt id itself expands back to the parent.
    out, _ = _one(resolver.drugs_to_chembl, ['CHEMBL1200561'])
    assert out == {'CHEMBL707', 'CHEMBL1200561'}


def test_unichem_fills_the_drugbank_gap(resolver):
    # Varenicline has no DrugBank xref in Open Targets; UniChem supplies it.
    out, st = _one(resolver.drugs_to_chembl, ['DB01273'])
    assert out == {'CHEMBL1076903'} and st.n_resolved == 1


def test_go_and_uberon_accept_primekg_bare_numbers(resolver):
    # PrimeKG writes GO:0008150 as '8150' and UBERON:0000002 as '2'.
    out, st = _one(resolver.bp_terms, ['8150', 'GO:0008150', 'bioprocess::GO:0008150'])
    assert out == {'GO:0008150'} and st.n_resolved == 3
    assert resolver.bp_terms(['5575'])[1].n_resolved == 0       # a CC term is not a BP term
    assert resolver.cc_terms(['5575'])[0] == {'GO:0005575'}
    out, st = _one(resolver.anatomy_to_uberon, ['2', 'UBERON:0000002', 'Anatomy::UBERON:0000002'])
    assert out == {'UBERON:0000002'} and st.n_resolved == 3


def test_hpo_accepts_primekg_bare_numbers(resolver):
    out, st = _one(resolver.phenotypes_to_hpo, ['118', 'HP:0000118', 'effect/phenotype::HP:0000118'])
    assert out == {'HP:0000118'} and st.n_resolved == 3
    assert resolver.phenotypes_to_hpo(['1'])[1].n_resolved == 0     # obsolete, no successor
    assert resolver.phenotypes_to_hpo(['100237'])[0] == {'HP:0000118'}   # obsolete, merged
    assert 'Phenotype' in resolver.RESOLVERS


def test_name_bridge_reaches_drugs_without_a_structure_mapping(resolver):
    out, st = _one(resolver.drugs_to_chembl, ['DB00047'])
    assert out == {'CHEMBL1201497'} and st.n_resolved == 1
    assert resolver.name_bridged == 1
    # A structure-mapped accession is not touched by the name bridge, and an
    # ambiguous name (Lantus) is never used.
    assert resolver.drugs_to_chembl(['DB00590'])[0] == {'CHEMBL707', 'CHEMBL1200561'}
    assert resolver.drugs_to_chembl(['DB99001'])[1].n_resolved == 0
