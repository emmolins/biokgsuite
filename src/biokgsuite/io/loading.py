"""Knowledge graph loaders.

Each loader returns a standardized tuple:
    (kg_df, nodes_df)

where:
  - kg_df has columns: relation, x_index, x_id, x_type, x_name, y_index, y_id, y_type, y_name
    (plus optional: display_relation, x_source, y_source)
  - nodes_df has columns: idx, id, type, name
"""

import gzip
import re
import pandas as pd
from pathlib import Path


# find_config() and load_config() moved to biokgsuite.io.config, which
# understands the config/ directory layout (kgs.yaml + params.yaml) and can
# locate the repo root from any depth. Import them from biokgsuite.io.


def load_primekg(path):
    """Load PrimeKG from its CSV file.

    Parameters
    ----------
    path : str or Path
        Path to primekg.csv.

    Returns
    -------
    kg : pd.DataFrame
        Edge table with standardized columns.
    nodes : pd.DataFrame
        Node table with columns [idx, id, type, name].
    """
    cols = {
        'relation', 'display_relation',
        'x_index', 'x_id', 'x_type', 'x_name', 'x_source',
        'y_index', 'y_id', 'y_type', 'y_name', 'y_source',
    }
    # Chunked read prevents peak-memory spikes on the ~8 M-row CSV.
    # pandas still parses full rows before filtering, so low_memory=False on
    # the full file routinely exceeds 3 GB; 200 k-row chunks stay comfortably
    # below that ceiling.
    chunks = []
    for chunk in pd.read_csv(path, usecols=lambda c: c in cols,
                              dtype={'x_id': str, 'y_id': str},
                              low_memory=True, chunksize=200_000):
        chunks.append(chunk)
    kg = pd.concat(chunks, ignore_index=True)
    del chunks

    nodes = _build_node_table(kg)
    return kg, nodes


def load_hetionet(nodes_path, edges_path):
    """Load Hetionet from TSV files.

    Parameters
    ----------
    nodes_path : str or Path
        Path to hetionet nodes TSV.
    edges_path : str or Path
        Path to hetionet edges TSV.

    Returns
    -------
    kg : pd.DataFrame   (standardized columns)
    nodes : pd.DataFrame (columns: idx, id, type, name)
    """
    raw_nodes = pd.read_csv(nodes_path, sep='\t')
    raw_edges = pd.read_csv(edges_path, sep='\t')

    # Hetionet nodes: id, name, kind
    raw_nodes = raw_nodes.rename(columns={'kind': 'type'})
    raw_nodes['idx'] = range(len(raw_nodes))
    id_to_idx = dict(zip(raw_nodes['id'], raw_nodes['idx']))

    # Hetionet edges: source, metaedge, target
    kg = pd.DataFrame({
        'relation': raw_edges['metaedge'],
        'x_index': raw_edges['source'].map(id_to_idx),
        'x_id': raw_edges['source'].astype(str),
        'y_index': raw_edges['target'].map(id_to_idx),
        'y_id': raw_edges['target'].astype(str),
    })

    # Merge node metadata
    node_meta = dict(zip(raw_nodes['id'], zip(raw_nodes['type'], raw_nodes['name'])))
    kg['x_type'] = raw_edges['source'].map(lambda x: node_meta.get(x, ('', ''))[0])
    kg['x_name'] = raw_edges['source'].map(lambda x: node_meta.get(x, ('', ''))[1])
    kg['y_type'] = raw_edges['target'].map(lambda x: node_meta.get(x, ('', ''))[0])
    kg['y_name'] = raw_edges['target'].map(lambda x: node_meta.get(x, ('', ''))[1])

    # Per-edge properties, present when the edge table was built from Hetionet's
    # hetnet JSON by tools/convert_hetionet_json.py and absent from the plain TSV
    # export. `edge_source` is the database the edge was gathered from, which is
    # what the trustworthiness dimension scores as per-edge provenance; reading
    # the TSV alone made Hetionet appear to have none. Carried through only when
    # present, so both export formats still load.
    for _col in ('edge_source', 'unbiased'):
        if _col in raw_edges.columns:
            kg[_col] = raw_edges[_col]

    kg = kg.dropna(subset=['x_index', 'y_index'])
    kg['x_index'] = kg['x_index'].astype(int)
    kg['y_index'] = kg['y_index'].astype(int)

    nodes_df = raw_nodes[['idx', 'id', 'type', 'name']].copy()
    return kg, nodes_df


def load_drkg(path):
    """Load DRKG from its TSV file.

    DRKG format: head \\t relation \\t tail  (no header).
    Entity IDs encode their type, e.g. "Gene::1234", "Compound::DB00001".

    Parameters
    ----------
    path : str or Path

    Returns
    -------
    kg : pd.DataFrame   (standardized columns)
    nodes : pd.DataFrame (columns: idx, id, type, name)
    """
    raw = pd.read_csv(path, sep='\t', header=None, names=['head', 'relation', 'tail'])

    # Extract entity type from ID prefix
    def _parse_entity(eid):
        parts = eid.split('::', 1)
        return (parts[0], parts[1]) if len(parts) == 2 else ('unknown', eid)

    all_entities = set(raw['head']) | set(raw['tail'])
    entity_list = sorted(all_entities)
    id_to_idx = {eid: i for i, eid in enumerate(entity_list)}

    nodes_records = []
    for eid in entity_list:
        etype, ename = _parse_entity(eid)
        nodes_records.append({
            'idx': id_to_idx[eid],
            'id': eid,
            'type': etype,
            'name': ename,
        })
    nodes_df = pd.DataFrame(nodes_records)
    type_map = dict(zip(nodes_df['id'], nodes_df['type']))
    name_map = dict(zip(nodes_df['id'], nodes_df['name']))

    kg = pd.DataFrame({
        'relation': raw['relation'],
        'x_index': raw['head'].map(id_to_idx),
        'x_id': raw['head'],
        'x_type': raw['head'].map(type_map),
        'x_name': raw['head'].map(name_map),
        'y_index': raw['tail'].map(id_to_idx),
        'y_id': raw['tail'],
        'y_type': raw['tail'].map(type_map),
        'y_name': raw['tail'].map(name_map),
    })

    return kg, nodes_df


def _openbilink_hq_mask(score_series, source_series):
    """Return boolean mask selecting high-quality OpenBioLink edges.

    Quality criteria (per the OpenBioLink benchmark specification):
    * STRING / STITCH (numeric confidence scores): keep score >= 700
    * Bgee (expression quality labels): keep 'gold quality' or 'high quality'
    * All other sources (curated): keep unconditionally

    This reduces the ~29 M-row ALL_DIR dump to ~5 M HQ edges while
    retaining *all* entity types and curated relation classes.
    """
    score  = score_series.fillna('')
    source = source_series.fillna('')
    # Numeric-score sources: STRING / STITCH
    numeric_src = source.isin(['STRING', 'STITCH'])
    score_num   = pd.to_numeric(score, errors='coerce')
    num_pass    = numeric_src & (score_num >= 700)
    # Bgee quality labels
    bgee_src    = source == 'Bgee'
    bgee_pass   = bgee_src & score.isin(['gold quality', 'high quality'])
    # All other sources pass unconditionally
    other_src   = ~numeric_src & ~bgee_src
    return num_pass | bgee_pass | other_src


def load_openbilink(path, pubchem_to_drugbank_path=None, hq_filter=True):
    """Load OpenBioLink from its TSV edges file.

    OpenBioLink format (TSV, no header, 5 columns):
        subject \\t relation \\t object \\t score \\t source

    Entity types are inferred from ID prefixes:
        NCBIGENE:*          → Gene
        PUBCHEM.COMPOUND:*  → Drug
        DOID:*              → Disease
        GO:*                → GO
        UBERON:*            → Anatomy
        HP:*                → Phenotype
        CL:*                → Cell
        REACTOME:*          → Pathway
        KEGG:*              → Pathway

    Parameters
    ----------
    path : str or Path
        Path to edges.csv (tab-separated, no header).
    pubchem_to_drugbank_path : str or Path, optional
        Path to pubchem_to_drugbank.csv mapping file.  When provided the
        mapping is attached to the nodes DataFrame as a 'drugbank_id' column
        on Drug nodes (empty string for unmapped nodes).
    hq_filter : bool, default True
        Apply the OpenBioLink high-quality filter (STRING/STITCH >= 700,
        Bgee gold/high quality, all curated kept).  Reduces ~29 M edges
        to ~5 M.  Set False only if sufficient RAM is available.

    Returns
    -------
    kg : pd.DataFrame   (standardized columns)
    nodes : pd.DataFrame (columns: idx, id, type, name, drugbank_id)
    """
    import gc

    # ── ID-prefix → canonical type mapping ─────────────────────────────────────
    _PREFIX_TO_TYPE = {
        'NCBIGENE':          'Gene',
        'PUBCHEM.COMPOUND':  'Drug',
        'DOID':              'Disease',
        'GO':                'GO',
        'UBERON':            'Anatomy',
        'HP':                'Phenotype',
        'CL':                'Cell',
        'REACTOME':          'Pathway',
        'KEGG':              'Pathway',
    }

    def _type_from_id(eid):
        colon_pos = str(eid).find(':')
        if colon_pos == -1:
            return 'unknown'
        return _PREFIX_TO_TYPE.get(str(eid)[:colon_pos], str(eid)[:colon_pos])

    def _name_from_id(eid):
        s = str(eid)
        colon_pos = s.find(':')
        return s[colon_pos + 1:] if colon_pos != -1 else s

    # ── Read, apply HQ filter, and collect entity IDs in one pass ─────────────
    all_entities = set()
    filtered_chunks = []
    for chunk in pd.read_csv(path, sep='\t', header=None,
                             names=['head', 'relation', 'tail', 'score', 'source'],
                             dtype=str, chunksize=2_000_000,
                             encoding='utf-8', lineterminator='\n'):
        # Strip Windows \r from all string columns
        for col in chunk.columns:
            chunk[col] = chunk[col].str.rstrip('\r')
        chunk = chunk.fillna('')

        if hq_filter:
            mask = _openbilink_hq_mask(chunk['score'], chunk['source'])
            chunk = chunk[mask].copy()

        if chunk.empty:
            continue

        all_entities.update(chunk['head'])
        all_entities.update(chunk['tail'])
        # Keep only columns needed for the edge table
        filtered_chunks.append(
            chunk[['head', 'relation', 'tail', 'source']].copy()
        )
        gc.collect()

    # ── Build node table ───────────────────────────────────────────────────────
    entity_list = sorted(all_entities)
    del all_entities
    id_to_idx = {eid: i for i, eid in enumerate(entity_list)}

    ids_series = pd.Series(entity_list, dtype='object')
    nodes_df = pd.DataFrame({
        'idx':  range(len(entity_list)),
        'id':   ids_series,
        'type': ids_series.map(_type_from_id),
        'name': ids_series.map(_name_from_id),
    })
    del ids_series, entity_list

    # ── Attach DrugBank IDs via PubChem→DrugBank mapping ──────────────────────
    nodes_df['drugbank_id'] = ''
    if pubchem_to_drugbank_path is not None:
        pcdb = pd.read_csv(pubchem_to_drugbank_path, dtype=str)
        pc_map = dict(zip(
            'PUBCHEM.COMPOUND:' + pcdb['pubchem_cid'].astype(str),
            pcdb['drugbank_id'].astype(str),
        ))
        drug_mask = nodes_df['type'] == 'Drug'
        nodes_df.loc[drug_mask, 'drugbank_id'] = (
            nodes_df.loc[drug_mask, 'id'].map(pc_map).fillna('')
        )

    # ── Build fast lookup arrays ───────────────────────────────────────────────
    _type_arr = nodes_df['type'].values
    _name_arr = nodes_df['name'].values

    # ── Build edge table from filtered chunks ─────────────────────────────────
    kg_parts = []
    for fchunk in filtered_chunks:
        xi = fchunk['head'].map(id_to_idx)
        yi = fchunk['tail'].map(id_to_idx)
        valid = xi.notna() & yi.notna()
        xi = xi[valid].astype('int32')
        yi = yi[valid].astype('int32')

        part = pd.DataFrame({
            'relation': fchunk.loc[valid, 'relation'].values,
            'x_index':  xi.values,
            'x_id':     fchunk.loc[valid, 'head'].values,
            'x_type':   _type_arr[xi.values],
            'x_name':   _name_arr[xi.values],
            'y_index':  yi.values,
            'y_id':     fchunk.loc[valid, 'tail'].values,
            'y_type':   _type_arr[yi.values],
            'y_name':   _name_arr[yi.values],
            'x_source': fchunk.loc[valid, 'source'].values,
            'y_source': fchunk.loc[valid, 'source'].values,
        })
        kg_parts.append(part)

    del filtered_chunks, id_to_idx, _type_arr, _name_arr
    gc.collect()

    kg = pd.concat(kg_parts, ignore_index=True)
    del kg_parts
    gc.collect()

    return kg, nodes_df


def load_biokg(links_path, metadata_dir=None):
    """Load BioKG from its consolidated links TSV file.

    BioKG format: subject \\t predicate \\t object  (no header, 3 columns).
    Entity types are inferred from the predicate (e.g. PPI → Gene/Protein ×
    Gene/Protein, DPI → Drug × Gene/Protein).  Entity names are resolved from
    optional metadata TSV files when *metadata_dir* is given.

    Parameters
    ----------
    links_path : str or Path
        Path to ``biokg.links.tsv`` (or ``.tsv.gz``).
    metadata_dir : str or Path, optional
        Directory containing ``biokg.metadata.*.tsv`` (or ``.tsv.gz``) files
        used for entity-name lookups.  If *None*, entity IDs are used as names.

    Returns
    -------
    kg : pd.DataFrame   (standardised columns)
    nodes : pd.DataFrame (columns: idx, id, type, name)
    """
    links_path = Path(links_path)

    # ── Predicate → (subject_type, object_type) ────────────────────────────
    _PRED_TYPES = {
        'PPI':                         ('Gene/Protein', 'Gene/Protein'),
        'DPI':                         ('Drug',         'Gene/Protein'),
        'DDI':                         ('Drug',         'Drug'),
        'DRUG_TARGET':                 ('Drug',         'Gene/Protein'),
        'DRUG_CARRIER':                ('Drug',         'Gene/Protein'),
        'DRUG_ENZYME':                 ('Drug',         'Gene/Protein'),
        'DRUG_TRANSPORTER':            ('Drug',         'Gene/Protein'),
        'PROTEIN_DISEASE_ASSOCIATION': ('Gene/Protein', 'Disease'),
        'DRUG_DISEASE_ASSOCIATION':    ('Drug',         'Disease'),
        'PROTEIN_PATHWAY_ASSOCIATION': ('Gene/Protein', 'Pathway'),
        'DRUG_PATHWAY_ASSOCIATION':    ('Drug',         'Pathway'),
        'DISEASE_PATHWAY_ASSOCIATION': ('Disease',      'Pathway'),
        'MEMBER_OF_COMPLEX':           ('Gene/Protein', 'Complex'),
        'MEMBER_OF_PATHWAY':           ('Complex',      'Pathway'),
        'MEMBER_OF_TOP_LEVEL_PATHWAY': ('Complex',      'Pathway'),
        'COMPLEX_IN_PATHWAY':          ('Complex',      'Pathway'),
        'COMPLEX_TOP_LEVEL_PATHWAY':   ('Complex',      'Pathway'),
        'HAS_PARENT_PATHWAY':          ('Pathway',      'Pathway'),
        'DISEASE_GENETIC_DISORDER':    ('Disease',      'GeneticDisorder'),
        'RELATED_GENETIC_DISORDER':    ('Gene/Protein', 'GeneticDisorder'),
        'DISEASE_SUPERGRP':            ('Disease',      'DiseaseCategory'),
        'DRUG_SIDEEFFECT_ASSOCIATION': ('Drug',         'SideEffect'),
        'DRUG_INDICATION_ASSOCIATION': ('Drug',         'Disease'),
        'DRUG_ATC_CODE':               ('Drug',         'ATC'),
        'PROTEIN_EXPRESSED_IN':        ('Gene/Protein', 'Tissue'),
        'PART_OF_TISSUE':              ('Cell',         'Tissue'),
    }

    # Type priority: higher-priority types override lower-priority ones when
    # the same entity is seen in multiple predicates (e.g. an R-HSA-* node
    # first seen as Complex target in MEMBER_OF_COMPLEX, later as Pathway
    # target in PROTEIN_PATHWAY_ASSOCIATION).
    _TYPE_PRIORITY = {
        'Drug': 5, 'Disease': 5, 'Gene/Protein': 5, 'Pathway': 5,
        'GeneticDisorder': 4, 'SideEffect': 4, 'ATC': 4,
        'DiseaseCategory': 3, 'Tissue': 3, 'Cell': 3,
        'Complex': 2,
        'unknown': 0,
    }

    # ── Read triples ────────────────────────────────────────────────────────
    open_fn = gzip.open if str(links_path).endswith('.gz') else open
    with open_fn(links_path, 'rt') as fh:
        raw = pd.read_csv(fh, sep='\t', header=None,
                          names=['head', 'relation', 'tail'], dtype=str)
    raw = raw.dropna(subset=['head', 'relation', 'tail'])

    # ── Assign entity types from predicate ──────────────────────────────────
    x_types = raw['relation'].map(lambda r: _PRED_TYPES.get(r, ('unknown', 'unknown'))[0])
    y_types = raw['relation'].map(lambda r: _PRED_TYPES.get(r, ('unknown', 'unknown'))[1])

    # ── Build entity → type mapping (highest-priority wins) ─────────────────
    entity_type = {}
    for eid, etype in zip(raw['head'], x_types):
        if etype == 'unknown':
            continue
        cur = entity_type.get(eid, 'unknown')
        if _TYPE_PRIORITY.get(etype, 0) > _TYPE_PRIORITY.get(cur, 0):
            entity_type[eid] = etype
    for eid, etype in zip(raw['tail'], y_types):
        if etype == 'unknown':
            continue
        cur = entity_type.get(eid, 'unknown')
        if _TYPE_PRIORITY.get(etype, 0) > _TYPE_PRIORITY.get(cur, 0):
            entity_type[eid] = etype

    # ── Load optional metadata for entity names ─────────────────────────────
    id_to_name = {}
    if metadata_dir is not None:
        metadata_dir = Path(metadata_dir)
        for meta_file in sorted(metadata_dir.glob('biokg.metadata.*.tsv*')):
            mopen = gzip.open if str(meta_file).endswith('.gz') else open
            with mopen(meta_file, 'rt') as fh:
                for line in fh:
                    parts = line.rstrip('\n').split('\t')
                    if len(parts) >= 3 and parts[1] == 'NAME':
                        eid, _, name = parts[0], parts[1], parts[2]
                        if eid not in id_to_name:
                            id_to_name[eid] = name

    # ── Build node table ────────────────────────────────────────────────────
    # BioKG uses bare IDs (e.g. "DB00001", "D000006", "P21917").  While
    # current ID namespaces don't collide in practice, we add a type-
    # qualified ``uid`` column for safe lookups, consistent with PrimeKG.
    all_entities = sorted(set(raw['head']) | set(raw['tail']))
    id_to_idx = {eid: i for i, eid in enumerate(all_entities)}

    _types_list = [entity_type.get(eid, 'unknown') for eid in all_entities]
    nodes_df = pd.DataFrame({
        'idx':  range(len(all_entities)),
        'id':   all_entities,
        'type': _types_list,
        'name': [id_to_name.get(eid, eid) for eid in all_entities],
    })
    nodes_df['uid'] = nodes_df['type'].astype(str) + '/' + nodes_df['id'].astype(str)

    # ── Build edge table ────────────────────────────────────────────────────
    type_arr = nodes_df['type'].values
    name_arr = nodes_df['name'].values

    xi = raw['head'].map(id_to_idx)
    yi = raw['tail'].map(id_to_idx)

    kg = pd.DataFrame({
        'relation': raw['relation'].values,
        'x_index':  xi.values,
        'x_id':     raw['head'].values,
        'x_type':   type_arr[xi.values.astype(int)],
        'x_name':   name_arr[xi.values.astype(int)],
        'y_index':  yi.values,
        'y_id':     raw['tail'].values,
        'y_type':   type_arr[yi.values.astype(int)],
        'y_name':   name_arr[yi.values.astype(int)],
    })

    kg['x_index'] = kg['x_index'].astype(int)
    kg['y_index'] = kg['y_index'].astype(int)

    return kg, nodes_df


def load_kg(kg_name, config):
    """Load a knowledge graph by name using config settings.

    Parameters
    ----------
    kg_name : str
        One of 'primekg', 'hetionet', 'drkg', 'matrix'.
    config : dict
        Parsed config (from biokgsuite.io.load_config).

    Returns
    -------
    kg : pd.DataFrame
    nodes : pd.DataFrame
    """
    cfg = config['knowledge_graphs'][kg_name]
    base = config.get('_base_dir', Path('.'))

    def _resolve(p):
        """Resolve a config path relative to the repo root."""
        p = Path(p)
        return p if p.is_absolute() else base / p

    if kg_name == 'primekg':
        return load_primekg(_resolve(cfg['path']))
    elif kg_name == 'hetionet':
        return load_hetionet(_resolve(cfg['nodes_path']), _resolve(cfg['edges_path']))
    elif kg_name == 'drkg':
        return load_drkg(_resolve(cfg['path']))
    elif kg_name == 'openbilink':
        pc_path = cfg.get('pubchem_to_drugbank_path')
        if pc_path:
            pc_path = _resolve(pc_path)
        hq = cfg.get('hq_filter', True)
        return load_openbilink(_resolve(cfg['path']), pc_path, hq_filter=hq)
    elif kg_name == 'matrix':
        kw = {}
        if 'keep_categories' in cfg:
            kc = cfg['keep_categories']
            # YAML 'all' (string) means keep every category. Pass through verbatim.
            if isinstance(kc, str):
                kw['keep_categories'] = kc
            else:
                kw['keep_categories'] = tuple(kc)
        if cfg.get('doid_to_mondo_path'):
            kw['doid_to_mondo_path'] = _resolve(cfg['doid_to_mondo_path'])
        if cfg.get('mesh_to_doid_path'):
            kw['mesh_to_doid_path'] = _resolve(cfg['mesh_to_doid_path'])
        if cfg.get('mondo_sssom_path'):
            sssom = _resolve(cfg['mondo_sssom_path'])
            if Path(sssom).exists():
                kw['mondo_sssom_path'] = sssom
        if cfg.get('pubchem_to_drugbank_path'):
            kw['pubchem_to_drugbank_path'] = _resolve(cfg['pubchem_to_drugbank_path'])
        if cfg.get('uniprot_to_entrez_path'):
            kw['uniprot_to_entrez_path'] = _resolve(cfg['uniprot_to_entrez_path'])
        if cfg.get('drugbank_xref_path'):
            xref = _resolve(cfg['drugbank_xref_path'])
            if Path(xref).exists():
                kw['drugbank_xref_path'] = xref
        return load_matrix(_resolve(cfg['nodes_path']),
                           _resolve(cfg['edges_path']), **kw)
    elif kg_name == 'biokg':
        meta_dir = cfg.get('metadata_dir')
        if meta_dir:
            meta_dir = _resolve(meta_dir)
        return load_biokg(_resolve(cfg['path']), metadata_dir=meta_dir)
    else:
        raise ValueError(
            f"Unknown KG: {kg_name!r}. "
            f"Valid names: 'primekg', 'hetionet', 'drkg', 'openbilink', 'matrix', 'biokg'"
        )


# Private helpers 

def _build_node_table(kg):
    """Build a deduplicated node table from PrimeKG edge table.

    .. warning::

       PrimeKG external IDs (the ``id`` column) are **not unique across node
       types**.  For example, ``"9891"`` is both gene NUAK1 (NCBI) and disease
       "acquired polycythemia vera" (MONDO).  Any lookup dict built from
       ``id`` alone will silently collide.  Use the ``uid`` column (which
       prefixes the type) or the helper :func:`graph_utils.node_id_lookup`
       for safe lookups.
    """
    nodes = pd.concat([
        kg[['x_index', 'x_id', 'x_type', 'x_name']].rename(
            columns={'x_index': 'idx', 'x_id': 'id', 'x_type': 'type', 'x_name': 'name'}),
        kg[['y_index', 'y_id', 'y_type', 'y_name']].rename(
            columns={'y_index': 'idx', 'y_id': 'id', 'y_type': 'type', 'y_name': 'name'}),
    ]).drop_duplicates(subset=['idx'])
    nodes = nodes.reset_index(drop=True)
    # Type-qualified unique ID — safe for dict lookups across node types.
    nodes['uid'] = nodes['type'].astype(str) + '/' + nodes['id'].astype(str)
    return nodes


def load_matrix(nodes_path, edges_path,
                keep_categories=('drug', 'disease', 'gene/protein',
                                 'pathway', 'effect/phenotype'),
                doid_to_mondo_path=None,
                mesh_to_doid_path=None,
                mondo_sssom_path=None,
                pubchem_to_drugbank_path=None,
                uniprot_to_entrez_path=None,
                drugbank_xref_path=None,
                trim_isolated_nodes=True,
                use_cache=True,
                cache_dir=None,
                node_chunksize=500_000, edge_chunksize=2_000_000,
                verbose=True):
    """Load Every Cure MATRIX knowledge graph from TSV files (streaming).

    The MATRIX KG is large (~5 GB nodes.tsv, ~14 GB edges.tsv at full scale,
    ~44 M nodes / ~87 M edges). This loader streams both files in pandas
    chunks and **filters to the canonical entity-type subset** declared in
    ``keep_categories`` so the resulting DataFrame fits in memory and is
    apples-to-apples with the other 5 BioKGSuite KGs (which only model
    drug/disease/gene/pathway/phenotype).

    The loader also **canonicalises primary IDs** so the standard
    ``entity_set_from_kg`` / ``extract_pairs`` pipeline joins against gold
    standards without Matrix-specific code paths:

    * Drug   — bare DrugBank accession (e.g. ``'DB00001'``) extracted from
               ``id`` or ``equivalent_identifiers``; falls back to original
               CURIE when no DrugBank mapping exists.
    * Disease — bare MONDO numeric (e.g. ``'5812'`` from ``'MONDO:0005812'``)
               extracted from ``id`` or ``equivalent_identifiers``. When a
               MONDO ID isn't directly available, the loader falls back to
               DOID → MONDO and MESH → DOID → MONDO bridges (if the
               corresponding gold-standard CSVs are passed in via
               ``doid_to_mondo_path`` / ``mesh_to_doid_path``). Disease
               nodes that still can't be canonicalised after both bridges
               are dropped (consistent with ``disease_id_scheme: mondo``).
    * Gene/Protein — bare NCBI gene id (e.g. ``'1234'``) from ``id`` /
               ``equivalent_identifiers``; falls back to original CURIE.
    * Pathway — Reactome ID (e.g. ``'R-HSA-1234'``) from ``id`` /
               ``equivalent_identifiers``; falls back to original CURIE.
    * Effect/Phenotype — original CURIE preserved (HP / MP / etc.).

    Auxiliary CURIE-keyed columns (``drugbank_id``, ``ncbigene_id``,
    ``reactome_id``, ``doid_id``) are retained on ``nodes_df`` for
    downstream notebooks that need richer crosswalks.

    Expected file formats
    ---------------------
    nodes_path (TSV) — columns from Every Cure ``kg-nodes``:
        id                      : CURIE-style identifier, e.g. 'DRUGBANK:DB00001'
        name                    : human-readable label
        category                : Biolink category, e.g. 'biolink:SmallMolecule'
        equivalent_identifiers  : repr-style list of alternate CURIEs
        upstream_data_source    : originating database(s) (optional)

    edges_path (TSV) — columns from Every Cure ``kg-edges``:
        subject                 : source node CURIE (matches nodes 'id')
        predicate               : Biolink predicate, 'biolink:' prefix stripped
        object                  : target node CURIE
        primary_knowledge_source / aggregator_knowledge_source (optional)

    Biolink category → canonical type mapping (stored in nodes_df 'type' column):
        biolink:SmallMolecule / Drug / ChemicalEntity  → 'drug'
        biolink:Disease / DiseaseOrPhenotypicFeature   → 'disease'
        biolink:PhenotypicFeature                      → 'effect/phenotype'
        biolink:Gene / Protein / GeneOrGeneProduct     → 'gene/protein'
        biolink:GenomicEntity / NucleicAcidEntity      → 'gene/protein'
        biolink:Pathway                                → 'pathway'
        (other categories are dropped unless added to keep_categories)

    Relation predicates
    -------------------
    The 'biolink:' prefix is stripped so predicates are stored as bare
    strings ('treats_or_applied_or_studied_to_treat',
    'directly_physically_interacts_with', etc.).

    Parameters
    ----------
    nodes_path : str or Path
    edges_path : str or Path
    keep_categories : iterable of str | str | None, default
        ('drug', 'disease', 'gene/protein', 'pathway', 'effect/phenotype')
        Canonical type values to retain. Edges are filtered to those whose
        endpoints both map to a kept category — drops the long tail of
        OrganismTaxon, Transcript, Procedure, Publication, etc. that has no
        analog in the other BioKGSuite KGs. Pass ``'all'`` (string) or
        ``None`` to disable filtering and retain every typed node and edge.
        Warning: at full scale the unfiltered graph is ~44 M nodes / 87 M
        edges and requires ~10 GB RAM for the DataFrames alone.
    doid_to_mondo_path : str or Path, optional
        Path to ``do_diseases.csv`` (gold-standard table with ``mondo_id``
        and ``doid`` columns). When provided, disease nodes lacking a
        direct MONDO ID are bridged via DOID → MONDO so they survive
        ``disease_id_scheme: mondo`` filtering.
    mesh_to_doid_path : str or Path, optional
        Path to ``mesh_to_doid.csv``. When provided alongside
        ``doid_to_mondo_path``, MESH-only disease nodes are bridged via
        MESH → DOID → MONDO. Has no effect without ``doid_to_mondo_path``.
    mondo_sssom_path : str or Path, optional
        Path to MONDO's SSSOM crosswalk (``mondo.sssom.tsv`` from the MONDO
        GitHub repo, fetched via ``scripts/download_mondo_sssom.sh``). When
        provided, disease nodes whose primary CURIE is UMLS / OMIM /
        Orphanet / ICD9 / NCIT / etc. are bridged to MONDO via this table —
        recovers the long tail of non-MONDO diseases that the other two
        bridges can't reach.
    pubchem_to_drugbank_path : str or Path, optional
        Path to a CSV with columns ``drugbank_id`` / ``pubchem_cid`` (the
        same file used by load_openbilink). When provided, drug nodes that
        have no DrugBank crosswalk in equivalent_identifiers but DO have a
        PUBCHEM.COMPOUND ID get canonicalised to DrugBank via this table.
        Significantly raises drug-coverage joins for Matrix.
    drugbank_xref_path : str or Path, optional
        Path to a multi-namespace DrugBank xref CSV with columns
        ``drugbank_id`` / ``namespace`` / ``external_id``. Supported
        namespaces: ``UNII`` / ``RxCUI`` / ``RXCUI`` / ``ATC`` /
        ``KEGG.DRUG`` / ``ChEBI`` / ``CHEBI`` / ``ChEMBL`` /
        ``CHEMBL.COMPOUND`` / ``CAS``. Generated by
        ``scripts/build_drugbank_xref.py`` from DrugCentral and other
        public sources. Recovers DrugBank-mappable drugs that don't have
        a PubChem CID in the smaller pubchem_to_drugbank.csv.
    uniprot_to_entrez_path : str or Path, optional
        Path to ``uniprot_genesproteins.csv`` (gold standard). When provided,
        gene/protein nodes that lack an NCBIGene crosswalk but have a
        UniProtKB ID are bridged to Entrez gene IDs. Significantly raises
        drug-target coverage joins for Matrix.
    trim_isolated_nodes : bool, default True
        Drop nodes from nodes_df that don't participate in any retained
        edge. Matches the implicit behaviour of the other KG loaders and
        prevents Matrix's huge typed-but-edgeless node tail (~3 M nodes at
        full scale) from dominating topology metrics in nb04 (LCC fraction,
        component count, clustering coefficient). Set False to keep every
        canonical-type node — useful for entity-coverage diagnostics where
        edgeless nodes still count.
    use_cache : bool, default True
        Cache the (kg, nodes_df) result as parquet under ``cache_dir`` so
        repeat loads (e.g. across notebooks 01–08) take ~30 sec instead of
        15-20 minutes. Cache is invalidated when (a) source files change
        mtime, or (b) any of the loader parameters affecting output change.
    cache_dir : str or Path, optional
        Directory for cache files. Defaults to ``Path(nodes_path).parent``
        (i.e. data/matrix/). Cache filenames embed a parameter-hash so
        different configurations don't clobber each other.
    node_chunksize, edge_chunksize : int
        pandas read_csv chunk sizes. Defaults are tuned for ~16 GB RAM hosts.
    verbose : bool
        Print progress every chunk (default True).

    Returns
    -------
    kg : pd.DataFrame   (standardised BioKGSuite edge schema)
    nodes : pd.DataFrame (columns: idx, id, type, name + crosswalk columns)
    """
    # ── Biolink category → BioKGSuite canonical type ──────────────────────────
    _BIOLINK_TO_CANONICAL = {
        # ── Drug / Chemical ───────────────────────────────────────────────
        'SmallMolecule':             'drug',
        'Drug':                      'drug',
        'ChemicalEntity':            'drug',
        'ChemicalMixture':           'drug',
        'MolecularMixture':          'drug',
        'ComplexMolecularMixture':   'drug',
        'ProcessedMaterial':         'drug',
        'ChemicalOrDrugOrTreatment': 'drug',
        'Food':                      'drug',       # food compounds (e.g. vitamins, nutrients)
        # ── Disease / Phenotype ────────────────────────────────────────────
        'Disease':                   'disease',
        'DiseaseOrPhenotypicFeature':'disease',
        'PathologicalProcess':       'disease',    # pathological processes treated as diseases
        'PhenotypicFeature':         'effect/phenotype',
        'BehavioralFeature':         'effect/phenotype',
        # ── Gene / Protein ─────────────────────────────────────────────────
        'Gene':                      'gene/protein',
        'Protein':                   'gene/protein',
        'Polypeptide':               'gene/protein',
        'GeneOrGeneProduct':         'gene/protein',
        'GeneFamily':                'gene/protein',  # gene family nodes (Entrez via equiv IDs)
        'GenomicEntity':             'gene/protein',
        'NucleicAcidEntity':         'gene/protein',
        'Transcript':                'gene/protein',
        'Exon':                      'gene/protein',
        'RNAProduct':                'gene/protein',
        'NoncodingRNAProduct':       'gene/protein',
        'MicroRNA':                  'gene/protein',
        # ── Other biological entities ──────────────────────────────────────
        'BiologicalProcess':         'biological_process',
        'PhysiologicalProcess':      'biological_process',
        'MolecularActivity':         'molecular_function',
        'CellularComponent':         'cellular_component',
        'Pathway':                   'pathway',
        'AnatomicalEntity':          'anatomy',
        'GrossAnatomicalStructure':  'anatomy',
        'Cell':                      'anatomy',
        'CellLine':                  'anatomy',
        'OrganismTaxon':             'organism',
        'SequenceVariant':           'variant',
        'Haplotype':                 'variant',
        'Procedure':                 'procedure',
        'MolecularEntity':           'molecular_entity',
    }

    def _map_biolink(cat: str) -> str:
        """Strip 'biolink:' prefix and map to canonical type.

        Handles both simple strings ('biolink:Disease') and multi-value
        array-like strings produced by numpy repr, e.g.:
            "['biolink:Disease' 'biolink:DiseaseOrPhenotypicFeature' 'biolink:NamedThing']"
        For multi-value strings the first CamelCase token that has a mapping
        in _BIOLINK_TO_CANONICAL is used, so the most-specific declared type
        wins (callers should list the most-specific category first).
        """
        bare = str(cat).replace('biolink:', '').strip()
        # Fast path: clean single-category string
        result = _BIOLINK_TO_CANONICAL.get(bare)
        if result:
            return result
        # Slow path: array-like multi-value string, e.g. "['Disease' 'NamedThing']"
        # Extract capitalised CamelCase tokens and return the first recognised one.
        for token in re.findall(r'[A-Z][a-zA-Z0-9/]*', bare):
            result = _BIOLINK_TO_CANONICAL.get(token)
            if result:
                return result
        return bare.lower()

    # ── Cache check (fast path for repeat loads) ─────────────────────────────
    if use_cache:
        import hashlib
        # keep_categories may be None / 'all' / empty here (means keep-all);
        # normalise to a hashable token for the cache key.
        if (keep_categories is None or
            (isinstance(keep_categories, str) and keep_categories.lower() == 'all') or
            (hasattr(keep_categories, '__len__') and len(keep_categories) == 0)):
            keep_token = '__ALL__'
        else:
            keep_token = sorted(keep_categories)
        param_blob = repr((
            keep_token,
            str(doid_to_mondo_path) if doid_to_mondo_path else '',
            str(mesh_to_doid_path) if mesh_to_doid_path else '',
            str(mondo_sssom_path) if mondo_sssom_path else '',
            str(pubchem_to_drugbank_path) if pubchem_to_drugbank_path else '',
            str(uniprot_to_entrez_path) if uniprot_to_entrez_path else '',
            str(drugbank_xref_path) if drugbank_xref_path else '',
            bool(trim_isolated_nodes),
        )).encode()
        param_hash = hashlib.md5(param_blob).hexdigest()[:8]
        _cache_dir = Path(cache_dir) if cache_dir else Path(nodes_path).parent
        kg_cache    = _cache_dir / f"_cache_matrix_{param_hash}_kg.parquet"
        nodes_cache = _cache_dir / f"_cache_matrix_{param_hash}_nodes.parquet"

        if kg_cache.exists() and nodes_cache.exists():
            cache_mtime  = min(kg_cache.stat().st_mtime, nodes_cache.stat().st_mtime)
            source_mtime = max(Path(nodes_path).stat().st_mtime,
                               Path(edges_path).stat().st_mtime)
            if cache_mtime > source_mtime:
                if verbose:
                    print(f"  matrix: loading from cache "
                          f"{kg_cache.name} / {nodes_cache.name}")
                kg       = pd.read_parquet(kg_cache)
                nodes_df = pd.read_parquet(nodes_cache)
                if verbose:
                    print(f"  matrix: cached load — "
                          f"{len(nodes_df):,} nodes, {len(kg):,} edges")
                return kg, nodes_df
            elif verbose:
                print("  matrix: cache stale (source files newer), rebuilding")

    # ── Pre-compiled regexes for ID canonicalisation ─────────────────────────
    DRUGBANK_RE = re.compile(r'\b(DB\d{5,})\b')
    MONDO_RE    = re.compile(r'MONDO:(\d+)')
    DOID_RE     = re.compile(r'DOID:(\d+)')
    NCBIGENE_RE = re.compile(r'NCBIGene:(\d+)')
    PUBCHEM_RE  = re.compile(r'PUBCHEM\.COMPOUND:(\d+)')
    UNIPROT_RE  = re.compile(r'UniProtKB:([A-Z0-9\-]+)')
    REACTOME_RE = re.compile(r'(R-HSA-\d+)')
    MESH_RE     = re.compile(r'MESH:([CD]\w*)')

    # 'all' / None / empty → keep every category (no node-type filter).
    keep_all = (keep_categories is None or
                (isinstance(keep_categories, str) and keep_categories.lower() == 'all') or
                (hasattr(keep_categories, '__len__') and len(keep_categories) == 0))
    keep_categories = None if keep_all else set(keep_categories)
    if verbose and keep_all:
        print("  matrix: keep_categories=all, retaining every typed node and edge")

    # ── Optional drug / gene crosswalks for canonicalising non-DrugBank /
    #   non-Entrez node IDs (raises gold-standard join rates dramatically) ──
    pubchem_to_drugbank: dict[str, str] = {}
    if pubchem_to_drugbank_path:
        _pc = pd.read_csv(pubchem_to_drugbank_path).dropna(subset=['drugbank_id', 'pubchem_cid'])
        for r in _pc[['pubchem_cid', 'drugbank_id']].itertuples(index=False):
            pubchem_to_drugbank.setdefault(str(r.pubchem_cid).strip(), str(r.drugbank_id).strip())
        if verbose:
            print(f"  matrix: loaded {len(pubchem_to_drugbank):,} PubChem→DrugBank bridges")

    uniprot_to_entrez: dict[str, str] = {}
    if uniprot_to_entrez_path:
        # Two file shapes are in circulation: UniProt's own export, whose
        # columns are 'Entry' and 'GeneID', and the normalised crosswalk this
        # repository builds, whose columns are 'uniprot_id', 'gene_symbol' and
        # 'ensembl_id'. Only the first carries an Entrez id at all. A file with
        # neither shape used to raise a KeyError from inside dropna, which took
        # the whole MATRIX load down over an optional gene-side bridge; it is now
        # reported and skipped, like the DrugBank xref above it.
        _up = pd.read_csv(uniprot_to_entrez_path)
        _entry_col = next((c for c in ('Entry', 'uniprot_id', 'UniProt', 'uniprot')
                           if c in _up.columns), None)
        _gene_col = next((c for c in ('GeneID', 'entrez_id', 'EntrezID', 'gene_id')
                          if c in _up.columns), None)
        if _entry_col and _gene_col:
            _up = _up.dropna(subset=[_entry_col, _gene_col])
            for r in _up[[_entry_col, _gene_col]].itertuples(index=False):
                entry = str(r[0]).strip()
                # GeneID can be multi-value semicolon-separated; take first numeric.
                ids = re.findall(r'\d+', str(r[1]))
                if ids:
                    uniprot_to_entrez.setdefault(entry, ids[0])
            if verbose:
                print(f"  matrix: loaded {len(uniprot_to_entrez):,} UniProt→Entrez bridges")
        elif verbose:
            print(f"  matrix: {Path(uniprot_to_entrez_path).name} carries "
                  f"{list(_up.columns)} and no Entrez gene id, so UniProt→Entrez "
                  "bridging is unavailable; UniProt-only genes keep their original "
                  "CURIE. Drug and disease resolution are unaffected.")

    # Multi-namespace DrugBank xref: {namespace: {external_id: drugbank_id}}.
    # Catches DrugBank drugs that don't have a PubChem CID in our PubChem
    # crosswalk but DO have UNII / RxCUI / ATC / KEGG.DRUG / ChEBI / ChEMBL.
    db_xref: dict[str, dict[str, str]] = {}
    if drugbank_xref_path:
        _xr = pd.read_csv(drugbank_xref_path).dropna(
            subset=['drugbank_id', 'namespace', 'external_id'])
        for r in _xr[['namespace', 'external_id', 'drugbank_id']].itertuples(index=False):
            ns = str(r.namespace).strip().upper()
            db_xref.setdefault(ns, {}).setdefault(
                str(r.external_id).strip(), str(r.drugbank_id).strip())
        if verbose:
            n_total = sum(len(v) for v in db_xref.values())
            print(f"  matrix: loaded {n_total:,} DrugBank xrefs across "
                  f"{len(db_xref)} namespaces ({sorted(db_xref.keys())})")

    # ── Optional disease ID crosswalks (DOID→MONDO and MESH→DOID→MONDO) ──────
    doid_to_mondo: dict[str, str] = {}
    if doid_to_mondo_path:
        _do = pd.read_csv(doid_to_mondo_path).dropna(subset=['mondo_id', 'doid'])
        _do['mondo_num'] = (
            _do['mondo_id'].astype(str)
            .str.replace('MONDO:', '', regex=False).str.lstrip('0').replace('', '0')
        )
        _do['doid_num'] = (
            _do['doid'].astype(str)
            .str.replace('DOID:', '', regex=False).str.lstrip('0').replace('', '0')
        )
        # First-wins: many DOIDs map to multiple MONDOs; the gold-standard pipeline
        # uses the first match for compat with how doid_to_mondo is collapsed in nb01.
        for r in _do[['doid_num', 'mondo_num']].itertuples(index=False):
            doid_to_mondo.setdefault(r.doid_num, r.mondo_num)
        if verbose:
            print(f"  matrix: loaded {len(doid_to_mondo):,} DOID→MONDO bridges")

    mesh_to_mondo: dict[str, str] = {}
    if mesh_to_doid_path and doid_to_mondo:
        # mesh_to_doid.csv contains rows with unquoted commas inside annotation
        # braces (e.g. ``{source="X", source="Y"}``). Skip those bad rows — same
        # behaviour as nb01's gold-standard load.
        _mh = (
            pd.read_csv(mesh_to_doid_path, on_bad_lines='skip')
            .dropna(subset=['mesh_id', 'doid'])
        )
        _mh['mesh_clean'] = (
            _mh['mesh_id'].astype(str)
            .str.replace('MESH:', '', regex=False).str.strip()
        )
        _mh['doid_num'] = (
            _mh['doid'].astype(str)
            .str.replace('DOID:', '', regex=False).str.lstrip('0').replace('', '0')
        )
        for r in _mh[['mesh_clean', 'doid_num']].itertuples(index=False):
            m = doid_to_mondo.get(r.doid_num)
            if m and r.mesh_clean not in mesh_to_mondo:
                mesh_to_mondo[r.mesh_clean] = m
        if verbose:
            print(f"  matrix: loaded {len(mesh_to_mondo):,} MESH→MONDO bridges")

    # MONDO SSSOM crosswalk: external CURIE → MONDO numeric.
    # Recovers UMLS / OMIM / Orphanet / ICD9 / NCIT-only disease nodes.
    sssom_to_mondo: dict[str, str] = {}
    if mondo_sssom_path:
        # SSSOM TSVs start with a YAML metadata header (lines beginning with '#').
        # Use comment='#' so pandas skips it; the table has columns
        # subject_id / predicate_id / object_id / mapping_justification / ...
        _ss = pd.read_csv(mondo_sssom_path, sep='\t', comment='#', dtype=str,
                          low_memory=False).fillna('')
        # Keep only equivalence-class mappings (close & exact). 'broadMatch' / 'narrowMatch'
        # are intentionally excluded — they could collapse non-equivalent diseases.
        equiv_predicates = {'skos:exactMatch', 'skos:closeMatch',
                            'owl:equivalentClass'}
        _ss = _ss[_ss['predicate_id'].isin(equiv_predicates)]
        # Determine which side carries MONDO and pivot accordingly.
        for r in _ss[['subject_id', 'object_id']].itertuples(index=False):
            s, o = r.subject_id, r.object_id
            if s.startswith('MONDO:') and not o.startswith('MONDO:'):
                m = s.replace('MONDO:', '').lstrip('0') or '0'
                sssom_to_mondo.setdefault(o, m)
            elif o.startswith('MONDO:') and not s.startswith('MONDO:'):
                m = o.replace('MONDO:', '').lstrip('0') or '0'
                sssom_to_mondo.setdefault(s, m)
        if verbose:
            print(f"  matrix: loaded {len(sssom_to_mondo):,} SSSOM ext→MONDO bridges")

    # ── Pass 1: stream nodes, filter by category, canonicalise IDs ───────────
    if verbose:
        print(f"  matrix: indexing nodes from {nodes_path} ...")
    node_records: dict[str, tuple] = {}   # CURIE → (idx, canonical_id, type, name, drugbank_id, ncbigene_id, reactome_id, doid_id)
    next_idx = 0
    n_nodes_seen = 0
    n_disease_no_mondo = 0
    n_disease_doid_bridge = 0
    n_disease_mesh_bridge = 0
    n_disease_sssom_bridge = 0

    node_cols = ['id', 'name', 'category', 'equivalent_identifiers',
                 'upstream_data_source']
    for chunk in pd.read_csv(nodes_path, sep='\t', dtype=str,
                             usecols=lambda c: c in node_cols,
                             chunksize=node_chunksize):
        chunk = chunk.fillna('')
        chunk['type'] = chunk['category'].map(_map_biolink)
        if keep_categories is not None:
            chunk = chunk[chunk['type'].isin(keep_categories)]
        n_nodes_seen += len(chunk)

        for row in chunk.itertuples(index=False):
            curie    = row.id
            ctype    = row.type
            name     = row.name
            haystack = curie + ' ' + (row.equivalent_identifiers or '')
            # Per-NODE upstream pipeline (rtxkg2 / robokop / primekg / ...).
            # This is true node-level provenance — distinct from the edge-level
            # primary_knowledge_source carried by edge_source.
            raw_node_upstream = getattr(row, 'upstream_data_source', '') or ''
            node_upstream = '|'.join(re.findall(r"[a-zA-Z][\w-]*", raw_node_upstream))

            # Crosswalk extraction (always recorded for downstream nb usage)
            mb = DRUGBANK_RE.search(haystack)
            mn = NCBIGENE_RE.search(haystack)
            mr = REACTOME_RE.search(haystack)
            md = DOID_RE.search(haystack)
            mm = MONDO_RE.search(haystack)
            drugbank_id = mb.group(1) if mb else ''
            ncbigene_id = mn.group(1) if mn else ''
            reactome_id = mr.group(1) if mr else ''
            doid_id     = md.group(1) if md else ''
            mondo_num   = str(int(mm.group(1))) if mm else ''   # strip leading zeros

            # Canonicalise the primary ID for joinability with gold standards
            if ctype == 'drug':
                if drugbank_id:
                    canonical_id = drugbank_id
                else:
                    # Cascade: PubChem first (most populated), then multi-namespace xref
                    bridged = ''
                    if pubchem_to_drugbank:
                        pc_match = PUBCHEM_RE.search(haystack)
                        if pc_match and pc_match.group(1) in pubchem_to_drugbank:
                            bridged = pubchem_to_drugbank[pc_match.group(1)]
                    if not bridged and db_xref:
                        # Walk every CURIE in haystack, look up by namespace.
                        for tok in re.findall(r'[A-Za-z][\w.]*:[\w.\-]+', haystack):
                            ns, _, ext = tok.partition(':')
                            ns_up = ns.strip().upper()
                            ns_dict = db_xref.get(ns_up)
                            if ns_dict and ext in ns_dict:
                                bridged = ns_dict[ext]
                                break
                    if bridged:
                        drugbank_id = bridged
                        canonical_id = bridged
                    else:
                        canonical_id = curie
            elif ctype == 'disease':
                # Cascade: MONDO direct → DOID bridge → MESH bridge → SSSOM bridge
                bridged = ''
                if mondo_num:
                    bridged = mondo_num
                elif doid_id and doid_id in doid_to_mondo:
                    bridged = doid_to_mondo[doid_id]
                    n_disease_doid_bridge += 1
                else:
                    mesh_match = MESH_RE.search(haystack)
                    mesh_key = mesh_match.group(1) if mesh_match else ''
                    if mesh_key and mesh_key in mesh_to_mondo:
                        bridged = mesh_to_mondo[mesh_key]
                        n_disease_mesh_bridge += 1
                    elif sssom_to_mondo:
                        # Try SSSOM on primary CURIE first, then any CURIE in
                        # equivalent_identifiers. First match wins.
                        sssom_hit = sssom_to_mondo.get(curie, '')
                        if not sssom_hit:
                            for tok in re.findall(r'[A-Za-z][\w.]*:[\w.\-]+', haystack):
                                hit = sssom_to_mondo.get(tok)
                                if hit:
                                    sssom_hit = hit
                                    break
                        if sssom_hit:
                            bridged = sssom_hit
                            n_disease_sssom_bridge += 1
                if bridged:
                    canonical_id = bridged
                elif keep_all:
                    # keep_categories=all: retain disease even without MONDO
                    # bridge — it just won't join MONDO-keyed gold standards.
                    canonical_id = curie
                else:
                    n_disease_no_mondo += 1
                    continue   # disease_id_scheme: mondo requires a MONDO mapping
            elif ctype == 'gene/protein':
                if ncbigene_id:
                    canonical_id = ncbigene_id
                elif uniprot_to_entrez:
                    # UniProt → Entrez: try primary id then equivalent_identifiers
                    up_match = UNIPROT_RE.search(haystack)
                    # Strip isoform suffix (P12345-2 → P12345) for the lookup
                    if up_match:
                        up_acc = up_match.group(1).split('-')[0]
                        if up_acc in uniprot_to_entrez:
                            ncbigene_id = uniprot_to_entrez[up_acc]
                            canonical_id = ncbigene_id
                        else:
                            canonical_id = curie
                    else:
                        canonical_id = curie
                else:
                    canonical_id = curie
            elif ctype == 'pathway':
                canonical_id = reactome_id or curie
            else:   # 'effect/phenotype' and any other kept category
                canonical_id = curie

            node_records[curie] = (
                next_idx, canonical_id, ctype, name,
                drugbank_id, ncbigene_id, reactome_id, doid_id,
                node_upstream,
            )
            next_idx += 1

        if verbose:
            print(f"    ... {n_nodes_seen:,} typed candidates seen, "
                  f"{len(node_records):,} kept")
    if verbose:
        if n_disease_doid_bridge or n_disease_mesh_bridge or n_disease_sssom_bridge:
            print(f"  matrix: bridged {n_disease_doid_bridge:,} disease nodes via DOID→MONDO, "
                  f"{n_disease_mesh_bridge:,} via MESH→DOID→MONDO, "
                  f"{n_disease_sssom_bridge:,} via MONDO SSSOM")
        if n_disease_no_mondo:
            print(f"  matrix: dropped {n_disease_no_mondo:,} disease nodes "
                  f"with no MONDO crosswalk (disease_id_scheme=mondo)")

    # ── Pass 2: stream edges, filter to those joining two kept nodes ─────────
    if verbose:
        print(f"  matrix: streaming edges from {edges_path} ...")
    node_id_set = set(node_records.keys())
    kg_chunks: list[pd.DataFrame] = []
    n_edges_seen = 0
    n_edges_kept = 0
    # Biolink knowledge_level → confidence-score ordinal mapping (used by nb03's
    # uncertainty quantification). Higher = more trustworthy.
    KL_TO_SCORE = {
        'knowledge_assertion':     1.00,
        'logical_entailment':      0.90,
        'statistical_association': 0.60,
        'observation':             0.50,
        'prediction':              0.30,
        'not_provided':            float('nan'),
        '':                        float('nan'),
    }
    edge_cols = ['subject', 'predicate', 'object',
                 'primary_knowledge_source', 'aggregator_knowledge_source',
                 'upstream_data_source',
                 'knowledge_level', 'agent_type', 'num_references']
    for chunk in pd.read_csv(edges_path, sep='\t', dtype=str,
                             usecols=lambda c: c in edge_cols,
                             chunksize=edge_chunksize):
        chunk = chunk.fillna('')
        n_edges_seen += len(chunk)
        mask = (chunk['subject'].isin(node_id_set) &
                chunk['object'].isin(node_id_set))
        if not mask.any():
            if verbose:
                print(f"    ... {n_edges_seen:,} edges seen, {n_edges_kept:,} kept")
            continue
        sub = chunk.loc[mask].copy()

        # Vectorised lookup of canonical fields by mapping CURIE → record
        sub['_s'] = sub['subject'].map(node_records)
        sub['_t'] = sub['object'].map(node_records)
        out = pd.DataFrame({
            'relation':  sub['predicate'].str.replace('biolink:', '', regex=False),
            'x_index':   sub['_s'].map(lambda r: r[0]),
            'x_id':      sub['_s'].map(lambda r: r[1]),
            'x_type':    sub['_s'].map(lambda r: r[2]),
            'x_name':    sub['_s'].map(lambda r: r[3]),
            # Per-NODE upstream pipeline (slot 8 in node_records). True
            # node-level provenance — analogous to PrimeKG's x_source / y_source.
            'x_source':  sub['_s'].map(lambda r: r[8]),
            'y_index':   sub['_t'].map(lambda r: r[0]),
            'y_id':      sub['_t'].map(lambda r: r[1]),
            'y_type':    sub['_t'].map(lambda r: r[2]),
            'y_name':    sub['_t'].map(lambda r: r[3]),
            'y_source':  sub['_t'].map(lambda r: r[8]),
        })
        # ── Edge-level provenance (per-edge, NOT per-endpoint) ──────────────
        if 'primary_knowledge_source' in sub.columns:
            # Asserting database for the relationship — Biolink infores: CURIE.
            # Lives in its own column so x_source/y_source remain pure node-level.
            out['edge_source'] = sub['primary_knowledge_source'].values
        if 'aggregator_knowledge_source' in sub.columns:
            out['aggregator_source'] = sub['aggregator_knowledge_source'].values
        if 'upstream_data_source' in sub.columns:
            # Edge-level upstream (which Translator pipeline ingested the edge).
            # Distinct from node-level x_source/y_source (which pipeline ingested
            # each endpoint).
            out['edge_upstream'] = sub['upstream_data_source'].values

        # Provenance / uncertainty columns (Biolink-native, used by nb03).
        if 'knowledge_level' in sub.columns:
            out['knowledge_level'] = sub['knowledge_level'].values
            # Numeric confidence_score derived from knowledge_level so nb03's
            # keyword scan ('score') picks it up without Matrix-specific code.
            out['confidence_score'] = (
                sub['knowledge_level'].map(KL_TO_SCORE).astype(float).values
            )
        if 'agent_type' in sub.columns:
            out['agent_type'] = sub['agent_type'].values
        if 'num_references' in sub.columns:
            # Numeric — count of supporting publications.
            out['num_references'] = pd.to_numeric(
                sub['num_references'], errors='coerce'
            ).values

        kg_chunks.append(out)
        n_edges_kept += len(out)
        if verbose:
            print(f"    ... {n_edges_seen:,} edges seen, {n_edges_kept:,} kept")

    if not kg_chunks:
        kg = pd.DataFrame(columns=[
            'relation','x_index','x_id','x_type','x_name',
            'y_index','y_id','y_type','y_name'])
    else:
        kg = pd.concat(kg_chunks, ignore_index=True)
        del kg_chunks
    kg['x_index'] = kg['x_index'].astype(int)
    kg['y_index'] = kg['y_index'].astype(int)

    # ── Build node DataFrame from the dict ───────────────────────────────────
    nodes_df = pd.DataFrame.from_records(
        list(node_records.values()),
        columns=['idx', 'id', 'type', 'name',
                 'drugbank_id', 'ncbigene_id', 'reactome_id', 'doid_id',
                 'upstream_pipeline'],
    )
    # Backward-compat with older notebook code that referenced drugbank_ids_all.
    nodes_df['drugbank_ids_all'] = nodes_df['drugbank_id']

    # ── Trim isolated nodes (default) ────────────────────────────────────────
    # The other BioKGSuite loaders emit only nodes participating in edges —
    # so for apples-to-apples topology / clustering / LCC metrics, drop Matrix
    # nodes that ended up with zero edges in the retained edge set. Keeps
    # x_index / y_index integer values intact (no re-indexing needed).
    if trim_isolated_nodes:
        n_before = len(nodes_df)
        edge_node_set = set(kg['x_index'].astype(int)) | set(kg['y_index'].astype(int))
        nodes_df = nodes_df[nodes_df['idx'].isin(edge_node_set)].reset_index(drop=True)
        n_dropped = n_before - len(nodes_df)
        if verbose and n_dropped:
            print(f"  matrix: trimmed {n_dropped:,} isolated nodes "
                  f"(no edges in kept-edge subset)")

    if verbose:
        print(f"  matrix: kept {len(nodes_df):,} nodes, {len(kg):,} edges")

    # ── Write cache for next time ────────────────────────────────────────────
    if use_cache:
        try:
            _cache_dir.mkdir(parents=True, exist_ok=True)
            kg.to_parquet(kg_cache, index=False, compression='snappy')
            nodes_df.to_parquet(nodes_cache, index=False, compression='snappy')
            if verbose:
                size_mb = (kg_cache.stat().st_size + nodes_cache.stat().st_size) / 1e6
                print(f"  matrix: cached to {_cache_dir.name}/ "
                      f"({size_mb:.0f} MB total)")
        except Exception as e:
            if verbose:
                print(f"  matrix: cache write failed ({e}); "
                      f"next load will re-stream from source")

    return kg, nodes_df
