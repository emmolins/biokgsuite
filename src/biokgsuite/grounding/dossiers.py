"""Knowledge-graph evidence blocks (dossiers) for the grounding experiment.

Entity-type policy
------------------
Instead of hand-listing each KG's relation names (DRKG alone has around 26
drug-to-gene relations), edges are kept by the ENTITY TYPE of the neighbour,
resolved per KG from the config entity_types. From a drug node, gene
neighbours render as "targets" and pathway neighbours as "modulates pathway";
from a disease node, gene as "associated gene", pathway as "associated
pathway" and disease as "related disease". The drug-disease bridge is excluded
structurally (a disease neighbour of a drug, or vice versa, is never rendered)
so every KG's leak relation is handled without naming it. The only per-KG list
is PK_EXCLUDE: ADME relations whose gene neighbours are not mechanistic targets.
"""
from collections import defaultdict
from pathlib import Path
import re

import pandas as pd

from biokgsuite.graph.graph_utils import build_lookup_maps
from biokgsuite.io import load_kg

PK_EXCLUDE = {
    'primekg':    set(),
    'drkg':       {'DRUGBANK::enzyme::Compound:Gene', 'DRUGBANK::carrier::Compound:Gene',
                   'DRUGBANK::transporter::Compound:Gene'},
    'biokg':      {'DRUG_ENZYME', 'DRUG_TRANSPORTER', 'DRUG_CARRIER'},
    'openbilink': set(),
    'matrix':     set(),
    'hetionet':   set(),
}
DRUG_PHRASE = {'gene': 'targets', 'pathway': 'modulates pathway'}
DISEASE_PHRASE = {'gene': 'is associated with gene', 'pathway': 'is associated with pathway',
                  'disease': 'is related to'}

# Canonical config roles, in the order the dossier policy refers to them.
_ROLE_OF_CANON = {'Drug': 'drug', 'Disease': 'disease', 'Gene/Protein': 'gene', 'Pathway': 'pathway'}

_BARE_MESH = re.compile(r'^[DC]\d+$')

# A gene name that is only digits and separators is an unmapped Entrez ID (a
# non-coding or predicted locus absent from the symbol map). It is noise to
# the model, so such gene neighbours are dropped from dossiers.
_UNREADABLE_GENE = re.compile(r'^[\d;,/|\s]+$')


def gene_verb(role, rels):
    """Interaction or regulation verb for a drug/disease-to-gene edge, read
    from the KG's relation name(s). Falls back to the generic phrase when the
    KG encodes no direction (PrimeKG 'drug_protein', BioKG 'DRUG_TARGET')."""
    s = ' '.join(str(r).lower() for r in rels)
    if role == 'drug':
        if any(k in s for k in ('inhibitor', 'antagonist', 'blocker', 'suppressor', '::n::', 'inhibits')):
            return 'inhibits'
        if any(k in s for k in ('agonist', 'activator', '::a+::', 'activates')):
            return 'activates'
        if any(k in s for k in ('cdg', 'downregulat', '::e-::')):
            return 'downregulates'
        if any(k in s for k in ('cug', 'upregulat', '::e+::', 'inducer', 'induces')):
            return 'upregulates'
        if any(k in s for k in ('binder', 'binds', 'cbg', '::b::')):
            return 'binds'
        if 'modulator' in s:
            return 'modulates'
        return DRUG_PHRASE['gene']
    if any(k in s for k in ('dug', 'upregulat', '::e+::')):
        return 'upregulates'
    if any(k in s for k in ('ddg', 'downregulat', '::e-::')):
        return 'downregulates'
    return DISEASE_PHRASE['gene']


def _id_variants(cands):
    """Every spelling a KG might use for each identifier: colon or underscore
    namespace separator, either namespace case, the bare local part, and the
    local part without leading zeros."""
    out = set(cands)
    for c in list(cands):
        if not isinstance(c, str):
            continue
        if ':' in c:
            ns, _, rest = c.partition(':')
            out.update({c.replace(':', '_', 1), f'{ns.upper()}:{rest}', f'{ns.lower()}:{rest}', rest})
            try:
                stripped = str(int(rest))
                if stripped != rest:
                    out.add(stripped)
            except (ValueError, TypeError):
                pass
        if '_' in c and ':' not in c:
            out.add(c.replace('_', ':', 1))
    return out


def _strip_prefix(value, prefix):
    """Numeric part of a CURIE without leading zeros ('0' when all zeros)."""
    return value.replace(prefix, '', 1).lstrip('0') or '0'


def load_crosswalks(base):
    """Disease and drug identifier bridges from data/gold_standards/bridges.

    do_diseases.csv, mesh_to_doid.csv and mondo.sssom.tsv are required.
    drugbank_xref.csv is optional: without it ``drug_equiv`` is empty and
    drugs resolve by their DrugBank accession alone.
    """
    bridges = Path(base) / 'data' / 'gold_standards' / 'bridges'
    do = pd.read_csv(bridges / 'do_diseases.csv')
    mesh_doid = pd.read_csv(bridges / 'mesh_to_doid.csv', on_bad_lines='skip')
    sssom_path = bridges / 'mondo.sssom.tsv'
    if not sssom_path.exists():
        raise FileNotFoundError(
            f'{sssom_path} is missing; run scripts/download_mondo_sssom.sh to fetch it. '
            'Without it cross-vocabulary disease equivalences would be silently dropped.')
    sssom = pd.read_csv(sssom_path, sep='\t', comment='#')
    sssom = sssom[sssom['predicate_id'] == 'skos:exactMatch']
    xref_path = bridges / 'drugbank_xref.csv'
    xref = pd.read_csv(xref_path) if xref_path.exists() else None

    do = do.dropna(subset=['mondo_id', 'doid']).copy()
    do['m'] = do['mondo_id'].map(lambda v: _strip_prefix(str(v), 'MONDO:'))
    do['d'] = do['doid'].map(lambda v: _strip_prefix(str(v), 'DOID:'))
    mondo_to_doid_all, doid_to_mondo_all = defaultdict(set), defaultdict(set)
    for m, d in zip(do['m'], do['d']):
        mondo_to_doid_all[m].add(d)
        doid_to_mondo_all[d].add(m)

    # Some MeSH cells carry a trailing "{...}" annotation; strip it. When one
    # MeSH term maps to several DOIDs, prefer DO's own mapping over MONDO's
    # and either over an unlabelled source.
    mesh_doid['mesh_clean'] = mesh_doid['mesh_id'].str.replace(r'\s*\{.*\}', '', regex=True).str.strip()
    if 'source' in mesh_doid.columns:
        src_rank = mesh_doid['source'].map({'DO': 0, 'MONDO': 1}).fillna(2)
    else:
        src_rank = 2
    mh = (mesh_doid.assign(src_rank=src_rank)
                   .sort_values('src_rank', kind='stable').drop_duplicates('mesh_clean', keep='first'))
    mesh_to_doid = {r['mesh_clean']: _strip_prefix(r['doid'], 'DOID:')
                    for _, r in mh.iterrows() if pd.notna(r['doid'])}
    doid_to_mesh = defaultdict(set)
    for m, d in mesh_to_doid.items():
        if m.startswith('MESH:'):
            doid_to_mesh[d].add(m)

    drug_equiv = defaultdict(set)
    if xref is not None:
        for db, ns, eid in zip(xref['drugbank_id'].astype(str), xref['namespace'].astype(str),
                               xref['external_id'].astype(str)):
            drug_equiv[db].update({f'{ns}:{eid}', eid})
            if ns == 'PUBCHEM.COMPOUND':
                drug_equiv[db].update({f'CID:{eid}', f'pubchem.compound:{eid}'})

    # Exact matches are symmetric and, to one hop, transitive: an ID is
    # equivalent to its direct matches and to their matches.
    eq = defaultdict(set)
    for s, o in zip(sssom['subject_id'], sssom['object_id']):
        eq[s].add(o)
        eq[o].add(s)
    disease_equiv = defaultdict(set)
    for k, vs in eq.items():
        s = set(vs)
        for v in vs:
            s |= eq.get(v, set())
        s.discard(k)
        disease_equiv[k] = s

    return {'mondo_to_doid_all': mondo_to_doid_all, 'doid_to_mondo_all': doid_to_mondo_all,
            'mesh_to_doid': mesh_to_doid, 'doid_to_mesh': doid_to_mesh,
            'drug_equiv': drug_equiv, 'disease_equiv': disease_equiv}


def load_entrez_to_symbol(base):
    """Entrez gene id to HGNC symbol, from data/gold_standards/bridges/entrez_to_symbol.tsv
    (built by tools/build_entrez_symbol_map.py from the vendored HGNC set).

    Required: DRKG and OpenBioLink key genes by Entrez id, and a dossier that
    reads "targets 3716" tells the model nothing, so an absent map would
    silently empty every dossier for those graphs.
    """
    path = Path(base) / 'data' / 'gold_standards' / 'bridges' / 'entrez_to_symbol.tsv'
    if not path.exists():
        raise FileNotFoundError(f'{path} is missing; run tools/build_entrez_symbol_map.py')
    out = {}
    for line in path.read_text().splitlines():
        p = line.split('\t')
        if len(p) == 2:
            out[p[0]] = p[1]
    return out


def make_resolvers(config, kg_name, xw):
    """(to_kg_drug_ids, to_kg_disease_ids): each maps a gold-standard ID to the
    set of identifier spellings under which the KG might key that entity.

    Disease IDs are first normalised to a bare MONDO number, then re-expressed
    in the KG's ``disease_id_scheme`` (mondo, doid, doid_mesh or mesh) and
    widened by SSSOM exact matches.
    """
    drug_equiv = xw['drug_equiv']
    mondo_to_doid_all = xw['mondo_to_doid_all']
    doid_to_mondo_all = xw['doid_to_mondo_all']
    mesh_to_doid = xw['mesh_to_doid']
    doid_to_mesh = xw['doid_to_mesh']
    disease_equiv = xw['disease_equiv']
    scheme = config['knowledge_graphs'][kg_name]['disease_id_scheme']

    def to_kg_drug_ids(db_id):
        s = str(db_id).strip()
        return _id_variants({s} | drug_equiv.get(s, set()))

    def to_kg_disease_ids(disease_id):
        s = str(disease_id).strip()
        if s.startswith('MONDO:'):
            mondo = _strip_prefix(s, 'MONDO:')
        elif s.startswith('DOID:'):
            mondo = sorted(doid_to_mondo_all.get(_strip_prefix(s, 'DOID:'), ['?']))[0]
        elif _BARE_MESH.match(s):
            d = mesh_to_doid.get(f'MESH:{s}')
            mondo = sorted(doid_to_mondo_all.get(d, ['?']))[0] if d else '?'
        else:
            mondo = s.lstrip('0') or '0'
        cands = set()
        if scheme == 'mondo':
            cands.update({mondo, mondo.zfill(7), f'MONDO:{mondo}'})
        elif scheme == 'doid':
            for d in mondo_to_doid_all.get(mondo, set()):
                cands.update({d, f'DOID:{d}'})
        elif scheme == 'doid_mesh':
            for d in mondo_to_doid_all.get(mondo, set()):
                cands.update({d, f'DOID:{d}'})
                cands |= doid_to_mesh.get(d, set())
        elif scheme == 'mesh':
            for d in mondo_to_doid_all.get(mondo, set()):
                for m in doid_to_mesh.get(d, set()):
                    cands.update({m, m.split(':', 1)[-1]})
        for c in list(cands) + [f'MONDO:{mondo}']:
            cands |= disease_equiv.get(c, set())
        return _id_variants(cands)

    return to_kg_drug_ids, to_kg_disease_ids


def build_edge_index(kg_df):
    """Undirected adjacency: node index -> [(neighbour index, relation), ...]."""
    ix = defaultdict(list)
    for x, y, r in zip(kg_df['x_index'].values, kg_df['y_index'].values, kg_df['relation'].values):
        u, v, rs = int(x), int(y), str(r)
        ix[u].append((v, rs))
        ix[v].append((u, rs))
    return ix


def _entrez_syms(raw_id, entrez_to_symbol):
    """Symbol(s) for a gene-node ID; handles ';'-separated multi-gene nodes
    ('Gene::3716;3717' -> 'JAK1;JAK2'). None when nothing maps."""
    tail = str(raw_id).split('::')[-1].split(':')[-1]
    syms = [entrez_to_symbol[x] for x in re.split(r'[;,/|]', tail)
            if x.isdigit() and x in entrez_to_symbol]
    return ';'.join(dict.fromkeys(syms)) if syms else None


class KGDossiers:
    """Drug and disease evidence blocks from one knowledge graph.

    The graph is loaded once at construction; the resolvers, identifier index
    and adjacency are built then and reused for every query, so a cap sweep
    or a rerun of one model never re-reads the graph.

    Parameters
    ----------
    cap : int
        Neighbours rendered per phrase group, lowest degree first, so hub
        genes and catch-all pathways do not crowd out specific ones.
    anonymize_genes : bool
        Replace gene names with per-disease codes ("Gene-1"), for a control
        arm where the model cannot recognise the targets by name.
    bridge_mode : str
        'full' renders the disease profile; 'mechanism_only' leaves it empty
        so the model sees drug mechanisms but no disease-side evidence.
    entrez_to_symbol : dict or None
        Entrez ID -> HGNC symbol; renders numeric gene IDs (DRKG 'Gene::3716',
        OpenBioLink 'NCBIGENE:3716') as symbols the model can reason over.
    """

    def __init__(self, kg_name, config, base, *, cap, anonymize_genes=False,
                 bridge_mode='full', entrez_to_symbol=None):
        self.kg_name = kg_name
        self.cap = cap
        self.anonymize_genes = anonymize_genes
        self.bridge_mode = bridge_mode
        self.pk_exclude = PK_EXCLUDE.get(kg_name, set())
        spec = config['knowledge_graphs'][kg_name]
        et = spec['entity_types']

        kg_df, nodes_df = load_kg(kg_name, config)
        self.name_map = build_lookup_maps(nodes_df)['node_name_map']
        self.idx_to_type = dict(zip(nodes_df['idx'].astype(int), nodes_df['type'].astype(str)))

        # entity_types values are tuples after load_config; every label in a
        # tuple plays the same role (BioKG: Disease and GeneticDisorder).
        self.drug_types = tuple(et.get('Drug', ()))
        self.disease_types = tuple(et.get('Disease', ()))
        self.gene_types = tuple(et.get('Gene/Protein', ()))
        self.canon = {}
        for canon_label, role in _ROLE_OF_CANON.items():
            for kg_t in et.get(canon_label, ()):
                self.canon[str(kg_t)] = role

        self.id_to_idx = self._index_ids(nodes_df, spec.get('disease_id_separator'))
        if entrez_to_symbol and self.gene_types:
            self._symbolise_genes(nodes_df, entrez_to_symbol)

        self.to_drug_ids, self.to_disease_ids = make_resolvers(config, kg_name, load_crosswalks(base))
        self.edge_index = build_edge_index(kg_df)
        self.degree = {u: len(adj) for u, adj in self.edge_index.items()}
        self._gene_code_maps = defaultdict(dict)

    def _index_ids(self, nodes_df, sep):
        """Alias key -> [(idx, type)]: the raw ID, its local part, the parts of
        a compound disease ID, and any DrugBank columns the loader carried."""
        id_to_idx = defaultdict(list)
        ids = nodes_df['id'].astype(str).tolist()
        idxs = nodes_df['idx'].astype(int).tolist()
        # Some KGs (OpenBioLink, MATRIX) key drug nodes by PubChem but carry a
        # DrugBank column; index those as alias keys so DrugBank-based gold IDs resolve.
        db = nodes_df['drugbank_id'].astype(str).tolist() if 'drugbank_id' in nodes_df.columns else None
        dball = (nodes_df['drugbank_ids_all'].astype(str).tolist()
                 if 'drugbank_ids_all' in nodes_df.columns else None)
        for i, (raw, idx) in enumerate(zip(ids, idxs)):
            nt = self.idx_to_type.get(idx)
            keys = {raw, raw.split('::', 1)[-1] if '::' in raw else raw}
            if sep and nt in self.disease_types and sep in raw:
                keys.update(p for p in raw.split(sep) if p)
            if db is not None:
                v = db[i]
                if v and v != 'nan':
                    keys.add(v)
            if dball is not None:
                v = dball[i]
                if v and v != 'nan':
                    keys.update(t.strip() for t in re.split(r'[|,;]', v) if t.strip())
            for k in keys:
                id_to_idx[k].append((idx, nt))
        return id_to_idx

    def _symbolise_genes(self, nodes_df, entrez_to_symbol):
        for raw, idx in zip(nodes_df['id'].astype(str), nodes_df['idx'].astype(int)):
            if self.idx_to_type.get(idx) in self.gene_types:
                sym = _entrez_syms(raw, entrez_to_symbol)
                if sym:
                    self.name_map[idx] = sym

    def _resolve(self, cands, want_types):
        # Candidates are a set; iterating it in hash order made the chosen node
        # depend on the process when several spellings hit different nodes.
        for cand in sorted(cands):
            for idx, t in self.id_to_idx.get(cand, []):
                if t in want_types:
                    return idx
        return None

    def resolves(self, drug_id):
        return self._resolve(self.to_drug_ids(drug_id), self.drug_types) is not None

    def resolves_disease(self, disease_id):
        return self._resolve(self.to_disease_ids(disease_id), self.disease_types) is not None

    def _gene_coder(self, disease_id):
        # Codes are assigned per disease so the same gene gets the same code
        # across every candidate of one query, but reveals nothing across queries.
        m = self._gene_code_maps[disease_id]

        def code(name):
            if name not in m:
                m[name] = f"Gene-{len(m) + 1}"
            return m[name]
        return code

    def _triples(self, node_idx, subject, disease_id, role, exclude_idx=()):
        """One "<subject> <phrase> <object>" per kept neighbour. Gene edges
        render the KG's interaction or regulation modality where available,
        else the generic phrase; pathway and disease use fixed phrases.
        Drug-disease neighbours are excluded structurally (leakage)."""
        allowed = ('gene', 'pathway') if role == 'drug' else ('gene', 'pathway', 'disease')
        nb_rels, nb_role_of = defaultdict(set), {}
        for v, r in self.edge_index.get(node_idx, []):
            nb = self.canon.get(self.idx_to_type.get(v))
            if nb is None or nb not in allowed:
                continue
            # Enzymes, carriers and transporters are ADME proteins, not targets.
            if role == 'drug' and nb == 'gene' and str(r) in self.pk_exclude:
                continue
            nb_rels[v].add(r)
            nb_role_of[v] = nb
        groups = defaultdict(list)
        phrases = DRUG_PHRASE if role == 'drug' else DISEASE_PHRASE
        for v, rels in nb_rels.items():
            nb = nb_role_of[v]
            phrase = gene_verb(role, rels) if nb == 'gene' else phrases[nb]
            groups[phrase].append(v)
        coder = self._gene_coder(disease_id) if self.anonymize_genes else None
        out = []
        for phrase in sorted(groups):
            seen, seen_names, nbrs = set(exclude_idx), set(), []
            for v in sorted(groups[phrase], key=lambda n: self.degree.get(n, 0)):
                if v in seen:
                    continue
                seen.add(v)
                name = str(self.name_map.get(v, v))
                if self.idx_to_type.get(v) in self.gene_types and _UNREADABLE_GENE.match(name):
                    continue
                # Some graphs carry several nodes for one gene (BioKG lists a
                # protein per isoform); one line per rendered name is what the
                # model should see, and the cap counts names, not nodes.
                if name in seen_names:
                    continue
                seen_names.add(name)
                nbrs.append(v)
                if len(nbrs) >= self.cap:
                    break
            for v in nbrs:
                nm = self.name_map.get(v, str(v))
                if coder and self.idx_to_type.get(v) in self.gene_types:
                    nm = coder(nm)
                out.append(f"{subject} {phrase} {nm}")
        return "; ".join(out)

    def drug_dossier(self, drug_id, drug_name, disease_id):
        """Evidence block for one candidate, or None when the drug is not in
        the graph. The target disease is masked even though disease neighbours
        are never rendered for a drug, as a second guard against leakage."""
        d_idx = self._resolve(self.to_drug_ids(drug_id), self.drug_types)
        if d_idx is None:
            return None
        x_idx = self._resolve(self.to_disease_ids(disease_id), self.disease_types)
        excl = {x_idx} if x_idx is not None else set()
        return self._triples(d_idx, drug_name, disease_id, 'drug', excl) or None

    def disease_profile(self, disease_id, disease_name=None):
        """Disease-side evidence shown once in the prompt header; '' when the
        disease is not in the graph or bridge_mode is 'mechanism_only'."""
        if self.bridge_mode == 'mechanism_only':
            return ''
        x_idx = self._resolve(self.to_disease_ids(disease_id), self.disease_types)
        if x_idx is None:
            return ''
        subject = disease_name or self.name_map.get(x_idx, 'disease')
        return self._triples(x_idx, subject, disease_id, 'disease')
