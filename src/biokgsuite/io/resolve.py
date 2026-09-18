"""One identifier-resolution pipeline, applied identically to every KG.

Why this exists
---------------
Coverage was previously computed after a per-KG normalisation pipeline whose
depth varied by graph: MATRIX was wired to six crosswalks in config, BioKG to
one deliberately lossy bridge, Hetionet to none. Coverage then measured
bridging effort as much as graph content.

This module resolves by SNIFFING THE NAMESPACE OF THE IDENTIFIER, not by
looking up which KG it came from, so every graph gets every bridge. It also
returns, for every call, how many inputs resolved and which namespaces failed,
so resolvability can be reported alongside coverage rather than folded into it.

Canonical target spaces — these match the reference sets in
data/reference_sets/, which are built from Open Targets and MONDO:

    disease  -> bare MONDO numeric      (e.g. '7254')
    drug     -> ChEMBL id               (e.g. 'CHEMBL941')
    gene     -> Ensembl gene id         (e.g. 'ENSG00000097007')
    pathway  -> Reactome stable id      (e.g. 'R-HSA-109582')

Mappings are many-to-many throughout, for every KG.
"""
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
import re
import pandas as pd

_PREFIX_RE = re.compile(r'^\s*([A-Za-z0-9._]+)\s*[:_]\s*(.+?)\s*$')


def split_curie(raw: str):
    """('DOID:1612') -> ('DOID','1612'); ('Compound::DB00001') -> (None,'DB00001')."""
    s = str(raw).strip()
    if '::' in s:
        s = s.split('::')[-1]
    s = re.sub(r'\s*\{.*\}\s*$', '', s)
    m = _PREFIX_RE.match(s)
    return (m.group(1).upper(), m.group(2)) if m else (None, s)


def _num(x):
    x = str(x).strip()
    return x.lstrip('0') or '0'


@dataclass
class ResolutionStats:
    n_input: int = 0
    n_resolved: int = 0
    n_output: int = 0
    by_namespace: dict = field(default_factory=lambda: defaultdict(int))
    failed_by_namespace: dict = field(default_factory=lambda: defaultdict(int))

    @property
    def rate(self):
        return self.n_resolved / self.n_input if self.n_input else float('nan')

    def as_row(self, **extra):
        top = sorted(self.failed_by_namespace.items(), key=lambda kv: -kv[1])[:3]
        row = {'raw': self.n_input, 'resolved': self.n_resolved,
               'resolvability': round(self.rate, 4) if self.n_input else float('nan'),
               'canonical_ids': self.n_output,
               'top_unresolved': '; '.join(f'{k}:{v:,}' for k, v in top)}
        row.update(extra)
        return row


class IdResolver:
    """Loads every crosswalk once; resolves any KG's ids the same way."""

    # Registries for entity types outside drug, disease, gene and pathway.
    # Off by default: nb01 scores only the types it holds reference sets for,
    # and nb02 reports registry agreement only for drug, disease and gene, so
    # loading these costs ~52 MB and several seconds to produce numbers nothing
    # reads. The files stay on disk, dated and checksummed, and one flag brings
    # them back if a later dimension needs them.
    AUX_REGISTRIES = ('Biological Process', 'Molecular Function',
                      'Cellular Component', 'Anatomy', 'Side Effect')

    def __init__(self, crosswalks, verbose=True, aux_registries=False):
        g = Path(crosswalks)
        self._v = verbose
        self._aux = aux_registries
        self.x_to_mondo = defaultdict(set)      # (NS, id) -> MONDO numeric
        self.drug_to_chembl = defaultdict(set)  # (NS, id) -> CHEMBLxxx
        self.gene_to_ensg = defaultdict(set)    # (NS, id) -> ENSGxxx
        self.ensg_all, self.chembl_all = set(), set()
        self.go_aspect = {}                     # GO id -> canonical type
        self.uberon = set()                     # UBERON ids
        self.obsolete = set()                   # retired GO / UBERON ids
        self.wikipathways = set()               # valid WP accessions
        self.reactome = set()                   # Reactome PATHWAY stable ids
        self.umls_se = set()                    # UMLS CUIs naming a side effect
        self._load_disease(g)
        self._load_reactome(g)
        self._load_gene(g)
        self._load_drug(g)
        if aux_registries:
            self._load_obo_types(g)
        if verbose:
            if aux_registries:
                print(f'IdResolver: GO {len(self.go_aspect):,} | UBERON {len(self.uberon):,} | '
                      f'WikiPathways {len(self.wikipathways):,} | MedDRA {len(self.umls_se):,}')
            print(f'IdResolver: disease {len(self.x_to_mondo):,} ids -> MONDO | '
                  f'drug {len(self.drug_to_chembl):,} -> ChEMBL | '
                  f'gene {len(self.gene_to_ensg):,} -> Ensembl')

    # ------------------------------------------------------------- loading
    def _load_reactome(self, g):
        """Reactome stable ids that name a PATHWAY, from Open Targets' hierarchy.

        Needed because a Reactome stable id does not say what kind of object it
        names: R-HSA-nnnnnn is used for pathways, reactions, complexes and
        physical entities alike (reactome.org/dev/graph-database). Testing the
        `R-HSA-` prefix therefore accepts any Reactome object as a pathway; on
        BioKG that admitted 12,687 protein complexes. Only membership of an
        actual pathway list can separate them.
        """
        d = g / 'ot_reactome'
        files = sorted(d.glob('*.parquet')) if d.exists() else []
        if not files:
            if self._v:
                print('  ot_reactome/ missing — pathway validity falls back to '
                      'the R-HSA- prefix, which accepts non-pathway Reactome ids')
            return
        for f in files:
            self.reactome |= set(pd.read_parquet(f, columns=['id'])['id'].astype(str))

    def _load_disease(self, g):
        p = g / 'mondo.sssom.tsv'
        if not p.exists():
            if self._v:
                print('  mondo.sssom.tsv missing — disease resolution will be poor')
            return
        ss = pd.read_csv(p, sep='\t', comment='#', low_memory=False)
        ss = ss[ss['predicate_id'] == 'skos:exactMatch']
        for subj, obj in zip(ss['subject_id'], ss['object_id']):
            if not isinstance(subj, str) or not subj.startswith('MONDO:'):
                continue
            m = _num(subj.split(':', 1)[1])
            if isinstance(obj, str) and ':' in obj:
                ns, lid = obj.split(':', 1)
                ns = ns.upper()
                self.x_to_mondo[(ns, lid)].add(m)
                if ns in ('DOID', 'MESH', 'OMIM', 'UMLS', 'NCIT', 'EFO', 'ORPHANET'):
                    self.x_to_mondo[(ns, _num(lid))].add(m)
                # MIM and OMIM name one registry. MONDO writes OMIM; BioKG and
                # others write MIM. Indexed under both so the prefix spelling
                # does not decide whether an identifier resolves.
                if ns == 'OMIM':
                    self.x_to_mondo[('MIM', lid)].add(m)
                    self.x_to_mondo[('MIM', _num(lid))].add(m)

    def _load_gene(self, g):
        p = g / 'hgnc_complete_set.txt'
        if p.exists():
            h = pd.read_csv(p, sep='\t', low_memory=False,
                            usecols=lambda c: c in ('hgnc_id', 'symbol', 'entrez_id',
                                                    'ensembl_gene_id', 'uniprot_ids',
                                                    'alias_symbol', 'prev_symbol'))
            for _, r in h.iterrows():
                ensg = r.get('ensembl_gene_id')
                if not isinstance(ensg, str) or not ensg.startswith('ENSG'):
                    continue
                self.ensg_all.add(ensg)
                ez = r.get('entrez_id')
                if pd.notna(ez):
                    self.gene_to_ensg[('ENTREZ', str(int(float(ez))))].add(ensg)
                sym = r.get('symbol')
                if isinstance(sym, str):
                    self.gene_to_ensg[('SYMBOL', sym.upper())].add(ensg)
                hg = r.get('hgnc_id')
                if isinstance(hg, str):
                    self.gene_to_ensg[('HGNC', hg.replace('HGNC:', ''))].add(ensg)
                up = r.get('uniprot_ids')
                if isinstance(up, str):
                    for u in up.split('|'):
                        if u.strip():
                            self.gene_to_ensg[('UNIPROTKB', u.strip())].add(ensg)
        elif self._v:
            print('  hgnc_complete_set.txt missing — Entrez/symbol genes will not resolve')

        # Open Targets target: UniProt and symbol -> Ensembl, as a second source
        d = g / 'ot_target'
        parts = sorted(d.glob('*.parquet')) if d.exists() else []
        for f in parts:
            t = pd.read_parquet(f, columns=['id', 'approvedSymbol', 'proteinIds'])
            for gid, sym, prots in zip(t['id'], t['approvedSymbol'], t['proteinIds']):
                if not isinstance(gid, str) or not gid.startswith('ENSG'):
                    continue
                self.ensg_all.add(gid)
                if isinstance(sym, str):
                    self.gene_to_ensg[('SYMBOL', sym.upper())].add(gid)
                if prots is not None:
                    for pr in prots:
                        if 'uniprot' in str(pr.get('source', '')):
                            self.gene_to_ensg[('UNIPROTKB', pr['id'])].add(gid)

    def _load_drug(self, g):
        d = g / 'ot_drug_molecule'
        parts = sorted(d.glob('*.parquet')) if d.exists() else []
        if not parts and self._v:
            print('  ot_drug_molecule missing — drugs will resolve only if already ChEMBL')
        for f in parts:
            dm = pd.read_parquet(f, columns=['id', 'crossReferences'])
            for cid, xrefs in zip(dm['id'], dm['crossReferences']):
                if not isinstance(cid, str):
                    continue
                self.chembl_all.add(cid)
                if xrefs is None:
                    continue
                for x in xrefs:
                    src = str(x.get('source', '')).upper()
                    # `ids` is a numpy array; `or []` raises on truth-testing it
                    ids = x.get('ids')
                    if ids is None:
                        continue
                    for ext in ids:
                        self.drug_to_chembl[(src, str(ext))].add(cid)
                        if src == 'DRUGBANK':
                            self.drug_to_chembl[(None, str(ext))].add(cid)

        # PubChem CID -> ChEMBL, from UniChem src1src22.
        #
        # Open Targets' crossReferences carries DrugBank, Probes&Drugs,
        # DailyMed, USAN, INN and EMA, but NOT PubChem. Without this file a
        # graph that keys drugs by PubChem CID (OpenBioLink does) resolves
        # nothing at all and is scored as holding no approved drugs, which is a
        # namespace gap being reported as a content gap.
        #
        # Stored as int -> int in its own dict rather than as (namespace, str)
        # keys in drug_to_chembl. The table is 2.8M rows, of which only ~0.7%
        # reach a molecule Open Targets knows; keeping it as string keys under
        # three namespace spellings cost 3.5 GB and 12 s per resolver, and a
        # resolver is built in every notebook. Both ids are numeric by
        # construction ('CHEMBL941' -> 941, CID -> int), so the compact form is
        # lossless: no mapping is dropped, and resolvability still counts a CID
        # that reaches ANY ChEMBL id, not only an approved one.
        self.pubchem_to_chembl = {}
        uc = g / 'unichem_pubchem_chembl.tsv'
        if uc.exists():
            multi = 0
            with uc.open() as fh:
                for line in fh:
                    a, _, b = line.partition('\t')
                    a = a.strip().upper()
                    b = b.strip()
                    if not a.startswith('CHEMBL') or not b.isdigit():
                        continue
                    num = a[6:]
                    if not num.isdigit():
                        continue
                    cid, ch = int(b), int(num)
                    prev = self.pubchem_to_chembl.get(cid)
                    if prev is None:
                        self.pubchem_to_chembl[cid] = ch
                    elif prev != ch:
                        multi += 1
                        self.pubchem_to_chembl[cid] = (
                            prev + (ch,) if isinstance(prev, tuple) else (prev, ch))
            if self._v:
                print(f'  unichem PubChem->ChEMBL: {len(self.pubchem_to_chembl):,} CIDs'
                      + (f' ({multi:,} multi-mapped)' if multi else ''))
        elif self._v:
            print('  unichem_pubchem_chembl.tsv MISSING — PubChem CIDs will not '
                  'resolve, so any KG keyed by PubChem (OpenBioLink) will score '
                  'near zero on drugs. Run tools/fetch_unichem.sh')

    _PUBCHEM_NS = ('PUBCHEM.COMPOUND', 'PUBCHEM', 'PUBCHEM_COMPOUND', 'CID')

    def _pubchem(self, local):
        """PubChem CID (any spelling of the namespace) -> {CHEMBLxxx}."""
        if not str(local).isdigit():
            return set()
        v = self.pubchem_to_chembl.get(int(local))
        if v is None:
            return set()
        return ({f'CHEMBL{x}' for x in v} if isinstance(v, tuple)
                else {f'CHEMBL{v}'})

    _GO_ASPECT = {'biological_process': 'Biological Process',
                  'molecular_function': 'Molecular Function',
                  'cellular_component': 'Cellular Component'}

    @staticmethod
    def _iter_obo(path, prefix):
        """(id, namespace, is_obsolete) for every term with the given prefix."""
        cur = ns = None
        obsolete = False
        with open(path, encoding='utf-8', errors='replace') as fh:
            for line in fh:
                line = line.rstrip('\n')
                if line == '[Term]':
                    if cur:
                        yield cur, ns, obsolete
                    cur, ns, obsolete = None, None, False
                elif line.startswith('id: ') and line[4:].startswith(prefix):
                    cur = line[4:].strip()
                elif line.startswith('namespace: '):
                    ns = line[11:].strip()
                elif line.startswith('is_obsolete: true'):
                    obsolete = True
        if cur:
            yield cur, ns, obsolete

    def _load_obo_types(self, g):
        """GO and UBERON, the entity types the four canonical namespaces miss.

        GO carries the aspect on each term, so a node a graph calls a
        Biological Process can be checked against whether its GO term really is
        one. That distinction is invisible to any pattern over the identifier,
        since all three aspects share the form GO:nnnnnnn.

        Obsolete terms are skipped, so a retired identifier reads as unresolved
        rather than as valid.
        """
        go = g / 'go.obo'
        if go.exists():
            for tid, ns, dead in self._iter_obo(go, 'GO:'):
                if dead:
                    self.obsolete.add(tid)
                    continue
                canon = self._GO_ASPECT.get(ns)
                if canon:
                    self.go_aspect[tid] = canon
        elif self._v:
            print('  go.obo missing — GO-typed nodes will not resolve; '
                  'run python tools/build_crosswalks.py --only go')

        md = g / 'meddra.tsv'
        if md.exists():
            with md.open(encoding='utf-8', errors='replace') as fh:
                self.umls_se = {ln.split('\t', 1)[0].strip() for ln in fh if ln.strip()}
        elif self._v:
            print('  meddra.tsv missing — Side Effect nodes will not resolve; '
                  'download from https://sideeffects.embl.de/download/')

        wp = g / 'wikipathways_ids.txt'
        if wp.exists():
            self.wikipathways = {l.strip() for l in wp.open() if l.strip()}

        ub = g / 'uberon.obo'
        if ub.exists():
            for tid, _, dead in self._iter_obo(ub, 'UBERON:'):
                (self.obsolete if dead else self.uberon).add(tid)
        elif self._v:
            print('  uberon.obo missing — Anatomy nodes will not resolve; '
                  'run python tools/build_crosswalks.py --only uberon')

    # ------------------------------------------------------------ resolving
    def _resolve(self, raw_ids, fn):
        out, st = set(), ResolutionStats()
        for raw in raw_ids:
            st.n_input += 1
            ns, local = split_curie(raw)
            hits = fn(ns, local, raw)
            if hits:
                st.n_resolved += 1
                st.by_namespace[ns or 'bare'] += 1
                out |= hits
            else:
                st.failed_by_namespace[ns or 'bare'] += 1
        st.n_output = len(out)
        return out, st

    def diseases_to_mondo(self, raw_ids):
        def f(ns, local, raw):
            if ns == 'MONDO':
                return {_num(local)}
            if ns is None:
                hits = set(self.x_to_mondo.get(('MESH', local), set()))
                return hits or set(self.x_to_mondo.get(('MONDO', local), set()))
            key = _num(local) if ns in ('DOID', 'OMIM', 'MIM', 'ORPHANET') else local
            return (set(self.x_to_mondo.get((ns, key), set()))
                    or set(self.x_to_mondo.get((ns, local), set())))
        return self._resolve(raw_ids, f)

    def drugs_to_chembl(self, raw_ids):
        def f(ns, local, raw):
            s = str(raw).strip().split('::')[-1]
            if s.upper().startswith('CHEMBL'):
                return {s.upper()}
            if ns == 'CHEMBL':
                return {('CHEMBL' + local).upper()}
            if re.fullmatch(r'DB\d+', s):
                return set(self.drug_to_chembl.get(('DRUGBANK', s), set()))
            if ns in self._PUBCHEM_NS:
                return self._pubchem(local)
            if ns:
                return set(self.drug_to_chembl.get((ns, local), set()))
            return set(self.drug_to_chembl.get((None, s), set()))
        return self._resolve(raw_ids, f)

    def genes_to_ensembl(self, raw_ids):
        def f(ns, local, raw):
            s = str(raw).strip().split('::')[-1]
            if s.startswith('ENSG'):
                return {s}
            if ns in ('NCBIGENE', 'ENTREZ', 'GENEID'):
                return set(self.gene_to_ensg.get(('ENTREZ', _num(local)), set()))
            if ns == 'UNIPROTKB':
                return set(self.gene_to_ensg.get(('UNIPROTKB', local), set()))
            if ns == 'HGNC':
                return set(self.gene_to_ensg.get(('HGNC', _num(local)), set()))
            if s.isdigit():                       # bare Entrez, e.g. Hetionet
                return set(self.gene_to_ensg.get(('ENTREZ', _num(s)), set()))
            return (set(self.gene_to_ensg.get(('UNIPROTKB', s), set()))
                    or set(self.gene_to_ensg.get(('SYMBOL', s.upper()), set())))
        return self._resolve(raw_ids, f)

    def pathways_to_reactome(self, raw_ids):
        """Reactome ids that are pathways, checked against the pathway list.

        Falls back to the prefix test only when the list could not be loaded,
        so a missing crosswalk degrades loudly rather than silently widening
        what counts as a pathway.
        """
        def f(ns, local, raw):
            s = str(raw).strip().split('::')[-1]
            if ns == 'REACTOME' and str(local).startswith('R-'):
                s = str(local)
            if not s.startswith('R-'):
                return set()
            if not self.reactome:
                return {s} if s.startswith('R-HSA-') else set()
            return {s} if s in self.reactome else set()
        return self._resolve(raw_ids, f)

    def _go(self, raw_ids, want):
        def f(ns, local, raw):
            tid = raw if str(raw).startswith('GO:') else (f'GO:{local}' if ns == 'GO' else None)
            return {tid} if tid and self.go_aspect.get(tid) == want else set()
        return self._resolve(raw_ids, f)

    def bp_terms(self, raw_ids):
        return self._go(raw_ids, 'Biological Process')

    def mf_terms(self, raw_ids):
        return self._go(raw_ids, 'Molecular Function')

    def cc_terms(self, raw_ids):
        return self._go(raw_ids, 'Cellular Component')

    def side_effects_to_umls(self, raw_ids):
        """UMLS CUI naming a MedDRA side-effect concept (SIDER 4.1 registry)."""
        def f(ns, local, raw):
            cui = str(raw).split('::')[-1].strip()
            if ns == 'UMLS':
                cui = local
            return {cui} if cui in self.umls_se else set()
        return self._resolve(raw_ids, f)

    def anatomy_to_uberon(self, raw_ids):
        def f(ns, local, raw):
            tid = raw if str(raw).startswith('UBERON:') else (
                f'UBERON:{local}' if ns == 'UBERON' else None)
            return {tid} if tid in self.uberon else set()
        return self._resolve(raw_ids, f)

    _WP_REV = re.compile(r'^(WP\d+)(?:_r\d+)?$')

    def is_valid_pathway(self, raw):
        """Is this a real pathway accession, in any vocabulary we can check?

        Distinct from ``pathways_to_reactome``, which asks whether an
        identifier reaches the Reactome reference set. A WikiPathways accession
        is a perfectly valid pathway that maps to no Reactome id, so coverage
        (nb01) and type validity (nb02) must not share one predicate.

        Hetionet pins WikiPathways by revision (``WP100_r74146``); the revision
        is stripped, since it names a version of the pathway, not a different
        pathway.
        """
        s = str(raw).split('::')[-1].strip()
        if s.startswith('R-') or s.upper().startswith(('REACTOME:', 'REACT:')):
            return bool(self.pathways_to_reactome([raw])[0])
        m = self._WP_REV.match(s)
        return bool(m and m.group(1) in self.wikipathways)

    # Check 1 (which identifier shapes each entity type accepts) moved to
    # config/id_schemes.yaml and biokgsuite.io.schemes.IdSchemes. The
    # namespace-prefix version that lived here could not see a bare accession
    # such as DB12494, which carries no prefix at all.

    _ALL_RESOLVERS = {'Drug': 'drugs_to_chembl', 'Disease': 'diseases_to_mondo',
                      'Gene/Protein': 'genes_to_ensembl', 'Pathway': 'pathways_to_reactome',
                      'Biological Process': 'bp_terms', 'Molecular Function': 'mf_terms',
                      'Cellular Component': 'cc_terms', 'Anatomy': 'anatomy_to_uberon',
                      'Side Effect': 'side_effects_to_umls'}

    @property
    def RESOLVERS(self):
        """Types this instance can resolve; the auxiliary ones only if loaded."""
        if self._aux:
            return self._ALL_RESOLVERS
        return {k: v for k, v in self._ALL_RESOLVERS.items()
                if k not in self.AUX_REGISTRIES}

    def resolve(self, entity_type, raw_ids):
        return getattr(self, self.RESOLVERS[entity_type])(raw_ids)
