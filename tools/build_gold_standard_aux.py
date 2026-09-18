#!/usr/bin/env python3
"""Rebuild the derivable files in data/gold_standards/bridges/ from crosswalks already here.

Notebooks 06 and 07 read nine files from data/gold_standards/bridges/ (and ctd/). Most of them were
never external downloads: they are projections of ontologies and Open Targets
tables that data/crosswalks/ already carries, so they can be rebuilt offline and
carry the same version stamp as their source.

    mesh_to_doid.csv          MESH:* -> DOID:*        from doid.obo xrefs
    do_diseases.csv           DOID:* -> MONDO id/name from mondo.sssom.tsv + mondo.obo
    drugbank_drugs.csv        DB* -> drug name        from ot_drug_molecule
    uniprot_genesproteins.csv UniProt accessions      from ot_target
    mondo.sssom.tsv           copied from data/crosswalks/ (it was duplicated)

NOT rebuilt here, and why:
    CTD_chem_gene_ixns.tsv.gz        external, see tools/fetch_ctd.sh
    CTD_curated_genes_diseases.tsv.gz  external, see tools/fetch_ctd.sh
    drugbank_xref.csv                see tools/build_drugbank_xref.py
    fda_drugs_id_mapping.csv         hand-curated (Combo_Flag, a manual exclusion
                                     for entry 203); no download reproduces it

Key formats matter and are not obvious: notebook 07's ontology classifier looks
up `mesh_to_doid[node_id]` where node_id is `MESH:D006394`, and notebook 06 maps
CTD's `disease_id`, also `MESH:D*`. So mesh_id carries the `MESH:` prefix and
doid carries `DOID:`. Writing them bare would map nothing and fail silently.

    python tools/build_gold_standard_aux.py
"""
from __future__ import annotations

import glob
import shutil
import sys
from pathlib import Path

import pandas as pd

ROOT = next(p for p in Path(__file__).resolve().parents
            if (p / 'config' / 'kgs.yaml').exists())
XW = ROOT / 'data' / 'crosswalks'
GS = ROOT / 'data' / 'gold_standards' / 'bridges'


def _need(path: Path) -> Path:
    if not path.exists():
        sys.exit(f'missing input: {path.relative_to(ROOT)}\n'
                 f'  run tools/build_crosswalks.py / fetch_opentargets.sh first')
    return path


def mesh_to_doid() -> pd.DataFrame:
    """MESH:* -> DOID:*, from the Disease Ontology's own xrefs."""
    rows, cur = [], None
    for line in _need(XW / 'doid.obo').read_text(errors='replace').splitlines():
        line = line.strip()
        if line == '[Term]':
            cur = None
        elif line.startswith('id: DOID:'):
            cur = line[4:].strip()
        elif line.startswith('xref: MESH:') and cur:
            rows.append({'mesh_id': line[6:].strip(), 'doid': cur})
    df = pd.DataFrame(rows).drop_duplicates(subset='mesh_id')
    return df.sort_values('mesh_id')


def do_diseases() -> pd.DataFrame:
    """DOID:* -> MONDO id and label, from the MONDO SSSOM plus MONDO labels."""
    ss = pd.read_csv(_need(XW / 'mondo.sssom.tsv'), sep='\t', comment='#',
                     low_memory=False)
    ss = ss[(ss['predicate_id'] == 'skos:exactMatch')
            & ss['object_id'].astype(str).str.startswith('DOID:')
            & ss['subject_id'].astype(str).str.startswith('MONDO:')]

    names = {}
    cur = None
    for line in _need(XW / 'mondo.obo').read_text(errors='replace').splitlines():
        line = line.strip()
        if line == '[Term]':
            cur = None
        elif line.startswith('id: MONDO:'):
            cur = line[4:].strip()
        elif line.startswith('name:') and cur:
            names[cur] = line[5:].strip()

    df = pd.DataFrame({'doid': ss['object_id'].astype(str),
                       'mondo_id': ss['subject_id'].astype(str)})
    df['mondo_name'] = df['mondo_id'].map(names)
    df = df.dropna(subset=['mondo_name']).drop_duplicates(subset='doid')
    return df.sort_values('doid')[['doid', 'mondo_id', 'mondo_name']]


def drugbank_drugs() -> pd.DataFrame:
    """DrugBank id -> drug name, from Open Targets crossReferences.

    Needs no DrugBank account: Open Targets carries the DrugBank ids alongside
    its own names, and is already version-pinned in data/REFERENCE_SETS.md.
    """
    files = sorted(glob.glob(str(XW / 'ot_drug_molecule' / '*.parquet')))
    if not files:
        sys.exit('missing input: data/crosswalks/ot_drug_molecule/*.parquet')
    rows = []
    for f in files:
        d = pd.read_parquet(f, columns=['id', 'name', 'crossReferences'])
        for name, xrefs in zip(d['name'], d['crossReferences']):
            if xrefs is None or (hasattr(xrefs, '__len__') and len(xrefs) == 0):
                continue
            for entry in xrefs:
                if str(entry.get('source', '')).lower() != 'drugbank':
                    continue
                # `ids` is a numpy array here, so `x or []` raises on its
                # truthiness rather than defaulting. Check for None explicitly.
                ids = entry.get('ids')
                if ids is None:
                    continue
                for db in ids:
                    rows.append({'drugbank_id': str(db), 'name': name})
    return pd.DataFrame(rows).drop_duplicates(subset='drugbank_id').sort_values('drugbank_id')


def uniprot_genesproteins() -> pd.DataFrame:
    """UniProt accessions with their gene symbol and Ensembl id, from ot_target."""
    files = sorted(glob.glob(str(XW / 'ot_target' / '*.parquet')))
    if not files:
        sys.exit('missing input: data/crosswalks/ot_target/*.parquet')
    rows = []
    for f in files:
        d = pd.read_parquet(f, columns=['id', 'approvedSymbol', 'proteinIds'])
        for ensg, sym, prots in zip(d['id'], d['approvedSymbol'], d['proteinIds']):
            if prots is None:
                continue
            for p in prots:
                if str(p.get('source', '')).startswith('uniprot'):
                    rows.append({'uniprot_id': str(p['id']), 'gene_symbol': sym,
                                 'ensembl_id': ensg})
    return pd.DataFrame(rows).drop_duplicates(subset='uniprot_id').sort_values('uniprot_id')


def main() -> None:
    GS.mkdir(parents=True, exist_ok=True)
    for name, fn in [('mesh_to_doid.csv', mesh_to_doid),
                     ('do_diseases.csv', do_diseases),
                     ('drugbank_drugs.csv', drugbank_drugs),
                     ('uniprot_genesproteins.csv', uniprot_genesproteins)]:
        df = fn()
        df.to_csv(GS / name, index=False)
        print(f'  {name:28s} {len(df):>7,} rows   {list(df.columns)}')

    src = XW / 'mondo.sssom.tsv'
    if src.exists():
        shutil.copyfile(src, GS / 'mondo.sssom.tsv')
        print(f'  {"mondo.sssom.tsv":28s} copied from data/crosswalks/')

    print('\nStill needed, not derivable here:')
    print('  CTD_chem_gene_ixns.tsv.gz         tools/fetch_ctd.sh')
    print('  CTD_curated_genes_diseases.tsv.gz tools/fetch_ctd.sh')
    print('  drugbank_xref.csv                 tools/build_drugbank_xref.py')
    print('  fda_drugs_id_mapping.csv          hand-curated; no download reproduces it')


if __name__ == '__main__':
    main()
