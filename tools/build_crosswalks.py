#!/usr/bin/env python3
"""Fetch every identifier crosswalk from a pinned, citable source.

Run:  python tools/build_crosswalks.py [--outdir data/crosswalks] [--only NAME]

Writes each crosswalk to --outdir and appends a row to data/CROSSWALKS.md
recording source URL, pinned version, retrieval date, SHA-256 and row count,
so every number in the paper traces to a specific release.

Network note: MONDO and the Disease Ontology are served from GitHub. UniChem,
UniProt and Reactome are on ebi.ac.uk / uniprot.org / reactome.org, which some
networks block. Sources that cannot be reached are reported and skipped, not
silently omitted; rerun from a network that can reach them.

Nothing here is licence-gated. DrugBank's own external_links.csv requires an
account and is deliberately NOT a dependency: UniChem supplies the same
DrugBank cross-references under CC-BY-SA.
"""
from __future__ import annotations
import argparse, datetime as dt, gzip, hashlib, io, json, sys, urllib.request
from pathlib import Path

UA = {'User-Agent': 'biokgsuite-crosswalk-builder/1.0'}
TIMEOUT = 180

# Pinned releases. Bump deliberately, never silently.
MONDO_TAG = 'v2026-08-04'
DOID_TAG  = 'v2026-08-31'
OT_RELEASE = '26.06'          # Open Targets Platform, released 2026-06-24
OT_FTP = f'https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/{OT_RELEASE}/output'
CHEMBL_RELEASE = 'chembl_36'  # verify against ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest

# GO and UBERON supply the entity types the four canonical namespaces do not
# cover. They are pinned like everything else here: UBERON by release tag, GO by
# commit SHA, because the GO repository publishes no releases.
UBERON_TAG = 'v2026-06-23'
GO_COMMIT  = 'f90ad02eb5c6072eb90bf73b84123eda087e16c7'   # go-edit.obo, 2026-09-14
# HPO is the registry behind PrimeKG's effect/phenotype nodes, which are written
# as bare numbers ('118' for HP:0000118) and would otherwise be unverifiable.
HPO_TAG = 'v2026-09-01'
# WikiPathways publishes one directory per pathway; the git tree is the registry.
# PathwayCommons publishes no comparable registry, so PC7 accessions are not verified.
WIKIPATHWAYS_TREE = '970b71ee0bedabb5025c2f3b4cc6dfcb2e475bb5'   # pathways/, 2026-09

SOURCES = {
    'mondo_sssom': dict(
        url=f'https://raw.githubusercontent.com/monarch-initiative/mondo/{MONDO_TAG}'
            '/src/ontology/mappings/mondo.sssom.tsv',
        version=MONDO_TAG, out='mondo.sssom.tsv',
        maps='MONDO <-> DOID, UMLS, MEDGEN, OMIM, Orphanet, SCTID, MeSH, NCIT, ICD, EFO',
        licence='CC-BY 4.0'),
    'doid': dict(
        url=f'https://raw.githubusercontent.com/DiseaseOntology/HumanDiseaseOntology/{DOID_TAG}'
            '/src/ontology/doid.obo',
        version=DOID_TAG, out='doid.obo',
        maps='DOID terms and their xrefs (MeSH, UMLS, NCI, ICD, OMIM)',
        licence='CC0 1.0'),
    'go': dict(
        url=f'https://raw.githubusercontent.com/geneontology/go-ontology/{GO_COMMIT}'
            '/src/ontology/go-edit.obo',
        version=GO_COMMIT[:12], out='go.obo',
        maps='GO terms -> aspect (biological_process / molecular_function / '
             'cellular_component); the aspect is what lets a GO-typed node be '
             'checked against the type its graph declares',
        licence='CC-BY 4.0'),
    'uberon': dict(
        url=f'https://github.com/obophenotype/uberon/releases/download/{UBERON_TAG}'
            '/uberon-basic.obo',
        version=UBERON_TAG, out='uberon.obo',
        maps='UBERON anatomy terms',
        licence='CC-BY 3.0'),
    'hpo': dict(
        url=f'https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/{HPO_TAG}'
            '/hp.obo',
        version=HPO_TAG, out='hp.obo',
        maps='HPO phenotype terms, live and obsolete; confirms that a bare '
             'phenotype number names a real term',
        licence='HPO licence (free for academic use, attribution required)'),
    'wikipathways': dict(
        url='https://api.github.com/repos/wikipathways/wikipathways-database'
            f'/git/trees/{WIKIPATHWAYS_TREE}',
        version=WIKIPATHWAYS_TREE[:12], out='wikipathways_ids.txt', wp_tree=True,
        maps='valid WikiPathways accessions (WP*), for pathway type validity',
        licence='CC0 1.0'),
    'reactome': dict(
        url='https://reactome.org/download/current/ReactomePathways.txt',
        version='current', out='reactome_pathways.txt',
        maps='Reactome stable IDs -> pathway name, species',
        licence='CC-BY 4.0'),
    'ot_clinical_indication': dict(
        url=f'{OT_FTP}/clinical_indication/',
        version=OT_RELEASE, out='ot_drug_indication/', parquet=True,
        maps='drug (ChEMBL id) -> disease (EFO/MONDO), by clinical phase',
        licence='CC0 1.0'),
    'ot_molecule': dict(
        url=f'{OT_FTP}/drug_molecule/',
        version=OT_RELEASE, out='ot_drug_molecule/', parquet=True,
        maps='ChEMBL drug id -> name, type, maximumClinicalTrialPhase, isApproved',
        licence='CC0 1.0'),
    'ot_mechanism_of_action': dict(
        url=f'{OT_FTP}/drug_mechanism_of_action/',
        version=OT_RELEASE, out='ot_drug_mechanism_of_action/', parquet=True,
        maps='drug -> target, the mechanistic bridge (there is no known_drug '
             'dataset in 26.06; names changed from earlier releases)',
        licence='CC0 1.0'),
    'ot_disease': dict(
        url=f'{OT_FTP}/disease/', version=OT_RELEASE, out='ot_disease/', parquet=True,
        maps='Open Targets disease list with EFO/MONDO ids', licence='CC0 1.0'),
    'ot_target': dict(
        url=f'{OT_FTP}/target/', version=OT_RELEASE, out='ot_target/', parquet=True,
        maps='gene list: Ensembl id, UniProt accessions, HGNC symbol', licence='CC0 1.0'),
    'ot_reactome': dict(
        url=f'{OT_FTP}/reactome/', version=OT_RELEASE, out='ot_reactome/', parquet=True,
        maps='Reactome pathway list as harmonised by Open Targets', licence='CC0 1.0'),
    'hgnc': dict(
        url='https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/'
            'hgnc_complete_set.txt',
        version='current', out='hgnc_complete_set.txt',
        maps='HGNC id <-> symbol <-> Entrez <-> Ensembl gene id. The ONLY source '
             'here for Entrez->Ensembl, which Hetionet, DRKG and PrimeKG need '
             '(Open Targets `target` carries no NCBI/Entrez xref).',
        licence='CC0 1.0'),
    'uniprot_human': dict(
        url='https://rest.uniprot.org/uniprotkb/stream?format=tsv&compressed=true'
            '&query=%28reviewed%3Atrue%29%20AND%20%28model_organism%3A9606%29'
            '&fields=accession%2Cid%2Cxref_geneid%2Cgene_names%2Cprotein_name',
        version='current', out='uniprot_human_reviewed.tsv', gz=True,
        maps='UniProt accession -> Entrez GeneID (human, reviewed only)',
        licence='CC-BY 4.0'),
}

UNICHEM_API = 'https://www.ebi.ac.uk/unichem/api/v1/sources'
UNICHEM_FTP = 'https://ftp.ebi.ac.uk/pub/databases/chembl/UniChem/data/wholeSourceMapping'
# Names resolved to numeric src_ids at runtime; the IDs are not hardcoded
# because they are a property of the live UniChem release.
UNICHEM_WANT = ['drugbank', 'chembl', 'pubchem', 'chebi', 'kegg', 'drugcentral']


def fetch(url, gz=False):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
        raw = r.read()
    return gzip.decompress(raw) if gz else raw


def sha256(b): return hashlib.sha256(b).hexdigest()


def count_records(fname, data):
    """Records, not lines. An OBO file's unit is a [Term] stanza, not a line."""
    if fname.endswith('.obo'):
        return data.count(b'\n[Term]')
    return sum(1 for ln in data.split(b'\n') if ln and not ln.startswith(b'#')) - 1


def unichem_pairs(outdir, rows):
    """Resolve source names -> src_id from the live API, then pull DrugBank<->X."""
    try:
        srcs = json.loads(fetch(UNICHEM_API).decode())
    except Exception as e:
        print(f'  unichem      SKIP ({type(e).__name__}: {str(e)[:60]})')
        print(f'               source IDs must come from {UNICHEM_API}; not hardcoded')
        return
    items = srcs.get('sources', srcs) if isinstance(srcs, dict) else srcs
    idx = {}
    for s in items:
        name = str(s.get('name') or s.get('nameLabel') or '').lower()
        sid = s.get('sourceID') or s.get('src_id') or s.get('id')
        for want in UNICHEM_WANT:
            if want in name and want not in idx:
                idx[want] = int(sid)
    print(f'  unichem      resolved src_ids: {idx}')
    if 'drugbank' not in idx:
        print('  unichem      DrugBank src_id not found; aborting UniChem step')
        return
    db = idx['drugbank']
    for want, sid in idx.items():
        if want == 'drugbank':
            continue
        lo, hi = min(db, sid), max(db, sid)
        url = f'{UNICHEM_FTP}/src_id{lo}/src{lo}src{hi}.txt.gz'
        try:
            data = fetch(url, gz=True)
        except Exception as e:
            print(f'  unichem      {want}: SKIP ({str(e)[:50]})')
            continue
        name = f'unichem_drugbank_{want}.tsv'
        (outdir / name).write_bytes(data)
        n = data.count(b'\n') - 1
        rows.append(dict(name=f'unichem:drugbank-{want}', out=name, url=url,
                         version='current', rows=n, sha=sha256(data),
                         maps=f'DrugBank <-> {want}', licence='CC-BY-SA 3.0'))
        print(f'  unichem      {want}: {n:,} rows')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--outdir', default='data/crosswalks')
    ap.add_argument('--only', help='comma-separated source names; '
                                   'default is every source')
    a = ap.parse_args()

    # An unknown name used to select nothing and exit silently, which looks
    # exactly like a successful run that fetched nothing. Fail loudly instead.
    only = {x.strip() for x in a.only.split(',')} if a.only else None
    if only:
        unknown = only - set(SOURCES)
        if unknown:
            ap.error('unknown source(s): ' + ', '.join(sorted(unknown))
                     + '\navailable: ' + ', '.join(sorted(SOURCES)))

    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    today = dt.date.today().isoformat()
    rows = []

    for name, spec in SOURCES.items():
        if only and name not in only:
            continue
        if spec.get('wp_tree'):
            import re as _re
            tree = json.loads(fetch(spec['url']).decode())
            ids = sorted({x['path'] for x in tree.get('tree', [])
                          if _re.fullmatch(r'WP\d+', x['path'])},
                         key=lambda x: int(x[2:]))
            data = ('\n'.join(ids) + '\n').encode()
            (out / spec['out']).write_bytes(data)
            rows.append(dict(name=name, out=spec['out'], url=spec['url'],
                             version=spec['version'], rows=len(ids), sha=sha256(data),
                             maps=spec['maps'], licence=spec['licence']))
            print(f"  {name:15s} {len(ids):,} rows  {spec['version']}")
            continue
        if spec.get('parquet'):
            # Open Targets datasets are parquet DIRECTORIES; they are fetched by
            # tools/fetch_opentargets.sh, which handles multi-part listings.
            continue
        dest = out / spec['out']
        if dest.is_file() and dest.stat().st_size > 0:
            data = dest.read_bytes()
            n = count_records(spec['out'], data)
            rows.append(dict(name=name, out=spec['out'], url=spec['url'],
                             version=spec['version'], rows=n, sha=sha256(data),
                             maps=spec['maps'], licence=spec['licence']))
            print(f'  {name:12s} {n:>9,} rows  {spec["version"]}  (already present)')
            continue
        try:
            data = fetch(spec['url'], gz=spec.get('gz', False))
        except Exception as e:
            print(f'  {name:12s} SKIP ({type(e).__name__}: {str(e)[:60]})')
            continue
        dest.write_bytes(data)
        n = count_records(spec['out'], data)
        rows.append(dict(name=name, out=spec['out'], url=spec['url'],
                         version=spec['version'], rows=n, sha=sha256(data),
                         maps=spec['maps'], licence=spec['licence']))
        print(f'  {name:12s} {n:>9,} rows  {spec["version"]}')

    if not a.only or a.only == 'unichem':
        unichem_pairs(out, rows)

    if a.only:
        print('\n--only run: leaving data/CROSSWALKS.md untouched '
              '(a partial run must not overwrite the full manifest).')
        return 0

    man = Path('data/CROSSWALKS.md')
    got = {r['name'] for r in rows}
    missing = [(k, v) for k, v in SOURCES.items() if k not in got]

    out_md = ['# Identifier crosswalks', '',
              'Generated by `tools/build_crosswalks.py` on ' + today + '.',
              'Every row pins a release, so any reported coverage figure traces to a',
              'specific version of each reference. Re-running with a different pinned',
              'version changes the SHA-256 and is therefore visible in the diff.', '',
              '## Retrieved', '',
              '| Crosswalk | File | Version | Retrieved | Records | Maps | Licence | SHA-256 (first 16) |',
              '|---|---|---|---|---|---|---|---|']
    for r in rows:
        out_md.append(f"| {r['name']} | `{r['out']}` | {r['version']} | {today} | "
                      f"{r['rows']:,} | {r['maps']} | {r['licence']} | `{r['sha'][:16]}` |")
    out_md += ['', '### Source URLs', '']
    for r in rows:
        out_md.append(f"- **{r['name']}** — {r['url']}")

    if missing:
        out_md += ['', '## Not retrieved on this run', '',
                   'These hosts were unreachable from the machine that ran the script',
                   '(commonly a proxy or egress allowlist). They are NOT optional; rerun',
                   'from a network that can reach them before reporting final numbers.', '',
                   '| Crosswalk | Maps | Licence | URL |', '|---|---|---|---|']
        for k, v in missing:
            out_md.append(f"| {k} | {v['maps']} | {v['licence']} | {v['url']} |")
        out_md += ['',
                   '| unichem | DrugBank <-> ChEMBL, PubChem, ChEBI, KEGG, DrugCentral | '
                   f'CC-BY-SA 3.0 | {UNICHEM_FTP} |', '',
                   'UniChem numeric source IDs are resolved at runtime from',
                   f'`{UNICHEM_API}` rather than hardcoded, because they are a property',
                   'of the live release.']

    out_md += ['', '## Why these sources', '',
        'DrugBank\'s own `external_links.csv` requires a licensed account and cannot be',
        'redistributed or re-fetched by a reader, so it is deliberately not a dependency.',
        'UniChem supplies the same DrugBank cross-references under CC-BY-SA and adds',
        'ChEMBL, which `drugbank_xref.csv` does not cover.', '',
        '## Citations', '',
        '- **MONDO** — Vasilevsky NA, et al. Mondo: Unifying diseases for the world, by the',
        '  world. medRxiv 2022. doi:10.1101/2022.04.13.22273750',
        '- **Human Disease Ontology** — Schriml LM, et al. The Human Disease Ontology 2022',
        '  update. Nucleic Acids Research 2022;50(D1):D1255-D1261. doi:10.1093/nar/gkab1063',
        '- **UniChem** — Chambers J, et al. UniChem: a unified chemical structure',
        '  cross-referencing and identifier tracking system. J Cheminform 2013;5:3.',
        '  doi:10.1186/1758-2946-5-3',
        '- **HGNC** — Seal RL, et al. Genenames.org: the HGNC resources.',
        '  Nucleic Acids Research (current release). NOTE: HGNC moved off the EBI',
        '  FTP to Google Cloud Storage; the old',
        '  ftp.ebi.ac.uk/pub/databases/genenames/ path now 404s. For a stable pin',
        '  use a dated monthly archive under',
        '  storage.googleapis.com/public-download-files/hgnc/archive/archive/monthly/tsv/',
        '  instead of the rolling current file.',
        '- **UniProt** — The UniProt Consortium. UniProt: the Universal Protein',
        '  Knowledgebase. Nucleic Acids Research (current release).',
        '- **Reactome** — Milacic M, et al. The Reactome Pathway Knowledgebase.',
        '  Nucleic Acids Research (current release).', '',
        'Verify the UniProt and Reactome citations against the release actually',
        'downloaded before submission; both publish a new database issue each year.', '']
    man.write_text('\n'.join(out_md))
    print(f'\nwrote {man} ({len(rows)} retrieved, {len(missing) + 1} pending)')

if __name__ == '__main__':
    sys.exit(main())
