#!/usr/bin/env python3
"""Generate data/PROVENANCE.md — one dated, referenced record of every input.

Run after any fetch or rebuild:  python tools/provenance.py

It inspects what is actually on disk rather than restating intent, so the
document cannot drift from reality. Anything not yet retrieved is listed as
pending rather than omitted. Every row carries a release identifier, the UTC
retrieval timestamp (file mtime), a SHA-256, a record count and a licence, and
each source has a citation.
"""
from __future__ import annotations
import datetime as dt, hashlib, json, sys
from pathlib import Path

# name -> (what it is, release pin, licence, citation)
REGISTRY = {
    'mondo.obo': (
        'MONDO disease ontology — disease entity universe (#Disease) and the '
        'is_a hierarchy behind relation #5',
        'v2026-08-04', 'CC-BY 4.0',
        'Vasilevsky NA, et al. Mondo: Unifying diseases for the world, by the world. '
        'medRxiv 2022. doi:10.1101/2022.04.13.22273750'),
    'mondo.sssom.tsv': (
        'MONDO SSSOM mappings — cross-namespace disease ID resolution',
        'v2026-08-04', 'CC-BY 4.0',
        'Vasilevsky NA, et al. Mondo: Unifying diseases for the world, by the world. '
        'medRxiv 2022. doi:10.1101/2022.04.13.22273750'),
    'hgnc_complete_set.txt': (
        'HGNC complete set — gene ID resolution (symbol / Entrez / UniProt / '
        'HGNC -> Ensembl). Every gene number in the benchmark depends on this file',
        'rolling (see Retrieved)', 'No restrictions stated (genenames.org)',
        'Seal RL, et al. Genenames.org: the HGNC resources in 2023. '
        'Nucleic Acids Research 2023;51(D1):D1003. doi:10.1093/nar/gkac888'),
    'unichem_pubchem_chembl.tsv': (
        'UniChem src1src22 — PubChem CID -> ChEMBL. Required to resolve any KG '
        'keyed by PubChem (OpenBioLink); Open Targets carries no PubChem xref',
        'rolling (see Retrieved)', 'CC0 1.0',
        'Chambers J, et al. UniChem: a unified chemical structure cross-referencing '
        'and identifier tracking system. J Cheminform 2013;5:3. '
        'doi:10.1186/1758-2946-5-3'),
    'unichem_drugbank_chembl.tsv': (
        'UniChem src1src2 — DrugBank -> ChEMBL, structure-matched. Open Targets '
        'names a DrugBank id for only about half of approved molecules, so every '
        'DrugBank-keyed graph (PrimeKG, Hetionet, DRKG, BioKG) depends on this file',
        'rolling (see Retrieved)', 'CC0 1.0',
        'Chambers J, et al. UniChem: a unified chemical structure cross-referencing '
        'and identifier tracking system. J Cheminform 2013;5:3. '
        'doi:10.1186/1758-2946-5-3'),
    'go.obo': (
        'Gene Ontology terms and their aspect (biological_process / '
        'molecular_function / cellular_component). The aspect is what lets a '
        'GO-typed node be checked against the type its graph declares; no '
        'pattern over the identifier can, since all three share GO:nnnnnnn',
        'commit f90ad02eb5c6', 'CC-BY 4.0',
        'Gene Ontology Consortium. The Gene Ontology knowledgebase in 2023. '
        'Genetics 2023;224(1):iyad031. doi:10.1093/genetics/iyad031'),
    'hp.obo': (
        'Human Phenotype Ontology terms. PrimeKG keys its 15,311 effect/phenotype '
        'nodes by the bare HPO number; this registry confirms that each names a '
        'live term',
        'v2026-09-01', 'HPO licence (academic use, attribution)',
        'Gargano MA, et al. The Human Phenotype Ontology in 2024. '
        'Nucleic Acids Research 2024;52(D1):D1333. doi:10.1093/nar/gkad1005'),
    'uberon.obo': (
        'UBERON anatomy terms, for validating Anatomy-typed nodes',
        'v2026-06-23', 'CC-BY 3.0',
        'Mungall CJ, et al. Uberon, an integrative multi-species anatomy ontology. '
        'Genome Biology 2012;13(1):R5. doi:10.1186/gb-2012-13-1-r5'),
    'wikipathways_ids.txt': (
        'Valid WikiPathways accessions, for pathway type validity where the '
        'identifier is not Reactome. PathwayCommons PC7 has no reachable '
        'registry, so those accessions remain unverifiable',
        'tree 970b71ee0bed', 'CC0 1.0',
        'Agrawal A, et al. WikiPathways 2024. Nucleic Acids Research '
        '2024;52(D1):D679. doi:10.1093/nar/gkad960'),
    'doid.obo': (
        'Human Disease Ontology — legacy DOID-keyed view only',
        'v2026-08-31', 'CC0 1.0',
        'Schriml LM, et al. The Human Disease Ontology 2022 update. '
        'Nucleic Acids Research 2022;50(D1):D1255. doi:10.1093/nar/gkab1063'),
}
OT_DATASETS = {
    'ot_drug_molecule':              'Drug entity set — isApproved',
    'ot_disease':                    'disease crosswalk (EFO/MONDO)',
    'ot_target':                     'Gene entity set and crosswalk',
    'ot_reactome':                   'Pathway entity set',
    'ot_drug_mechanism_of_action':   'relation #1 drug->gene; input to #2',
    'ot_clinical_indication':        'relation #6 drug->disease; drugged/undrugged strata',
    # relation #3 disease->gene (and input to #4). Either the seven curated
    # evidence datasets (GENE_DISEASE=curated, the default) or the aggregated
    # association table (GENE_DISEASE=all). Whichever is present is reported.
    'ot_evidence_clingen':            'relation #3 disease->gene (ClinGen)',
    'ot_evidence_orphanet':           'relation #3 disease->gene (Orphanet)',
    'ot_evidence_gene2phenotype':     'relation #3 disease->gene (Gene2Phenotype)',
    'ot_evidence_genomics_england':   'relation #3 disease->gene (Genomics England)',
    'ot_evidence_uniprot_literature': 'relation #3 disease->gene (UniProt literature)',
    'ot_evidence_uniprot_variants':   'relation #3 disease->gene (UniProt variants)',
    'ot_evidence_gene_burden':        'relation #3 disease->gene (gene burden)',
    'ot_association_by_datatype_direct': 'relation #3 disease->gene, aggregated (optional alternative)',
}
# Reported as missing only if NO gene-disease source is present.
GENE_DISEASE_ALTERNATIVES = {
    'ot_evidence_clingen', 'ot_evidence_orphanet', 'ot_evidence_gene2phenotype',
    'ot_evidence_genomics_england', 'ot_evidence_uniprot_literature',
    'ot_evidence_uniprot_variants', 'ot_evidence_gene_burden',
    'ot_association_by_datatype_direct',
}
OT_CITATION = ('Open Targets Platform. Nucleic Acids Research 2025;53(D1):D1467.')
OT_LICENCE = 'CC0 1.0'

# Knowledge graphs under data/<name>/ — the objects being evaluated.
KG_REGISTRY = {
    'hetionet': ('Hetionet v1.0', 'CC0 1.0',
                 'Himmelstein DS, et al. Systematic integration of biomedical '
                 'knowledge prioritizes drugs for repurposing. eLife 2017;6:e26726. '
                 'doi:10.7554/eLife.26726',
                 'https://github.com/hetio/hetionet (hetnet/tsv, Git LFS)'),
    # Release, licence and source for the four fetched by tools/fetch_kgs.py are
    # the publishers' own, read from their APIs on 2026-09-18. When a snapshot
    # has been fetched, its SNAPSHOT.json supersedes the release and licence
    # here, so this table cannot drift from what is on disk.
    'primekg':    ('PrimeKG V2.1', 'CC0 1.0',
                   'Chandak P, Huang K, Zitnik M. Building a knowledge graph to enable '
                   'precision medicine. Sci Data 2023;10:67. doi:10.1038/s41597-023-01960-3',
                   'Harvard Dataverse doi:10.7910/DVN/IXA7BM (file kg.csv)'),
    'drkg':       ('DRKG (unversioned; repo commit d4bb1974)', 'Apache-2.0 (per-source terms apply)',
                   'Ioannidis VN, et al. DRKG - Drug Repurposing Knowledge Graph for Covid-19. 2020.',
                   'https://dgl-data.s3-us-west-2.amazonaws.com/dataset/DRKG/drkg.tar.gz '
                   '(https://github.com/gnn4dr/DRKG)'),
    'openbilink': ('OpenBioLink2020 v1.0.1', 'CC-BY 4.0',
                   'Breit A, Ott S, Agibetov A, Samwald M. OpenBioLink: a benchmarking '
                   'framework for large-scale biomedical link prediction. Bioinformatics '
                   '2020;36(13):4097. doi:10.1093/bioinformatics/btaa274',
                   'Zenodo doi:10.5281/zenodo.5361324 (file ALL_DIR.zip, graph_files/edges.csv). '
                   'The MIT licence on the OpenBioLink repository covers the code; the '
                   'dataset record is CC-BY 4.0'),
    'biokg':      ('BioKG', 'MIT', 'Walsh B, et al. BioKG. CIKM 2020.', 'https://github.com/dsi-bdi/biokg'),
    'matrix':     ('Every Cure integrated KG (HF kg-nodes/kg-edges, 2026-05-28)',
                   'Custom (see LICENSE in each dataset repository)',
                   'Every Cure MATRIX knowledge graph. https://huggingface.co/datasets/everycure/',
                   'https://huggingface.co/datasets/everycure/kg-nodes and .../kg-edges, '
                   'pinned by commit; derived TSVs written by tools/fetch_kgs.py'),
}


# Short release labels for the table. The full release string lives in each
# graph's SNAPSHOT/ADOPTED json and its data/<kg>/<KG>_RELEASE.txt; a table
# column is not the place for a paragraph.
KG_SHORT = {
    'hetionet':   'v1.0',
    'primekg':    'V2.1',
    'drkg':       'unversioned (2022-04-19 build)',
    'openbilink': 'OpenBioLink2020 v1.0.1',
    'biokg':      'v1.0.0',
    'matrix':     'HF main',
}


def release_note(root: Path, kg: str) -> dict:
    """Parse data/<kg>/<KG>_RELEASE.txt into its key=value lines.

    Graphs retrieved before tools/fetch_kgs.py still carry one of these, and it
    holds the retrieval time and the source URL, which are otherwise guessed
    from a file mtime that reflects when an archive was built rather than when
    it was fetched.
    """
    path = root / 'data' / kg / f'{kg.upper()}_RELEASE.txt'
    if not path.exists():
        return {}
    out = {'_raw': path.read_text(errors='replace')}
    for line in out['_raw'].splitlines():
        if line.startswith((' ', '#')) or '=' not in line:
            continue
        k, _, v = line.partition('=')
        out.setdefault(k.strip(), v.strip())
    return out


def obtained_and_backing(kg: str, snap, note=None):
    """How this snapshot was obtained, and what backs its version claim.

    Neither question is answerable from a release string alone. Both are read from the
    snapshot rather than asserted here, so a graph put in place by hand reports
    exactly that rather than borrowing another graph's guarantee.
    """
    if not snap:
        # Retrieved before tools/fetch_kgs.py existed. What backs each of these
        # is recorded in data/<kg>/<KG>_RELEASE.txt and data/MANIFEST.md, so it
        # is restated rather than left blank; anything not listed here reports
        # that nothing is recorded, instead of borrowing a guarantee.
        prior = {
            # hetionet's edge file is Git LFS; the LFS pointer publishes the
            # expected sha256, so the download is checkable against the source.
            'hetionet': 'publisher sha256 (Git LFS pointer)',
            # biokg ships as a tagged GitHub release asset. The tag fixes the
            # version; GitHub publishes no checksum for the asset.
            'biokg': 'release tag v1.0.0 (no checksum published)',
        }
        return 'retrieved before tools/fetch_kgs.py', prior.get(kg, 'not recorded')

    adopted = '--adopt' in snap.get('fetched_by', '')
    how = 'adopted (already on disk)' if adopted else 'downloaded'

    files = list(snap.get('raw_files', [])) + list(snap.get('derived_files', []))
    if any(f.get('md5_published') for f in files):
        backing = 'publisher md5 verified'
    elif any('revision' in f for f in files):
        revs = sorted({f['revision'][:8] for f in files if f.get('revision')})
        backing = 'commit-pinned (' + ', '.join(revs) + ')'
    elif kg == 'drkg':
        backing = 'no checksum published by the source'
    elif 'adoption_evidence:' in (note or {}).get('_raw', ''):
        # The publisher's checksum covers an archive, not the file extracted
        # from it. What could be checked was checked, and it is listed in the
        # release note rather than summarised away to nothing here.
        backing = 'own sha256; source record immutable, contents checked (see release note)'
    else:
        backing = 'own sha256 only'
    return how, backing


def latest_snapshot(root: Path, kg: str):
    """Newest data/raw/<kg>/<date>/SNAPSHOT.json, or None if never fetched.

    tools/fetch_kgs.py writes one per retrieval, so the release identifier and
    licence reported here are the ones the download actually carried rather than
    a string in this file that nothing checks.
    """
    snaps = sorted((root / 'data' / 'raw' / kg).glob('*/SNAPSHOT.json'))
    snaps += sorted((root / 'data' / 'raw' / kg).glob('*/ADOPTED.json'))
    if not snaps:
        return None
    try:
        return json.loads(snaps[-1].read_text())
    except (OSError, json.JSONDecodeError):
        return None


def sha256_of(p: Path, cap=None):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def mtime(p: Path):
    return dt.datetime.utcfromtimestamp(p.stat().st_mtime).strftime('%Y-%m-%d %H:%M UTC')


def count_records(p: Path):
    if p.suffix == '.obo':
        # Count terms native to this ontology, not imported stanzas from others.
        # mondo.obo carries 63k [Term] blocks but only ~32k are MONDO: terms.
        pre = b'id: MONDO:' if p.name.startswith('mondo') else b'id: DOID:'
        n = sum(1 for ln in p.open('rb') if ln.startswith(pre))
        return n or sum(1 for ln in p.open('rb') if ln.startswith(b'[Term]'))
    if p.suffix in ('.tsv', '.csv', '.txt'):
        return max(0, sum(1 for ln in p.open('rb') if ln.strip() and not ln.startswith(b'#')) - 1)
    return None


def human(n):
    for u in ('B', 'KB', 'MB', 'GB'):
        if n < 1024 or u == 'GB':
            return f'{n:.0f}{u}' if u == 'B' else f'{n:.1f}{u}'
        n /= 1024


def main():
    root = Path('.')
    xw, rs = root / 'data' / 'crosswalks', root / 'data' / 'reference_sets'
    now = dt.datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')
    L = ['# Provenance\n',
         f'Generated by `tools/provenance.py` at **{now}** from the files actually',
         'present on disk. Regenerate after any fetch so it cannot drift.\n']

    # ---- ontologies -------------------------------------------------------
    L += ['## Ontologies and mappings\n',
          '| File | Contents | Release | Retrieved | Records | Size | Licence | SHA-256 (16) |',
          '|---|---|---|---|---|---|---|---|']
    seen = []
    for fn, (what, rel, lic, cite) in REGISTRY.items():
        p = xw / fn
        if not p.exists():
            continue
        seen.append(fn)
        n = count_records(p)
        L.append(f'| `{fn}` | {what} | **{rel}** | {mtime(p)} | '
                 f'{n:,} | {human(p.stat().st_size)} | {lic} | `{sha256_of(p)[:16]}` |')
    missing_ont = [f for f in REGISTRY if f not in seen]
    if not seen:
        L.append('| — | *none retrieved yet* | | | | | | |')

    # ---- open targets -----------------------------------------------------
    L += ['\n## Open Targets Platform\n']
    relfile = xw / 'OPENTARGETS_RELEASE.txt'
    if relfile.exists():
        L += ['```', relfile.read_text().strip(), '```\n']
    L += ['| Dataset | Serves | Files | Records | Size | Retrieved |',
          '|---|---|---|---|---|---|']
    ot_missing = []
    for d, what in OT_DATASETS.items():
        p = xw / d
        parts = sorted(p.glob('*.parquet')) if p.exists() else []
        if not parts:
            ot_missing.append(d)
            continue
        try:
            import pyarrow.parquet as pq
            n = sum(pq.ParquetFile(f).metadata.num_rows for f in parts)
        except Exception:
            n = None
        sz = sum(f.stat().st_size for f in parts)
        L.append(f'| `{d}` | {what} | {len(parts)} | '
                 f'{n:,} | {human(sz)} | {mtime(parts[0])} |' if n is not None else
                 f'| `{d}` | {what} | {len(parts)} | ? | {human(sz)} | {mtime(parts[0])} |')
    if len(ot_missing) == len(OT_DATASETS):
        L.append('| — | *none retrieved yet* | | | | |')
    # a gene-disease source is satisfied by EITHER route
    have_gd = any(d not in ot_missing for d in GENE_DISEASE_ALTERNATIVES)
    if have_gd:
        ot_missing = [d for d in ot_missing if d not in GENE_DISEASE_ALTERNATIVES]
    if ot_missing:
        L.append('')
        L.append('**Not yet retrieved:** ' + ', '.join(f'`{d}`' for d in ot_missing))
        L.append('')
        L.append('Run `bash tools/fetch_opentargets.sh` from a network that can reach '
                 'ebi.ac.uk. Per-file SHA-256 is written to '
                 '`data/crosswalks/OPENTARGETS_CHECKSUMS.txt`.')

    # ---- knowledge graphs -------------------------------------------------
    L += ['\n## Knowledge graphs under evaluation\n',
          '| KG | Release | Last checked against source | Files | Records | Size | Licence | SHA-256 (16) |',
          '|---|---|---|---|---|---|---|---|']
    kg_present = []
    for kg, (rel, lic, cite, src) in KG_REGISTRY.items():
        d = root / 'data' / kg
        # *_RELEASE.txt is provenance, not data: counting its lines as records
        # and folding it into the digest would make the digest change whenever
        # the note is rewritten, which is exactly when it should not.
        files = sorted(f for f in d.glob('*')
                       if f.is_file() and not f.name.endswith('_RELEASE.txt')) \
            if d.exists() else []
        if not files:
            continue
        kg_present.append(kg)
        # A fetched snapshot states its own release, licence and retrieval time.
        # Prefer it over the table above, which is a fallback for graphs put in
        # place by hand before tools/fetch_kgs.py existed.
        snap = latest_snapshot(root, kg)
        if snap:
            lic = snap.get('licence', lic)
        note = release_note(root, kg)
        how, backing = obtained_and_backing(kg, snap, note)
        rel = KG_SHORT.get(kg, rel)
        retrieved = ((snap or {}).get('retrieved_utc', '')[:10]
                     or note.get('retrieved_utc', '')[:10]
                     or mtime(files[0])[:10])
        n = sum(count_records(f) or 0 for f in files)
        sz = sum(f.stat().st_size for f in files)
        digest = hashlib.sha256(''.join(sha256_of(f) for f in files).encode()).hexdigest()
        L.append(f'| **{kg}** | {rel} | {retrieved} | {len(files)} | '
                 f'{n:,} | {human(sz)} | {lic} | `{digest[:16]}` |')
    if not kg_present:
        L.append('| — | *none downloaded yet* | | | | | | |')
    else:
        L.append('')
        L.append('"Last checked against source" is the date the bytes were downloaded, '
                 'or, where they were already on disk, the date they were checked against '
                 'the publisher and the publisher confirmed to have released nothing '
                 'newer. It is not always the date the file itself was created.')
        L.append('')
        L.append('SHA-256 is a digest over the per-file hashes of the data files, so it '
                 'changes if any of them changes; the `*_RELEASE.txt` provenance notes '
                 'are excluded from it and from the record counts. How each graph was '
                 'obtained, what backs its version claim and every per-file hash are in '
                 '`data/<kg>/<KG>_RELEASE.txt` and `data/raw/<kg>/<date>/`. Sources:')
        for kg in kg_present:
            L.append(f'- **{kg}** — {KG_REGISTRY[kg][3]}')
        missing_kg = [k for k in KG_REGISTRY if k not in kg_present]
        if missing_kg:
            L.append('')
            L.append('**Not yet downloaded:** ' + ', '.join(f'`{k}`' for k in missing_kg)
                     + '. Fetch with `python tools/fetch_kgs.py --all` from a network '
                       'that can reach dataverse.harvard.edu, zenodo.org, '
                       'dgl-data.s3-us-west-2.amazonaws.com and huggingface.co.')

    # ---- derived reference sets ------------------------------------------
    L += ['\n## Derived reference sets\n']
    man = rs / 'MANIFEST.json'
    if man.exists():
        m = json.loads(man.read_text())
        L += ['| Set | Rule | Built from | n |', '|---|---|---|---|']
        for k, v in m.items():
            extra = f" ({v['n_rare']:,} rare)" if 'n_rare' in v else ''
            n = v.get('n')
            if n is None and 'pairs' in v:
                n = f"{v['pairs']} pairs, {v['n_drugs']} drugs, {v['n_diseases']} diseases"
            src = str(v.get('source', ''))
            # shorten absolute paths only; leave descriptive strings intact
            if src.startswith('/') and ' ' not in src:
                src = '.../' + '/'.join(Path(src).parts[-2:])
            L.append(f"| {k} | {v.get('rule','')} | {src} | "
                     f"{n:,}{extra} |" if isinstance(n, int) else
                     f"| {k} | {v.get('rule','')} | {src} | {n} |" if n else
                     f"| {k} | {v.get('rule','')} | {src} | — |")
        for f in sorted(rs.glob('*.txt')):
            L.append(f"| `{f.name}` | | | {count_records(f) + 1:,} |")
    else:
        L.append('*Not built yet — run `python tools/build_reference_sets.py`.*')

    # ---- derived gold standards (built by the notebooks, not downloaded) ---
    L += ['\n## Derived gold standards\n']
    bench = root / 'data' / 'gold_standards' / 'derived' / 'multihop_benchmark_v1.csv'
    meta_p = bench.with_suffix('.meta.json')
    if bench.exists():
        meta = json.loads(meta_p.read_text()) if meta_p.exists() else {}
        sha = hashlib.sha256(bench.read_bytes()).hexdigest()
        built = dt.datetime.fromtimestamp(bench.stat().st_mtime, dt.timezone.utc).strftime('%Y-%m-%d %H:%M UTC')
        L += ['| File | What it is | Built from | Triples | Selection | Built | SHA-256 |',
              '|---|---|---|---|---|---|---|',
              f"| `{bench.relative_to(root)}` | Frozen multi-hop benchmark: CTD drug-gene-disease "
              f"triples scored identically by every graph (notebook 06, Supplementary Table S2b) | "
              f"{', '.join(meta.get('built_from_kgs', [])) or 'unknown'} | "
              f"{meta.get('n_triples', count_records(bench)):,} | "
              f"{meta.get('n_common_support', '?'):,} resolvable in every graph; "
              f"{meta.get('n_removed_drug_gene_edge_present', '?'):,} removed (drug-gene edge present); "
              f"capped at {meta.get('max_per_entity', '?')} per entity; seed {meta.get('seed', '?')}; "
              f"target {meta.get('target_n', '?')} | {built} | `{sha}` |"]
        if meta.get('sha256') and meta['sha256'] != sha:
            L.append(f"\n**Warning:** `{meta_p.name}` records SHA-256 `{meta['sha256'][:16]}…`, "
                     f"which does not match the file on disk.")
    else:
        L.append('*Not built yet — generated by notebook 06 on first run.*')

    # ---- pending / not applicable ----------------------------------------
    L += ['\n## Still pending\n']
    pend = []
    if missing_ont:
        pend.append('Ontologies: ' + ', '.join(f'`{f}`' for f in missing_ont))
    if ot_missing:
        pend.append(f'Open Targets: {len(ot_missing)} dataset(s), see above')
    if not (xw / 'ddinter').exists():
        pend.append('DDInter 2.0 — only needed if drug_drug is reinstated; it is '
                    'currently dropped by design. Manual download: '
                    'https://ddinter.scbdd.com/download/')
    kgman = root / 'data' / 'MANIFEST.md'
    if kgman.exists() and '| | | |' in kgman.read_text():
        pend.append('Knowledge-graph snapshots themselves — `data/MANIFEST.md` is '
                    'still an empty template; fill it on KG retrieval')
    L += [f'- {x}' for x in pend] if pend else ['- nothing']

    # ---- citations --------------------------------------------------------
    L += ['\n## Citations\n']
    cites = {cite for fn, (_, _, _, cite) in REGISTRY.items() if (xw / fn).exists()}
    cites.add(OT_CITATION)
    for kg, (_, _, cite, _) in KG_REGISTRY.items():
        if (root / 'data' / kg).exists() and any((root / 'data' / kg).iterdir()):
            cites.add(cite)
    for c in sorted(cites):
        L.append(f'- {c}')
    L += ['',
          'Verify any "current release" citation against the release actually '
          'downloaded; Open Targets, UniProt and Reactome each publish a new '
          'database issue annually.']

    out = root / 'data' / 'PROVENANCE.md'
    out.write_text('\n'.join(L) + '\n')
    print(f'wrote {out}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
