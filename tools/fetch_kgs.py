#!/usr/bin/env python3
"""Fetch the four knowledge graphs that are not yet in data/, at pinned releases.

Hetionet and BioKG were retrieved on 2026-09-14 and are recorded in
data/PROVENANCE.md; this script covers the remaining four:

    primekg     Harvard Dataverse   doi:10.7910/DVN/IXA7BM V2.1, md5 verified
    drkg        AWS S3              unversioned, no checksum published anywhere
    openbilink  Zenodo              doi:10.5281/zenodo.5361324 v1.0.1, md5 verified
    matrix      Hugging Face        everycure/kg-{nodes,edges}, tracked at main

The three GitHub repositories (mims-harvard/PrimeKG, gnn4dr/DRKG,
OpenBioLink/OpenBioLink) carry the build code, not the graphs; the URLs above
are where each project's own README sends you for the data. Before each
download the script asks the publisher whether the pinned release is still the
latest and says so either way, so a pin cannot go stale unnoticed.

Each snapshot lands immutably under

    data/raw/<kg>/<YYYY-MM-DD>/

together with a SNAPSHOT.json recording the release identifier, the source URL,
the UTC retrieval time, the publisher's own checksum where one exists, our
sha256, the byte size and the derived files. The files config/kgs.yaml points at
are then written under data/<kg>/, so the raw download is never edited in place
and a rerun can always be checked against the record.

Publisher checksums are verified before anything is extracted. Where a publisher
publishes none (DRKG), the script records its own sha256 and says so rather than
implying a verification that did not happen.

Usage
-----
    python tools/fetch_kgs.py --list               # what would be fetched
    python tools/fetch_kgs.py --kg primekg         # one graph
    python tools/fetch_kgs.py --all                # all four
    python tools/fetch_kgs.py --all --verify-only  # recheck what is on disk

Network: the sources live on dataverse.harvard.edu, dgl-data.s3-us-west-2
.amazonaws.com, zenodo.org and huggingface.co. Run this where those are
reachable. Disk: about 25 GB for MATRIX, 2 GB for the other three together.
Afterwards run `python tools/provenance.py` to regenerate data/PROVENANCE.md.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import shutil
import sys
import tarfile
import time
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

CHUNK = 1 << 20


# --------------------------------------------------------------------------- specs
# Every field here is a fact taken from the publisher's own API on 2026-09-18,
# not from a README that may have drifted from what is served.
SPECS = {
    'primekg': {
        'release': 'PrimeKG V2.1 (Harvard Dataverse, published 2022-05-02)',
        'citation': ('Chandak P, Huang K, Zitnik M. Building a knowledge graph to enable '
                     'precision medicine. Sci Data 2023;10:67. doi:10.1038/s41597-023-01960-3'),
        'licence': 'CC0 1.0',
        'doi': 'doi:10.7910/DVN/IXA7BM',
        'downloads': [{
            # File id 6180620 is kg.csv in version 2.1. The dataset also carries
            # kg_raw / kg_giant / kg_grouped; kg.csv is the released graph the
            # PrimeKG paper describes and the one config/kgs.yaml expects.
            'url': 'https://dataverse.harvard.edu/api/access/datafile/6180620',
            'name': 'kg.csv',
            'md5': 'aac8191d4fbc5bf09cdf8c3c78b4e75f',
            'size': 981751236,
        }],
        'targets': ['data/primekg/primekg.csv'],
    },
    'drkg': {
        'release': 'DRKG, no versioned release; repository gnn4dr/DRKG at commit '
                   'd4bb1974312013c4bd79a13e42c1d9492033f8c7 (2022-04-19)',
        'citation': ('Ioannidis VN, et al. DRKG - Drug Repurposing Knowledge Graph for '
                     'Covid-19. 2020. https://github.com/gnn4dr/DRKG/'),
        'licence': 'Apache-2.0 (per-source licences apply; see the DRKG licences table)',
        'doi': '',
        'downloads': [{
            'url': 'https://dgl-data.s3-us-west-2.amazonaws.com/dataset/DRKG/drkg.tar.gz',
            'name': 'drkg.tar.gz',
            'md5': '',          # DRKG publishes no checksum; ours is recorded instead
            'size': 0,          # and no size either
        }],
        'targets': ['data/drkg/drkg.tsv'],
    },
    'openbilink': {
        'release': 'OpenBioLink2020 v1.0.1 (Zenodo record 5361324, 2021-09-01)',
        'citation': ('Breit A, Ott S, Agibetov A, Samwald M. OpenBioLink: a benchmarking '
                     'framework for large-scale biomedical link prediction. Bioinformatics '
                     '2020;36(13):4097. doi:10.1093/bioinformatics/btaa274'),
        'licence': 'CC-BY 4.0',
        'doi': '10.5281/zenodo.5361324',
        'downloads': [{
            # ALL_DIR, not HQ_DIR: config/kgs.yaml sets hq_filter: true and the
            # loader applies the quality cutoff itself (STRING/STITCH >= 700,
            # Bgee gold/high), which needs the uncut graph as input.
            'url': 'https://zenodo.org/api/records/5361324/files/ALL_DIR.zip/content',
            # Zenodo throttles hard and closes the connection every few MB. These
            # are the same file served elsewhere, tried in order when one stalls;
            # the md5 below is checked whichever one answers, so a mirror can
            # only ever supply the identical bytes or be rejected.
            'mirrors': [
                'https://zenodo.org/records/5361324/files/ALL_DIR.zip?download=1',
                'https://zenodo.org/record/3834052/files/ALL_DIR.zip?download=1',
                'https://samwald.info/res/OpenBioLink_2020/ALL_DIR.zip',
            ],
            'name': 'ALL_DIR.zip',
            'md5': '9c150bdd971c786e7d25c8db12dcc21f',
            'size': 700146163,
        }],
        'targets': ['data/openbilink/edges.csv'],
        # Recorded when this graph is adopted rather than downloaded. The zip's
        # md5 cannot cover a file extracted from it, so these are the checks
        # that can be made on the extracted file itself.
        'adopt_evidence': [
            'Zenodo record 5361324 confirmed latest on 2026-09-18; the record is '
            'immutable and OpenBioLink has published nothing since 2021-09-01',
            'edges.csv holds 28,973,555 rows in 5 tab-separated columns, the '
            'ALL_DIR (no quality cutoff) variant the loader expects, not HQ_DIR',
            'the loader quality cutoff (STRING/STITCH >= 700, Bgee gold/high, '
            'curated sources unconditional) yields 4,986,521 edges, matching the '
            'benchmark size OpenBioLink documents',
            'line endings are CRLF; the loader strips \\r before comparing the '
            'source column, so the cutoff is applied as documented',
        ],
    },
    'matrix': {
        'release': 'Every Cure integrated KG, Hugging Face everycure/kg-nodes and '
                   'everycure/kg-edges, tracked at main (the commit each resolved to '
                   'is appended by the fetch)',
        'citation': 'Every Cure MATRIX. https://huggingface.co/datasets/everycure/',
        'licence': 'CC BY 4.0 per the Every Cure public data release notes (Hugging Face tag: custom); attribution to Every Cure and to upstream sources including RTX-KG2 and ROBOKOP',
        'doi': '',
        # Tracked at main, not pinned: MATRIX is the one graph here that is still
        # being republished, so a fetch should get what Every Cure serves on the
        # day. The commit main resolves to is written into SNAPSHOT.json, which
        # is what makes the run reproducible afterwards; pass --revision to
        # reproduce an earlier one. Observed on 2026-09-18: kg-nodes
        # c9c60b2c55e72ce, kg-edges fcd0585957ab513, both published 2026-05-28.
        'hf': [
            {'repo': 'everycure/kg-nodes', 'revision': 'main', 'pattern': 'data/nodes/*'},
            {'repo': 'everycure/kg-edges', 'revision': 'main', 'pattern': 'data/edges/*'},
        ],
        'downloads': [],
        'targets': ['data/matrix/nodes.tsv', 'data/matrix/edges.tsv'],
    },
}

# The columns biokgsuite.io.loading.load_matrix reads. Writing only these keeps
# the derived TSVs to a fraction of the 28 GB the full parquet schema would make,
# and a column the loader does not read cannot silently change a result.
MATRIX_NODE_COLS = ['id', 'name', 'category', 'equivalent_identifiers',
                    'upstream_data_source']
MATRIX_EDGE_COLS = ['subject', 'predicate', 'object',
                    'primary_knowledge_source', 'aggregator_knowledge_source',
                    'upstream_data_source', 'knowledge_level', 'agent_type',
                    'num_references']


# --------------------------------------------------------------------------- helpers
def _repo_root() -> Path:
    for p in [Path.cwd().resolve(), *Path.cwd().resolve().parents]:
        if (p / 'config' / 'kgs.yaml').exists():
            return p
    sys.exit('Run this from inside the repository (config/kgs.yaml not found).')


def _digests(path: Path) -> tuple[str, str, int]:
    """Return (sha256, md5, size) in one pass over the file."""
    sha, md5, size = hashlib.sha256(), hashlib.md5(), 0
    with open(path, 'rb') as fh:
        while chunk := fh.read(CHUNK):
            sha.update(chunk); md5.update(chunk); size += len(chunk)
    return sha.hexdigest(), md5.hexdigest(), size


def _download(url: str, dest: Path, expect: int = 0, stall_limit: int = 8) -> None:
    """Stream a URL to dest, resuming a partial transfer rather than restarting.

    These are 0.7 to 4 GB files over links that drop, and Zenodo in particular
    closes the connection every few MB. What counts against the retry limit is
    therefore a *stall*, an attempt that moved no bytes at all: an attempt that
    fetched even one byte resets the counter, so a server dribbling the file out
    over a hundred connections still finishes, while one that is truly refusing
    stops quickly. The partial file is kept as <name>.part and renamed into
    place only once the byte count is what the publisher said it would be.
    """
    tmp = dest.with_suffix(dest.suffix + '.part')
    print(f'    GET {url}')
    stalls = 0
    while True:
        have = tmp.stat().st_size if tmp.exists() else 0
        req = urllib.request.Request(url)
        if have:
            req.add_header('Range', f'bytes={have}-')
        try:
            with urllib.request.urlopen(req) as resp:
                resumed = resp.status == 206
                if have and not resumed:
                    # The server ignored the range: start over rather than
                    # append a second copy of the file onto the first.
                    print('    server does not support resume, restarting')
                    have = 0
                total = int(resp.headers.get('Content-Length') or 0) + (have if resumed else 0)
                mode = 'ab' if (have and resumed) else 'wb'
                done = have if resumed else 0
                with open(tmp, mode) as out:
                    while chunk := resp.read(CHUNK):
                        out.write(chunk); done += len(chunk)
                        if total:
                            print(f'\r    {done / 1e6:,.0f} / {total / 1e6:,.0f} MB '
                                  f'({100 * done / total:4.1f}%)', end='', flush=True)
                print()
        except Exception as exc:                     # noqa: BLE001 - any transport fault
            # A proxy policy denial is a decision, not a transport fault, and
            # retrying it just spends time. Report the blocked host and stop.
            if '403' in str(exc) or '407' in str(exc):
                raise SystemExit(
                    f'    {urllib.parse.urlsplit(url).netloc} is refused by the network '
                    f'this is running on ({exc}).\n'
                    '    Run tools/fetch_kgs.py from a machine whose network allows '
                    'dataverse.harvard.edu, zenodo.org,\n'
                    '    dgl-data.s3-us-west-2.amazonaws.com and huggingface.co.')
            print(f'\n    connection dropped ({exc})')

        got = tmp.stat().st_size if tmp.exists() else 0
        if expect and got == expect:
            tmp.rename(dest)
            return
        if not expect and got and got == have:
            # Nothing more is arriving and the publisher declared no size to
            # check against, so this is as complete as it can be shown to be.
            tmp.rename(dest)
            return
        if got > have:
            stalls = 0
        else:
            stalls += 1
            if stalls >= stall_limit:
                raise RuntimeError(
                    f'{dest.name}: stalled at {got:,} bytes after {stalls} attempts that '
                    f'moved nothing'
                    + (f' (publisher says {expect:,})' if expect else '')
                    + f'. The partial file is kept at {tmp.name}; rerun to resume.')
            time.sleep(min(30, 2 ** stalls))
        if expect:
            print(f'    have {got / 1e6:,.0f} of {expect / 1e6:,.0f} MB; resuming')


def _download_any(sources: list, dest: Path, expect: int = 0) -> None:
    """Try each source in turn, keeping the partial file across the switch.

    A mirror is only ever allowed to supply bytes the publisher's own checksum
    then accepts, so falling back cannot quietly change what is benchmarked. The
    .part file carries over, so a mirror resumes where the primary stalled.
    """
    for i, src in enumerate(sources):
        try:
            _download(src, dest, expect=expect)
            return
        except RuntimeError as exc:
            if i == len(sources) - 1:
                raise
            print(f'    {exc}')
            print(f'    trying source {i + 2} of {len(sources)}')


def _check_latest(kg: str) -> None:
    """Say whether the publisher has moved past the release pinned above.

    Informational only. The pinned release is still what gets downloaded, so a
    newer one cannot change the benchmark without somebody deciding to let it;
    what this stops is a pin quietly going stale while the paper claims the
    current release.
    """
    import json as _json
    try:
        if kg == 'primekg':
            with urllib.request.urlopen(
                    'https://dataverse.harvard.edu/api/datasets/:persistentId/'
                    '?persistentId=doi:10.7910/DVN/IXA7BM', timeout=30) as r:
                v = _json.load(r)['data']['latestVersion']
                latest = f"V{v['versionNumber']}.{v['versionMinorNumber']}"
                print(f'    publisher latest: {latest}, released '
                      f"{v['releaseTime'][:10]}"
                      + ('  (matches the pin)' if latest == 'V2.1'
                         else '  <-- NEWER THAN THE PIN, review before citing'))
        elif kg == 'openbilink':
            with urllib.request.urlopen(
                    'https://zenodo.org/api/records/3826918/versions/latest',
                    timeout=30) as r:
                d = _json.load(r)
                rec, ver = str(d['id']), d['metadata'].get('version', '?')
                print(f'    publisher latest: Zenodo record {rec}, version {ver}'
                      + ('  (matches the pin)' if rec == '5361324'
                         else '  <-- NEWER THAN THE PIN, review before citing'))
        elif kg == 'drkg':
            print('    DRKG publishes no version and no checksum; the tarball at the '
                  'URL above is whatever is served today')
    except Exception as exc:                          # noqa: BLE001 - never blocks a fetch
        print(f'    could not check for a newer release ({exc})')


def _member(archive_names, *, endswith: str, must_contain: str = '') -> str:
    """Find exactly one archive member by file name, or fail with the listing.

    ``endswith`` is matched against the member's own file name, not its path,
    and the whole name must match: OpenBioLink's zip holds both
    ``graph_files/edges.csv`` and ``graph_files/TN_edges.csv`` (the true
    negatives), and a suffix match on ``edges.csv`` finds both. The archives
    are looked into rather than assumed: a publisher that moves a file inside
    its own zip should make this stop, not make it read the wrong file.
    """
    # Skip macOS packaging debris. DRKG's published tarball was rolled on a Mac
    # and carries an AppleDouble sidecar for every file ('._drkg.tsv'), which
    # matches any suffix its partner does; OpenBioLink's zip can carry the
    # __MACOSX/ form of the same thing. Neither is data.
    def _is_debris(n: str) -> bool:
        parts = n.split('/')
        return parts[0] == '__MACOSX' or any(q.startswith('._') for q in parts)

    hits = [n for n in archive_names
            if n.rsplit('/', 1)[-1] == endswith and must_contain in n and not _is_debris(n)]
    if len(hits) != 1:
        listing = '\n      '.join(sorted(archive_names)[:40])
        raise RuntimeError(
            f'expected exactly one member named {endswith!r}'
            f'{" under " + must_contain if must_contain else ""}, found {len(hits)}: {hits}\n'
            f'    archive contains:\n      {listing}')
    return hits[0]


# --------------------------------------------------------------------------- derive
def _derive_primekg(raw_dir: Path, root: Path) -> list[Path]:
    dest = root / 'data' / 'primekg' / 'primekg.csv'
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(raw_dir / 'kg.csv', dest)
    return [dest]


def _derive_drkg(raw_dir: Path, root: Path) -> list[Path]:
    dest = root / 'data' / 'drkg' / 'drkg.tsv'
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(raw_dir / 'drkg.tar.gz', 'r:gz') as tf:
        name = _member(tf.getnames(), endswith='drkg.tsv')
        src = tf.extractfile(name)
        if src is None:
            raise RuntimeError(f'{name} is not a regular file in the tarball')
        with open(dest, 'wb') as out:
            shutil.copyfileobj(src, out, CHUNK)
    return [dest]


def _derive_openbilink(raw_dir: Path, root: Path) -> list[Path]:
    dest = root / 'data' / 'openbilink' / 'edges.csv'
    dest.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(raw_dir / 'ALL_DIR.zip') as zf:
        name = _member(zf.namelist(), endswith='edges.csv', must_contain='graph_files')
        with zf.open(name) as src, open(dest, 'wb') as out:
            shutil.copyfileobj(src, out, CHUNK)
    return [dest]


def _is_seq(v) -> bool:
    """True for a list-valued cell, whatever container pyarrow chose for it."""
    return hasattr(v, '__len__') and not isinstance(v, (str, bytes))


def _parquet_to_tsv(parquet_dir: Path, dest: Path, columns: list[str]) -> int:
    """Stream every parquet shard in a directory into one TSV, columns only.

    Row groups are written one at a time so peak memory stays in the hundreds of
    MB rather than the tens of GB the full table would need. List-valued columns
    are written as their Python repr, which is the form load_matrix's regexes
    already read.
    """
    import pyarrow.parquet as pq

    shards = sorted(p for p in parquet_dir.rglob('*.parquet'))
    if not shards:
        raise RuntimeError(f'no .parquet files under {parquet_dir}')
    present = [c for c in columns if c in pq.ParquetFile(shards[0]).schema_arrow.names]
    missing = [c for c in columns if c not in present]
    if missing:
        print(f'    note: columns absent from this release and left out: {missing}')

    dest.parent.mkdir(parents=True, exist_ok=True)
    rows = 0
    with open(dest, 'w', encoding='utf-8', newline='') as out:
        out.write('\t'.join(present) + '\n')
        for shard in shards:
            pf = pq.ParquetFile(shard)
            for batch in pf.iter_batches(batch_size=100_000, columns=present):
                df = batch.to_pandas()
                for col in present:
                    # pyarrow hands list columns back as numpy arrays, whose str()
                    # drops the commas ("['a' 'b']"). Normalise to a real Python
                    # repr so the TSV round-trips to the form load_matrix reads.
                    if df[col].map(_is_seq).any():
                        df[col] = df[col].map(
                            lambda v: repr([str(x) for x in v]) if _is_seq(v) else v)
                df.to_csv(out, sep='\t', header=False, index=False,
                          na_rep='', lineterminator='\n')
                rows += len(df)
            print(f'\r    {dest.name}: {rows:,} rows', end='', flush=True)
    print()
    return rows


def _derive_matrix(raw_dir: Path, root: Path) -> list[Path]:
    out = []
    for sub, dest_name, cols in [('kg-nodes', 'nodes.tsv', MATRIX_NODE_COLS),
                                 ('kg-edges', 'edges.tsv', MATRIX_EDGE_COLS)]:
        dest = root / 'data' / 'matrix' / dest_name
        _parquet_to_tsv(raw_dir / sub, dest, cols)
        out.append(dest)
    return out


DERIVE = {'primekg': _derive_primekg, 'drkg': _derive_drkg,
          'openbilink': _derive_openbilink, 'matrix': _derive_matrix}


# --------------------------------------------------------------------------- fetch
def fetch(kg: str, root: Path, verify_only: bool = False,
          revision_override: str = '') -> None:
    spec = SPECS[kg]
    print(f'\n{kg}  {spec["release"]}')
    if verify_only:
        # Check the newest snapshot that exists, not a directory named after
        # today: verifying is for a download that already happened.
        snaps = sorted(d for d in (root / 'data' / 'raw' / kg).glob('*') if d.is_dir())
        if not snaps:
            print('    no snapshot on disk')
            return
        raw_dir = snaps[-1]
        print(f'    checking {raw_dir.relative_to(root)}')
    else:
        raw_dir = root / 'data' / 'raw' / kg / dt.date.today().isoformat()
        raw_dir.mkdir(parents=True, exist_ok=True)

    if not verify_only:
        _check_latest(kg)

    files = []
    resolved_revisions = {}
    for d in spec['downloads']:
        path = raw_dir / d['name']
        if not path.exists():
            if verify_only:
                print(f'    MISSING {path.relative_to(root)}')
                continue
            _download_any([d['url'], *d.get('mirrors', [])], path, expect=d['size'])
        sha, md5, size = _digests(path)
        if d['md5']:
            if md5 != d['md5']:
                path.unlink()
                raise SystemExit(
                    f'    CHECKSUM MISMATCH for {d["name"]}\n'
                    f'      published md5 {d["md5"]}\n'
                    f'      downloaded    {md5}\n'
                    '    The file has been removed. Do not use it.')
            print(f'    md5 verified against the publisher: {md5}')
        else:
            print(f'    no published checksum; recording ours: sha256 {sha[:16]}')
        if d['size'] and size != d['size']:
            print(f'    warning: size {size:,} != published {d["size"]:,}')
        files.append({'name': d['name'], 'url': d['url'], 'bytes': size,
                      'sha256': sha, 'md5': md5,
                      'md5_published': d['md5'] or None})

    for hf in spec.get('hf', []):
        sub = raw_dir / hf['repo'].split('/')[-1]
        if verify_only:
            if not sub.exists():
                print(f'    MISSING {sub.relative_to(root)}')
                continue
        else:
            try:
                from huggingface_hub import snapshot_download
            except ImportError:
                raise SystemExit(
                    '    huggingface_hub is needed for MATRIX: pip install huggingface_hub')
            revision = revision_override or hf['revision']
            try:
                from huggingface_hub import HfApi
                resolved = HfApi().dataset_info(hf['repo'], revision=revision).sha
            except Exception:                        # noqa: BLE001 - informational only
                resolved = revision
            print(f'    snapshot_download {hf["repo"]} @ {revision} '
                  f'-> commit {resolved[:12]}')
            snapshot_download(repo_id=hf['repo'], repo_type='dataset',
                              revision=revision, local_dir=str(sub),
                              allow_patterns=[hf['pattern']])
            resolved_revisions[hf['repo']] = resolved
        for p in sorted(sub.rglob('*')):
            if p.is_file():
                sha, md5, size = _digests(p)
                files.append({'name': str(p.relative_to(raw_dir)), 'url': hf['repo'],
                              'revision': resolved_revisions.get(hf['repo'], hf['revision']),
                              'bytes': size, 'sha256': sha, 'md5': md5,
                              'md5_published': None})

    if verify_only:
        print(f'    {len(files)} file(s) present, {sum(f["bytes"] for f in files) / 1e9:.2f} GB')
        return

    print('    deriving the files config/kgs.yaml points at')
    derived = DERIVE[kg](raw_dir, root)
    derived_rec = []
    for p in derived:
        sha, _, size = _digests(p)
        derived_rec.append({'path': str(p.relative_to(root)), 'bytes': size, 'sha256': sha})
        print(f'    {p.relative_to(root)}  {size / 1e6:,.0f} MB  sha256 {sha[:16]}')

    snapshot = {
        'kg': kg,
        'release': (spec['release'] + '; commits '
                    + ', '.join(f'{r}@{s[:12]}' for r, s in resolved_revisions.items())
                    if resolved_revisions else spec['release']),
        'doi': spec['doi'],
        'licence': spec['licence'],
        'citation': spec['citation'],
        'retrieved_utc': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds'),
        'fetched_by': 'tools/fetch_kgs.py',
        'raw_files': files,
        'derived_files': derived_rec,
    }
    (raw_dir / 'SNAPSHOT.json').write_text(json.dumps(snapshot, indent=2) + '\n')
    note = _write_release_note(kg, root, snapshot)
    print(f'    recorded {raw_dir.relative_to(root)}/SNAPSHOT.json '
          f'and {note.relative_to(root)}')


# For --adopt: what each file config/kgs.yaml points at actually is. A file that
# IS the publisher's artifact can be checked against the publisher's own
# checksum; one extracted or converted from an archive cannot, and saying so is
# the point of this table.
ADOPTABLE = {
    'primekg':    {'data/primekg/primekg.csv': 'kg.csv'},
    'drkg':       {'data/drkg/drkg.tsv': None},
    'openbilink': {'data/openbilink/edges.csv': None},
    'matrix':     {'data/matrix/nodes.tsv': None, 'data/matrix/edges.tsv': None},
}


def adopt(kg: str, root: Path) -> None:
    """Record provenance for files already sitting in data/<kg>/.

    A copy that verifies byte-for-byte against the publisher's checksum is as
    citable as one downloaded a minute ago, and more citable than a fresh
    download of a resource that publishes no checksum at all. This mode
    establishes which of the two a file on disk is, and never implies the
    stronger claim.
    """
    spec = SPECS[kg]
    published = {d['name']: d for d in spec['downloads']}
    print(f'\n{kg}  {spec["release"]}')
    records, verified, unverifiable = [], [], []
    for rel, artifact in ADOPTABLE[kg].items():
        path = root / rel
        if not path.exists():
            print(f'    MISSING {rel}')
            return
        sha, md5, size = _digests(path)
        pub = published.get(artifact) if artifact else None
        if pub and pub['md5']:
            ok = (md5 == pub['md5'])
            state = ('verified against the publisher'
                     if ok else 'DOES NOT MATCH the published checksum')
            (verified if ok else unverifiable).append(rel)
            print(f'    {rel}  md5 {md5}  {state}')
            if not ok:
                raise SystemExit(
                    f'    {rel} does not match the md5 published for {artifact}.\n'
                    f'      published  {pub["md5"]}\n      on disk    {md5}\n'
                    '    Do not adopt this file; fetch it instead.')
        else:
            state = ('derived from the publisher archive; no publisher checksum '
                     'covers this file')
            unverifiable.append(rel)
            print(f'    {rel}  sha256 {sha[:16]}  {state}')
        records.append({'path': rel, 'bytes': size, 'sha256': sha, 'md5': md5,
                        'publisher_artifact': artifact,
                        'md5_published': (pub or {}).get('md5') or None,
                        'verification': state})

    snapshot = {
        'kg': kg,
        'release': spec['release'],
        'doi': spec['doi'],
        'licence': spec['licence'],
        'citation': spec['citation'],
        'retrieved_utc': dt.datetime.fromtimestamp(
            (root / next(iter(ADOPTABLE[kg]))).stat().st_mtime,
            dt.timezone.utc).isoformat(timespec='seconds'),
        'fetched_by': 'tools/fetch_kgs.py --adopt (file already on disk, not downloaded here)',
        'raw_files': [],
        'derived_files': records,
    }
    raw_dir = root / 'data' / 'raw' / kg / dt.date.today().isoformat()
    raw_dir.mkdir(parents=True, exist_ok=True)
    (raw_dir / 'ADOPTED.json').write_text(json.dumps(snapshot, indent=2) + '\n')
    note = _write_release_note(kg, root, snapshot, adopted=True)
    print(f'    recorded {(raw_dir / "ADOPTED.json").relative_to(root)} '
          f'and {note.relative_to(root)}')
    if unverifiable:
        print('    note: ' + ', '.join(unverifiable) + ' rest on our own hash only. '
              'Run without --adopt to fetch and verify against the publisher.')


def _write_release_note(kg: str, root: Path, snapshot: dict, adopted: bool = False) -> Path:
    """Drop data/<kg>/<KG>_RELEASE.txt beside the derived files.

    Same key=value shape as the notes hetionet and biokg already carry, so all
    six graphs announce their release the same way to anyone reading data/.
    """
    dest = root / 'data' / kg / f'{kg.upper()}_RELEASE.txt'
    lines = [f'dataset={kg}',
             f'release={snapshot["release"]}']
    if adopted:
        lines.append('provenance=adopted: these files were already on disk and were '
                     'recorded, not downloaded, by tools/fetch_kgs.py --adopt')
    if snapshot['doi']:
        lines.append(f'doi={snapshot["doi"]}')
    lines += [f'licence={snapshot["licence"]}',
              f'retrieved_utc={snapshot["retrieved_utc"]}',
              f'fetched_by={snapshot["fetched_by"]}',
              f'citation={snapshot["citation"]}',
              'raw_files:']
    for f in snapshot['raw_files']:
        pub = 'md5 verified' if f.get('md5_published') else 'no published checksum'
        lines.append(f'  {f["name"]:<44s} {f["bytes"]:>14,} bytes  '
                     f'sha256 {f["sha256"][:16]}  ({pub})')
    if adopted and SPECS[kg].get('adopt_evidence'):
        lines.append('adoption_evidence:')
        lines += [f'  - {e}' for e in SPECS[kg]['adopt_evidence']]
    lines.append('derived_files:')
    for f in snapshot['derived_files']:
        note = f'  ({f["verification"]})' if 'verification' in f else ''
        lines.append(f'  {f["path"]:<44s} {f["bytes"]:>14,} bytes  '
                     f'sha256 {f["sha256"][:16]}{note}')
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text('\n'.join(lines) + '\n')
    return dest


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--kg', choices=sorted(SPECS), action='append', default=[])
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--list', action='store_true')
    ap.add_argument('--verify-only', action='store_true',
                    help='rehash a snapshot already on disk without downloading')
    ap.add_argument('--revision', default='',
                    help='for MATRIX only: fetch this Hugging Face revision instead of '
                         'main, to reproduce an earlier snapshot exactly')
    ap.add_argument('--adopt', action='store_true',
                    help='record provenance for files already in data/<kg>/, checking '
                         'them against the publisher checksum where one applies')
    args = ap.parse_args()

    if args.list:
        for name, spec in SPECS.items():
            print(f'{name:12s} {spec["release"]}')
            print(f'{"":12s} licence {spec["licence"]}')
            for d in spec['downloads']:
                size = f'{d["size"] / 1e6:,.0f} MB' if d['size'] else 'size unpublished'
                pub = f'md5 {d["md5"]}' if d['md5'] else 'no published checksum'
                print(f'{"":12s} {d["name"]:14s} {size:>18s}  {pub}')
            for hf in spec.get('hf', []):
                print(f'{"":12s} {hf["repo"]:14s} @ {hf["revision"][:12]}  {hf["pattern"]}')
            print(f'{"":12s} -> {", ".join(spec["targets"])}')
        return

    kgs = sorted(SPECS) if args.all else args.kg
    if not kgs:
        ap.error('pass --kg NAME (repeatable), --all, or --list')
    root = _repo_root()
    for kg in kgs:
        if args.adopt:
            adopt(kg, root)
        else:
            fetch(kg, root, verify_only=args.verify_only,
                  revision_override=args.revision)
    print('\nNow run:  python tools/provenance.py')


if __name__ == '__main__':
    main()
