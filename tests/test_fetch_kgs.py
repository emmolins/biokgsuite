"""Guards on tools/fetch_kgs.py.

The script runs once over about 25 GB on a machine with the right network, so
the cheap failures have to be caught here rather than an hour into a download.
Each test covers a way the derived MATRIX TSVs could silently stop matching what
biokgsuite.io.loading.load_matrix reads.
"""
import importlib.util
import re
import sys
from pathlib import Path

import pandas as pd
import pytest

BASE = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location('fetch_kgs', BASE / 'tools' / 'fetch_kgs.py')
fk = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fk)

pa = pytest.importorskip('pyarrow')
pq = pytest.importorskip('pyarrow.parquet')

# The regexes load_matrix runs over `id + ' ' + equivalent_identifiers`.
DRUGBANK_RE = re.compile(r'\b(DB\d{5,})\b')
MONDO_RE = re.compile(r'MONDO:(\d+)')
NCBIGENE_RE = re.compile(r'NCBIGene:(\d+)')


@pytest.fixture
def nodes_tsv(tmp_path):
    table = pa.table({
        'id': ['DRUGBANK:DB00001', 'MONDO:0005812', 'NCBIGene:1234'],
        'name': ['lepirudin', 'tab\there', 'newline\nhere'],
        'category': ['biolink:SmallMolecule', 'biolink:Disease', 'biolink:Gene'],
        # A column the loader does not read, to prove it is dropped.
        'description': ['x', 'y', 'z'],
        'equivalent_identifiers': [['DRUGBANK:DB00001', 'CHEBI:1'],
                                   ['MONDO:0005812', 'DOID:8469'],
                                   ['NCBIGene:1234']],
        'upstream_data_source': [['rtxkg2', 'robokop'], ['robokop'], []],
    })
    src = tmp_path / 'kg-nodes' / 'data' / 'nodes'
    src.mkdir(parents=True)
    pq.write_table(table, src / 'train-00000.parquet')
    dest = tmp_path / 'nodes.tsv'
    fk._parquet_to_tsv(tmp_path / 'kg-nodes', dest, fk.MATRIX_NODE_COLS)
    return dest


def _read(dest, cols):
    return pd.read_csv(dest, sep='\t', dtype=str,
                       usecols=lambda c: c in cols).fillna('')


def test_only_the_columns_the_loader_reads_are_written(nodes_tsv):
    got = _read(nodes_tsv, fk.MATRIX_NODE_COLS)
    assert list(got.columns) == fk.MATRIX_NODE_COLS
    assert 'description' not in got.columns


def test_list_columns_round_trip_as_python_repr(nodes_tsv):
    got = _read(nodes_tsv, fk.MATRIX_NODE_COLS)
    # pyarrow hands these back as numpy arrays, whose str() drops the commas.
    assert got.loc[0, 'equivalent_identifiers'] == "['DRUGBANK:DB00001', 'CHEBI:1']"
    assert got.loc[2, 'equivalent_identifiers'] == "['NCBIGene:1234']"


def test_embedded_tabs_and_newlines_survive(nodes_tsv):
    got = _read(nodes_tsv, fk.MATRIX_NODE_COLS)
    assert got.loc[1, 'name'] == 'tab\there'
    assert got.loc[2, 'name'] == 'newline\nhere'


def test_curie_regexes_still_match_the_written_form(nodes_tsv):
    got = _read(nodes_tsv, fk.MATRIX_NODE_COLS)
    for i, rx in [(0, DRUGBANK_RE), (1, MONDO_RE), (2, NCBIGENE_RE)]:
        haystack = got.loc[i, 'id'] + ' ' + got.loc[i, 'equivalent_identifiers']
        assert rx.search(haystack), haystack


def test_upstream_data_source_parses_to_pipe_form(nodes_tsv):
    got = _read(nodes_tsv, fk.MATRIX_NODE_COLS)
    parsed = '|'.join(re.findall(r"[a-zA-Z][\w-]*", got.loc[0, 'upstream_data_source']))
    assert parsed == 'rtxkg2|robokop'


def test_edge_columns(tmp_path):
    table = pa.table({
        'subject': ['DRUGBANK:DB00001'], 'predicate': ['biolink:treats'],
        'object': ['MONDO:0005812'], 'primary_knowledge_source': ['infores:drugbank'],
        'aggregator_knowledge_source': [['infores:rtx-kg2']],
        'upstream_data_source': [['rtxkg2']], 'knowledge_level': ['knowledge_assertion'],
        'agent_type': ['manual_agent'], 'num_references': [3], 'affinity': [0.5],
    })
    src = tmp_path / 'kg-edges' / 'data' / 'edges'
    src.mkdir(parents=True)
    pq.write_table(table, src / 'train-00000.parquet')
    dest = tmp_path / 'edges.tsv'
    fk._parquet_to_tsv(tmp_path / 'kg-edges', dest, fk.MATRIX_EDGE_COLS)
    got = _read(dest, fk.MATRIX_EDGE_COLS)
    assert list(got.columns) == fk.MATRIX_EDGE_COLS
    assert 'affinity' not in got.columns
    # The 'biolink:' prefix is left on; load_matrix strips it itself.
    assert got.loc[0, 'predicate'] == 'biolink:treats'


def test_missing_parquet_is_an_error_not_an_empty_file(tmp_path):
    (tmp_path / 'empty').mkdir()
    with pytest.raises(RuntimeError, match='no .parquet'):
        fk._parquet_to_tsv(tmp_path / 'empty', tmp_path / 'out.tsv', ['id'])


def test_archive_member_is_found_not_assumed():
    names = ['ALL_DIR/graph_files/edges.csv', 'ALL_DIR/graph_files/nodes.csv',
             'ALL_DIR/train_test_data/test_edges.csv']
    assert fk._member(names, endswith='edges.csv',
                      must_contain='graph_files') == 'ALL_DIR/graph_files/edges.csv'


def test_ambiguous_archive_member_raises():
    with pytest.raises(RuntimeError):
        fk._member(['a/edges.csv', 'b/edges.csv'], endswith='edges.csv')


def test_absent_archive_member_raises():
    with pytest.raises(RuntimeError):
        fk._member(['a/nodes.csv'], endswith='edges.csv')


def test_every_spec_declares_a_release_licence_and_target():
    for kg, spec in fk.SPECS.items():
        assert spec['release'], kg
        assert spec['licence'], kg
        assert spec['targets'], kg
        assert spec['downloads'] or spec.get('hf'), kg


def test_published_checksums_are_recorded_where_they_exist():
    # DRKG publishes none; the other two file-based sources do, and a silent
    # loss of one would turn a verified download into an unverified one.
    assert fk.SPECS['primekg']['downloads'][0]['md5']
    assert fk.SPECS['openbilink']['downloads'][0]['md5']
    assert fk.SPECS['drkg']['downloads'][0]['md5'] == ''


def test_derive_drkg_extracts_the_tsv(tmp_path):
    import io, tarfile
    raw = tmp_path / 'raw'
    raw.mkdir()
    payload = b'Compound::DB00001\tDRUGBANK::treats\tDisease::DOID:8469\n'
    with tarfile.open(raw / 'drkg.tar.gz', 'w:gz') as tf:
        for name, data in [('drkg.tsv', payload), ('embed/entities.tsv', b'x\n')]:
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tf.addfile(info, io.BytesIO(data))
    (out,) = fk._derive_drkg(raw, tmp_path)
    assert out == tmp_path / 'data' / 'drkg' / 'drkg.tsv'
    assert out.read_bytes() == payload


def test_derive_openbilink_takes_the_graph_files_edges(tmp_path):
    import zipfile
    raw = tmp_path / 'raw'
    raw.mkdir()
    wanted = b'NCBIGENE:1\tGENE_GENE\tNCBIGENE:2\t900\tSTRING\n'
    with zipfile.ZipFile(raw / 'ALL_DIR.zip', 'w') as zf:
        zf.writestr('ALL_DIR/graph_files/edges.csv', wanted)
        zf.writestr('ALL_DIR/graph_files/nodes.csv', b'NCBIGENE:1\tGENE\n')
        # The split directory carries files with confusable names.
        zf.writestr('ALL_DIR/train_test_data/train_sample.csv', b'x\n')
    (out,) = fk._derive_openbilink(raw, tmp_path)
    assert out == tmp_path / 'data' / 'openbilink' / 'edges.csv'
    assert out.read_bytes() == wanted


def test_release_note_names_every_raw_and_derived_file(tmp_path):
    snap = {'release': 'R', 'doi': 'd', 'licence': 'L', 'citation': 'C',
            'retrieved_utc': '2026-09-18T10:00:00+00:00', 'fetched_by': 'tools/fetch_kgs.py',
            'raw_files': [{'name': 'kg.csv', 'bytes': 10, 'sha256': 'a' * 64,
                           'md5_published': 'x'}],
            'derived_files': [{'path': 'data/primekg/primekg.csv', 'bytes': 10,
                               'sha256': 'b' * 64}]}
    note = fk._write_release_note('primekg', tmp_path, snap)
    text = note.read_text()
    assert note.name == 'PRIMEKG_RELEASE.txt'
    for key in ('dataset=primekg', 'release=R', 'licence=L', 'retrieved_utc='):
        assert key in text
    assert 'kg.csv' in text and 'md5 verified' in text
    assert 'data/primekg/primekg.csv' in text


def _fake_primekg(tmp_path, content: bytes):
    (tmp_path / 'data' / 'primekg').mkdir(parents=True)
    (tmp_path / 'data' / 'primekg' / 'primekg.csv').write_bytes(content)


def test_adopt_accepts_a_file_matching_the_published_checksum(tmp_path, monkeypatch, capsys):
    import hashlib
    payload = b'relation,x_id,y_id\n'
    _fake_primekg(tmp_path, payload)
    monkeypatch.setitem(fk.SPECS['primekg']['downloads'][0], 'md5',
                        hashlib.md5(payload).hexdigest())
    fk.adopt('primekg', tmp_path)
    out = capsys.readouterr().out
    assert 'verified against the publisher' in out
    note = (tmp_path / 'data' / 'primekg' / 'PRIMEKG_RELEASE.txt').read_text()
    assert 'provenance=adopted' in note
    assert (tmp_path / 'data' / 'raw' / 'primekg').exists()


def test_adopt_refuses_a_file_that_does_not_match(tmp_path, monkeypatch):
    _fake_primekg(tmp_path, b'wrong bytes\n')
    monkeypatch.setitem(fk.SPECS['primekg']['downloads'][0], 'md5', '0' * 32)
    with pytest.raises(SystemExit, match='does not match'):
        fk.adopt('primekg', tmp_path)


def test_adopt_never_claims_verification_for_a_derived_file(tmp_path, capsys):
    (tmp_path / 'data' / 'drkg').mkdir(parents=True)
    (tmp_path / 'data' / 'drkg' / 'drkg.tsv').write_bytes(b'a\tb\tc\n')
    fk.adopt('drkg', tmp_path)
    out = capsys.readouterr().out
    assert 'no publisher checksum covers this file' in out
    assert 'verified against the publisher' not in out


def test_every_adoptable_target_is_a_declared_target():
    for kg, files in fk.ADOPTABLE.items():
        assert set(files) == set(fk.SPECS[kg]['targets']), kg


# ---------------------------------------------------------------- downloads
# These run against a local server rather than the real hosts. The resume path
# is the one that must not be got wrong: a server that ignores a Range header
# and is appended to anyway produces a file that is the right size for a while
# and wrong for ever, and the natural response to a download that keeps failing
# is to switch the checksum off.
import http.server  # noqa: E402
import socketserver  # noqa: E402
import threading  # noqa: E402

PAYLOAD = bytes(range(256)) * 4000          # 1,024,000 bytes


@pytest.fixture
def server(monkeypatch):
    mode = {'ignore_range': False, 'cut_after': 0}

    class Handler(http.server.BaseHTTPRequestHandler):
        def log_message(self, *a):
            pass

        def do_GET(self):
            rng = self.headers.get('Range')
            start = 0
            if rng and not mode['ignore_range']:
                start = int(rng.split('=')[1].split('-')[0])
                self.send_response(206)
                self.send_header('Content-Range',
                                 f'bytes {start}-{len(PAYLOAD) - 1}/{len(PAYLOAD)}')
            else:
                self.send_response(200)
            body = PAYLOAD[start:]
            if mode['cut_after'] > 0:
                body = body[:mode['cut_after']]
            elif mode['cut_after'] < 0:
                body = b''                       # a server that moves nothing
            self.send_header('Content-Length', str(len(PAYLOAD) - start))
            self.end_headers()
            try:
                self.wfile.write(body)
            except OSError:
                pass

    srv = socketserver.TCPServer(('127.0.0.1', 0), Handler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    for var in ('HTTPS_PROXY', 'https_proxy', 'HTTP_PROXY', 'http_proxy'):
        monkeypatch.delenv(var, raising=False)
    yield f'http://127.0.0.1:{srv.server_address[1]}/file.bin', mode
    srv.shutdown()


def test_download_writes_the_whole_file(server, tmp_path):
    url, _ = server
    fk._download(url, tmp_path / 'a.bin', expect=len(PAYLOAD))
    assert (tmp_path / 'a.bin').read_bytes() == PAYLOAD


def test_download_resumes_a_partial_transfer(server, tmp_path):
    url, _ = server
    (tmp_path / 'b.bin.part').write_bytes(PAYLOAD[:300_000])
    fk._download(url, tmp_path / 'b.bin', expect=len(PAYLOAD))
    assert (tmp_path / 'b.bin').read_bytes() == PAYLOAD


def test_download_restarts_when_the_server_ignores_range(server, tmp_path):
    url, mode = server
    mode['ignore_range'] = True
    (tmp_path / 'c.bin.part').write_bytes(PAYLOAD[:300_000])
    fk._download(url, tmp_path / 'c.bin', expect=len(PAYLOAD))
    # Appending would have produced 1,324,000 bytes of corrupt file.
    assert (tmp_path / 'c.bin').read_bytes() == PAYLOAD


def test_download_completes_across_repeated_cutoffs(server, tmp_path):
    """Zenodo closes the connection every few MB; the file still has to arrive.

    Counting attempts rather than stalls is what broke here: each reconnection
    made progress but was charged against the retry budget, so a 700 MB file
    gave up at about 30 MB.
    """
    url, mode = server
    mode['cut_after'] = 400_000                  # forces three round trips
    fk._download(url, tmp_path / 'd.bin', expect=len(PAYLOAD))
    assert (tmp_path / 'd.bin').read_bytes() == PAYLOAD


def test_download_gives_up_on_a_server_that_moves_nothing(server, tmp_path):
    url, mode = server
    mode['cut_after'] = -1                       # sends an empty body every time
    with pytest.raises(RuntimeError, match='moved nothing'):
        fk._download(url, tmp_path / 'e.bin', expect=len(PAYLOAD), stall_limit=3)
    assert not (tmp_path / 'e.bin').exists()
    assert (tmp_path / 'e.bin.part').exists(), 'the partial file must survive for a resume'


def test_a_short_download_is_never_published(server, tmp_path):
    url, mode = server
    mode['cut_after'] = -1
    with pytest.raises(RuntimeError):
        fk._download(url, tmp_path / 'f.bin', expect=len(PAYLOAD), stall_limit=2)
    assert not (tmp_path / 'f.bin').exists()
    mode['cut_after'] = 0
    fk._download(url, tmp_path / 'f.bin', expect=len(PAYLOAD))
    assert (tmp_path / 'f.bin').read_bytes() == PAYLOAD


def test_applellouble_sidecars_are_not_mistaken_for_data():
    """DRKG's published tarball was rolled on a Mac and ships '._drkg.tsv'."""
    names = ['._drkg.tsv', 'drkg.tsv', '._embed', 'embed/entities.tsv',
             'embed/._entities.tsv', 'entity2src.tsv']
    assert fk._member(names, endswith='drkg.tsv') == 'drkg.tsv'


def test_macosx_folder_is_not_mistaken_for_data():
    names = ['__MACOSX/ALL_DIR/graph_files/edges.csv', 'ALL_DIR/graph_files/edges.csv']
    assert (fk._member(names, endswith='edges.csv', must_contain='graph_files')
            == 'ALL_DIR/graph_files/edges.csv')


def test_a_mirror_takes_over_when_the_primary_stalls(server, tmp_path, monkeypatch):
    """Zenodo throttles; the same file is served elsewhere and must be tried."""
    good, mode = server
    mode['cut_after'] = 0
    calls = []

    real = fk._download

    def flaky(url, dest, expect=0, stall_limit=8):
        calls.append(url)
        if url == 'https://primary.invalid/ALL_DIR.zip':
            raise RuntimeError('ALL_DIR.zip: stalled at 0 bytes, moved nothing')
        return real(url, dest, expect=expect, stall_limit=stall_limit)

    monkeypatch.setattr(fk, '_download', flaky)
    fk._download_any(['https://primary.invalid/ALL_DIR.zip', good],
                     tmp_path / 'z.bin', expect=len(PAYLOAD))
    assert len(calls) == 2
    assert (tmp_path / 'z.bin').read_bytes() == PAYLOAD


def test_all_sources_failing_still_raises(tmp_path, monkeypatch):
    monkeypatch.setattr(fk, '_download',
                        lambda *a, **k: (_ for _ in ()).throw(RuntimeError('stalled')))
    with pytest.raises(RuntimeError, match='stalled'):
        fk._download_any(['https://a.invalid/x', 'https://b.invalid/x'], tmp_path / 'x')


def test_openbiolink_declares_mirrors_and_one_checksum_for_all_of_them():
    d = fk.SPECS['openbilink']['downloads'][0]
    assert len(d['mirrors']) >= 2
    assert all(m.startswith('https://') for m in d['mirrors'])
    assert d['md5'] == '9c150bdd971c786e7d25c8db12dcc21f'
