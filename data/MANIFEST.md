# Data acquisition

**The record of what is actually on disk is `data/PROVENANCE.md`**, regenerated
from the files themselves by `python tools/provenance.py`. It carries, per
graph, the release identifier, whether the bytes were downloaded or recorded in
place, what backs the version claim, the retrieval date, record counts and a
SHA-256 digest. This file is the how-to; that one is the record. Nothing is
transcribed between them, so they cannot disagree.

Per-graph detail sits beside the data: `data/<kg>/<KG>_RELEASE.txt` for a
readable note, `data/raw/<kg>/<YYYY-MM-DD>/SNAPSHOT.json` (or `ADOPTED.json`)
for the machine-readable one with every per-file hash.

## Layout

    data/raw/<kg>/<YYYY-MM-DD>/     the download exactly as published, never edited (primekg, drkg, openbilink, matrix)
    data/<kg>/                      the files config/kgs.yaml reads, with <KG>_RELEASE.txt beside them
    data/crosswalks/                ontologies, ID maps and the Open Targets datasets (CROSSWALKS.md, PROVENANCE.md)
    data/reference_sets/            the coverage dimension's reference sets, built from crosswalks/ (REFERENCE_SETS.md)
    data/gold_standards/            the sets notebooks 06 and 07 read; the post-cutoff set lives in post_cutoff_approvals/ (its README.md, POST_CUTOFF_APPROVALS.md)

Hetionet and BioKG were recorded in place rather than fetched by
`tools/fetch_kgs.py`, so they have no `raw/` folder; their release notes and
per-file hashes are in `data/hetionet/HETIONET_RELEASE.txt` and
`data/biokg/BIOKG_RELEASE.txt`. `data/matrix/_cache_matrix_*.parquet` are
loader caches written by `biokgsuite.io.loading` on first use; they are
derived from `edges.tsv` and `nodes.tsv`, excluded from the provenance digest,
and safe to delete.

The documents in this folder, and what each is for:

| file | what it records |
|---|---|
| `PROVENANCE.md` | generated: release, retrieval date, counts and SHA-256 for every file on disk |
| `MANIFEST.md` | this file: how the data were obtained and what backs each version claim |
| `CROSSWALKS.md` | the identifier crosswalks and the ontology releases they pin |
| `SOURCE_COUNTS.md` | the per-graph source counts behind the trustworthiness dimension |
| `REFERENCE_SETS.md` | the coverage dimension's reference sets: rule, source and independence per set |
| `POST_CUTOFF_APPROVALS.md` | the 128-pair post-cutoff set: source, funnel, checks, files, version log |

`config/kgs.yaml` points at the active snapshot under `data/<kg>/`.

## Fetching

PrimeKG, DRKG, OpenBioLink and MATRIX are handled by `tools/fetch_kgs.py`. It
pins each release, asks the publisher whether anything newer exists, verifies
the publisher's own checksum before unpacking where one is published, and
writes the snapshot and release note:

    python tools/fetch_kgs.py --list          # releases, sizes, checksums
    python tools/fetch_kgs.py --all           # fetch and derive
    python tools/fetch_kgs.py --adopt --all   # record a copy already on disk
    python tools/provenance.py                # regenerate data/PROVENANCE.md

Downloads resume, so a dropped connection costs nothing; rerun the same command.
Sources live on `dataverse.harvard.edu`, `zenodo.org`,
`dgl-data.s3-us-west-2.amazonaws.com` and `huggingface.co`, so run it from a
network that allows those.

### Adopting rather than downloading

`--adopt` records a file already on disk. It checks the file against the
publisher's checksum where the file **is** the published artifact, and says
plainly that it cannot where the file was extracted or converted from an
archive, listing instead what could be checked. A copy that matches the
publisher's checksum is as citable as one downloaded a minute ago.

## What backs each version claim

The four graphs differ in what is actually checkable, and the record says so per
graph rather than implying one standard across all six.

- **PrimeKG** publishes an md5 per file on Dataverse, so the file is verified
  outright. Dataverse also reports the latest version, and V2.1 (2022-05-02)
  is still it.
- **OpenBioLink** publishes an md5 for the zip, which cannot cover a file
  extracted from it. Zenodo records are immutable and record 5361324 is
  confirmed latest, so the version is certain; the extracted file is checked on
  its own terms instead (row count, column layout, and the quality cutoff
  yielding the documented benchmark size). Those checks are listed in
  `data/openbilink/OPENBILINK_RELEASE.txt`.
- **DRKG** publishes no checksum and no version of any kind. A fresh download
  is therefore no more verifiable than a copy, and the record says so rather
  than implying otherwise. Two independent retrievals six months apart returned
  byte-identical files, which is the strongest statement available for it.
- **MATRIX** is the only graph still being republished. It is tracked at `main`
  and identified by the commit each dataset resolved to, recorded at fetch time;
  `--revision <sha>` reproduces an earlier pull.

## Hetionet

Hetionet's edge file is stored in Git LFS, so `raw.githubusercontent.com`
returns a pointer rather than the data. Fetch it from
`media.githubusercontent.com/media/...`. The LFS pointer publishes the expected
SHA-256, so the download is verifiable against the source without trusting the
transfer.
