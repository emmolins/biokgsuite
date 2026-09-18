#!/usr/bin/env python3
"""Convert Hetionet's hetnet JSON into the flat edge/node tables the pipeline reads.

WHY. Hetionet is distributed in four formats. The TSV export carries only
(source, metaedge, target); the JSON and Neo4j formats additionally carry node
and edge PROPERTIES, and Hetionet's own README names JSON and Neo4j as the
recommended formats for that reason. The dropped property matters here:
describe/definitions.json defines

    properties.edges = {"source": "The database from which the edge and its
                        properties were gathered"}

which is the per-edge provenance the trustworthiness dimension scores. Reading
the TSV made Hetionet look like it had none.

The JSON holds 2.25M edges, which will not fit in memory as parsed dicts on a
small machine, so this streams with ijson and writes a flat TSV once. The
notebooks then read that, cheaply, as before.
"""
from __future__ import annotations

import argparse
import bz2
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import ijson


def sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open('rb') as fh:
        for block in iter(lambda: fh.read(chunk), b''):
            h.update(block)
    return h.hexdigest()


def metaedge_abbrev(k2a: dict, source_kind: str, edge_kind: str,
                    target_kind: str, direction: str) -> str:
    """Rebuild the metaedge abbreviation the TSV export uses, e.g. `GpBP`, `Gr>G`.

    The JSON stores an edge's `kind` (`regulates`) but not its metaedge
    abbreviation, and `kind` alone is ambiguous: Anatomy, Compound and Disease
    all have `upregulates` edges to Gene, which the TSV distinguishes as AuG,
    CuG and DuG. The abbreviation is source-abbrev + kind-abbrev + direction
    marker + target-abbrev, verified to reproduce all 24 of the TSV's metaedges
    exactly.
    """
    arrow = {'forward': '>', 'backward': '<'}.get(direction, '')
    return f'{k2a[source_kind]}{k2a[edge_kind]}{arrow}{k2a[target_kind]}'


def convert(src: Path, out_dir: Path) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    stats = {'nodes': 0, 'edges': 0, 'edges_with_source': 0,
             'edges_unbiased': 0, 'sources': {}}

    with bz2.open(src, 'rb') as fh:
        k2a = next(ijson.items(fh, 'kind_to_abbrev'))

    # ---- nodes: id, name, kind  (same columns as the TSV export)
    with bz2.open(src, 'rb') as fh, (out_dir / 'nodes.tsv').open('w', newline='') as o:
        w = csv.writer(o, delimiter='\t', lineterminator='\n')
        w.writerow(['id', 'name', 'kind'])
        for n in ijson.items(fh, 'nodes.item'):
            w.writerow([f"{n['kind']}::{n['identifier']}", n.get('name', ''), n['kind']])
            stats['nodes'] += 1

    # ---- edges: the TSV columns PLUS the per-edge properties the TSV drops
    with bz2.open(src, 'rb') as fh, (out_dir / 'edges.tsv').open('w', newline='') as o:
        w = csv.writer(o, delimiter='\t', lineterminator='\n')
        w.writerow(['source', 'metaedge', 'target', 'edge_source', 'unbiased'])
        for e in ijson.items(fh, 'edges.item'):
            sk, si = e['source_id']
            tk, ti = e['target_id']
            data = e.get('data') or {}
            es = data.get('source') or ''
            metaedge = metaedge_abbrev(k2a, sk, e['kind'], tk, e.get('direction', 'both'))
            w.writerow([f'{sk}::{si}', metaedge, f'{tk}::{ti}', es,
                        '' if data.get('unbiased') is None else int(bool(data['unbiased']))])
            stats['edges'] += 1
            if es:
                stats['edges_with_source'] += 1
                stats['sources'][es] = stats['sources'].get(es, 0) + 1
            if data.get('unbiased'):
                stats['edges_unbiased'] += 1
    return stats


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--src', default='data/hetionet/hetionet-v1.0.json.bz2')
    ap.add_argument('--out', default='data/hetionet')
    a = ap.parse_args()

    src, out_dir = Path(a.src), Path(a.out)
    if not src.exists():
        raise SystemExit(
            f'{src} not found. Download it with:\n'
            f'  curl -sL "https://media.githubusercontent.com/media/hetio/hetionet/'
            f'main/hetnet/json/hetionet-v1.0.json.bz2" -o {src}')

    digest = sha256(src)
    stats = convert(src, out_dir)
    pct = 100 * stats['edges_with_source'] / stats['edges'] if stats['edges'] else 0

    (out_dir / 'HETIONET_RELEASE.txt').write_text(f"""dataset=Hetionet
release_tag=v1.0
format=hetnet JSON, converted to TSV by tools/convert_hetionet_json.py
source_url=https://media.githubusercontent.com/media/hetio/hetionet/main/hetnet/json/hetionet-v1.0.json.bz2
retrieved_utc={datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')}
sha256={digest}
bytes={src.stat().st_size}
nodes={stats['nodes']}
edges={stats['edges']}
edges_with_edge_source={stats['edges_with_source']} ({pct:.1f}%)
edges_flagged_unbiased={stats['edges_unbiased']}
distinct_edge_sources={len(stats['sources'])}
note=The previously used TSV export carries only (source, metaedge, target). The
  JSON additionally carries per-edge properties, which is why Hetionet's README
  names JSON and Neo4j the recommended formats. describe/definitions.json defines
  properties.edges as {{"source": "The database from which the edge and its
  properties were gathered"}}. That field is the per-edge provenance scored by
  the trustworthiness dimension; reading the TSV made Hetionet score zero for it.
""")
    (out_dir / 'edge_sources.json').write_text(
        json.dumps(dict(sorted(stats['sources'].items(), key=lambda kv: -kv[1])), indent=1))
    print(f"nodes {stats['nodes']:,} | edges {stats['edges']:,} | "
          f"with edge_source {stats['edges_with_source']:,} ({pct:.1f}%) | "
          f"{len(stats['sources'])} distinct sources")


if __name__ == '__main__':
    main()
