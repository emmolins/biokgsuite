#!/usr/bin/env python3
"""Download the Every Cure MATRIX KG and convert to the TSV format
the BioKGSuite loader expects.

The MATRIX dataset on HuggingFace was restructured in 2025 from a single
`everycure/matrix-kg` repo (with nodes.tsv + edges.tsv) into two separate
sharded-parquet datasets (`everycure/kg-nodes`, `everycure/kg-edges`).
This script downloads both, concatenates the shards, and writes them as
the TSV files at `data/matrix/{nodes,edges}.tsv` that
`src/loading.py :: load_matrix` reads.

Usage (from repo root, with venv activated):
    python scripts/hpc/download_matrix.py

Wall time on HPC fast network: 10-20 min. RAM peak: ~30 GB during the
edge conversion (streaming to keep it bounded).
"""
import sys
from pathlib import Path

from huggingface_hub import snapshot_download
import pandas as pd
import pyarrow.parquet as pq

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
MATRIX_DIR = REPO_ROOT / 'data' / 'matrix'
MATRIX_DIR.mkdir(parents=True, exist_ok=True)


def main():
    print('=' * 60)
    print('Downloading MATRIX from HuggingFace and converting to TSV')
    print(f'  target dir: {MATRIX_DIR}')
    print('=' * 60)

    # 1. Download parquet shards
    print('\n[1/3] Downloading nodes (11 shards)...')
    snapshot_download(
        repo_id='everycure/kg-nodes', repo_type='dataset',
        local_dir=str(MATRIX_DIR / '_dl_nodes'),
        allow_patterns=['data/nodes/*.parquet'],
    )

    print('\n[2/3] Downloading edges (36 shards)...')
    snapshot_download(
        repo_id='everycure/kg-edges', repo_type='dataset',
        local_dir=str(MATRIX_DIR / '_dl_edges'),
        allow_patterns=['data/edges/*.parquet'],
    )

    # 2. Nodes: small enough to do in one shot
    print('\n[3/3] Converting nodes -> nodes.tsv ...')
    node_files = sorted((MATRIX_DIR / '_dl_nodes' / 'data' / 'nodes').glob('*.parquet'))
    nodes_df = pd.concat([pd.read_parquet(f) for f in node_files], ignore_index=True)
    print(f'  {len(nodes_df):,} nodes')
    print(f'  columns: {list(nodes_df.columns)}')
    nodes_path = MATRIX_DIR / 'nodes.tsv'
    nodes_df.to_csv(nodes_path, sep='\t', index=False)
    print(f'  wrote {nodes_path} ({nodes_path.stat().st_size / 1e9:.2f} GB)')

    # 3. Edges: stream shard-by-shard so we don't blow 125 GB of RAM
    print('\nConverting edges -> edges.tsv (streaming)...')
    edge_files = sorted((MATRIX_DIR / '_dl_edges' / 'data' / 'edges').glob('*.parquet'))
    edges_path = MATRIX_DIR / 'edges.tsv'
    header_written = False
    total_rows = 0
    with open(edges_path, 'w') as out:
        for i, f in enumerate(edge_files, 1):
            df = pq.read_table(f).to_pandas()
            df.to_csv(out, sep='\t', index=False,
                      header=not header_written, mode='a')
            header_written = True
            total_rows += len(df)
            print(f'  shard {i:2d}/{len(edge_files)}: '
                  f'+{len(df):>10,} rows  (running total {total_rows:,})',
                  flush=True)
    print(f'  wrote {edges_path} ({edges_path.stat().st_size / 1e9:.2f} GB)')

    print('\n' + '=' * 60)
    print('Done. Inspect the schemas with:')
    print(f'  head -1 {nodes_path}')
    print(f'  head -1 {edges_path}')
    print('\nDelete the intermediate download dirs to save space:')
    print(f'  rm -rf {MATRIX_DIR}/_dl_nodes {MATRIX_DIR}/_dl_edges')
    print('=' * 60)


if __name__ == '__main__':
    sys.exit(main())
