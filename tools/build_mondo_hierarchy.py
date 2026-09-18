#!/usr/bin/env python3
"""Rebuild data/gold_standards/bridges/mondo_hierarchy.csv from mondo.obo.

The file is a MONDO `is_a` edge list (child_id, parent_id) used by
src/biokgsuite/gold/ontology_classifier.py to walk a
disease up to a therapeutic-domain ancestor. It was absent from
data/gold_standards/, so the cross-domain section of notebook 07 could not run.

Unlike fda_drugs_id_mapping.csv this needs no hand curation: it is a mechanical
projection of the ontology, so it is rebuilt here rather than reconstructed.
Only `is_a` is followed. Obsolete terms are skipped, since classifying a
disease through a deprecated ancestor would assign a domain the ontology no
longer asserts.

Writes the CSV and a provenance sidecar recording the source file's bytes and
SHA-256 and the MONDO release stated in the obo header.

    python tools/build_mondo_hierarchy.py
"""
import csv
import hashlib
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / 'data' / 'crosswalks' / 'mondo.obo'
OUT = ROOT / 'data' / 'gold_standards' / 'bridges' / 'mondo_hierarchy.csv'
PROV = ROOT / 'data' / 'gold_standards' / 'bridges' / 'MONDO_HIERARCHY_PROVENANCE.txt'


def main():
    if not SRC.exists():
        sys.exit(f'Source ontology not found: {SRC}')

    edges, n_terms, n_obsolete = [], 0, 0
    data_version = ''
    cur_id, cur_parents, cur_obsolete = None, [], False

    def flush():
        nonlocal n_terms, n_obsolete
        if cur_id is None:
            return
        n_terms += 1
        if cur_obsolete:
            n_obsolete += 1
            return
        for p in cur_parents:
            edges.append((cur_id, p))

    with open(SRC, encoding='utf-8') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('data-version:'):
                data_version = line.split(':', 1)[1].strip()
            elif line == '[Term]':
                flush()
                cur_id, cur_parents, cur_obsolete = None, [], False
            elif line.startswith('id: '):
                cur_id = line[4:].strip()
            elif line.startswith('is_a: '):
                # 'is_a: MONDO:0003274 {source="DOID:1612"} ! thoracic cancer'
                # The trailing comment AND the optional {..} qualifier block
                # both have to go. Keeping the qualifier leaves parent IDs that
                # match no child key, so every is_a chain dies after one hop and
                # no term ever reaches a domain ancestor.
                parent = line[6:].split('!', 1)[0].split('{', 1)[0].strip()
                if parent:
                    cur_parents.append(parent)
            elif line.startswith('is_obsolete: true'):
                cur_obsolete = True
    flush()

    edges = sorted(set(edges))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    tmp = OUT.with_suffix('.csv.tmp')
    with open(tmp, 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh)
        w.writerow(['child_id', 'parent_id'])
        w.writerows(edges)
    tmp.replace(OUT)

    digest = hashlib.sha256(SRC.read_bytes()).hexdigest()
    PROV.write_text(
        '# Provenance for mondo_hierarchy.csv\n'
        f'# built:        {datetime.now(timezone.utc).isoformat(timespec="seconds")}\n'
        f'# built_by:     tools/build_mondo_hierarchy.py\n'
        f'# source:       data/crosswalks/mondo.obo\n'
        f'# source_bytes: {SRC.stat().st_size}\n'
        f'# source_sha256:{digest}\n'
        f'# mondo_release:{data_version or "(not stated in obo header)"}\n'
        f'# relation:     is_a only\n'
        f'# terms_seen:   {n_terms}\n'
        f'# obsolete_skipped: {n_obsolete}\n'
        f'# edges:        {len(edges)}\n'
        '# Cite the MONDO release above, not the build date.\n')

    print(f'terms seen        : {n_terms:,}')
    print(f'obsolete skipped  : {n_obsolete:,}')
    print(f'is_a edges written: {len(edges):,} -> {OUT.relative_to(ROOT)}')
    print(f'MONDO release     : {data_version or "(not stated)"}')


if __name__ == '__main__':
    main()
