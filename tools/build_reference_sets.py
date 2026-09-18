#!/usr/bin/env python3
"""Build the task-derived entity and relation reference sets.

Definitions are asymmetric by design (see data/REFERENCE_SETS.md):
  DISEASES  every disease, including undrugged ones — repurposing exists to
            serve them. Rule: MONDO non-obsolete, non-grouping.
  DRUGS     approved only — you can only repurpose a drug that exists.
            Rule: Open Targets maximumClinicalStage == 'APPROVAL'.
  GENES     targets of the approved drug set.
  PATHWAYS  Reactome pathways of those targets.

All sources are Open Targets (ChEMBL/DailyMed derived) or MONDO, so none is
circular with the DrugBank-ingesting graphs being scored.

Run: python tools/build_reference_sets.py
"""
from __future__ import annotations
import argparse, json, sys
from collections import defaultdict
from pathlib import Path
import re
import pandas as pd


def load_parquet(d: Path, columns=None):
    parts = sorted(d.glob('*.parquet'))
    if not parts:
        return None
    return pd.concat([pd.read_parquet(p, columns=columns) for p in parts],
                     ignore_index=True)


def mondo_universe(obo: Path):
    """Non-obsolete, non-grouping MONDO terms, plus rare flag and is_a edges."""
    live, rare, group, isa, cur, parents = set(), set(), set(), set(), None, []
    for line in obo.open(encoding='utf-8'):
        line = line.rstrip('\n')
        if line == '[Term]':
            cur, parents = {'id': None, 'obs': False, 'rare': False, 'grp': False}, []
            continue
        if cur is None:
            continue
        if line.startswith('id: MONDO:'):
            mm = re.match(r'id: MONDO:(\d+)', line)
            if mm:
                cur['id'] = str(int(mm.group(1)))
        elif line.startswith('is_obsolete: true'):
            cur['obs'] = True
        elif line.startswith('subset: rare'):
            cur['rare'] = True
        elif line.startswith('subset: disease_grouping'):
            cur['grp'] = True
        elif line.startswith('is_a: MONDO:'):
            # is_a lines may carry a {source=...} qualifier and a ! comment:
            #   is_a: MONDO:0002816 {source="DOID:10493"} ! carcinoma
            mm = re.match(r'is_a: MONDO:(\d+)', line)
            if mm:
                parents.append(str(int(mm.group(1))))
        elif line == '' and cur.get('id') and not cur['obs']:
            live.add(cur['id'])
            if cur['rare']:
                rare.add(cur['id'])
            if cur['grp']:
                group.add(cur['id'])
            for p in parents:
                isa.add((cur['id'], p))
            cur = None
    uni = live - group
    return uni, rare & uni, {(c, p) for c, p in isa if c in uni and p in uni}


def to_mondo_map(xw: Path):
    """OT disease id (EFO_/MONDO_/HP_/...) -> bare MONDO numeric.

    OT ids are a mix. MONDO_ ones are direct; the rest go through the pinned
    MONDO SSSOM bridge, which carries EFO among its namespaces.
    """
    m = {}
    sssom = xw / 'mondo.sssom.tsv'
    if sssom.exists():
        ss = pd.read_csv(sssom, sep='\t', comment='#', low_memory=False)
        ss = ss[ss['predicate_id'] == 'skos:exactMatch']
        for subj, obj in zip(ss['subject_id'], ss['object_id']):
            if not isinstance(subj, str) or not subj.startswith('MONDO:'):
                continue
            mid = str(int(subj.split(':')[1]))
            if isinstance(obj, str) and ':' in obj:
                ns, lid = obj.split(':', 1)
                m[f'{ns.upper()}_{lid}'] = mid
    return m


def ot_to_mondo(series, bridge):
    out = []
    for v in series:
        if not isinstance(v, str):
            out.append(None); continue
        if v.startswith('MONDO_'):
            out.append(str(int(v.split('_')[1])))
        else:
            ns, _, lid = v.partition('_')
            out.append(bridge.get(f'{ns.upper()}_{lid}'))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--crosswalks', default='data/crosswalks')
    ap.add_argument('--out', default='data/reference_sets')
    ap.add_argument('--benchmark', default='data/gold_standards/legacy/gold_standard_v4_bigset.tsv')
    a = ap.parse_args()
    xw, out = Path(a.crosswalks), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    man, rel_counts = {}, {}

    def write(name, rows, header):
        p = out / name
        pd.DataFrame(rows, columns=header).drop_duplicates().to_csv(p, sep='\t', index=False)
        return sum(1 for _ in p.open()) - 1

    # ---------------- entities: disease ------------------------------------
    obo = xw / 'mondo.obo'
    if not obo.exists():
        print('  FATAL: mondo.obo missing. Run tools/build_crosswalks.py'); return 1
    universe, rare, isa = mondo_universe(obo)
    (out / 'entity_disease.txt').write_text('\n'.join(sorted(universe, key=int)) + '\n')
    (out / 'entity_disease_rare.txt').write_text('\n'.join(sorted(rare, key=int)) + '\n')
    man['entity_disease'] = {'rule': 'MONDO non-obsolete, non-grouping',
                             'source': 'mondo.obo', 'n': len(universe), 'n_rare': len(rare)}
    print(f'  entity  disease        {len(universe):>8,}  ({len(rare):,} rare)')

    bridge = to_mondo_map(xw)

    # ---------------- entities: drug ---------------------------------------
    dm = load_parquet(xw / 'ot_drug_molecule', ['id', 'name', 'maximumClinicalStage'])
    if dm is None:
        print('  FATAL: ot_drug_molecule missing'); return 1
    appr = dm[dm['maximumClinicalStage'] == 'APPROVAL']
    drugs = set(appr['id'])
    (out / 'entity_drug.txt').write_text('\n'.join(sorted(drugs)) + '\n')
    man['entity_drug'] = {'rule': "maximumClinicalStage == 'APPROVAL' (there is no "
                                  "isApproved column in 26.06)",
                          'source': 'ot_drug_molecule', 'n': len(drugs)}
    print(f'  entity  drug           {len(drugs):>8,}  (approved, of {len(dm):,} molecules)')

    # ---------------- relation #1 drug -> gene ------------------------------
    moa = load_parquet(xw / 'ot_drug_mechanism_of_action', ['chemblIds', 'targets'])
    dg = set()
    if moa is not None:
        for cids, tgts in zip(moa['chemblIds'], moa['targets']):
            if cids is None or tgts is None:
                continue
            for c in cids:
                if c in drugs:
                    for t in tgts:
                        dg.add((c, t))
    genes = {t for _, t in dg}
    (out / 'entity_gene.txt').write_text('\n'.join(sorted(genes)) + '\n')
    rel_counts['drug_gene'] = write('relation_drug_gene.tsv', sorted(dg), ['drug', 'gene'])
    man['entity_gene'] = {'rule': 'targets of the approved drug set',
                          'source': 'ot_drug_mechanism_of_action', 'n': len(genes)}
    print(f'  entity  gene           {len(genes):>8,}  (targets of approved drugs)')
    print(f'  rel #1  drug->gene     {len(dg):>8,}')

    # ---------------- gene -> pathway (building block) ----------------------
    tg = load_parquet(xw / 'ot_target', ['id', 'pathways'])
    gp = set()
    pw_names = {}
    if tg is not None:
        for gid, pws in zip(tg['id'], tg['pathways']):
            if pws is None:
                continue
            for pw in pws:
                gp.add((gid, pw['pathwayId']))
                if pw.get('pathway'):
                    pw_names[pw['pathwayId']] = str(pw['pathway']).strip()
    pathways = {p for g, p in gp if g in genes}
    (out / 'entity_pathway.txt').write_text('\n'.join(sorted(pathways)) + '\n')

    # Reactome pathway NAMES, over every pathway OT mentions rather than only
    # those of the target gene set. Used by the nb01 pathway diagnostic to
    # measure name agreement between a KG's pathway vocabulary and Reactome's,
    # which is what distinguishes "different accessions, same biology" from
    # "no pathway biology". Not a scored reference set.
    with (out / 'entity_pathway_names.tsv').open('w') as fh:
        fh.write('reactome_id\tname\n')
        for pid, nm in sorted(pw_names.items()):
            fh.write(f'{pid}\t{nm}\n')
    man['entity_pathway_names'] = {
        'rule': 'Reactome id -> name, all pathways mentioned by ot_target',
        'source': 'ot_target.pathways', 'n': len(pw_names),
        'scored': False}
    print(f'  aux     pathway names  {len(pw_names):>8,}  (for the nb01 pathway diagnostic)')
    man['entity_pathway'] = {'rule': 'Reactome pathways of the target gene set',
                             'source': 'ot_target.pathways', 'n': len(pathways)}
    print(f'  entity  pathway        {len(pathways):>8,}  (of target genes)')

    g2p = defaultdict(set)
    for g, p in gp:
        g2p[g].add(p)

    # ---------------- relation #2 drug -> pathway ---------------------------
    dp = {(d, p) for d, g in dg for p in g2p.get(g, ())}
    rel_counts['drug_pathway'] = write('relation_drug_pathway.tsv', sorted(dp), ['drug', 'pathway'])
    print(f'  rel #2  drug->pathway  {len(dp):>8,}  (derived: #1 o gene->pathway)')

    # ---------------- relation #3 disease -> gene ---------------------------
    ev_dirs = sorted(xw.glob('ot_evidence_*'))
    dgene, by_src = set(), {}
    for d in ev_dirs:
        df = load_parquet(d, ['targetId', 'diseaseId', 'datasourceId'])
        if df is None:
            continue
        df = df.dropna(subset=['targetId', 'diseaseId'])
        df['m'] = ot_to_mondo(df['diseaseId'], bridge)
        df = df.dropna(subset=['m'])
        pairs = set(zip(df['m'], df['targetId']))
        by_src[d.name.replace('ot_evidence_', '')] = len(pairs)
        dgene |= pairs
    rel_counts['disease_gene'] = write('relation_disease_gene.tsv', sorted(dgene), ['disease', 'gene'])
    print(f'  rel #3  disease->gene  {len(dgene):>8,}  from {len(by_src)} curated sources')
    for k, v in sorted(by_src.items(), key=lambda kv: -kv[1]):
        print(f'             {k:28s} {v:>8,}')

    # ---------------- relation #4 disease -> pathway ------------------------
    dpw = {(dz, p) for dz, g in dgene for p in g2p.get(g, ())}
    rel_counts['disease_pathway'] = write('relation_disease_pathway.tsv', sorted(dpw), ['disease', 'pathway'])
    print(f'  rel #4  dis->pathway   {len(dpw):>8,}  (derived: #3 o gene->pathway)')

    # ---------------- relation #5 disease -> disease ------------------------
    rel_counts['disease_disease'] = write('relation_disease_disease.tsv', sorted(isa), ['child', 'parent'])
    print(f'  rel #5  dis->dis       {len(isa):>8,}  (MONDO is_a, within the universe)')

    # ---------------- relation #6 drug -> disease ---------------------------
    ci = load_parquet(xw / 'ot_clinical_indication', ['drugId', 'diseaseId', 'maxClinicalStage'])
    dd, drugged = set(), set()
    if ci is not None:
        ci = ci.dropna(subset=['drugId', 'diseaseId'])
        ci['m'] = ot_to_mondo(ci['diseaseId'], bridge)
        ci = ci.dropna(subset=['m'])
        appr_ind = ci[(ci['maxClinicalStage'] == 'APPROVAL') & (ci['drugId'].isin(drugs))]
        dd = set(zip(appr_ind['drugId'], appr_ind['m']))
        drugged = set(appr_ind['m']) & universe
    rel_counts['drug_disease'] = write('relation_drug_disease.tsv', sorted(dd), ['drug', 'disease'])
    print(f'  rel #6  drug->disease  {len(dd):>8,}  (approved indications only)')

    # ---------------- disease strata ---------------------------------------
    undrugged = universe - drugged
    (out / 'entity_disease_drugged.txt').write_text('\n'.join(sorted(drugged, key=int)) + '\n')
    (out / 'entity_disease_undrugged.txt').write_text('\n'.join(sorted(undrugged, key=int)) + '\n')
    man['disease_strata'] = {'rule': 'drugged = has >=1 APPROVAL-stage indication',
                             'source': 'ot_clinical_indication',
                             'n_drugged': len(drugged), 'n_undrugged': len(undrugged),
                             'n_rare_undrugged': len(rare & undrugged)}
    print(f'\n  strata  drugged        {len(drugged):>8,}  ({100*len(drugged)/len(universe):.1f}%)')
    print(f'          undrugged      {len(undrugged):>8,}  ({100*len(undrugged)/len(universe):.1f}%)')
    print(f'          rare+undrugged {len(rare & undrugged):>8,}')

    # ---------------- tier 1 -----------------------------------------------
    b = Path(a.benchmark)
    if b.exists():
        import csv
        rows = list(csv.DictReader(b.open(), delimiter='\t'))
        man['tier1'] = {'rule': 'entities of the prospective benchmark', 'source': b.name,
                        'pairs': len(rows),
                        'n_drugs': len({r['drug_id'] for r in rows if r.get('drug_id')}),
                        'n_diseases': len({r['disease_id'] for r in rows if r.get('disease_id')})}

    man['relations'] = rel_counts
    (out / 'MANIFEST.json').write_text(json.dumps(man, indent=2))
    print(f'\nwrote {out}/MANIFEST.json')
    return 0


if __name__ == '__main__':
    sys.exit(main())
