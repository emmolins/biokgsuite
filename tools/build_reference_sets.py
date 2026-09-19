#!/usr/bin/env python3
"""Build the task-derived entity and relation reference sets.

Definitions are asymmetric by design (see data/REFERENCE_SETS.md):
  DISEASES  every human disease, including undrugged ones: repurposing exists
            to serve them. Rule: MONDO non-obsolete, non-grouping, descendant
            of MONDO:0700096 (human disease or disorder).
  DRUGS     approved therapeutics, one row per active molecule: you can only
            repurpose a drug that exists and is used for something. Rule:
            Open Targets maximumClinicalStage == 'APPROVAL', with at least one
            APPROVAL-stage indication (on the molecule or its ChEMBL parent),
            not a vaccine component, collapsed to the ChEMBL parent so a salt
            and its free base are one drug.
  GENES     the genes the task touches: targets of the approved drugs plus
            genes with curated disease evidence.
  PATHWAYS  Reactome pathways of that gene set.

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
    return uni, rare & uni, {(c, p) for c, p in isa if c in uni and p in uni}, isa


HUMAN_DISEASE_ROOT = '700096'     # MONDO:0700096 human disease or disorder


def descendants(root: str, isa: set) -> set:
    """Every term below ``root`` in the is_a graph (root included)."""
    kids = defaultdict(set)
    for c, p in isa:
        kids[p].add(c)
    seen, stack = {root}, [root]
    while stack:
        for k in kids.get(stack.pop(), ()):
            if k not in seen:
                seen.add(k); stack.append(k)
    return seen


def approved_therapeutics(dm: pd.DataFrame, ci: pd.DataFrame):
    """(set of parent ChEMBL ids, child->parent map) for the drug reference set.

    Rules, each stated so it can be argued with:
      1. maximumClinicalStage == 'APPROVAL' (Open Targets, ChEMBL-derived).
      2. At least one APPROVAL-stage indication in clinical_indication, on the
         molecule or its parent. Drops imaging agents, excipients and
         hormones that carry an approval but no therapeutic use.
      3. Not a vaccine component: the post-cutoff approval set excludes
         vaccines, and one rule should decide what a drug is in both parts.
      4. Collapsed to the ChEMBL parent, so doxazosin and doxazosin mesylate
         are one drug and coverage counts active molecules, not salt forms.
    """
    ap = dm[dm['maximumClinicalStage'] == 'APPROVAL'].copy()
    ap['parent'] = ap['parentId'].fillna(ap['id'])
    parent_of = dict(zip(dm['id'], dm['parentId'].fillna(dm['id'])))
    indicated = set(ci.loc[ci['maxClinicalStage'] == 'APPROVAL', 'drugId'])
    indicated |= {parent_of.get(d, d) for d in indicated}
    keep = ap[(ap['id'].isin(indicated) | ap['parent'].isin(indicated))
              & (ap['drugType'] != 'Vaccine component')]
    return set(keep['parent']), parent_of


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
    ap.add_argument('--benchmark', default='data/gold_standards/post_cutoff_approvals/post_cutoff_approvals.csv')
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
    universe, rare, isa, isa_all = mondo_universe(obo)
    human = descendants(HUMAN_DISEASE_ROOT, isa_all)
    dropped = universe - human
    universe, rare = universe & human, rare & human
    isa = {(c, p) for c, p in isa if c in universe and p in universe}
    (out / 'entity_disease.txt').write_text('\n'.join(sorted(universe, key=int)) + '\n')
    (out / 'entity_disease_rare.txt').write_text('\n'.join(sorted(rare, key=int)) + '\n')
    man['entity_disease'] = {'rule': 'MONDO non-obsolete, non-grouping, descendant of '
                                     'MONDO:0700096 (human disease or disorder)',
                             'source': 'mondo.obo', 'n': len(universe), 'n_rare': len(rare),
                             'n_outside_human_disease': len(dropped)}
    print(f'  entity  disease        {len(universe):>8,}  ({len(rare):,} rare; '
          f'{len(dropped)} live terms outside human disease dropped)')

    bridge = to_mondo_map(xw)

    # ---------------- entities: drug ---------------------------------------
    dm = load_parquet(xw / 'ot_drug_molecule', ['id', 'name', 'drugType', 'parentId',
                                                 'maximumClinicalStage'])
    ci = load_parquet(xw / 'ot_clinical_indication', ['drugId', 'diseaseId', 'maxClinicalStage'])
    if dm is None or ci is None:
        print('  FATAL: ot_drug_molecule or ot_clinical_indication missing'); return 1
    drugs, parent_of = approved_therapeutics(dm, ci)
    n_approved = int((dm['maximumClinicalStage'] == 'APPROVAL').sum())
    (out / 'entity_drug.txt').write_text('\n'.join(sorted(drugs)) + '\n')
    man['entity_drug'] = {'rule': "maximumClinicalStage == 'APPROVAL', with an APPROVAL-stage "
                                  "indication on the molecule or its parent, not a vaccine "
                                  "component, collapsed to the ChEMBL parent",
                          'source': 'ot_drug_molecule + ot_clinical_indication',
                          'n': len(drugs), 'n_approved_molecules': n_approved}
    print(f'  entity  drug           {len(drugs):>8,}  (approved therapeutics, parent molecules; '
          f'{n_approved:,} approved molecule rows)')

    # ---------------- relation #1 drug -> gene ------------------------------
    moa = load_parquet(xw / 'ot_drug_mechanism_of_action', ['chemblIds', 'targets'])
    dg = set()
    if moa is not None:
        for cids, tgts in zip(moa['chemblIds'], moa['targets']):
            if cids is None or tgts is None:
                continue
            for c in cids:
                c = parent_of.get(c, c)          # a salt's mechanism is its parent's
                if c in drugs:
                    for t in tgts:
                        dg.add((c, t))
    targets = {t for _, t in dg}
    rel_counts['drug_gene'] = write('relation_drug_gene.tsv', sorted(dg), ['drug', 'gene'])
    print(f'  rel #1  drug->gene     {len(dg):>8,}  ({len(targets):,} targets)')

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
    g2p = defaultdict(set)
    for g, p in gp:
        g2p[g].add(p)

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

    # ---------------- entities: gene and pathway ----------------------------
    # The genes the task touches: what approved drugs act on, and what curated
    # evidence implicates in disease. Targets alone (about 1,100 genes) are
    # held by every gene-centric graph, so the set would not discriminate.
    disease_genes = {g for _, g in dgene}
    genes = targets | disease_genes
    (out / 'entity_gene.txt').write_text('\n'.join(sorted(genes)) + '\n')
    man['entity_gene'] = {'rule': 'targets of the approved drugs, plus genes with curated '
                                  'disease evidence (relation #3)',
                          'source': 'ot_drug_mechanism_of_action + ot_evidence_*',
                          'n': len(genes), 'n_targets': len(targets),
                          'n_disease_genes': len(disease_genes)}
    print(f'  entity  gene           {len(genes):>8,}  ({len(targets):,} targets + '
          f'{len(disease_genes):,} disease genes)')
    pathways = {p for g, p in gp if g in genes}
    (out / 'entity_pathway.txt').write_text('\n'.join(sorted(pathways)) + '\n')
    man['entity_pathway'] = {'rule': 'Reactome pathways of the gene reference set',
                             'source': 'ot_target.pathways', 'n': len(pathways)}
    print(f'  entity  pathway        {len(pathways):>8,}  (of the gene set)')

    # ---------------- relation #4 disease -> pathway ------------------------
    dpw = {(dz, p) for dz, g in dgene for p in g2p.get(g, ())}
    rel_counts['disease_pathway'] = write('relation_disease_pathway.tsv', sorted(dpw), ['disease', 'pathway'])
    print(f'  rel #4  dis->pathway   {len(dpw):>8,}  (derived: #3 o gene->pathway)')

    # ---------------- relation #5 disease -> disease ------------------------
    rel_counts['disease_disease'] = write('relation_disease_disease.tsv', sorted(isa), ['child', 'parent'])
    print(f'  rel #5  dis->dis       {len(isa):>8,}  (MONDO is_a, within the universe)')

    # ---------------- relation #6 drug -> disease ---------------------------
    dd, drugged = set(), set()
    cim = ci.dropna(subset=['drugId', 'diseaseId']).copy()
    cim['m'] = ot_to_mondo(cim['diseaseId'], bridge)
    cim = cim.dropna(subset=['m'])
    cim['drugId'] = cim['drugId'].map(lambda d: parent_of.get(d, d))
    appr_ind = cim[(cim['maxClinicalStage'] == 'APPROVAL') & (cim['drugId'].isin(drugs))
                   & (cim['m'].isin(universe))]
    dd = set(zip(appr_ind['drugId'], appr_ind['m']))
    drugged = set(appr_ind['m'])
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
        bench = pd.read_csv(b, dtype=str)
        dcol = next((c for c in bench.columns if 'drugbank' in c.lower()), None)
        zcol = next((c for c in bench.columns if 'mondo' in c.lower()), None)
        man['post_cutoff_set'] = {'rule': 'entities of the post-cutoff approval set (notebook 07, Part II)',
                                  'source': b.name, 'pairs': len(bench),
                                  'n_drugs': int(bench[dcol].nunique()) if dcol else None,
                                  'n_diseases': int(bench[zcol].nunique()) if zcol else None}

    man['relations'] = rel_counts
    (out / 'MANIFEST.json').write_text(json.dumps(man, indent=2))
    print(f'\nwrote {out}/MANIFEST.json')
    return 0


if __name__ == '__main__':
    sys.exit(main())
