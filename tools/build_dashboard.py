"""Build the interactive results dashboard, docs/index.html.

The query explorer draws the evidence neighbourhood each prompt carried, read
from results/tables/G02_run/evidence_graphs.json (tools/extract_evidence_graphs.py).

Reads the result tables under results/tables, assembles one JSON payload
and injects it into tools/dashboard_template.html. The output is a single
self-contained page (Chart.js from cdnjs is the only external file), served
by GitHub Pages from docs/.

Usage (from the repo root):
    python tools/build_dashboard.py
"""
import glob
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
T = ROOT / 'results' / 'tables'
TEMPLATE = ROOT / 'tools' / 'dashboard_template.html'
OUT = ROOT / 'docs' / 'index.html'

KG_ORDER = ['primekg', 'hetionet', 'drkg', 'openbilink', 'biokg', 'matrix']
KG_NAMES = {'primekg': 'PrimeKG', 'hetionet': 'Hetionet', 'drkg': 'DRKG',
            'openbilink': 'OpenBioLink', 'biokg': 'BioKG', 'matrix': 'MATRIX'}
KG_PALETTE = {'primekg': '#0072B2', 'hetionet': '#D55E00', 'drkg': '#009E73',
              'openbilink': '#E69F00', 'biokg': '#CC79A7', 'matrix': '#56B4E9'}
DIMS = ['Coverage', 'Annotation accuracy', 'Trustworthiness', 'Topology',
        'Stability', 'Task performance', 'Generalisation']
METRICS = ['Entity coverage', 'Relation coverage', 'Entity validity', 'Relational consistency',
           'Edge traceability', 'Source diversity', 'Uncertainty quantification',
           'Graph cohesion', 'Clustering coefficient', 'Graph reachability', 'Community purity',
           'Random-dropout stability', 'Periphery-dropout stability',
           'Link prediction', 'Neighbourhood retrieval', 'Multi-hop reasoning',
           'Data-sparse generalisation', 'Cross-domain generalisation', 'Prospective generalisation']
METRIC_DIM = dict(zip(METRICS, ['Coverage'] * 2 + ['Annotation accuracy'] * 2 + ['Trustworthiness'] * 3
                      + ['Topology'] * 4 + ['Stability'] * 2 + ['Task performance'] * 3 + ['Generalisation'] * 3))
RELEASES = {
    'primekg': 'V2.1, May 2022', 'hetionet': 'v1.0, June 2017', 'drkg': 'May 2020 build',
    'openbilink': '2020 v1.0.1', 'biokg': 'v1.0.0, June 2020', 'matrix': 'May 2026 snapshot',
}


def csv(rel, **kw):
    return pd.read_csv(T / rel, **kw)


def records(df, cols=None, r=4):
    d = df if cols is None else df[cols]
    d = d.round(r)
    return json.loads(d.to_json(orient='records'))


def kg_key(s):
    return str(s).lower().replace('openbiolink', 'openbilink')


def facts():
    coh = csv('04_topology/graph_cohesion.csv').set_index('KG')
    rel = csv('02_annotation_accuracy/relational_consistency_by_kg.csv').set_index('KG')
    sch = csv('01_coverage/schema_coverage.csv').set_index('KG')
    src = csv('03_trustworthiness/source_counts.csv').set_index('KG')
    out = []
    for k in KG_ORDER:
        edges = int(rel.loc[k, 'Edges'])
        if k == 'primekg':
            edges //= 2  # PrimeKG lists each edge in both directions; report undirected edges as in Table 2
        out.append(dict(key=k, name=KG_NAMES[k], colour=KG_PALETTE[k], release=RELEASES[k],
                        nodes=int(coh.loc[k, 'Nodes']), edges=edges,
                        node_types=int(sch.loc[k, 'Node types']), sources=int(src.loc[k, 'Sources'])))
    return out


def summary():
    bs = csv('00_summary/benchmark_summary.csv').set_index('KG')
    rows = []
    for k in KG_ORDER:
        r = bs.loc[k]
        rows.append(dict(kg=k, metrics={m: (None if pd.isna(r[m]) else round(float(r[m]), 4)) for m in METRICS},
                         dims={d: round(float(r[f'{d} [dim]']), 4) for d in DIMS},
                         overall=round(float(r['Overall Mean']), 4)))
    cat = csv('00_summary/category_scores.csv').set_index('KG')
    for row in rows:
        row['categories'] = {c: round(float(cat.loc[row['kg'], c]), 4) for c in ['Content', 'Structure', 'Inference']}
    disp = csv('00_summary/rank_dispersion.csv')
    disp['KG'] = disp['KG'].map(kg_key)
    return dict(rows=rows, metric_dim=METRIC_DIM, dims=DIMS, metrics=METRICS,
                rank_dispersion=records(disp))


def rankings():
    diss = csv('00_summary/dimension_dissociations.csv')
    sens = csv('00_summary/reversal_sensitivity.csv')
    conc = csv('00_summary/concordance_specifications.csv')
    return dict(dissociations=records(diss), sensitivity=records(sens), concordance=records(conc.fillna('')))


def coverage():
    ec = csv('01_coverage/entity_coverage.csv'); ec['KG'] = ec['KG'].map(kg_key)
    rc = csv('01_coverage/relation_coverage.csv'); rc['KG'] = rc['KG'].map(kg_key)
    res = csv('01_coverage/resolvability.csv')
    return dict(entity=records(ec, ['KG', 'Entity Type', 'KG Entities', 'Gold Standard', 'Overlap', 'Coverage (%)']),
                relation=records(rc, ['KG', 'Relation', 'KG pairs', 'Overlap', 'Gold', 'Coverage (%)']),
                resolvability=records(res, ['KG', 'Entity', 'raw', 'resolved', 'resolvability']))


def annotation():
    d = csv('02_annotation_accuracy/dimension_scores.csv')
    return dict(scores=records(d))


def trust():
    tr = csv('03_trustworthiness/edge_traceability.csv')
    sc = csv('03_trustworthiness/source_counts.csv')
    uq = csv('03_trustworthiness/uq_score_fields.csv')
    return dict(traceability=records(tr), sources=records(sc.fillna('')), uncertainty=records(uq.fillna('')))


def topology():
    coh = csv('04_topology/graph_cohesion.csv')
    cc = csv('04_topology/clustering_coefficient.csv')
    hop = csv('04_topology/hop_recovery.csv')
    cp = csv('04_topology/community_purity.csv')
    return dict(cohesion=records(coh.fillna(''), ['KG', 'Nodes', 'LCC nodes', 'LCC (%)', 'Components', 'Isolated']),
                clustering=records(cc, ['KG', 'C (empirical)', 'C_ER (random)', 'C / C_ER', 'SW score']),
                hops=records(hop), purity=records(cp, ['KG', 'Communities', 'AMI (scored)', 'NMI', 'Q (modularity)']))


def stability():
    st = csv('05_stability/stability.csv')
    st = st[st['Scorer'] == 'Adamic-Adar']
    st['Pair set'] = st['Pair set'].fillna('null')
    return dict(curves=records(st, ['KG', 'Strategy', 'Pair set', 'Dropout', 'Spearman r', 'Spearman SD']))


def tasks():
    lp = csv('06_task_performance/link_prediction.csv')
    lp = lp[lp['heuristic'] == 'Adamic-Adar']
    nr = csv('06_task_performance/neighbourhood_retrieval.csv')
    mh = csv('06_task_performance/multihop_reasoning.csv')
    tc = csv('06_task_performance/task_contrasts.csv')
    emb = csv('06_task_performance/link_prediction_embeddings.csv')
    return dict(link_prediction=records(lp, ['kg', 'neg_strategy', 'auroc', 'auroc_ci_lo', 'auroc_ci_hi', 'auprc', 'n_test_pos']),
                retrieval=records(nr, ['KG', 'recall@10', 'recall@50', 'recall@100', 'recall@frac', 'k_norm', 'n_evaluated', 'pool_size']),
                multihop=records(mh, ['KG', 'hits@10', 'hits@50', 'hits@100', 'hits@frac', 'k_norm', 'n_evaluated', 'pool_size']),
                contrasts=records(tc), embeddings=records(emb.fillna(''), ['kg', 'model', 'dim', 'test_auroc', 'adamic_adar_auroc']))


def generalisation():
    sp = csv('07_generalization/sparse_tiers.csv')
    cd = csv('07_generalization/cross_domain.csv')
    pr = csv('07_generalization/prospective.csv')
    pc = csv('07_generalization/prospective_coverage.csv')
    return dict(sparse=records(sp, ['KG', 'tier', 'auroc', 'auroc_ci_lo', 'auroc_ci_hi', 'n_pos']),
                areas=records(cd.fillna(''), ['KG', 'area', 'auroc', 'ci_lo', 'ci_hi', 'n_pos']),
                prospective=records(pr, ['KG', 'auroc', 'auroc_ci_lo', 'auroc_ci_hi', 'hits@10']),
                prospective_coverage=records(pc))


def grounding():
    out = {}
    for fam in ['main', 'covered', 'common']:
        a = csv(f'G03_analysis/arms_{fam}.csv')
        out[f'arms_{fam}'] = records(a, ['model', 'model_label', 'arm', 'arm_label', 'n_units', 'mrr', 'lo', 'hi', 'lift', 'lift_lo', 'lift_hi'])
    sig = csv('G03_analysis/significance_tests.csv')
    out['significance'] = records(sig, ['subset', 'model_label', 'contrast', 'n', 'diff', 'p_wilcoxon', 'q_bh', 'sig_fdr'])
    out['levers'] = records(csv('G03_analysis/lever_decomposition.csv'))
    out['absence'] = records(csv('G03_analysis/absence_penalty.csv'), ['model_label', 'arm', 'arm_label', 'positive', 'n_pairs', 'mrr', 'mrr_none', 'lift', 'lift_lo', 'lift_hi'])
    cap = csv('G04_scaling/capacity.csv')
    out['capacity'] = records(cap, ['model_label', 'arm', 'arm_label', 'mrr', 'lo', 'hi', 'lift', 'lift_lo', 'lift_hi', 'params_b'])
    rs = csv('G02_run/run_summary.csv')
    out['run_summary'] = records(rs, ['family', 'model_label', 'arm', 'arm_label', 'n', 'mrr', 'hit1', 'hit3', 'parse_rate', 'reliance', 'confidence'])
    return out


def explorer():
    b = csv('G01_benchmark/benchmark.csv')
    q = csv('G01_benchmark/queries_seed0.csv')
    q = q[q['family'] == 'main'].set_index('pair_id')
    runs = pd.concat([pd.read_csv(f) for f in glob.glob(str(T / 'G02_run' / 'runs' / 'main__*.csv'))], ignore_index=True)
    runs = runs[runs['pair_id'].isin(b['pair_id'])]
    agg = (runs.groupby(['pair_id', 'model', 'arm'])
           .agg(rr=('reciprocal_rank', 'mean'), rank=('rank', 'mean'), hit1=('hit1', 'mean'),
                evidence=('pos_has_evidence', 'max'), n_evidence=('n_with_evidence', 'mean'),
                parsed=('parsed', 'mean'), n=('rank', 'size')).reset_index())
    first = (runs.sort_values(['seed', 'shuffle']).groupby(['pair_id', 'model', 'arm'])
             .agg(rationale=('pos_rationale', 'first'), basis=('pos_basis', 'first'),
                  agreement=('pos_agreement', 'first'), confidence=('pos_confidence', 'first')).reset_index())
    agg = agg.merge(first, on=['pair_id', 'model', 'arm'])
    pairs = []
    for _, r in b.iterrows():
        pid = int(r['pair_id'])
        rows = agg[agg['pair_id'] == pid]
        res = {}
        for _, x in rows.iterrows():
            res.setdefault(x['model'], {})[x['arm']] = dict(
                rr=round(float(x['rr']), 3), rank=round(float(x['rank']), 2), hit1=round(float(x['hit1']), 2),
                evidence=bool(x['evidence']), n_evidence=round(float(x['n_evidence']), 1), n=int(x['n']),
                rationale=None if pd.isna(x['rationale']) else str(x['rationale'])[:300],
                basis=None if pd.isna(x['basis']) else x['basis'],
                agreement=None if pd.isna(x['agreement']) else x['agreement'],
                confidence=None if pd.isna(x['confidence']) else float(x['confidence']))
        pairs.append(dict(
            id=pid, drug=r['drug_name'], disease=r['disease_name'], drugbank=r['drugbank_id'], mondo=r['mondo_id'],
            date=r['approval_date'], agencies=r['agencies'], area=r['area'], modality=r['modality'],
            distractors=[s.strip() for s in str(q.loc[pid, 'distractors']).split(';')] if pid in q.index else [],
            evaluable={k: bool(r[f'evaluable_{k}']) for k in ['primekg', 'drkg', 'biokg']},
            shown={k: bool(r[f'shown_{k}']) for k in ['primekg', 'drkg', 'biokg']},
            common=bool(r['common']), results=res))
    ev_path = T / 'G02_run' / 'evidence_graphs.json'
    evidence = json.load(open(ev_path)) if ev_path.exists() else {}
    main_models = ['gpt-4.1-mini', 'gemini:gemini-3.5-flash-lite', 'llama3.3:70b']
    ladder = [m for m in sorted(runs['model'].unique().tolist()) if m not in main_models]
    return dict(pairs=pairs, models=main_models, ladder=ladder, arms=['none', 'primekg', 'drkg', 'biokg'], evidence=evidence)


def main():
    data = dict(
        kgs=facts(), summary=summary(), rankings=rankings(), coverage=coverage(),
        annotation=annotation(), trust=trust(), topology=topology(), stability=stability(),
        tasks=tasks(), generalisation=generalisation(), grounding=grounding(), explorer=explorer(),
        palette=KG_PALETTE, names=KG_NAMES, order=KG_ORDER,
    )
    payload = json.dumps(data, separators=(',', ':'), allow_nan=False)
    html = TEMPLATE.read_text()
    assert '/*__DATA__*/' in html
    html = html.replace('/*__DATA__*/', payload.replace('</', '<\\/'))
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(html)
    print(f'wrote {OUT} ({OUT.stat().st_size / 1e6:.1f} MB; {len(data["explorer"]["pairs"])} pairs in the explorer)')


if __name__ == '__main__':
    main()
