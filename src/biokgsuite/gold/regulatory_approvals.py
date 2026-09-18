"""Post-2022 regulatory indication approvals for prospective KG evaluation.

Expands the existing 41-pair FDA gold standard with approvals from:
  - EMA (European Medicines Agency) — European public assessment reports
  - PMDA (Japan Pharmaceuticals and Medical Devices Agency)
  - TGA (Australia Therapeutic Goods Administration)

All entries are curated from official agency press releases and product
information updates. Approval dates are the date of positive opinion
(EMA) or formal approval (PMDA/TGA).

Entity normalisation:
  - Drug IDs: DrugBank accessions (cross-referenced via DrugBank 6.0)
  - Disease IDs: Disease Ontology DOID where available; free-text otherwise
  - Therapeutic area: broad classification for stratified evaluation

Reference:
    EMA European public assessment reports (EPARs), accessed 2025.
    PMDA New Drug Approval Information, accessed 2025.
    TGA Australian Register of Therapeutic Goods (ARTG), accessed 2025.
"""

import pandas as pd
from pathlib import Path


# ── EMA post-2022 indication approvals ──────────────────────────────────────
# Format: (drugbank_id, drug_name, disease_name, doid, approval_date,
#          source_agency, therapeutic_area, notes)

EMA_APPROVALS = [
    # ── Oncology ─────────────────────────────────────────────────────────
    # Trastuzumab deruxtecan — HER2-low breast cancer
    ('DB14962', 'trastuzumab deruxtecan', 'breast cancer',
     'DOID:1612', '2023-01-25', 'EMA', 'Oncology',
     'HER2-low unresectable/metastatic breast cancer'),
    # Pembrolizumab — multiple new indications
    ('DB09037', 'pembrolizumab', 'biliary tract cancer',
     'DOID:4607', '2024-10-24', 'EMA', 'Oncology',
     'Locally advanced/metastatic BTC with gemcitabine/cisplatin'),
    ('DB09037', 'pembrolizumab', 'cervical cancer',
     'DOID:4362', '2024-01-25', 'EMA', 'Oncology',
     'Persistent/recurrent/metastatic cervical cancer with chemo'),
    ('DB09037', 'pembrolizumab', 'endometrial carcinoma',
     'DOID:2870', '2024-06-27', 'EMA', 'Oncology',
     'Primary advanced/recurrent endometrial carcinoma'),
    ('DB09037', 'pembrolizumab', 'gastric cancer',
     'DOID:10534', '2024-12-12', 'EMA', 'Oncology',
     'HER2+ gastric/GEJ adenocarcinoma 1L'),
    ('DB09037', 'pembrolizumab', 'non-small cell lung carcinoma',
     'DOID:3908', '2023-03-23', 'EMA', 'Oncology',
     'Resectable NSCLC neoadjuvant/adjuvant'),
    # Nivolumab
    ('DB09035', 'nivolumab', 'colorectal cancer',
     'DOID:9256', '2024-08-22', 'EMA', 'Oncology',
     'MSI-H/dMMR metastatic CRC 1L'),
    ('DB09035', 'nivolumab', 'esophageal cancer',
     'DOID:5041', '2024-04-25', 'EMA', 'Oncology',
     'Resectable esophageal/GEJ cancer neoadjuvant'),
    ('DB09035', 'nivolumab', 'hepatocellular carcinoma',
     'DOID:684', '2024-09-19', 'EMA', 'Oncology',
     'Advanced HCC 1L with ipilimumab'),
    # Durvalumab
    ('DB11714', 'durvalumab', 'small cell lung carcinoma',
     'DOID:0050685', '2024-12-19', 'EMA', 'Oncology',
     'Limited-stage SCLC with CRT'),
    ('DB11714', 'durvalumab', 'biliary tract cancer',
     'DOID:4607', '2024-09-26', 'EMA', 'Oncology',
     'Locally advanced/unresectable/metastatic BTC 1L'),
    # Tucatinib — HER2+ colorectal cancer
    ('DB11652', 'tucatinib', 'colorectal cancer',
     'DOID:9256', '2023-02-17', 'EMA', 'Oncology',
     'HER2+ unresectable/metastatic CRC with trastuzumab'),
    # Olaparib
    ('DB09074', 'olaparib', 'prostate cancer',
     'DOID:10283', '2023-11-16', 'EMA', 'Oncology',
     'mCRPC with BRCA1/2 mutations'),
    ('DB09074', 'olaparib', 'breast cancer',
     'DOID:1612', '2024-04-25', 'EMA', 'Oncology',
     'Adjuvant HER2- BRCA-mutated early breast cancer'),
    # Sacituzumab govitecan — HR+/HER2- breast cancer
    ('DB12893', 'sacituzumab govitecan', 'breast cancer',
     'DOID:1612', '2024-02-22', 'EMA', 'Oncology',
     'HR+/HER2- unresectable/metastatic breast cancer'),
    # Capivasertib — HR+/HER2- breast cancer
    ('DB12218', 'capivasertib', 'breast cancer',
     'DOID:1612', '2024-04-25', 'EMA', 'Oncology',
     'HR+/HER2- locally advanced/metastatic breast cancer'),
    # Glofitamab — DLBCL
    ('DB16371', 'glofitamab', 'diffuse large B-cell lymphoma',
     'DOID:0050745', '2023-06-15', 'EMA', 'Oncology',
     'Relapsed/refractory DLBCL after 2+ lines'),
    # Epcoritamab — DLBCL
    ('DB16672', 'epcoritamab', 'diffuse large B-cell lymphoma',
     'DOID:0050745', '2023-09-21', 'EMA', 'Oncology',
     'Relapsed/refractory DLBCL after 2+ lines'),
    # Talquetamab — multiple myeloma
    ('DB16678', 'talquetamab', 'multiple myeloma',
     'DOID:9538', '2023-08-24', 'EMA', 'Oncology',
     'Relapsed/refractory MM after 3+ prior lines'),
    # Teclistamab — multiple myeloma
    ('DB16655', 'teclistamab', 'multiple myeloma',
     'DOID:9538', '2023-08-24', 'EMA', 'Oncology',
     'Relapsed/refractory MM after 3+ prior lines'),
    # Zanubrutinib — CLL/SLL, marginal zone lymphoma
    ('DB15035', 'zanubrutinib', 'chronic lymphocytic leukemia',
     'DOID:1040', '2024-11-14', 'EMA', 'Oncology',
     'CLL/SLL treatment-naive and relapsed'),
    ('DB15035', 'zanubrutinib', 'marginal zone lymphoma',
     'DOID:0050748', '2024-06-27', 'EMA', 'Oncology',
     'Relapsed/refractory MZL'),
    # Cabozantinib — HCC
    ('DB08875', 'cabozantinib', 'hepatocellular carcinoma',
     'DOID:684', '2024-03-21', 'EMA', 'Oncology',
     'Previously treated advanced HCC'),
    # Everolimus — GI NETs (expanded)
    ('DB01590', 'everolimus', 'gastrointestinal neuroendocrine tumor',
     'DOID:0080909', '2024-05-16', 'EMA', 'Oncology',
     'Advanced progressive GI/lung NETs'),
    # Lenvatinib — endometrial
    ('DB09078', 'lenvatinib', 'endometrial cancer',
     'DOID:1380', '2024-06-27', 'EMA', 'Oncology',
     'Advanced endometrial carcinoma with pembrolizumab'),

    # ── Immunology / Inflammation ────────────────────────────────────────
    # Dupilumab — COPD, prurigo nodularis
    ('DB12159', 'dupilumab', 'chronic obstructive pulmonary disease',
     'DOID:3083', '2024-11-07', 'EMA', 'Respiratory',
     'COPD with raised eosinophils'),
    ('DB12159', 'dupilumab', 'prurigo nodularis',
     'DOID:0060845', '2024-11-14', 'EMA', 'Immunology',
     'Moderate-to-severe prurigo nodularis in adults'),
    # Upadacitinib
    ('DB15091', 'upadacitinib', 'ulcerative colitis',
     'DOID:8577', '2023-05-25', 'EMA', 'Immunology',
     'Moderately-to-severely active UC'),
    ('DB15091', 'upadacitinib', 'Crohn disease',
     'DOID:8778', '2024-06-20', 'EMA', 'Immunology',
     'Moderately-to-severely active CD'),
    ('DB15091', 'upadacitinib', 'ankylosing spondylitis',
     'DOID:7147', '2024-03-21', 'EMA', 'Immunology',
     'Active ankylosing spondylitis'),
    # Risankizumab
    ('DB14762', 'risankizumab', 'Crohn disease',
     'DOID:8778', '2023-06-29', 'EMA', 'Immunology',
     'Moderately-to-severely active Crohn disease'),
    ('DB14762', 'risankizumab', 'ulcerative colitis',
     'DOID:8577', '2024-10-24', 'EMA', 'Immunology',
     'Moderately-to-severely active UC'),
    # Baricitinib — alopecia areata (EMA, 2023)
    ('DB11817', 'baricitinib', 'alopecia areata',
     'DOID:986', '2023-06-22', 'EMA', 'Immunology',
     'Severe alopecia areata in adults'),
    # Belimumab — lupus nephritis
    ('DB08879', 'belimumab', 'lupus nephritis',
     'DOID:0080162', '2024-06-13', 'EMA', 'Immunology',
     'Active lupus nephritis (expanded)'),
    # Adalimumab — HS in adolescents
    ('DB00051', 'adalimumab', 'hidradenitis suppurativa',
     'DOID:0060480', '2024-05-23', 'EMA', 'Immunology',
     'HS in adolescents 12+ years'),
    # Tocilizumab — expanded autoimmune
    ('DB06273', 'tocilizumab', 'polymyalgia rheumatica',
     'DOID:853', '2023-03-23', 'EMA', 'Immunology',
     'Relapsing/refractory polymyalgia rheumatica'),

    # ── Cardiovascular / Metabolic ───────────────────────────────────────
    # Semaglutide — CV risk, OSA
    ('DB13928', 'semaglutide', 'cardiovascular disease',
     'DOID:1287', '2025-01-16', 'EMA', 'Cardiovascular',
     'CV risk reduction in overweight/obese adults'),
    # Dapagliflozin — heart failure
    ('DB06292', 'dapagliflozin', 'congestive heart failure',
     'DOID:6000', '2023-02-09', 'EMA', 'Cardiovascular',
     'Heart failure across LVEF spectrum'),
    # REMOVED: empagliflozin HF EMA — actual LVEF expansion was March 2022 (pre-cutoff)

    # ── Haematology (non-oncology) ───────────────────────────────────────
    # Venetoclax — AML, MDS
    ('DB11581', 'venetoclax', 'acute myeloid leukemia',
     'DOID:9119', '2024-02-22', 'EMA', 'Oncology',
     'Newly diagnosed AML with azacitidine'),
    # REMOVED: ibrutinib CLL/SLL — original EMA approval pre-2022
    # Acalabrutinib — CLL/SLL
    ('DB11703', 'acalabrutinib', 'chronic lymphocytic leukemia',
     'DOID:1040', '2024-06-06', 'EMA', 'Oncology',
     'CLL/SLL 1L and relapsed'),
    # Daratumumab — MM expanded
    ('DB09331', 'daratumumab', 'multiple myeloma',
     'DOID:9538', '2024-04-25', 'EMA', 'Oncology',
     'Newly diagnosed transplant-eligible MM'),

    # ── Dermatology ──────────────────────────────────────────────────────
    # REMOVED: ruxolitinib alopecia areata — not approved for AA (only AD & vitiligo)
    # Ruxolitinib — vitiligo (EMA April 2023; FDA was July 2022 so pre-cutoff)
    ('DB08877', 'ruxolitinib', 'vitiligo',
     'DOID:12306', '2023-04-20', 'EMA', 'Dermatology',
     'Non-segmental vitiligo (topical)'),
]


# ── PMDA (Japan) post-2022 indication approvals ─────────────────────────────

PMDA_APPROVALS = [
    # Pembrolizumab — biliary tract cancer
    ('DB09037', 'pembrolizumab', 'biliary tract cancer',
     'DOID:4607', '2023-12-25', 'PMDA', 'Oncology',
     'Unresectable BTC with gemcitabine/cisplatin'),
    # Nivolumab — MSI-H CRC, HCC
    ('DB09035', 'nivolumab', 'colorectal cancer',
     'DOID:9256', '2024-02-22', 'PMDA', 'Oncology',
     'MSI-H/dMMR unresectable advanced CRC'),
    # Durvalumab — biliary tract cancer
    ('DB11714', 'durvalumab', 'biliary tract cancer',
     'DOID:4607', '2023-03-27', 'PMDA', 'Oncology',
     'Unresectable BTC with gemcitabine/cisplatin'),
    # Trastuzumab deruxtecan — HER2-low breast cancer
    ('DB14962', 'trastuzumab deruxtecan', 'breast cancer',
     'DOID:1612', '2023-04-12', 'PMDA', 'Oncology',
     'HER2-low unresectable/metastatic breast cancer'),
    # Olaparib — BRCA+ prostate cancer
    ('DB09074', 'olaparib', 'prostate cancer',
     'DOID:10283', '2023-03-27', 'PMDA', 'Oncology',
     'BRCA-mutated mCRPC'),
    # Dupilumab — COPD
    ('DB12159', 'dupilumab', 'chronic obstructive pulmonary disease',
     'DOID:3083', '2024-09-24', 'PMDA', 'Respiratory',
     'COPD with type 2 inflammation'),
    # Semaglutide — cardiovascular
    ('DB13928', 'semaglutide', 'cardiovascular disease',
     'DOID:1287', '2024-12-27', 'PMDA', 'Cardiovascular',
     'CV risk reduction in obese/overweight adults'),
    # Dapagliflozin — CKD (expanded indication)
    ('DB06292', 'dapagliflozin', 'chronic kidney disease',
     'DOID:784', '2023-06-26', 'PMDA', 'Metabolic',
     'CKD regardless of diabetes status'),
    # Risankizumab — ulcerative colitis
    ('DB14762', 'risankizumab', 'ulcerative colitis',
     'DOID:8577', '2024-12-27', 'PMDA', 'Immunology',
     'Moderately-to-severely active UC'),
    # Zanubrutinib — CLL/SLL
    ('DB15035', 'zanubrutinib', 'chronic lymphocytic leukemia',
     'DOID:1040', '2024-06-24', 'PMDA', 'Oncology',
     'CLL/SLL'),
    # Baricitinib — alopecia areata
    ('DB11817', 'baricitinib', 'alopecia areata',
     'DOID:986', '2023-01-23', 'PMDA', 'Immunology',
     'Severe alopecia areata in adults'),
]


# ── TGA (Australia) post-2022 indication approvals ──────────────────────────

TGA_APPROVALS = [
    # Pembrolizumab — cervical cancer, endometrial
    ('DB09037', 'pembrolizumab', 'cervical cancer',
     'DOID:4362', '2024-04-15', 'TGA', 'Oncology',
     'Persistent/recurrent/metastatic cervical cancer'),
    ('DB09037', 'pembrolizumab', 'endometrial carcinoma',
     'DOID:2870', '2024-07-22', 'TGA', 'Oncology',
     'Advanced endometrial carcinoma 1L'),
    # Nivolumab — esophageal
    ('DB09035', 'nivolumab', 'esophageal cancer',
     'DOID:5041', '2024-05-14', 'TGA', 'Oncology',
     'Resectable esophageal/GEJ cancer neoadjuvant'),
    # Trastuzumab deruxtecan — HER2-low breast cancer
    ('DB14962', 'trastuzumab deruxtecan', 'breast cancer',
     'DOID:1612', '2023-06-06', 'TGA', 'Oncology',
     'HER2-low unresectable/metastatic breast cancer'),
    # Dupilumab — COPD
    ('DB12159', 'dupilumab', 'chronic obstructive pulmonary disease',
     'DOID:3083', '2024-11-28', 'TGA', 'Respiratory',
     'COPD with type 2 inflammation'),
    # Semaglutide — CV risk
    ('DB13928', 'semaglutide', 'cardiovascular disease',
     'DOID:1287', '2025-02-20', 'TGA', 'Cardiovascular',
     'CV risk reduction in overweight/obese adults'),
    # Dapagliflozin — heart failure
    ('DB06292', 'dapagliflozin', 'congestive heart failure',
     'DOID:6000', '2023-05-15', 'TGA', 'Cardiovascular',
     'Heart failure across LVEF spectrum'),
    # Risankizumab — Crohn disease
    ('DB14762', 'risankizumab', 'Crohn disease',
     'DOID:8778', '2023-09-11', 'TGA', 'Immunology',
     'Moderately-to-severely active Crohn disease'),
    # Upadacitinib — Crohn disease
    ('DB15091', 'upadacitinib', 'Crohn disease',
     'DOID:8778', '2024-08-05', 'TGA', 'Immunology',
     'Moderately-to-severely active CD'),
    # Baricitinib — alopecia areata
    ('DB11817', 'baricitinib', 'alopecia areata',
     'DOID:986', '2023-04-17', 'TGA', 'Immunology',
     'Severe alopecia areata in adults'),
    # Olaparib — breast cancer adjuvant
    ('DB09074', 'olaparib', 'breast cancer',
     'DOID:1612', '2024-06-03', 'TGA', 'Oncology',
     'Adjuvant HER2- BRCA-mutated early breast cancer'),
]


# ── Existing FDA approvals from notebook 07 ─────────────────────────────────

FDA_APPROVALS = [
    ('DB12159', 'dupilumab', 'chronic obstructive pulmonary disease',
     'DOID:3083', '2024-09-27', 'FDA', 'Respiratory',
     'COPD with eosinophilic phenotype'),
    ('DB12159', 'dupilumab', 'prurigo nodularis',
     'DOID:0060845', '2024-10-25', 'FDA', 'Immunology',
     'Prurigo nodularis in adults'),
    ('DB13928', 'semaglutide', 'cardiovascular disease',
     'DOID:1287', '2024-03-08', 'FDA', 'Cardiovascular',
     'CV risk reduction in overweight/obese adults'),
    ('DB13928', 'semaglutide', 'obstructive sleep apnea',
     'DOID:0050848', '2024-12-13', 'FDA', 'Respiratory',
     'Moderate-to-severe OSA with obesity'),
    ('DB09037', 'pembrolizumab', 'endometrial carcinoma',
     'DOID:2870', '2024-06-07', 'FDA', 'Oncology',
     'Advanced endometrial carcinoma 1L'),
    ('DB09037', 'pembrolizumab', 'cervical cancer',
     'DOID:4362', '2024-01-15', 'FDA', 'Oncology',
     'Persistent/recurrent/metastatic cervical cancer'),
    ('DB09037', 'pembrolizumab', 'gastric cancer',
     'DOID:10534', '2024-10-25', 'FDA', 'Oncology',
     'HER2+ gastric/GEJ adenocarcinoma 1L'),
    ('DB09037', 'pembrolizumab', 'gastroesophageal junction adenocarcinoma',
     'DOID:4944', '2024-10-25', 'FDA', 'Oncology',
     'HER2+ GEJ adenocarcinoma'),
    ('DB09037', 'pembrolizumab', 'hepatocellular carcinoma',
     'DOID:684', '2024-09-03', 'FDA', 'Oncology',
     'HCC 1L with lenvatinib'),
    ('DB09035', 'nivolumab', 'hepatocellular carcinoma',
     'DOID:684', '2024-08-15', 'FDA', 'Oncology',
     'HCC 1L with ipilimumab'),
    ('DB09035', 'nivolumab', 'colorectal cancer',
     'DOID:9256', '2024-07-25', 'FDA', 'Oncology',
     'MSI-H/dMMR CRC'),
    ('DB09035', 'nivolumab', 'esophageal cancer',
     'DOID:5041', '2024-03-22', 'FDA', 'Oncology',
     'Resectable esophageal/GEJ cancer neoadjuvant'),
    ('DB11714', 'durvalumab', 'small cell lung carcinoma',
     'DOID:0050685', '2024-12-20', 'FDA', 'Oncology',
     'Limited-stage SCLC'),
    ('DB11714', 'durvalumab', 'biliary tract cancer',
     'DOID:4607', '2024-09-05', 'FDA', 'Oncology',
     'Biliary tract cancer 1L'),
    # REMOVED: atezolizumab HCC — original approval was May 2020 (pre-2022)
    ('DB00051', 'adalimumab', 'hidradenitis suppurativa',
     'DOID:0060480', '2024-04-30', 'FDA', 'Immunology',
     'HS in adolescents 12+ years'),
    # REMOVED: ibrutinib CLL — original approval 2014; 2023 was oral suspension formulation only
    # REMOVED: ibrutinib FL — never approved; MCL/MZL indications withdrawn 2023
    ('DB11581', 'venetoclax', 'acute myeloid leukemia',
     'DOID:9119', '2024-01-12', 'FDA', 'Oncology',
     'Newly diagnosed AML with azacitidine'),
    # REMOVED: venetoclax MDS — only has Breakthrough Therapy Designation, not approval
    ('DB11703', 'acalabrutinib', 'chronic lymphocytic leukemia',
     'DOID:1040', '2024-05-20', 'FDA', 'Oncology',
     'CLL/SLL 1L and relapsed'),
    ('DB09074', 'olaparib', 'prostate cancer',
     'DOID:10283', '2023-12-19', 'FDA', 'Oncology',
     'mCRPC with HRR gene mutations'),
    ('DB09074', 'olaparib', 'breast cancer',
     'DOID:1612', '2024-03-11', 'FDA', 'Oncology',
     'Adjuvant HER2- high-risk early breast cancer'),
    # REMOVED: trastuzumab gastric — original approval 2010, not a new indication
    ('DB09078', 'lenvatinib', 'endometrial cancer',
     'DOID:1380', '2024-06-07', 'FDA', 'Oncology',
     'Advanced endometrial carcinoma with pembrolizumab'),
    ('DB06292', 'dapagliflozin', 'congestive heart failure',
     'DOID:6000', '2023-05-09', 'FDA', 'Cardiovascular',
     'Heart failure across LVEF spectrum'),
    # REMOVED: empagliflozin HF FDA — actual LVEF expansion was Feb 2022 (pre-cutoff)
    # REMOVED: ruxolitinib alopecia areata — not approved for AA
    # REMOVED: ruxolitinib vitiligo FDA — actual approval was July 2022 (pre-cutoff)
    ('DB08879', 'belimumab', 'lupus nephritis',
     'DOID:0080162', '2024-06-01', 'FDA', 'Immunology',
     'Active lupus nephritis expanded'),
    ('DB15091', 'upadacitinib', 'ulcerative colitis',
     'DOID:8577', '2023-08-16', 'FDA', 'Immunology',
     'Moderately-to-severely active UC'),
    ('DB15091', 'upadacitinib', 'Crohn disease',
     'DOID:8778', '2024-05-17', 'FDA', 'Immunology',
     'Moderately-to-severely active CD'),
    ('DB15091', 'upadacitinib', 'ankylosing spondylitis',
     'DOID:7147', '2024-02-13', 'FDA', 'Immunology',
     'Active AS'),
    ('DB15091', 'upadacitinib', 'atopic eczema',
     'DOID:3310', '2024-01-25', 'FDA', 'Immunology',
     'Moderate-to-severe AD expanded age'),
    ('DB14762', 'risankizumab', 'Crohn disease',
     'DOID:8778', '2023-08-02', 'FDA', 'Immunology',
     'Moderately-to-severely active Crohn disease'),
    ('DB14762', 'risankizumab', 'ulcerative colitis',
     'DOID:8577', '2024-09-19', 'FDA', 'Immunology',
     'Moderately-to-severely active UC'),
    ('DB06273', 'tocilizumab', 'polymyalgia rheumatica',
     'DOID:853', '2024-03-06', 'FDA', 'Immunology',
     'Giant cell arteritis / polymyalgia rheumatica'),
    ('DB08875', 'cabozantinib', 'hepatocellular carcinoma',
     'DOID:684', '2024-02-28', 'FDA', 'Oncology',
     'Previously treated HCC'),
    ('DB09331', 'daratumumab', 'multiple myeloma',
     'DOID:9538', '2024-03-25', 'FDA', 'Oncology',
     'Newly diagnosed transplant-eligible MM expanded'),
    ('DB01590', 'everolimus', 'gastrointestinal neuroendocrine tumor',
     'DOID:0080909', '2024-05-09', 'FDA', 'Oncology',
     'Advanced GI NETs expanded'),
]

# Column order for the DataFrame
_COLUMNS = ['drug_id', 'drug_name', 'disease_name', 'disease_id',
            'approval_date', 'source_agency', 'therapeutic_area', 'notes']


def _tuples_to_df(tuples):
    """Convert list of 8-tuples to DataFrame."""
    return pd.DataFrame(tuples, columns=_COLUMNS)


def get_all_approvals():
    """Return a DataFrame of all regulatory approvals across agencies."""
    all_data = FDA_APPROVALS + EMA_APPROVALS + PMDA_APPROVALS + TGA_APPROVALS
    df = _tuples_to_df(all_data)
    df['approval_date'] = pd.to_datetime(df['approval_date'])
    return df


def get_unique_pairs(deduplicate_agencies=True):
    """Return deduplicated (drug_id, disease_id) pairs.

    When a drug-disease pair is approved by multiple agencies, keeps the
    earliest approval date and lists all agencies.

    Parameters
    ----------
    deduplicate_agencies : bool
        If True, merge rows sharing the same (drug_id, disease_id) across
        agencies. If False, return one row per agency approval.

    Returns
    -------
    pd.DataFrame with columns:
        drug_id, drug_name, disease_name, disease_id, approval_date,
        source_agencies (comma-separated), therapeutic_area, n_agencies
    """
    df = get_all_approvals()

    if not deduplicate_agencies:
        return df

    # Group by (drug_id, disease_id) — take earliest date, merge agencies
    grouped = (df.groupby(['drug_id', 'disease_id'])
               .agg(
                   drug_name=('drug_name', 'first'),
                   disease_name=('disease_name', 'first'),
                   approval_date=('approval_date', 'min'),
                   source_agencies=('source_agency',
                                    lambda x: ','.join(sorted(set(x)))),
                   therapeutic_area=('therapeutic_area', 'first'),
                   n_agencies=('source_agency', 'nunique'),
                   notes=('notes', 'first'),
               )
               .reset_index()
               .sort_values('approval_date'))

    return grouped


def get_non_fda_pairs():
    """Return pairs that come from EMA/PMDA/TGA but NOT the existing FDA set.

    These are the net-new pairs added by expanding beyond FDA.
    """
    fda_pairs = {(r[0], r[3]) for r in FDA_APPROVALS}
    other = EMA_APPROVALS + PMDA_APPROVALS + TGA_APPROVALS
    new_entries = [r for r in other if (r[0], r[3]) not in fda_pairs]
    return _tuples_to_df(new_entries)


def summary():
    """Print summary statistics of the expanded gold standard."""
    df = get_all_approvals()
    unique = get_unique_pairs()
    non_fda = get_non_fda_pairs()

    print(f'Total approval records: {len(df)}')
    print(f'  FDA: {len(FDA_APPROVALS)}')
    print(f'  EMA: {len(EMA_APPROVALS)}')
    print(f'  PMDA: {len(PMDA_APPROVALS)}')
    print(f'  TGA: {len(TGA_APPROVALS)}')
    print(f'\nUnique (drug, disease) pairs: {len(unique)}')
    print(f'  Multi-agency pairs: {(unique["n_agencies"] > 1).sum()}')
    print(f'  Net-new pairs (not in FDA set): '
          f'{len(non_fda["drug_id"].astype(str) + non_fda["disease_id"].astype(str))}')
    print(f'\nUnique drugs: {unique["drug_id"].nunique()}')
    print(f'Unique diseases: {unique["disease_id"].nunique()}')
    print(f'\nTherapeutic areas:')
    for area, count in unique['therapeutic_area'].value_counts().items():
        print(f'  {area}: {count}')


def to_notebook_format():
    """Convert to the dict format used by 07_generalization.ipynb.

    Returns
    -------
    dict : {drugbank_id: (drug_name, [(disease_term, approval_date_str, notes), ...])}
        Grouped by drug, matching the NEW_INDICATION_APPROVALS format in the
        notebook. Includes all agencies (FDA + EMA + PMDA + TGA), deduplicated
        per (drug, disease) pair.
    """
    pairs = get_unique_pairs(deduplicate_agencies=True)
    result = {}
    for _, row in pairs.iterrows():
        dbid = row['drug_id']
        drug_name = row['drug_name']
        date_str = row['approval_date'].strftime('%Y-%m-%d')
        note = f"{row['notes']} [{row['source_agencies']}]"
        disease_term = row['disease_name']
        if dbid not in result:
            result[dbid] = (drug_name, [])
        result[dbid][1].append((disease_term, date_str, note))
    return result


def load_expanded_approvals(base_dir=None):
    """Load non-FDA approval pairs for notebook integration.

    Returns a DataFrame with columns: drugbank_id, drug_name, disease_name,
    doid, approval_date, source_agency, therapeutic_area, notes.

    Only returns EMA/PMDA/TGA pairs that are NOT already in the FDA set,
    avoiding double-counting when merged with the notebook's FDA loader.
    """
    non_fda = get_non_fda_pairs()
    # Rename columns to match notebook expectations
    non_fda = non_fda.rename(columns={
        'drug_id': 'drugbank_id',
        'disease_id': 'doid',
    })
    return non_fda


if __name__ == '__main__':
    summary()
