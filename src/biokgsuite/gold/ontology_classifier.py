"""Ontology-based disease domain classification.

Two-tier approach:
  1. MONDO hierarchy traversal (gold-standard: uses actual ontology graph)
  2. MeSH C-category term matching (for diseases not in MONDO hierarchy)

Replaces the ad-hoc keyword classifier in notebook 07 with a more rigorous,
ontology-grounded approach.
"""
from __future__ import annotations
import pandas as pd
from collections import defaultdict
from pathlib import Path
from typing import Optional, Tuple, Dict, Set

# ── MONDO top-level → domain mapping ──────────────────────────────────
# Direct children of MONDO:0700096 (human disease)
MONDO_DOMAIN_MAP = {
    'MONDO:0005550': 'Infectious',       # infectious disease
    'MONDO:0045024': 'Oncology',         # cancer or benign tumor
    'MONDO:0004995': 'Cardiovascular',   # cardiovascular disorder
    'MONDO:0005071': 'Neurology',        # nervous system disorder
    'MONDO:0005066': 'Metabolic',        # metabolic disease
    'MONDO:0005151': 'Metabolic',        # endocrine system disorder
    'MONDO:0005046': 'Immunology',       # immune system disorder
    'MONDO:0021166': 'Immunology',       # inflammatory disease
    'MONDO:0005087': 'Respiratory',      # respiratory system disorder
    'MONDO:0002081': 'Musculoskeletal',  # musculoskeletal system disorder
    'MONDO:0005084': 'Psychiatry',       # mental disorder
}

# ── MeSH C-subcategory grounded term lists ─────────────────────────────
# Derived from MeSH Category C (Diseases) subcategory scope notes.
# Priority order follows MeSH primary tree assignment convention.
MESH_C_TERMS = [
    ('Oncology', [
        'neoplasm', 'cancer', 'carcinoma', 'tumor', 'tumour', 'lymphoma',
        'leukemia', 'leukaemia', 'melanoma', 'sarcoma', 'myeloma', 'glioma',
        'glioblastoma', 'blastoma', 'adenocarcinoma', 'mesothelioma',
        'malignant', 'myelodysplastic', 'polycythemia vera', 'myelofibrosis',
        'mastocytosis', 'neuroblastoma', 'retinoblastoma', 'medulloblastoma',
        'cholangiocarcinoma', 'hepatocellular carcinoma',
    ]),
    ('Infectious', [
        'infection', 'infectious', 'bacterial', 'viral', 'mycosis', 'mycoses',
        'parasitic', 'tuberculosis', 'hepatitis b', 'hepatitis c', 'hiv',
        'aids', 'malaria', 'sepsis', 'septicemia',
        'candidiasis', 'aspergillosis', 'herpes', 'influenza', 'covid',
        'dengue', 'ebola', 'zika', 'chlamydia', 'gonorrhea',
        'syphilis', 'trypanosomiasis', 'leishmaniasis', 'schistosomiasis',
        'pneumococcal', 'staphylococcal', 'streptococcal',
    ]),
    ('Psychiatry', [
        'mental disorder', 'psychiatric', 'schizophrenia', 'bipolar disorder',
        'major depressive', 'depressive disorder', 'anxiety disorder',
        'panic disorder', 'phobia', 'obsessive-compulsive', 'post-traumatic stress',
        'attention deficit', 'adhd', 'autism spectrum', 'eating disorder',
        'anorexia nervosa', 'bulimia', 'substance use', 'addiction',
        'alcohol dependence', 'opioid', 'psychosis', 'psychotic',
        'personality disorder', 'dissociative', 'conduct disorder',
    ]),
    ('Neurology', [
        'nervous system disease', 'neurodegenerative', 'alzheimer',
        'parkinson', 'epilepsy', 'seizure disorder', 'dementia',
        'multiple sclerosis', 'amyotrophic lateral', 'huntington',
        'peripheral neuropathy', 'migraine', 'cerebral palsy',
        'encephalopathy', 'myelitis', 'neuralgia', 'ataxia', 'dystonia',
        'chorea', 'myasthenia gravis', 'spinal muscular atrophy',
        'charcot-marie-tooth', 'guillain-barre',
    ]),
    ('Cardiovascular', [
        'cardiovascular', 'heart disease', 'heart failure', 'cardiac',
        'coronary', 'myocardial', 'atrial fibrillation', 'arrhythmia',
        'hypertension', 'hypotension', 'atherosclerosis', 'aneurysm',
        'thrombosis', 'pulmonary embolism', 'deep vein thrombosis',
        'stroke', 'cerebrovascular', 'cardiomyopathy', 'pericarditis',
        'aortic', 'peripheral arterial', 'angina', 'valvular heart',
    ]),
    ('Metabolic', [
        'metabolic disease', 'diabetes mellitus', 'diabetes type',
        'diabetic', 'obesity', 'hyperlipidemia', 'hypercholesterolemia',
        'gout', 'hyperuricemia', 'porphyria', 'glycogen storage',
        'lysosomal storage', 'phenylketonuria', 'galactosemia',
        'hypothyroidism', 'hyperthyroidism', 'thyroid disease',
        'cushing', 'addison disease', 'metabolic syndrome',
        'amyloidosis', 'wilson disease', 'hemochromatosis',
        'non-alcoholic steatohepatitis', 'fatty liver disease',
        'gaucher disease', 'fabry disease', 'mucopolysaccharidosis',
    ]),
    ('Immunology', [
        'autoimmune', 'immunodeficiency', 'allergy', 'allergic',
        'hypersensitivity', 'lupus erythematosus', 'rheumatoid arthritis',
        'systemic sclerosis', 'scleroderma', 'vasculitis',
        'psoriasis', 'atopic dermatitis', 'eczema', 'crohn disease',
        'ulcerative colitis', 'celiac disease', 'ankylosing spondylitis',
        'sarcoidosis', 'graft-versus-host', 'transplant rejection',
        'inflammatory bowel', 'sjogren', 'dermatomyositis', 'polymyositis',
        'immune thrombocytopeni', 'pemphigus',
    ]),
    ('Respiratory', [
        'respiratory disease', 'pulmonary disease', 'lung disease',
        'asthma', 'chronic obstructive pulmonary', 'copd', 'bronchitis',
        'emphysema', 'pulmonary fibrosis', 'idiopathic pulmonary fibrosis',
        'pulmonary hypertension', 'cystic fibrosis', 'bronchiectasis',
        'sleep apnea', 'interstitial lung', 'acute respiratory distress',
        'pneumothorax', 'pleural effusion',
    ]),
    ('Musculoskeletal', [
        'musculoskeletal', 'osteoarthritis', 'osteoporosis', 'bone disease',
        'skeletal dysplasia', 'muscular dystrophy', 'myopathy',
        'tendinitis', 'bursitis', 'fibromyalgia', 'osteogenesis imperfecta',
        'scoliosis', 'disc degeneration', 'spinal stenosis',
    ]),
]

# All valid domain labels (including the 9th: Musculoskeletal)
ONTOLOGY_DOMAINS = [d for d, _ in MESH_C_TERMS]


class OntologyClassifier:
    """Classify disease nodes into therapeutic domains using ontology data."""

    def __init__(self, base_dir: Path):
        self._load_reference_data(base_dir)

    def _load_reference_data(self, base_dir: Path):
        gs = base_dir / 'data' / 'gold_standards' / 'bridges'
        do_df = pd.read_csv(gs / 'do_diseases.csv')
        self.doid_to_name = dict(zip(do_df['doid'].astype(str),
                                     do_df['mondo_name'].astype(str)))
        self.doid_to_mondo = dict(zip(do_df['doid'].astype(str),
                                      do_df['mondo_id'].astype(str)))

        mesh_df = pd.read_csv(gs / 'mesh_to_doid.csv', on_bad_lines='skip')
        mesh_df['mesh_clean'] = mesh_df['mesh_id'].str.split(' ').str[0]
        self.mesh_to_doid = (mesh_df.drop_duplicates(subset='mesh_clean')
                             .set_index('mesh_clean')['doid'].to_dict())

        hier = pd.read_csv(gs / 'mondo_hierarchy.csv')
        self.child_to_parents: Dict[str, Set[str]] = defaultdict(set)
        for _, row in hier.iterrows():
            self.child_to_parents[str(row['child_id'])].add(str(row['parent_id']))

    # ── MONDO hierarchy ────────────────────────────────────────────────
    def _ancestors(self, mondo_id: str, max_depth: int = 15) -> set:
        visited: set = set()
        frontier = {mondo_id}
        for _ in range(max_depth):
            nxt: set = set()
            for nid in frontier:
                if nid in visited:
                    continue
                visited.add(nid)
                for p in self.child_to_parents.get(nid, set()):
                    nxt.add(p)
            frontier = nxt
            if not frontier:
                break
        return visited

    def _classify_mondo(self, mondo_id: str) -> Optional[str]:
        anc = self._ancestors(mondo_id)
        for mid, domain in MONDO_DOMAIN_MAP.items():
            if mid in anc:
                return domain
        return None

    # ── MeSH C-category term matching ──────────────────────────────────
    @staticmethod
    def _classify_mesh_terms(disease_name: str) -> Optional[str]:
        name_lower = disease_name.lower()
        for domain, terms in MESH_C_TERMS:
            for term in terms:
                if term in name_lower:
                    return domain
        return None

    # ── Public API ─────────────────────────────────────────────────────
    def resolve_name(self, node_id: str) -> str:
        """Resolve a disease node ID to a human-readable name."""
        if node_id.startswith('DOID:'):
            return self.doid_to_name.get(node_id, '')
        if node_id.startswith('MESH:'):
            doid = self.mesh_to_doid.get(node_id)
            return self.doid_to_name.get(doid, '') if doid else ''
        cand = self.doid_to_name.get('DOID:' + node_id, '')
        if cand:
            return cand
        return node_id  # plain name (PrimeKG / Hetionet)

    def classify_mondo(self, mondo_id: str, disease_name: str = '') -> Tuple[Optional[str], str]:
        """Classify a disease already keyed by MONDO id.

        The post-cutoff approval set carries MONDO ids directly, so the DOID and
        MeSH resolution in ``classify`` would only lose terms that have no DOID.
        Same two tiers: MONDO hierarchy first, MeSH C-category terms on the name
        as the fallback.
        """
        r = self._classify_mondo(str(mondo_id)) if mondo_id else None
        if r:
            return r, 'mondo'
        r = self._classify_mesh_terms(disease_name)
        return (r, 'mesh_terms') if r else (None, 'unclassified')

    def classify(self, node_id: str, disease_name: str = '') -> Tuple[Optional[str], str]:
        """Classify a disease into a therapeutic domain.

        Returns (domain, method) where method is 'mondo', 'mesh_terms',
        or 'unclassified'.
        """
        # Resolve MONDO ID for hierarchy traversal
        mondo_id = None
        if node_id.startswith('DOID:'):
            mondo_id = self.doid_to_mondo.get(node_id)
        elif node_id.startswith('MESH:'):
            doid = self.mesh_to_doid.get(node_id)
            if doid:
                mondo_id = self.doid_to_mondo.get(doid)
        else:
            mondo_id = self.doid_to_mondo.get('DOID:' + node_id)

        if mondo_id:
            r = self._classify_mondo(mondo_id)
            if r:
                return r, 'mondo'

        # Fallback: MeSH C-category term matching
        name = disease_name or self.resolve_name(node_id)
        r = self._classify_mesh_terms(name)
        if r:
            return r, 'mesh_terms'

        return None, 'unclassified'
