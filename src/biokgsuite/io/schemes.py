"""Check 1: does an identifier's shape fit the entity type its graph claims?

A node labelled Drug whose id is ``DB12494`` is recognised by shape: DB plus
five digits is a DrugBank accession, so it is a valid drug identifier. This
rejects a DrugBank accession filed under Anatomy.

It cannot *confirm* a type. Some shapes are shared, and the table says which:
a bare integer is a valid Entrez gene and a valid bare PubChem CID; all three
GO aspects share ``GO:nnnnnnn``. Confirming is check 2's job (resolution
against a registry), which notebook 02 reports for Drug, Disease and
Gene/Protein and deliberately does NOT score -- see REGISTRY_CHECKED in
``biokgsuite.stats.annotation`` for why.

The schemes live in ``config/id_schemes.yaml`` rather than in code, so the list
can be read and argued with. See that file for the rule that keeps it from
growing into a list that admits every shape.
"""
from __future__ import annotations

import re
from pathlib import Path

import yaml

CONFIG = 'id_schemes.yaml'


def _find_config(name=CONFIG):
    for base in Path(__file__).resolve().parents:
        cand = base / 'config' / name
        if cand.exists():
            return cand
    raise FileNotFoundError(f'config/{name} not found above {Path(__file__)}')


class IdSchemes:
    """Identifier shapes per entity type, loaded once from config."""

    def __init__(self, path=None):
        self.path = Path(path) if path else _find_config()
        raw = yaml.safe_load(self.path.read_text())['id_schemes']
        self.schemes = {
            etype: [(e['scheme'], re.compile(e['pattern']), bool(e.get('ambiguous')))
                    for e in entries]
            for etype, entries in raw.items()}
        # Schemes that name concepts of every kind (UMLS, MeSH, EFO). Their
        # shape is recognised, but fitting one under another type is not
        # evidence of a mislabelling: see `types_evidencing`.
        self.agnostic = {
            etype: {e['scheme'] for e in entries if e.get('agnostic')}
            for etype, entries in raw.items()}
        self._raw = raw

    # ------------------------------------------------------------- matching
    @staticmethod
    def _bare(raw) -> str:
        """Strip a 'Type::' qualifier; the shape is matched on what remains."""
        return str(raw).split('::')[-1].strip()

    def match(self, raw, entity_type):
        """(scheme, ambiguous) for the first shape this type accepts, else (None, False)."""
        s = self._bare(raw)
        for scheme, rx, ambiguous in self.schemes.get(entity_type, ()):
            if rx.match(s):
                return scheme, ambiguous
        return None, False

    def types_matching(self, raw) -> set:
        """Every entity type whose declared shapes this identifier fits."""
        s = self._bare(raw)
        return {etype for etype, entries in self.schemes.items()
                if any(rx.match(s) for _, rx, _ in entries)}

    def types_evidencing(self, raw) -> set:
        """Types whose TYPE-SPECIFIC shapes this identifier fits.

        Used to tell a mislabelling from an identifier the check cannot judge:
        a DrugBank accession under Anatomy fits a drug-only shape and is a
        mislabelling; a UMLS CUI under Anatomy fits only schemes that name
        every kind of concept, and says nothing about the type.
        """
        s = self._bare(raw)
        return {etype for etype, entries in self.schemes.items()
                if any(rx.match(s) and scheme not in self.agnostic.get(etype, ())
                       for scheme, rx, _ in entries)}

    # ------------------------------------------------------------ reporting
    def table(self):
        """The declared schemes, for printing beside the results."""
        return [{'Entity type': etype, 'Scheme': e['scheme'], 'Pattern': e['pattern'],
                 'Example': e['example'], 'Shared shape': bool(e.get('ambiguous')),
                 'Type-agnostic': bool(e.get('agnostic')), 'Note': e.get('note', '')}
                for etype, entries in self._raw.items() for e in entries]

    def shared_shapes(self):
        """Shapes accepted by more than one entity type.

        The honest weakness of shape matching, reported rather than buried: a
        node passing check 1 on one of these has only been shown to look like
        something, not to be the thing its graph says it is.
        """
        by_pattern = {}
        for etype, entries in self._raw.items():
            for e in entries:
                by_pattern.setdefault(e['pattern'], set()).add(etype)
        return [{'Pattern': p, 'Entity types': ', '.join(sorted(t)), 'N types': len(t)}
                for p, t in sorted(by_pattern.items()) if len(t) > 1]

    def self_check(self):
        """Every declared example must match its own pattern."""
        return [{'Entity type': etype, 'Scheme': e['scheme'], 'Pattern': e['pattern'],
                 'Example': e['example']}
                for etype, entries in self._raw.items() for e in entries
                if not re.match(e['pattern'], str(e['example']))]
