"""The ranking loop: every query, in every arm, for every model, seed and shuffle.

One call to a model ranks one candidate pool for one disease. The only thing
that differs between the ``none`` arm and a KG arm is the evidence attached to
each candidate (and the disease profile in the prompt header), so a lift over
the ``none`` arm is attributable to that graph's evidence.

Results are written per (family, model, seed) so that re-running one model, or
adding a seed, never rewrites the others. Every row carries the pair id and
the arm, which is what lets the analysis pair queries across arms without
replaying the sampling.
"""
from __future__ import annotations

import re
import string
from pathlib import Path

import numpy as np
import pandas as pd

from .clients import rank_call
from .prompt import build_prompt, parse_response, rank_of

LETTERS = list(string.ascii_uppercase)
MOCK_EVIDENCE_P = {True: 0.9, False: 0.4}   # chance a positive / distractor gets a mock dossier
RUN_COLUMNS = ['mode', 'family', 'model', 'arm', 'condition', 'seed', 'shuffle', 'pair_id', 'drugbank_id',
               'mondo_id', 'disease', 'drug', 'pool_size', 'positive_letter', 'rank', 'reciprocal_rank',
               'hit1', 'hit3', 'hit5', 'parsed', 'error', 'pos_has_evidence', 'n_with_evidence',
               'pos_basis', 'pos_agreement', 'pos_confidence', 'pos_rationale', 'frac_basis_kg',
               'mean_confidence']


def model_slug(model: str) -> str:
    return re.sub(r'[^a-z0-9]+', '_', model.lower()).strip('_')


def run_path(out_dir, family: str, model: str, seed: int) -> Path:
    return Path(out_dir) / f'{family}__{model_slug(model)}__seed{seed}.csv'


def _shuffle_rng(seed: int, pair_id: int, shuffle: int, arm: str) -> np.random.Generator:
    """One stream per (seed, pair, shuffle, arm), so a candidate order depends on
    nothing that ran before it and a partial rerun reproduces the same prompts."""
    arm_code = sum(ord(c) for c in arm) % 97
    return np.random.default_rng(seed * 1_000_003 + pair_id * 1_009 + shuffle * 101 + arm_code)


def candidate_view(candidates, *, arm, dossiers, disease_id, rng, mock, anonymize):
    """Shuffle the pool, assign letters and attach evidence. Returns (view, positive letter)."""
    order = rng.permutation(len(candidates))
    code = None
    if anonymize:
        # Names are replaced by Drug-1, Drug-2, ... so a model cannot recall the
        # approval from the name alone; the numbering follows the shuffled order.
        code = {candidates[i]['drug_id']: f'Drug-{k + 1}' for k, i in enumerate(order)}
    view, pos_letter = [], None
    for k, i in enumerate(order):
        c = candidates[i]
        letter = LETTERS[k]
        name = code[c['drug_id']] if code else c['drug_name']
        evidence = None
        if arm != 'none':
            if mock:
                evidence = 'mechanistic link present' if rng.random() < MOCK_EVIDENCE_P[c['is_pos']] else None
            else:
                evidence = dossiers.drug_dossier(c['drug_id'], name, disease_id)
        if c['is_pos']:
            pos_letter = letter
        view.append({'letter': letter, 'drug_name': name, 'evidence': evidence})
    return view, pos_letter


# An error that no retry will cure: a wrong or inaccessible model, a missing
# or rejected key, a malformed request, exhausted prepaid credits (402) or an
# exhausted daily quota. One such
# reply means every further call in the arm would fail the same way.
FATAL_PREFIXES = ('400', '401', '402', '403', '404', 'no key', 'not an API model', '429 daily quota')


def is_fatal_error(error):
    return error is not None and error.startswith(FATAL_PREFIXES)


class FatalAPIError(RuntimeError):
    """Raised by ``run_family`` after saving what was scored, so a job with a
    dead model stops at once instead of filing a placeholder row per call."""


def run_arm(*, family, arm, queries, model, seed, shuffles, pool_size, dossiers, mock, anonymize,
            temperature, api_delay, example_dir=None, progress=None, abort_after=8, skip=frozenset()):
    """All queries in one arm for one model and seed. Returns (rows, stopped),
    where ``stopped`` is None, ``'rate_limited'`` or ``'fatal: <error>'``.
    ``skip`` holds (pair_id, shuffle) already scored, which are not rerun."""
    condition = 'no_kg' if arm == 'none' else 'kg'
    rows, rate_limited, example_written = [], 0, False
    for q in queries:
        for s in range(shuffles):
            if (q['pair_id'], s) in skip:
                continue
            rng = _shuffle_rng(seed, q['pair_id'], s, arm)
            view, pos_letter = candidate_view(q['candidates'], arm=arm, dossiers=dossiers,
                                              disease_id=q['mondo_id'], rng=rng, mock=mock, anonymize=anonymize)
            profile = None
            if arm != 'none' and not mock:
                profile = dossiers.disease_profile(q['mondo_id'], q['disease_name']) or None
            prompt = build_prompt(q['disease_name'], q['mondo_id'], view, disease_profile=profile)
            if example_dir is not None and not example_written:
                # The first prompt per arm is the worked example in the supplement.
                Path(example_dir).mkdir(parents=True, exist_ok=True)
                (Path(example_dir) / f'{family}__{arm}.txt').write_text(
                    f"# {q['disease_name']} ({q['mondo_id']})  arm={arm}  positive={pos_letter}\n\n{prompt}\n")
                example_written = True
            response = rank_call(model, prompt, view, seed * 100_000 + q['pair_id'] * 10 + s,
                                 mock=mock, temperature=temperature, api_delay=api_delay)
            error = (response[len('__ERROR__'):].strip()[:200]
                     if isinstance(response, str) and response.startswith('__ERROR__') else None)
            if is_fatal_error(error):
                return rows, f'fatal: {error}'
            rate_limited = rate_limited + 1 if (error and '429' in error) else 0
            ordered, fields = parse_response(response, [c['letter'] for c in view])
            r = rank_of(pos_letter, ordered, pool_size)
            bases = [f['basis'] for f in fields.values() if f.get('basis')]
            confs = [f['confidence'] for f in fields.values() if isinstance(f.get('confidence'), int)]
            pf = fields.get(pos_letter, {})
            rows.append({
                'mode': 'mock' if mock else 'real',
                'family': family, 'model': model, 'arm': arm, 'condition': condition, 'seed': seed,
                'shuffle': s, 'pair_id': q['pair_id'], 'drugbank_id': q['drugbank_id'],
                'mondo_id': q['mondo_id'], 'disease': q['disease_name'], 'drug': q['drug_name'],
                'pool_size': pool_size, 'positive_letter': pos_letter, 'rank': r,
                'reciprocal_rank': 1.0 / r, 'hit1': int(r <= 1), 'hit3': int(r <= 3), 'hit5': int(r <= 5),
                'parsed': int(bool(ordered)), 'error': error,
                'pos_has_evidence': bool(next(c['evidence'] for c in view if c['letter'] == pos_letter)),
                'n_with_evidence': sum(bool(c['evidence']) for c in view),
                'pos_basis': pf.get('basis'), 'pos_agreement': pf.get('evidence_agreement'),
                'pos_confidence': pf.get('confidence'), 'pos_rationale': pf.get('rationale'),
                'frac_basis_kg': (sum(b in ('kg', 'both') for b in bases) / len(bases)) if bases else None,
                'mean_confidence': (sum(confs) / len(confs)) if confs else None,
            })
            if progress is not None:
                progress.update(1)
            if rate_limited >= abort_after:
                # Every further call would fail the same way; keep what was scored
                # and let the caller decide whether to retry the arm later.
                return rows, 'rate_limited'
    return rows, None


def run_family(*, family, queries, models, kgs, dossiers, mock, seed, shuffles, pool_size, anonymize,
               temperature, api_delay, out_dir, example_dir=None, force=False, progress_factory=None):
    """Run every model over the ``none`` arm and each KG arm; one CSV per model.

    ``dossiers`` maps KG name to a ``KGDossiers`` (or to None in mock mode). An
    existing file for (family, model, seed) is loaded, not recomputed, unless
    ``force``. Only rows without an error count as scored: an arm that aborted
    on rate limits, or whose calls failed, is left incomplete in the file and
    reported; the next launch keeps every scored row and runs only the
    (pair, shuffle) calls that are missing or errored, so a rate-limited job
    is finished by resubmitting it, and ``force`` is only for redoing finished
    work. A fatal API error (wrong model, no key, daily quota) saves the file
    and raises ``FatalAPIError``.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    frames, aborted = [], []
    arms = ['none'] + list(kgs)
    per_arm = len(queries) * shuffles
    for model in models:
        path = run_path(out_dir, family, model, seed)
        rows, todo, done = [], arms, {}
        if path.exists() and not force:
            prev = pd.read_csv(path)
            ok = prev[prev['error'].isna()] if len(prev) else prev
            n_err = len(prev) - len(ok)
            counts = ok.groupby('arm').size() if len(ok) else {}
            complete = [a for a in arms if counts.get(a, 0) == per_arm]
            if len(complete) == len(arms):
                frames.append(ok.reset_index(drop=True))
                continue
            rows = ok.to_dict('records')
            done = {a: set(zip(g.pair_id, g.shuffle)) for a, g in ok.groupby('arm')}
            todo = [a for a in arms if a not in complete]
            print(f'{path.name}: {len(complete)} arms complete, completing {todo} '
                  f'({len(ok)} rows kept, {n_err} errored rows redone)')
        n_calls = sum(per_arm - len(done.get(a, ())) for a in todo)
        progress = progress_factory(n_calls, f'{family} {model} seed {seed}') if progress_factory else None
        fatal = None
        for arm in todo:
            arm_rows, stopped = run_arm(family=family, arm=arm, queries=queries, model=model, seed=seed,
                                        shuffles=shuffles, pool_size=pool_size,
                                        dossiers=None if arm == 'none' else dossiers.get(arm),
                                        mock=mock, anonymize=anonymize, temperature=temperature,
                                        api_delay=api_delay, example_dir=example_dir, progress=progress,
                                        skip=frozenset(done.get(arm, ())))
            rows += arm_rows
            if stopped:
                aborted.append((model, arm))
                if stopped.startswith('fatal'):
                    fatal = f'{model} on arm {arm}: {stopped}'
                break
        if progress is not None:
            progress.close()
        df = pd.DataFrame(rows, columns=RUN_COLUMNS)
        df = df.sort_values(['arm', 'pair_id', 'shuffle'], key=lambda c: c.map(arms.index) if c.name == 'arm' else c,
                            kind='stable').reset_index(drop=True)
        df.to_csv(path, index=False)
        frames.append(df)
        if fatal:
            raise FatalAPIError(f'{fatal}; {len(df)} scored rows saved to {path.name}')
    return pd.concat(frames, ignore_index=True), aborted
