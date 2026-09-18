#!/usr/bin/env python3
"""Execute a notebook's code cells headlessly, in order, in one namespace.

A minimal stand-in for `jupyter nbconvert --execute` for environments without
jupyter installed (the cluster, CI). Semantics match a fresh kernel: cells run
top to bottom sharing one global namespace, and the first failure stops the run
and reports the cell index and source.

    python tools/run_notebook.py biokgsuite_dimensions/01_coverage/notebook.ipynb

Outputs are not written back into the .ipynb; this is for verifying that a
notebook runs and for regenerating figures and checkpoints, not for producing
an executed copy.
"""
import argparse
import ast
import json
import os
import sys
import traceback
from pathlib import Path

os.environ.setdefault('MPLBACKEND', 'Agg')


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('notebook', type=Path)
    ap.add_argument('--stop-before', type=int, default=None,
                    help='stop before this code-cell index (for partial runs)')
    args = ap.parse_args()

    nb_path = args.notebook.resolve()
    nb = json.loads(nb_path.read_text())
    cells = [(i, ''.join(c['source']))
             for i, c in enumerate(nb['cells']) if c['cell_type'] == 'code']

    # Match Jupyter: cwd is the notebook's own directory.
    os.chdir(nb_path.parent)
    # Shim the handful of IPython builtins notebooks rely on, so cells written
    # for a kernel run unchanged here.
    def _display(*objs, **_kw):
        for o in objs:
            try:
                print(o.to_string() if hasattr(o, 'to_string') else o)
            except Exception:
                print(repr(o))

    ns = {
        '__name__': '__main__',
        '__file__': str(nb_path),
        'display': _display,
        'get_ipython': lambda: None,
    }

    for n, (idx, source) in enumerate(cells):
        if args.stop_before is not None and n >= args.stop_before:
            print(f'\n[run_notebook] stopping before code cell {n} (nb cell {idx})')
            break
        print(f'\n[run_notebook] === code cell {n} (nb cell {idx}) ===', flush=True)
        try:
            # Jupyter echoes a cell's trailing expression. Reproduce that, so a
            # headless run shows the same tables a kernel would and a cell whose
            # only output is a bare DataFrame is not silently invisible here.
            tree = ast.parse(source)
            if tree.body and isinstance(tree.body[-1], ast.Expr):
                last = ast.Expression(tree.body.pop().value)
                exec(compile(ast.Module(body=tree.body, type_ignores=[]),
                             f'<cell {idx}>', 'exec'), ns)
                value = eval(compile(last, f'<cell {idx}>', 'eval'), ns)
                if value is not None:
                    ns['display'](value)
            else:
                exec(compile(source, f'<cell {idx}>', 'exec'), ns)
        except Exception:
            print(f'\n[run_notebook] FAILED at code cell {n} (nb cell {idx})\n')
            print(source)
            print()
            traceback.print_exc()
            return 1

    print(f'\n[run_notebook] OK — {len(cells)} code cells')
    return 0


if __name__ == '__main__':
    sys.exit(main())
