# Notebook conventions

Every dimension notebook follows the same skeleton, so a reader who has read one
can navigate any of them. Applied to 00-07.

`00_summary` follows the same skeleton with three departures, because it scores
nothing and loads no graph: its title cell lists display items where a dimension
lists sub-metrics, it has no `## Dimension score` section, and its loading cell
reads the seven checkpoints rather than the KG files. It carries the panel guard
like the rest, plus a second gate on missing dimensions.

## Cell order

| # | Type | Contents |
|---|---|---|
| 0 | markdown | Title cell (below) |
| 1 | markdown | `## Set-up`, plus one paragraph on what loads and under what conditions |
| 2 | code | Bootstrap: autoreload, repo-root discovery, imports, `setup_style()` |
| 3 | code | Config, parameters, constants, the `Outputs` writer |
| 4 | markdown | `### Loading` |
| 5 | code | Free any prior KG state, then load |
| … | | Numbered sections, markdown heading first, then its code |
| n-2 | markdown | `## Dimension score` with the formula |
| n-1 | code | Score computation and checkpoint |
| n | code | Free KG state, print the per-KG scores |

## Title cell

Fixed skeleton, in this order:

1. `# 0N: Dimension name`
2. The dimension's question, in italics
3. A sub-metric table with the columns `Sub-metric | Definition | Quantification | Scoring`
4. One line stating how the sub-metrics combine
5. **Reported, not scored.** followed by what is excluded and why, in one paragraph
6. **Inputs:** and **Outputs:** as `·`-separated paths
7. A link to `README.md`

Methodological justifications go in the folder `README.md` under `## Design
notes`, not in the notebook. Result numbers go in the manuscript, not in either.

## Headings

Sentence case throughout: `## 3. Uncertainty quantification`, not `## 3.
Uncertainty Quantification`. Top-level sections are numbered and match the
sub-metrics in the title cell's table. Use "and", not "&".

## Comments

A comment explains **why**, not what. `# loop over diseases` is noise; `# All
diseases including undrugged ones: repurposing exists to serve them` is a
decision a reviewer can audit. Every non-obvious constant is justified where it
is defined. No working notes ("BOTH are replicated now"), no commented-out code.

## Imports

One group per source, in this order, alphabetical within each group:

```
stdlib
<blank>
third-party
<blank>
biokgsuite
```

All at the top of the bootstrap cell. Nothing imported mid-notebook.

## Parameters

Read with `require(params, key[, section])`, never `params.get(key, default)`.
A silent default lets a run use a value `config/params.yaml` never mentions.

## Outputs

Everything a notebook writes goes through `Outputs`, so figures, tables and
checkpoints land under `results/{figures,tables}/<notebook>/` and
`results/checkpoints/<notebook>.pkl`. Do not build result paths by hand.
Figures are drawn by functions under `src/biokgsuite/viz/`, never inline in a
notebook; how they look is fixed by `FIGURE_STYLE.md`.

## Memory

Each notebook frees the previous one's KG state before loading, and its own
before exiting, with one idiom:

```python
for _v in ('kgs', 'graphs'):
    globals().pop(_v, None)
gc.collect()
```

## Folder README

Same sections in every folder: title, **Computes.**, **Inputs.** /
**Outputs.** / **Display item.**, `## What is scored` (table plus formula),
`## Design notes`, `## Run`, and a closing pointer to the `src/` modules holding
the metric definitions and the rendering.
