#!/usr/bin/env python3
"""
Map the FDA indication-expansion census onto DrugBank IDs.

Kept separate from build_fda_census.py on purpose: the census is a statement
about the regulatory record and should not change when the crosswalk changes.
This step is the lossy one, so it is where the losses get counted and named.

Salt forms are handled by falling back to the parent ingredient name
(SILDENAFIL CITRATE -> SILDENAFIL), because DrugBank indexes the parent. Every
fallback is recorded in match_level so the cost of the heuristic is visible
rather than buried. Combination products (ingredients separated by ';') map to
all their components; a combination counts as mapped only if every component
resolves, and partial matches are labelled as such rather than silently kept.

USAGE
  python3 tools/map_census_drugbank.py \
    --census data/gold_standards/legacy/fda_indication_expansion_census.csv \
    --xref   data/gold_standards/bridges/drug_name_to_drugbank.csv
"""
import argparse, json, re
from pathlib import Path
import pandas as pd

# Salt and ester suffixes stripped only as a labelled fallback, never silently.
SALTS = [
    "hydrochloride", "hydrobromide", "sodium", "potassium", "calcium",
    "magnesium", "maleate", "mesylate", "besylate", "tosylate", "citrate",
    "tartrate", "succinate", "fumarate", "acetate", "sulfate", "phosphate",
    "nitrate", "bitartrate", "dihydrochloride", "hemifumarate", "malate",
    "lactate", "gluconate", "carbonate", "chloride", "bromide", "iodide",
    "pamoate", "palmitate", "propionate", "valerate", "furoate", "dipropionate",
    "monohydrate", "dihydrate", "anhydrous", "trihydrate", "disodium",
]


def norm(s: str) -> str:
    """Identical to tools/build_drug_name_xref.py so the keys line up."""
    return re.sub(r"\s+", " ", re.sub(r"[^a-z0-9 ]+", " ", str(s).lower())).strip()


BIOLOGIC_SUFFIX = re.compile(r"-[a-z]{4}$", re.I)


def strip_biologic_suffix(raw: str) -> str:
    """FDA appends a meaningless four-letter suffix to biologic nonproprietary
    names (trastuzumab deruxtecan-nxki). DrugBank indexes the name without it.
    Must run on the raw string: normalisation turns the hyphen into a space and
    the suffix stops being identifiable."""
    return BIOLOGIC_SUFFIX.sub("", str(raw).strip())


def strip_salt(n: str) -> str:
    toks = n.split()
    while toks and toks[-1] in SALTS:
        toks = toks[:-1]
    return " ".join(toks)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--census", required=True)
    ap.add_argument("--xref", required=True)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    cen = pd.read_csv(a.census, dtype=str)
    xr = pd.read_csv(a.xref, dtype=str)

    unambiguous = xr[xr["ambiguous"].fillna("0").isin(["0", "False", "false"])]
    lut = dict(zip(unambiguous["name_normalised"], unambiguous["drugbank_id"]))
    ambiguous_names = set(xr[~xr.index.isin(unambiguous.index)]["name_normalised"])

    def resolve(component: str):
        """Cascade of increasingly lossy attempts. Each returns its own label so
        the mapping's reliance on heuristics stays measurable."""
        raw = str(component).strip()
        if not raw:
            return None, "empty"

        n = norm(raw)
        if n in lut:
            return lut[n], "exact"
        if n in ambiguous_names:
            return None, "ambiguous"

        s = strip_salt(n)
        if s and s != n:
            if s in lut:
                return lut[s], "salt_stripped"
            if s in ambiguous_names:
                return None, "ambiguous_after_salt_strip"

        b = norm(strip_biologic_suffix(raw))
        if b and b != n:
            if b in lut:
                return lut[b], "biologic_suffix_stripped"
            if b in ambiguous_names:
                return None, "ambiguous_after_suffix_strip"
            # "fam-trastuzumab deruxtecan" is indexed as "trastuzumab deruxtecan"
            toks = b.split()
            if len(toks) > 1 and " ".join(toks[1:]) in lut:
                return lut[" ".join(toks[1:])], "suffix_and_prefix_stripped"

        return None, "unmatched"

    ids, levels, unmapped = [], [], []
    for raw in cen["ActiveIngredients"].fillna(""):
        # A single application can list several products; components within a
        # product are ';' separated, products are '|' separated.
        # Products are '|' separated, components within a product ';' separated,
        # and co-formulated biologics use the word 'and' (DARATUMUMAB AND
        # HYALURONIDASE-FIHJ), which is a combination, not one ingredient.
        comps = [c.strip()
                 for part in str(raw).split("|")
                 for piece in part.split(";")
                 for c in re.split(r"\s+\band\b\s+", piece, flags=re.I)
                 if c.strip()]
        comps = list(dict.fromkeys(comps))
        got, lev = [], []
        for c in comps:
            db, how = resolve(c)
            lev.append(how)
            if db:
                got.append(db)
            else:
                unmapped.append(c)
        ids.append("|".join(dict.fromkeys(got)))
        if not comps:
            levels.append("no_ingredient")
        elif len(got) == len(comps):
            if set(lev) == {"exact"}:
                levels.append("exact")
            else:
                heur = sorted({x for x in lev if x != "exact"})
                levels.append("heuristic:" + "+".join(heur))
        elif got:
            levels.append("incomplete_combination")
        else:
            levels.append("unmatched")

    cen["drugbank_ids"] = ids
    cen["drugbank_match_level"] = levels

    out = Path(a.out or a.census.replace(".csv", "_drugbank.csv"))
    cen.to_csv(out, index=False)

    counts = pd.Series(levels).value_counts().to_dict()
    top_unmapped = pd.Series(unmapped).value_counts().head(30).to_dict()
    summary = {
        "n_events": int(len(cen)),
        "n_events_fully_mapped": int(sum(
            1 for l in levels
            if l not in ("unmatched", "incomplete_combination", "no_ingredient"))),
        "n_events_exact": int(sum(1 for l in levels if l == "exact")),
        "match_levels": {k: int(v) for k, v in counts.items()},
        "n_distinct_unmapped_ingredients": int(len(set(unmapped))),
        "most_common_unmapped": {k: int(v) for k, v in top_unmapped.items()},
    }
    Path(str(out) + ".mapping.json").write_text(json.dumps(summary, indent=2))

    print(f"wrote {out}")
    for k, v in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"  {v:>5}  {k}")
    print(f"\nfully mapped: {summary['n_events_fully_mapped']}/{len(cen)} events "
          f"({100*summary['n_events_fully_mapped']/len(cen):.1f}%)")
    print(f"distinct unmapped ingredients: {summary['n_distinct_unmapped_ingredients']}")
    print("most common unmapped:")
    for k, v in list(top_unmapped.items())[:12]:
        print(f"    {v:>3}  {k}")


if __name__ == "__main__":
    main()
