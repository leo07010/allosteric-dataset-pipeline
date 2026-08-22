"""
filter_allobench.py
===================
Step 1 of the pipeline: turn the AlloBench release CSV into the entry list that the
download and build steps consume.

Two rows are dropped here, for reasons that matter to what the labels mean:

  site_overlap == Yes   AlloBench flags these itself: the annotated allosteric site
                        overlaps the active site. They are real biology, but they make
                        "find the allosteric site given the active site" degenerate --
                        the answer is partly the thing you were handed. Keeping them
                        would inflate any seeded model's score for the wrong reason.

  < MIN_SITE residues   A site of one or two residues is usually a truncated
                        annotation rather than a pocket, and it cannot support a 4A
                        Jaccard evaluation.

Everything else is kept, including rows whose structure later fails to parse -- that
filtering belongs to the build step, which can actually see the coordinates, and it is
counted there rather than silently here.

Output: metadata/allobench_filtered.json  {n, pdbs, entries[]}
"""

import ast
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_IN = os.environ.get("ALLOBENCH_CSV",
                        os.path.join(HERE, "..", "data", "AlloBench.csv"))
OUT = os.environ.get("FILTERED_JSON",
                     os.path.join(HERE, "..", "metadata", "allobench_filtered.json"))
MIN_SITE = 3


def parse_list(raw):
    """AlloBench stores residue lists as a python-literal string."""
    if not raw:
        return []
    try:
        v = ast.literal_eval(raw)
    except (ValueError, SyntaxError):
        return []
    return list(v) if isinstance(v, (list, tuple)) else [v]


def main():
    if not os.path.exists(CSV_IN):
        sys.exit("AlloBench.csv not found at {}\n"
                 "Download it from the AlloBench release (see README) and place it "
                 "there, or set ALLOBENCH_CSV.".format(CSV_IN))

    rows = list(csv.DictReader(open(CSV_IN)))
    entries, dropped = [], {"site_overlap": 0, "site_too_small": 0, "no_pdb": 0}

    for r in rows:
        if (r.get("site_overlap") or "").strip().lower() == "yes":
            dropped["site_overlap"] += 1
            continue
        pdb = (r.get("allosteric_pdb") or "").strip().upper()
        if not pdb:
            dropped["no_pdb"] += 1
            continue
        allo = parse_list(r.get("allosteric_site_residue"))
        if len(allo) < MIN_SITE:
            dropped["site_too_small"] += 1
            continue
        entries.append({
            "pdb": pdb,
            "uniprot": (r.get("pdb_uniprot") or "").strip(),
            "gene": (r.get("target_gene") or "").strip(),
            "organism": (r.get("organism") or "").strip(),
            "modulator": (r.get("modulator_alias") or "").strip(),
            "mod_chain": (r.get("modulator_chain") or "").strip(),
            "allo": allo,
            "active": parse_list(r.get("active_site_residue")),
            "sequence": (r.get("sequence") or "").strip(),
        })

    pdbs = sorted({e["pdb"] for e in entries})
    out = {"n": len(entries), "pdbs": pdbs, "to_download": pdbs, "entries": entries}
    os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)
    json.dump(out, open(OUT, "w"), indent=1)

    print("AlloBench rows read       : {}".format(len(rows)))
    for k, v in dropped.items():
        print("  dropped {:<16}: {}".format(k, v))
    print("entries kept              : {}".format(len(entries)))
    print("unique PDB ids            : {}".format(len(pdbs)))
    print("wrote {}".format(OUT))


if __name__ == "__main__":
    main()
