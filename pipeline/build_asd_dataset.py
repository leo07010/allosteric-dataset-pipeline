"""
build_asd_dataset.py
====================
Stage 1 of building our own training set: filter the full ASD release and report
what actually survives. No downloads, no structure parsing -- this decides whether
the download is worth doing before spending hours on it.

Why we need our own set at all: the current model trains on 123 ASD proteins, and we
measured that this is the binding constraint. Held-out enrichment tracks how similar
a protein is to the rest of the training set (spearman +0.53, p<1e-4); the
least-similar tertile enriches 1.40x against 7.04x for the most-similar. CardiacMyosin
sits in the least-similar tertile and its observed 1.49x is exactly what that predicts.
More coverage, not a bigger model, is what the evidence points at.

What gets filtered, and why each rule exists (all of these are mistakes we already
made or nearly made):
  * modulator_class == Ion        a lone Ca/Mg/Cl is usually structural or a
                                  crystallisation additive, not an allosteric effector
  * EXCLUDE_RESNAMES (164 names)  waters, buffers, cryoprotectants, modified residues,
                                  glycans -- reused from residual_rank.py rather than
                                  rewritten
  * DTT / DTU                     reducing agents. ASD lists them as modulators for six
                                  entries; labelling residues "within 8A of a reducing
                                  agent" is meaningless and we already dropped these
  * missing pdb / modulator code  cannot build a label without both

Nucleotides (AMP, ATP, ADP, ...) are deliberately NOT filtered. AMP is the textbook
allosteric activator of glycogen phosphorylase; excluding it by chemistry would throw
away real allostery. They are counted and reported separately so the decision stays
visible rather than buried in a filter.
"""

import collections
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from residual_rank import EXCLUDE_RESNAMES          # noqa: E402  (164 names)

ASD = os.path.join(HERE, "data", "ASD_Release_201909_AS.txt")
REDUCING = {"DTT", "DTU", "BME", "TCE"}
NUCLEOTIDES = {"AMP", "ADP", "ATP", "CMP", "CDP", "CTP", "GMP", "GDP", "GTP",
               "UMP", "UDP", "UTP", "TMP", "TDP", "TTP", "IMP", "NAD", "NAP",
               "NDP", "FAD", "FMN", "COA", "SAM", "SAH", "ANP", "AGS", "GNP"}


def main():
    rows = list(csv.DictReader(open(ASD), delimiter="\t"))
    print("ASD release: {} records".format(len(rows)))

    keep, dropped = [], collections.Counter()
    for r in rows:
        pdb = (r.get("allosteric_pdb") or "").strip().upper()
        mod = (r.get("modulator_alias") or "").strip().upper()
        cls = (r.get("modulator_class") or "").strip()
        if not pdb or len(pdb) != 4:
            dropped["no_pdb"] += 1; continue
        if not mod:
            dropped["no_modulator_code"] += 1; continue
        if cls == "Ion":
            dropped["modulator_class_Ion"] += 1; continue
        if mod in REDUCING:
            dropped["reducing_agent"] += 1; continue
        if mod in EXCLUDE_RESNAMES:
            dropped["buffer_cryo_water_glycan"] += 1; continue
        keep.append({"pdb": pdb, "modulator": mod, "class": cls,
                     "chain": (r.get("modulator_chain") or "").strip(),
                     "gene": (r.get("target_gene") or "").strip(),
                     "organism": (r.get("organism") or "").strip(),
                     "uniprot": (r.get("pdb_uniprot") or "").strip(),
                     "function": (r.get("function") or "").strip()})

    print("\n-- filtering --")
    for k, v in dropped.most_common():
        print("  dropped {:<28} {:>5}".format(k, v))
    print("  KEPT {:>31}".format(len(keep)))

    pdbs = {k["pdb"] for k in keep}
    mods = collections.Counter(k["modulator"] for k in keep)
    unis = {k["uniprot"] for k in keep if k["uniprot"]}
    orgs = collections.Counter(k["organism"] for k in keep)
    nuc = sum(v for m, v in mods.items() if m in NUCLEOTIDES)

    print("\n-- what survives --")
    print("  unique allosteric PDB entries : {}".format(len(pdbs)))
    print("  unique modulators             : {}".format(len(mods)))
    print("  unique UniProt (proxy for fold diversity) : {}".format(len(unis)))
    print("  organisms                     : {}".format(len(orgs)))
    print("  records whose modulator is a nucleotide/cofactor : {} ({:.0%} of kept)"
          .format(nuc, nuc / max(1, len(keep))))
    print("  top modulators: {}".format(
        ", ".join("{}({})".format(m, c) for m, c in mods.most_common(10))))

    have_pt = {os.path.basename(f)[:-3].upper()
               for f in os.listdir(os.path.join(HERE, "data", "asd_processed"))
               if f.endswith(".pt")}
    have_pdb = {os.path.basename(f)[:-4].upper()
                for f in os.listdir(os.path.join(HERE, "data", "asd_pdbs"))
                if f.endswith(".pdb")}
    print("\n-- against what we already have --")
    print("  already preprocessed (.pt)    : {}".format(len(have_pt)))
    print("  already downloaded  (.pdb)    : {}".format(len(have_pdb)))
    print("  kept entries already downloaded : {}".format(len(pdbs & have_pdb)))
    print("  NEW pdbs that would need downloading : {}".format(len(pdbs - have_pdb)))
    print("  -> potential training set = {} proteins ({:.1f}x the current 123)"
          .format(len(pdbs), len(pdbs) / 123.0))

    out = {"n_records_raw": len(rows), "n_records_kept": len(keep),
           "dropped": dict(dropped), "n_unique_pdb": len(pdbs),
           "n_unique_uniprot": len(unis), "n_organisms": len(orgs),
           "n_nucleotide_records": nuc,
           "already_downloaded": sorted(pdbs & have_pdb),
           "to_download": sorted(pdbs - have_pdb),
           "entries": keep}
    dst = os.path.join(HERE, "outputs", "dataset_v2", "asd_filtered.json")
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    json.dump(out, open(dst, "w"), indent=2)
    print("\nwrote {}".format(dst))


if __name__ == "__main__":
    main()
