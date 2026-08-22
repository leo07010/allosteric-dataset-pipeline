"""
build_dataset_v2.py
===================
Turn AlloBench annotations plus downloaded structures into a training set that
matches how we actually evaluate.

Three defects in the current 123-protein set that this fixes:

  labels    the old .pt files mark a residue positive if it lies within 8 A (some
            versions 10 A) of the modulator, yet every number we report is scored at
            4 A heavy-atom. Training on a looser target than you evaluate teaches the
            model to flag a whole region when the metric wants five specific
            residues -- a plausible reason Hit@5 never moved. AlloBench ships 4 A
            residue lists, so we adopt them directly.

  active    ASD has no active-site column; we previously recovered it from
            non-modulator ligands and it worked for 52 of 123 proteins. AlloBench
            annotates it for every entry, which makes the seeded formulation
            ("given the active site, find the allosteric site") trainable at scale
            for the first time.

  coverage  60 unique UniProts today. Held-out enrichment tracks similarity to the
            training set (spearman +0.53, p<1e-4; least-similar tertile 1.40x vs
            most-similar 7.04x), so coverage is the binding constraint. This reaches
            ~375.

Chain handling: AlloBench residue IDs carry their chain ("B-THR-7"), so the chain is
read from the label instead of assumed. Only that chain is kept -- the BCR-ABL1
collision earlier in this project came from matching residue numbers chain-agnostically
across two copies of the same domain.

A structure is dropped, and counted, if the labelled residues do not map, if fewer
than 3 allosteric or 1 active residue survive, or if the chain is too short. Residues
claimed by both sites are removed from the allosteric label rather than silently kept.
"""

import collections
import hashlib
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PDB_DIR = os.environ.get("PDB_DIR", os.path.join(HERE, "..", "data", "allobench_pdbs"))
OUT_DIR = os.environ.get("OUT_DIR", os.path.join(HERE, "..", "data", "processed"))
MIN_ALLO, MIN_ACT, MIN_RES = 3, 1, 40


def parse_resid(tag):
    """'B-THR-7' -> ('B', 7); a bare number -> (None, n)."""
    p = str(tag).split("-")
    if len(p) >= 3:
        try:
            return p[0].strip(), int(p[-1])
        except ValueError:
            return None
    try:
        return None, int(p[-1])
    except ValueError:
        return None


def main():
    from prody import parsePDB, confProDy
    confProDy(verbosity="none")

    data = json.load(open(os.environ.get(
        "FILTERED_JSON", os.path.join(HERE, "..", "metadata", "allobench_filtered.json"))))
    entries = data["entries"]
    os.makedirs(OUT_DIR, exist_ok=True)
    stats = collections.Counter()
    manifest = []
    seen_keys = set()

    for i, e in enumerate(entries):
        pdb_path = os.path.join(PDB_DIR, e["pdb"] + ".pdb")
        if not os.path.exists(pdb_path):
            stats["pdb_missing"] += 1
            continue
        allo = [x for x in (parse_resid(t) for t in e["allo"]) if x]
        act = [x for x in (parse_resid(t) for t in e["active"]) if x]
        if not allo:
            stats["unparsable_allo"] += 1
            continue
        chains = collections.Counter(c for c, _ in allo if c)
        chain = chains.most_common(1)[0][0] if chains else "A"

        try:
            st = parsePDB(pdb_path)
        except Exception:
            stats["parse_failed"] += 1
            continue
        if st is None:                      # ProDy returns None on some malformed files
            stats["parse_returned_none"] += 1
            continue
        sel = st.select("protein and name CA and chain {}".format(chain)) or \
            st.select("protein and name CA")
        if sel is None:
            stats["no_ca"] += 1
            continue
        resnums = [int(r) for r in sel.getResnums()]
        if len(resnums) < MIN_RES:
            stats["too_short"] += 1
            continue
        idx = {r: k for k, r in enumerate(resnums)}

        allo_r = sorted({n for c, n in allo if c in (None, chain) and n in idx})
        act_r = sorted({n for c, n in act if c in (None, chain) and n in idx})
        if len(allo_r) < MIN_ALLO:
            stats["allo_unmapped"] += 1
            continue
        if len(act_r) < MIN_ACT:
            stats["active_unmapped"] += 1
            continue
        both = set(allo_r) & set(act_r)
        if both:
            allo_r = [r for r in allo_r if r not in both]
            stats["overlap_trimmed"] += 1
            if len(allo_r) < MIN_ALLO:
                stats["allo_gone_after_trim"] += 1
                continue

        y = np.zeros(len(resnums), np.int8)
        for r in allo_r:
            y[idx[r]] = 1
        seed = np.zeros(len(resnums), np.int8)
        for r in act_r:
            seed[idx[r]] = 1

        meta = {"name": e["pdb"], "chain": chain, "uniprot": e["uniprot"],
                "gene": e["gene"], "organism": e["organism"],
                "modulator": e["modulator"], "n_residues": len(resnums),
                "n_allo": int(y.sum()), "n_active": int(seed.sum()),
                "label_definition": "AlloBench 4A heavy-atom to allosteric modulator"}
        # One (pdb, chain, modulator) can carry SEVERAL distinct allosteric sites --
        # 2JC9/ADN has two non-overlapping pockets, 4GQQ/0XR has three -- so keying on
        # pdb+chain+modulator still silently overwrote one site's labels while leaving
        # its manifest row behind. The key therefore ends in a hash of the site itself.
        # Content-addressed on purpose: it is stable across rebuilds (independent of
        # iteration order) and two rows describing the genuinely same site collapse to
        # one entry instead of becoming a phantom duplicate.
        mod_safe = "".join(ch if ch.isalnum() else "_" for ch in (e["modulator"] or "NA"))
        site_hash = hashlib.sha1(
            ",".join(str(r) for r in allo_r).encode()).hexdigest()[:6]
        key = "{}_{}_{}_{}".format(e["pdb"], chain, mod_safe, site_hash)
        if key in seen_keys:
            stats["duplicate_site_dropped"] += 1
            continue
        seen_keys.add(key)
        meta["key"] = key
        np.savez_compressed(
            os.path.join(OUT_DIR, key + ".npz"),
            coords=sel.getCoords().astype(np.float32),
            resnums=np.array(resnums, np.int32),
            allo_labels=y, active_site_mask=seed,
            meta=json.dumps(meta))
        manifest.append(meta)
        stats["ok"] += 1
        if (i + 1) % 300 == 0:
            print("  {}/{} scanned, {} kept".format(i + 1, len(entries), stats["ok"]),
                  flush=True)

    print("\n=== build summary ===")
    for k, v in stats.most_common():
        print("  {:<24}{:>6}".format(k, v))
    if not manifest:
        print("nothing built"); return

    unis = {m["uniprot"] for m in manifest if m["uniprot"]}
    print("\n  samples           : {}".format(len(manifest)))
    print("  unique UniProt    : {}".format(len(unis)))
    print("  median allo resid : {}".format(sorted(m["n_allo"] for m in manifest)[len(manifest)//2]))
    print("  median chain len  : {}".format(sorted(m["n_residues"] for m in manifest)[len(manifest)//2]))
    print("  positive rate     : {:.2%}".format(
        sum(m["n_allo"] for m in manifest) / sum(m["n_residues"] for m in manifest)))

    # folds grouped by UniProt so no protein straddles a split
    by_uni = collections.defaultdict(list)
    for m in manifest:
        by_uni[m["uniprot"] or m["name"]].append(m["key"])
    folds = [[] for _ in range(5)]
    for u in sorted(by_uni, key=lambda x: -len(by_uni[x])):
        folds[min(range(5), key=lambda k: len(folds[k]))].extend(by_uni[u])
    print("\n  5-fold split grouped by UniProt:")
    for k, f in enumerate(folds):
        print("    fold {}: {:>4}".format(k, len(f)))

    manifest_out = os.environ.get(
        "MANIFEST_OUT", os.path.join(HERE, "..", "metadata", "manifest.json"))
    on_disk = len([f for f in os.listdir(OUT_DIR) if f.endswith(".npz")])
    keys = {m["key"] for m in manifest}
    print("\n  files on disk     : {}".format(on_disk))
    print("  unique keys       : {}".format(len(keys)))
    if on_disk != len(manifest) or len(keys) != len(manifest):
        raise SystemExit(
            "FATAL: manifest rows {} vs unique keys {} vs files {} -- refusing to "
            "write a manifest whose rows do not each own a file".format(
                len(manifest), len(keys), on_disk))
    print("  OK: one file per manifest row, no collisions")

    json.dump({"manifest": manifest, "folds": {str(k): f for k, f in enumerate(folds)},
               "n_uniprot": len(unis), "stats": dict(stats),
               "output_dir": OUT_DIR},
              open(manifest_out, "w"),
              indent=2)
    print("\nwrote {}/*.npz".format(OUT_DIR))
    print("wrote {}".format(manifest_out))


if __name__ == "__main__":
    main()
