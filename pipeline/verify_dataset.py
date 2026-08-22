import json, os, collections, random
import numpy as np
import sys
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.environ.get("OUT_DIR", os.path.join(HERE, "..", "data", "processed"))
man = json.load(open(os.environ.get("MANIFEST_OUT", os.path.join(HERE, "..", "metadata", "manifest.json"))))
rows, folds = man["manifest"], man["folds"]
bykey = {r["key"]: r for r in rows}
fail = []

# 1. every fold key owns a file, and folds partition the manifest exactly
allf = [k for f in folds.values() for k in f]
if len(allf) != len(set(allf)): fail.append("fold keys contain duplicates")
if set(allf) != set(bykey): fail.append("folds != manifest keys")
missing = [k for k in allf if not os.path.exists(os.path.join(D, k + ".npz"))]
if missing: fail.append("{} fold keys have no file (e.g. {})".format(len(missing), missing[:3]))

# 2. no UniProt spans two folds
f_of = {k: fi for fi, ks in folds.items() for k in ks}
u2f = collections.defaultdict(set)
for r in rows: u2f[r["uniprot"] or r["name"]].add(f_of[r["key"]])
leak = [u for u, s in u2f.items() if len(s) > 1]
if leak: fail.append("{} UniProt span >1 fold: {}".format(len(leak), leak[:3]))

# 3. reload files, check contents match the manifest claim
random.seed(0)
bad = 0
for k in random.sample(list(bykey), 200):
    z = np.load(os.path.join(D, k + ".npz")); r = bykey[k]
    if (int(z["allo_labels"].sum()) != r["n_allo"] or
        int(z["active_site_mask"].sum()) != r["n_active"] or
        len(z["coords"]) != r["n_residues"] or
        len(z["resnums"]) != r["n_residues"] or
        int((z["allo_labels"] & z["active_site_mask"]).sum()) != 0):
        bad += 1
if bad: fail.append("{}/200 sampled files disagree with manifest".format(bad))

# 4. the 18 previously-colliding groups must now hold DISTINCT labels
grp = collections.defaultdict(list)
for r in rows: grp[(r["name"], r["chain"], r["modulator"])].append(r["key"])
multi = {g: ks for g, ks in grp.items() if len(ks) > 1}
ident = 0
for g, ks in multi.items():
    sets = [frozenset(np.load(os.path.join(D, k + ".npz"))["resnums"][
        np.load(os.path.join(D, k + ".npz"))["allo_labels"].astype(bool)].tolist()) for k in ks]
    if len(set(sets)) != len(sets): ident += 1
if ident: fail.append("{} multi-site groups still share identical labels".format(ident))

print("files on disk        :", len([f for f in os.listdir(D) if f.endswith('.npz')]))
print("manifest rows        :", len(rows))
print("fold keys (total)    :", len(allf))
print("multi-site groups    :", len(multi), "(all distinct labels:", ident == 0, ")")
print("uniprot leak across folds:", len(leak))
print("sampled files checked: 200, mismatches:", bad)
print()
print("VERDICT:", "PASS - dataset is internally consistent" if not fail else "FAIL")
for f in fail: print("  !!", f)
