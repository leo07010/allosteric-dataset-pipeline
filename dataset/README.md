# dataset/

The built dataset, 1439 samples. See [../DATASET.md](../DATASET.md) for provenance,
design rationale, composition and limitations.

```
samples/{pdb}_{chain}_{modulator}_{sitehash}.npz    coords, resnums, allo_labels,
                                                    active_site_mask, meta
cb/{same key}.npz                                   cb_coords, resnames, resnums
```

Both directories are keyed identically, so a sample and its Cβ sidecar are the same
filename in each. Fold assignments are in [../metadata/manifest.json](../metadata/manifest.json).

```python
import json, numpy as np

folds = json.load(open("metadata/manifest.json"))["folds"]
held_out = folds["0"]                       # grouped by UniProt; no protein spans folds

z  = np.load(f"dataset/samples/{held_out[0]}.npz")
zc = np.load(f"dataset/cb/{held_out[0]}.npz", allow_pickle=True)

coords = z["coords"]                         # (N,3) Calpha
cb     = zc["cb_coords"]                     # (N,3) Cbeta, Calpha for glycine
y      = z["allo_labels"].astype(bool)       # within 4 A heavy-atom of the modulator
anchor = z["active_site_mask"].astype(bool)  # annotated active site
meta   = json.loads(str(z["meta"]))
```

Coordinates originate from RCSB and are public domain. The residue-level annotations
derive from ASD via AlloBench and are **research use only**; cite AlloBench and ASD,
and do not redistribute them onward.
