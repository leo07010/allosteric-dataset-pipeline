# allosteric-dataset-pipeline

Builds a residue-level allosteric-site dataset from the [AlloBench](https://pmc.ncbi.nlm.nih.gov/articles/PMC12059942/) release, for training and evaluating allosteric pocket predictors.

**1439 samples · 327 unique UniProt · 1367 PDB structures · labels at 4 Å heavy-atom**

The repository ships the *pipeline*, not the data — see [Licensing](#licensing). Running four commands reproduces the dataset bit-for-bit.

## Why this exists

Work on allosteric pocket prediction kept running into the same three limits, and this dataset is built to remove them:

| Limit | What it did | Fix here |
|---|---|---|
| **Label/metric mismatch** | Training labels were drawn at 8 Å while every reported score was measured at 4 Å. A model trained that way learns to flag a whole region when the metric rewards a handful of specific residues. | Labels are AlloBench's own **4 Å heavy-atom** residue lists — the same definition the evaluation uses. |
| **Missing active sites** | ASD has no active-site column; recovering it from non-modulator ligands worked for 52 of 123 proteins, so the seeded formulation ("given the active site, find the allosteric site") could not be trained at scale. | AlloBench annotates the active site for **every** entry: 100 % coverage. |
| **Coverage** | 60 unique UniProt. Held-out enrichment tracked similarity to the training set (Spearman +0.53, p<1e-4; least-similar tertile 1.40× vs most-similar 7.04×), which makes coverage the binding constraint rather than architecture. | **327 unique UniProt**, a 5.5× increase. |

## Dataset

| | |
|---|---|
| Samples | 1439 |
| Unique UniProt | 327 |
| Unique PDB | 1367 |
| Multi-site groups | 18 (same PDB+chain+modulator, distinct pockets) |
| Median allosteric residues | 10 |
| Median chain length | 327 |
| Positive rate | 2.97 % |
| Label definition | 4 Å heavy-atom to the allosteric modulator |
| Cross-validation | 5 folds **grouped by UniProt** (288/288/288/288/287) |

Each sample is an `.npz`:

| Array | Shape | Meaning |
|---|---|---|
| `coords` | (N, 3) float32 | Cα coordinates, single chain |
| `resnums` | (N,) int32 | author residue numbers |
| `allo_labels` | (N,) int8 | 1 = within 4 Å of the allosteric modulator |
| `active_site_mask` | (N,) int8 | 1 = annotated active-site residue |
| `meta` | JSON string | PDB, chain, UniProt, gene, organism, modulator, counts |

Files are named `{pdb}_{chain}_{modulator}_{sitehash}.npz`. The trailing hash is a digest of the site's own residue list, which keeps two distinct pockets on the same chain and the same ligand from colliding — `2JC9`/ADN carries two non-overlapping sites, `4GQQ`/0XR carries three. Keying on PDB+chain+modulator alone silently overwrote one site's labels while leaving its manifest row in place, so the hash is load-bearing, not cosmetic.

### Design decisions worth knowing before you train on this

- **Chain-aware residue matching.** AlloBench residue IDs carry their chain (`B-THR-7`), and only that chain is kept. Matching residue numbers chain-agnostically merges repeated domains — in BCR-ABL1 (`1OPL`, two chains, resnums repeating 81–531) that collapsed a genuine 27.8 Å active↔allosteric separation to an apparent 4.7 Å.
- **Overlapping residues are removed from the allosteric label**, not silently kept in both (598 entries affected). A residue claimed by both sites is not evidence for the allosteric one.
- **`site_overlap = Yes` rows are dropped** (316). They are real biology, but they make the seeded task degenerate: part of the answer is the thing the model was handed.
- **Folds are grouped by UniProt**, so no protein appears on both sides of a split. With 1367 structures over 327 proteins, a random split would leak heavily.

### What is dropped, and why

`filter_allobench.py` drops 402 of 2257 CSV rows (316 `site_overlap`, 86 sites under 3 residues). `build_dataset_v2.py` then reads structures and drops 416 more:

| Reason | Count |
|---|---|
| kept | **1439** |
| active-site residues did not map to the structure | 331 |
| allosteric label fell below 3 residues after overlap removal | 51 |
| chain shorter than 40 residues | 18 |
| allosteric residues did not map | 14 |
| ProDy returned no structure | 2 |

The `active_unmapped` count is the largest single loss and is worth attention if you do not need the active-site channel: those 331 entries have usable allosteric labels and are excluded only because the seeded formulation requires a seed. Lowering `MIN_ACT` to 0 recovers them.

## Reproducing

Requires Python 3 with `numpy` and [`prody`](http://prody.csb.pitt.edu/); the download step needs `curl`.

```bash
# 1. obtain AlloBench.csv from the AlloBench release and place it at data/AlloBench.csv
# 2. filter the release into an entry list        (~1 s)
python pipeline/filter_allobench.py
# 3. fetch 1766 structures from RCSB              (~5 min, resumable, ~2 GB)
bash pipeline/download_allobench.sh
# 4. build labelled samples + UniProt-grouped folds  (~10 min)
python pipeline/build_dataset_v2.py
# 5. verify what landed on disk
python pipeline/verify_dataset.py
```

Every path is overridable: `ALLOBENCH_CSV`, `FILTERED_JSON`, `PDB_DIR`, `OUT_DIR`, `MANIFEST_OUT`. Put `OUT_DIR` on a filesystem with a few GB free.

`verify_dataset.py` is an independent check, not a re-run of the builder's own logic. It reloads every file from disk and asserts that each manifest row owns exactly one file, that fold keys partition the manifest, that no UniProt spans two folds, that sampled arrays agree with their manifest row, and that the 18 multi-site groups really do carry distinct labels. The build itself aborts rather than writing a manifest whose rows do not each own a file. Expected output:

```
files on disk        : 1439
manifest rows        : 1439
fold keys (total)    : 1439
multi-site groups    : 18 (all distinct labels: True )
uniprot leak across folds: 0
sampled files checked: 200, mismatches: 0

VERDICT: PASS - dataset is internally consistent
```

## Contents

```
pipeline/filter_allobench.py   AlloBench.csv -> entry list
pipeline/download_allobench.sh entry list -> PDB structures (8-way parallel, resumable)
pipeline/build_dataset_v2.py   structures + labels -> .npz + manifest + folds
pipeline/verify_dataset.py     independent consistency check
pipeline/build_asd_dataset.py  earlier ASD-release path, kept for provenance
metadata/manifest.json         per-sample metadata and fold assignment (no residue lists)
metadata/summary.json          aggregate statistics
```

`metadata/` holds identifiers, counts, and fold assignments so a rebuild can be checked against a known-good result. It contains no residue-level annotations.

## Known limitations

- **Cα only.** Side-chain geometry is discarded, though allosteric pockets are largely lined by side chains. The `.npz` schema has room for heavy-atom coordinates; the builder does not currently write them.
- **Apo/holo is not separated.** Structures are the allosteric-modulator-bound ones AlloBench annotates. Blind prediction from an unbound structure needs apo counterparts, which this pipeline does not pair.
- **Chain selection is by majority vote** over the labelled residues, falling back to all Cα atoms when the named chain is absent. For entries whose site spans an interface this keeps one side only.
- **Sequence redundancy is not clustered.** Folds are grouped by UniProt, which stops the same protein straddling a split but does not stop close homologues sitting in different folds. For a stricter split, cluster the sequences (MMseqs2 at 30 % identity) and group by cluster instead.
- **AlloBench's annotations are inherited as given**, including any errors in the underlying ASD records.

## Licensing

The dataset derives from AlloBench, itself derived from the Allosteric Database (ASD), which is released for **research use only and may not be redistributed to third parties**. This repository therefore contains no ASD- or AlloBench-derived residue annotations and no structures — only the code that rebuilds them, plus identifier-level metadata. Obtain `AlloBench.csv` from its own release under its own terms.

The code in `pipeline/` is MIT licensed.

## Citing

Cite AlloBench and ASD for the data. If the pipeline itself is useful, a link back to this repository is enough.
