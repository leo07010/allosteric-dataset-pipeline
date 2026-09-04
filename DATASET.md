# The dataset: what it is, where it came from, and why it is built this way

**1439 samples · 327 unique UniProt · 1367 PDB structures · labels at 4 Å heavy-atom · 5 folds grouped by UniProt**

This document describes the dataset itself. [README.md](README.md) covers running the pipeline.

---

## 1. Provenance

Three layers, each adding something and each inheriting the one below it:

| Layer | What it contributes | Terms |
|---|---|---|
| **RCSB PDB** | Atomic coordinates for 1766 structures, fetched by accession | Public domain |
| **ASD** (Allosteric Database) | The underlying curation: which proteins are allosteric, which ligand is the modulator, which residues line the site | **Research use only, not redistributable** |
| **AlloBench** ([PMC12059942](https://pmc.ncbi.nlm.nih.gov/articles/PMC12059942/)) | ASD re-annotated into a benchmark: residue lists computed at a stated 4 Å cutoff, active sites annotated for every entry, `site_overlap` flagged | Inherits ASD's terms |
| **This pipeline** | Structure parsing, chain resolution, label alignment, overlap handling, content-addressed site keys, UniProt-grouped folds | MIT |

The coordinates are public; the *annotations* are not. That distinction decides what this repository can contain — see [Redistribution](#7-redistribution).

## 2. What one sample is

One sample is **one allosteric site on one chain of one structure** — not one protein and not one PDB entry. A chain carrying two distinct pockets for the same ligand yields two samples.

| Array | Shape | Meaning |
|---|---|---|
| `coords` | (N, 3) float32 | Cα coordinates, single chain |
| `cb_coords` | (N, 3) float32 | Cβ, falling back to Cα for glycine *(sidecar)* |
| `resnames` | (N,) str | three-letter codes *(sidecar)* |
| `resnums` | (N,) int32 | author residue numbers |
| `allo_labels` | (N,) int8 | 1 = within 4 Å (heavy atom) of the allosteric modulator |
| `active_site_mask` | (N,) int8 | 1 = annotated active-site residue |
| `meta` | JSON | PDB, chain, UniProt, gene, organism, modulator, counts, key |

File name: `{pdb}_{chain}_{modulator}_{sitehash}.npz`.

## 3. Design rationale

Each decision below is a response to something that went wrong in prior work on this problem, not a preference.

### 3.1 Labels are drawn at the distance the evaluation measures

Earlier training labels marked a residue positive within 8 Å of the modulator while every reported score was computed at 4 Å. A model trained that way is rewarded for flagging a whole region when the metric wants a handful of specific residues, and the mismatch is invisible in the training curve. AlloBench publishes 4 Å heavy-atom residue lists, so they are adopted verbatim rather than recomputed.

### 3.2 Every sample has an active site

ASD has no active-site column. Recovering it from non-modulator ligands worked for 52 of 123 proteins, which left the seeded formulation — *given the active site, find the allosteric site* — untrainable at scale. AlloBench annotates the active site for every entry, so `active_site_mask` is populated for 100 % of samples.

This is what makes the dataset usable by methods that need an anchor: perturbation-response scanning, distance-conditional normalisation, seeded propagation, and any evaluation that restricts scoring to a distal pool.

### 3.3 Coverage was the binding constraint, so coverage is what grew

Held-out enrichment on the earlier 60-UniProt set tracked similarity to the training set (Spearman +0.53, p < 1e-4; least-similar tertile 1.40× versus most-similar 7.04×). Adding architecture to a model that has seen 60 proteins does not fix that. This dataset covers **327 UniProt accessions**, a 5.5× increase, which is at the scale the field's standard training sets occupy.

The effect is measurable. Under an identical training recipe, moving from 105 ASD proteins to ~1150 samples per fold raised held-out distance-stratified AUC from 0.5850 to 0.6563 (paired, p = 8.1e-11) and top-5 hit rate from 0.264 to 0.357 on the same 704 targets.

### 3.4 Sites are addressed by content, so distinct pockets cannot collide

One chain can host several distinct allosteric sites for the same ligand: `2JC9`/ADN has two non-overlapping pockets, `4GQQ`/0XR has three. Keying a sample on PDB + chain + modulator silently overwrote one site's labels with another's while leaving both rows in the manifest, so a loader trusting the manifest read the wrong labels for the shadowed rows.

The key therefore ends in a digest of the site's own residue list. It is content-addressed on purpose: stable across rebuilds regardless of iteration order, and two rows describing genuinely the same site collapse into one entry rather than becoming a phantom duplicate. 18 groups in the dataset are multi-site; all carry distinct labels.

### 3.5 Chains are resolved from the label, never assumed

AlloBench residue IDs carry their chain (`B-THR-7`), and only that chain is kept. Matching residue numbers chain-agnostically merges repeated domains: in BCR-ABL1 (`1OPL`, two chains with residue numbers repeating 81–531) that collapsed a genuine 27.8 Å active↔allosteric separation to an apparent 4.7 Å, which would have made a long-range case look like a short-range one.

### 3.6 Residues claimed by both sites are removed from the allosteric label

598 entries have residues appearing in both the allosteric and active-site lists. A residue the annotation assigns to both is not evidence for the allosteric one, so it is dropped from `allo_labels` rather than silently kept in both. Entries falling below 3 allosteric residues after this removal are excluded (51 of them).

### 3.7 `site_overlap = Yes` rows are dropped

AlloBench flags 316 entries whose allosteric site overlaps the active site. They are real biology, but they make the seeded task degenerate — part of the answer is the thing the model was handed — and they would inflate any anchored method's score for the wrong reason.

### 3.8 Folds are grouped by UniProt

1367 structures over 327 proteins means a random split puts different structures of the same protein on both sides. Folds are built by assigning whole UniProt groups, largest first, to the currently smallest fold: **288 / 288 / 288 / 288 / 287**, with zero UniProt spanning two folds.

This is not hypothetical. A model trained on the earlier 105-protein ASD set and evaluated on this data shares 11.8 % of UniProt accessions with its own training set; on that shared portion it scores 0.806 and on the unseen portion 0.585. Grouped folds make that failure impossible by construction rather than something to be caught afterwards.

## 4. Composition

| | |
|---|---|
| Samples | 1439 |
| Unique UniProt | 327 |
| Unique PDB | 1367 |
| Multi-site groups | 18 |
| Median allosteric residues per site | 10 |
| Median chain length | 327 |
| Positive rate | 2.97 % |
| Distinct modulators | 824 |
| Distinct genes | 263 |
| Organisms | 105 |

Structures per protein are very uneven: the median UniProt contributes 2 samples and
half contribute exactly one, while transthyretin alone contributes 114. See
[6.1](#61-the-327-proteins-are-not-evenly-represented).

## 5. What is excluded, and why

`filter_allobench.py` drops 402 of 2257 CSV rows; `build_dataset_v2.py` drops 416 more after reading coordinates.

| Stage | Reason | Count |
|---|---|---|
| filter | `site_overlap = Yes` | 316 |
| filter | site under 3 residues | 86 |
| build | **kept** | **1439** |
| build | active-site residues did not map to the structure | 331 |
| build | allosteric label fell under 3 residues after overlap removal | 51 |
| build | chain shorter than 40 residues | 18 |
| build | allosteric residues did not map | 14 |
| build | ProDy returned no structure | 2 |

`active_unmapped` is the largest single loss. Those 331 entries have usable allosteric labels and are excluded only because the seeded formulation requires an anchor; setting `MIN_ACT = 0` recovers them for unanchored use.

## 6. Known limitations

### 6.1 The 327 proteins are not evenly represented

Sample counts are dominated by a handful of well-crystallised targets:

| | |
|---|---|
| Top 10 proteins | 33 % of samples |
| Top 50 proteins | 66 % of samples |
| Proteins needed to cover half the samples | 25 |
| UniProt entries with exactly one sample | 163 (50 %) |

The largest groups are transthyretin (114), muscle glycogen phosphorylase (73), HIV-1
gag-pol (61), ER-alpha (46) and CRP (35).

This has a direct consequence for the folds. Grouping by UniProt keeps a protein out
of two folds, but a group that large lands whole in one of them, so folds are balanced
in sample count and not in composition:

| fold | samples | UniProt | largest single protein |
|---|---|---|---|
| 0 | 288 | 63 | TTR, 114 (**40 %**) |
| 1 | 288 | 65 | PYGM, 73 (25 %) |
| 2 | 288 | 66 | gag-pol, 61 (21 %) |
| 3 | 288 | 67 | ESR1, 46 (16 %) |
| 4 | 287 | 66 | CRP, 35 (12 %) |

Read per-fold numbers with that in mind: fold 0's held-out score is substantially a
statement about transthyretin. Report the pooled result across folds, or weight by
UniProt group rather than by sample, if per-protein generalisation is the claim.

### 6.2 Everything else

- **Cα and Cβ only.** No side-chain atoms, though allosteric pockets are largely lined by side chains. The schema has room; the builder does not write them.
- **No apo/holo pairing.** Structures are the modulator-bound ones AlloBench annotates. Blind prediction from an unbound structure needs apo counterparts this pipeline does not resolve.
- **Chain selection is by majority vote** over the labelled residues, falling back to all Cα when the named chain is absent. For a site spanning an interface this keeps one side.
- **Sequence redundancy is not clustered.** Grouping by UniProt stops the same protein straddling a split but not close homologues. For a stricter protocol, cluster sequences (MMseqs2 at 30 % identity) and group by cluster.
- **The annotations are inherited as given**, including any errors in the underlying ASD records. Nothing here re-curates the biology.
- **Positives are sparse and sites are small** (2.97 %, median 10 residues). Report AUPRC or a distance-controlled metric; plain accuracy is meaningless at this rate.

## 7. Redistribution

The coordinates come from RCSB and are public. The **residue-level annotations derive from ASD, which is research-use-only and may not be redistributed to third parties** — so the public pipeline repository carries no annotations and no structures, only the code that rebuilds them plus identifier-level metadata. Obtain `AlloBench.csv` from its own release under its own terms.

Anyone running the four commands in [README.md](README.md) reproduces this dataset bit-for-bit; the build has been run twice independently and matched on every key and every per-sample count.

## 8. Verification

`verify_dataset.py` re-reads every file from disk and checks, independently of the builder's own logic:

```
files on disk        : 1439
manifest rows        : 1439
fold keys (total)    : 1439
multi-site groups    : 18 (all distinct labels: True )
uniprot leak across folds: 0
sampled files checked: 200, mismatches: 0

VERDICT: PASS - dataset is internally consistent
```

The build aborts rather than writing a manifest whose rows do not each own a file.

## 9. Citing

Cite **AlloBench** and **ASD** for the data, and **RCSB PDB** for the structures. A link back to this repository is enough for the pipeline.
