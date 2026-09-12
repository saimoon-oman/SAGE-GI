# Method Overview

This page explains *why* SAGE-GI is built the way it is. For the code, see
**[[API Reference]]**; for how to run it, **[[Quick Start]]**.

---

## 1. The setting

A genomic island (GI) is a cluster of genes a bacterium acquired by horizontal gene
transfer. They matter because antibiotic-resistance determinants and virulence factors
concentrate in them.

Two families of method find them:

- **Comparative** methods (IslandPick, IslandViewer) align the query against close relatives
  and call the regions present in one and absent in others. Accurate — but they need those
  relatives to already exist in a database.
- **Composition-based** methods use only the query sequence, exploiting the fact that DNA
  from another organism carries a different compositional signature.

For a newly sequenced, unannotated, taxonomically isolated organism, only the second family
applies. That is the setting SAGE-GI targets.

## 2. The observation

Every composition-based detector describes a window with hand-designed statistics. Meanwhile
**DNABERT-S** is a genome foundation model whose training objective is, explicitly, that
*sequences from different species should separate in embedding space*.

A genomic island is, by definition, DNA that came from a different species.

That is an unusually direct match between what a model was trained to encode and what a task
needs to find — and DNABERT-S had only ever been applied *across* genomes, for metagenomic
binning. Nobody had asked what it does *within* one.

## 3. Why nobody had done it: three obstacles

### 3.1 Cost

Window-scanning detectors get their resolution from overlap. SSG-LUGIA slides a
**10 kb window in 100 bp steps**, so every base is covered by ~100 windows. Embedding each
window separately means pushing

```
L × w / Δw  =  L × 10000 / 100  =  100 L
```

bases through a 117M-parameter transformer. For the 485 Mbp IslandPick benchmark that is
**≈5 × 10¹⁰ bases** — hundreds of GPU-hours, tens of thousands of CPU-hours.

**Fix — tile-and-pool.** Split the chromosome into *non-overlapping* tiles of length
`T = 1 kb` and embed each once. Then the embedding of any window `[s, s+w)` is recovered
exactly as the overlap-weighted mean

```
ē(s) = Σᵢ oᵢ(s) · eᵢ  /  Σᵢ oᵢ(s) ,     oᵢ(s) = |[s, s+w) ∩ [iT, (i+1)T)|
```

computed in O(1) per window from a prefix sum over tiles plus fractional weights for the two
partial end tiles. Total transformer work drops from `100 L` to `L`.

Two details matter:

- Because the end weights are **fractional**, `ē(s)` changes smoothly as the window slides by
  100 bp, even though tiles are 1 kb apart. The window grid — and therefore the resolution of
  the anomaly signal — is *unchanged*.
- `T = 1 kb` is not arbitrary: ten tiles aggregate to the **10 kb sequence length DNABERT-S
  was contrastively trained on**, which is also SSG-LUGIA's window length.

The closed form is verified against a brute-force weighted mean in the test suite
(max error 3×10⁻⁸).

### 3.2 Conditioning the embedding to the host

Embeddings are **L2-normalised** (DNABERT-S optimises cosine similarity, so the unit sphere
is its natural geometry) and projected onto the leading `d = 24` principal directions **of
the chromosome being analysed**. That is not merely dimensionality reduction: a robust
covariance estimate in 768 dimensions is neither well-conditioned nor computationally
sensible, and the host should define its own frame of reference.

#### The part we got wrong: positional detrending

We expected a second step to be necessary, and it turned out not to be. The argument was
clean. A species-aware embedding is trained to separate *different* genomes, so inside **one**
chromosome its largest variation ought to be host structure rather than islands:

- the **GC-skew inversion** at the origin and terminus of replication,
- **replichore asymmetry** between leading and lagging strands,
- slow compositional **drift** along the chromosome.

The natural remedy is a high-pass filter along the genome coordinate — subtract a
long-wavelength running median from each coordinate:

```
z̃ᵢ = zᵢ − median_B(z)ᵢ ,        B = 250 kb / T
```

On a synthetic control this behaves exactly as designed: a genome-scale sinusoid is
attenuated from σ = 3.54 to σ = 0.03 while an inserted block of amplitude 8.0 survives
untouched.

**On real genomes it is harmful.** Swept over bandwidths of 100 kb, 250 kb and 1 Mbp and
dimensionalities 8/16/24/48, across all 118 chromosomes, it lowers island–backbone separation
at *every* setting, and in the full pipeline it costs **6.9 F1** (paired Wilcoxon
*p* = 5×10⁻¹²).

The reason is a collision of scales that no bandwidth escapes. Replichore structure spans
megabases; genomic islands reach 200 kb. A filter wide enough to model the trend treats an
island as part of the local baseline rather than as an outlier. Using a median instead of a
mean helps but cannot rescue it — a 134 kb island occupies more than half of a 250 kb window.

`detrend_bp` is kept in the code so the ablation stays reproducible, but it defaults to `0`.
We report this because the reasoning behind it is natural enough that others would otherwise
repeat the experiment.

### 3.3 Boundaries are not what F1 measures

A prediction can overlap an island generously — scoring well on nucleotide precision and
recall — while placing both edges tens of kilobases away. The metric the field reports simply
does not see this.

SAGE-GI reports, for every curated island a method overlaps at all, the absolute distance
between true and predicted start and end coordinates (**MABE**, median absolute boundary
error), and also the **signed** error, which separates bias from noise.

Doing so exposes something concrete: a window scanner that credits each window to its
**leading** `Δw` bases — as the reference implementation does — biases every predicted island
leftward by about `w/2`.

**Fix — centred projection plus CUSUM refinement.** Windows are credited to their centres.
Then each edge is placed at the maximiser of the normalised two-sample mean-shift statistic

```
k̂ = argmax_k  √( k(n−k)/n ) · | mean(a₁..a_k) − mean(a_{k+1}..a_n) |
```

over a local segment of the tile-resolution anomaly track `a` spanning 20 kb of putative host
outside the edge and 15 kb inside it. If refinement would invert an interval, the coarse
coordinates are kept.

## 4. The compositional channel (optional, and off by default)

The pooled, host-projected embedding block can be **concatenated** with the eleven
compositional features. No weighting constant is needed, because the detector scores windows by
**Mahalanobis distance**, which is invariant under any invertible affine transformation of the
feature space — the scale mismatch between the blocks is irrelevant by construction.

This was intended to be the method. It is not. The two configurations differ in exactly one
setting (`use_composition`), and over 118 genomes:

| operating point | fusion | embedding alone | Δ F1 | *p* |
|---|---|---|---|---|
| default | 59.63 | **60.32** | −0.69 | 0.18 |
| precision-oriented | 57.94 | **58.01** | −0.08 | 0.56 |
| recall-oriented | 48.09 | **49.21** | −1.12 | **0.0002** |

Precision, recall and boundary error agree, and so do the case studies (the fusion wins on
CT18, loses on *C. diphtheriae* and LESB58). Since the two are statistically tied at best, the
simpler method is the published one: `use_composition` defaults to `False`.

## 5. Detection and post-processing

Following SSG-LUGIA:

1. Fit an elliptic envelope (Minimum Covariance Determinant) with contamination `c₁` to all
   windows.
2. Re-fit with contamination `c₂` to the windows the first stage called native. Cascading
   recovers alien windows whose signal was masked by contamination of the first fit.
3. Median-smooth the score track over `f` windows.
4. Threshold at the detector's own decision boundary, project to nucleotides, discard runs
   shorter than 10 kb.

**One efficiency note.** Adjacent windows overlap by 99%, so the robust covariance can be
estimated from a decimated subsample. Fitting on every tenth window and scoring all of them
is **6.6× faster** (136 s → 21 s per chromosome) and agrees with the full fit on **99.84%**
of window calls (Pearson r = 0.99999). SAGE-GI uses this; the SSG-LUGIA baseline does not, so
that baseline is run exactly as published.

## 6. What is *not* claimed

- SAGE-GI is not a replacement for comparative or ensemble methods when close relatives and
  good annotation exist — those use strictly more information.
- The embedding does not beat well-tuned composition on *every* genome — it wins on 76 of 118.
  The contribution is that a frozen, task-agnostic model beats hand-engineered features on
  average without any training data, and that host referencing is what makes the learned
  channel usable at all.
- DNABERT-S was not trained for this task and is used entirely frozen. Fine-tuning it on
  curated islands is an obvious next step and is deliberately out of scope, since it would
  reintroduce the dependence on labelled data that motivates unsupervised detection.

---

**Next:** **[[Reproducing the Paper]]** · **[[API Reference]]** · **[[Benchmarks and Data]]**
