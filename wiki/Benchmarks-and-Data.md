# Benchmarks and Data

Every dataset used in this project is public, and every processing step is scripted. Nothing
was hand-curated by us.

---

## IslandPick (the main benchmark)

**Source.** Langille MGI, Hsiao WWL, Brinkman FSL. *Evaluation of genomic island predictors
using a comparative genomics approach.* **BMC Bioinformatics** 2008;**9**:329.
DOI [10.1186/1471-2105-9-329](https://doi.org/10.1186/1471-2105-9-329)

The benchmark is distributed only as Excel additional files.
`scripts/00_prepare_benchmarks.py` downloads and parses them:

| Additional file | What we take from it |
|---|---|
| **AF 2** | the positive set — islands called by comparative genomics |
| **AF 4** | the negative set — backbone regions conserved across all comparison genomes |
| **AF 6** | per-genome TP/FP/TN/FN for the six tools evaluated in the original study |

### Contents

| | |
|---|---|
| Chromosomes | **118** |
| Total sequence | **484,996,900 bp** (mean 4.11 Mbp, range 1.16–8.26 Mbp) |
| Curated islands (positive set) | **771**, 12.4 Mbp |
| Curated backbone regions (negative set) | **3,700**, 50.6 Mbp |

### The evaluation protocol — read this before comparing anything

Accuracy is **not** measured over the whole chromosome. Only the two curated region sets
count:

```
TP = |predicted ∩ positive|      FN = |positive \ predicted|
FP = |predicted ∩ negative|      TN = |negative \ predicted|
```

Sequence belonging to neither set — most of the genome — is **ignored**. This matters: a
method is never penalised for calling a region the benchmark has no opinion about.

Metrics are then averaged **over genomes**, not over pooled nucleotides. Pooling would let
the largest chromosomes and the largest islands dominate the mean; Ibtehaz *et al.* make the
same point and use the same convention, so our numbers are directly comparable to theirs.

### Baselines we quote rather than re-run

From Additional File 6, per genome: **SIGI-HMM**, **Centroid**, **IslandPath-DIMOB**,
**PAI-IDA**, **IslandPath-DINUC**, **AlienHunter**. Quoting them (as SSG-LUGIA also did)
means every method in the comparison is scored on identical chromosomes under an identical
protocol, with no re-implementation risk on our side.

**SSG-LUGIA** is the one baseline we *do* re-run, because we need its per-island boundary
coordinates — which were never published — to compute boundary error. Our re-implementation
is faithful to the released code, including two non-standard formula choices; see
**[[API Reference]]**, `sagegi.features`.

---

## Case-study genomes

| Accession | Organism | Why |
|---|---|---|
| `NC_003198.1` | *Salmonella enterica* serovar Typhi CT18 | 19 well-characterised islands, 10 of them pathogenicity islands; the standard qualitative case study |
| `NC_002935.2` | *Corynebacterium diphtheriae* NCTC13129 | gram-positive, used for horizontally transferred gene identification |
| `NC_011770.1` | *Pseudomonas aeruginosa* LESB58 | gram-negative, same |

---

## Sequences

All chromosomes are fetched from **NCBI Entrez `nuccore`** by accession
(`scripts/01_download_genomes.py`), throttled to three concurrent requests to stay inside the
3 req/s limit, with retries. Nothing is redistributed in the repository — only the accession
list and the interval tables.

---

## Models

| Model | Checkpoint | Role |
|---|---|---|
| **DNABERT-S** | [`zhihan1996/DNABERT-S`](https://huggingface.co/zhihan1996/DNABERT-S) | the species-aware embedding channel |
| **DNABERT-2** | [`zhihan1996/DNABERT-2-117M`](https://huggingface.co/zhihan1996/DNABERT-2-117M) | general-purpose control — same architecture, same size, no species-aware contrastive stage |

Both are used **frozen**, at inference only. DNABERT-2 is the right control precisely because
it differs from DNABERT-S in exactly one respect: the species-aware training. Any gap between
them is attributable to that and not to capacity, tokenisation or architecture.

A relevant fact from the DNABERT-S paper: it was trained on sequences of **10,000 bp** — the
same length as SSG-LUGIA's window. Our 1 kb tiles aggregate to exactly that.

---

## Synthetic insertions

Curated benchmarks inherit the boundary uncertainty of the comparative-genomics procedure
that produced them, so they cannot settle boundary questions. The synthetic benchmark can.

Donor blocks are spliced into host chromosomes at random, well-separated positions that avoid
the host's own annotated islands, sweeping insert length (5, 10, 20, 40, 80 kb) and
donor–host compositional distance.

**Why it is exact and cheap.** Insert offsets and lengths are constrained to be multiples of
the tile length, so the tiling of the chimera is an exact interleaving of host and donor
tiles — no tile ever straddles a junction. The chimera embeddings can therefore be assembled
by *concatenating precomputed matrices* rather than re-running the model, with no
approximation at all — provided host and donor tiles are full-dimensional embeddings in one
shared basis (per-genome PCA archives are refused, since their bases disagree).
Compositional features **are** recomputed on the real chimeric sequence, so its 10 kb
windows see the junctions exactly as they would in a real analysis.

Donor–host distance is quantified by Euclidean distance between tetranucleotide profiles and
by |ΔGC|.

---

## What we release

| Path | Contents |
|---|---|
| `data/benchmarks/*.csv` | tidy IslandPick tables (positive, negative, baseline accuracies) |
| `results/*_per_genome.csv` | per-genome metrics for every method and ablation |
| `results/*_intervals.csv` | the predicted island coordinates themselves |
| `results/*_boundaries.csv` | per-island boundary errors |
| `results/sweep_*.csv`, `tuning_*.json` | the full hyper-parameter grid and CV selections |
| `results/synthetic_insertions.csv` | the synthetic sweep |

---

**Next:** **[[Reproducing the Paper]]** · **[[Method Overview]]**
