<div align="center">

# SAGE-GI

**S**pecies-**A**ware **G**enomic **E**mbeddings for **G**enomic **I**sland prediction

*Unsupervised, annotation-free genomic island detection from a single bacterial chromosome,
built on a frozen genome foundation model.*

[![Website](https://img.shields.io/badge/website-saimoon--oman.github.io%2FSAGE--GI-1f6feb)](https://saimoon-oman.github.io/SAGE-GI/)
[![Wiki](https://img.shields.io/badge/docs-wiki-0a7)](https://github.com/saimoon-oman/SAGE-GI/wiki)
[![Python](https://img.shields.io/badge/python-3.10%2B-3776ab)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

</div>

---

## What this is

Genomic islands (GIs) are blocks of DNA a bacterium acquired by horizontal gene transfer.
They matter because antibiotic-resistance determinants and virulence factors are
concentrated in them. For a newly sequenced, unannotated genome with no close relative in
the databases, the only practical way to find them is to look for stretches whose
*sequence composition* does not look like the rest of the chromosome.

Every unsupervised detector built on that idea — SSG-LUGIA, AlienHunter, Centroid,
IslandPath-DINUC — describes each window of the genome with **hand-designed** statistics:
GC skew, dinucleotide bias, codon usage, k-mer frequencies. None of them brings any
*pretrained* knowledge of what bacterial DNA looks like.

Meanwhile **DNABERT-S** is a genome foundation model trained specifically so that its
embeddings separate sequences *by species of origin*. That is almost exactly the signal a
GI detector needs — a genomic island is, after all, a stretch of DNA that came from a
different species. But DNABERT-S has only ever been used *across* genomes, for
metagenomic binning. **SAGE-GI asks whether it works *within* one.**

## The three problems that had to be solved

Plugging a transformer into a sliding-window scanner does not work, for three reasons.
SAGE-GI's contribution is the fixes.

| Problem | Fix |
|---|---|
| **Cost.** SSG-LUGIA scans a 10 kb window every 100 bp, so every base is covered ~100x. Embedding each window separately costs ~100x the genome length in transformer compute — roughly 48,500 CPU-hours for the 485 Mbp IslandPick benchmark. | **Tile-and-pool.** Embed *non-overlapping* 1 kb tiles once; reconstruct any window embedding as the overlap-weighted mean of the tiles it covers. Compute becomes **linear** in genome length (a ~100x reduction) and the 100 bp window step is preserved exactly. |
| **Conditioning.** A 768-dimensional embedding cannot support a robust covariance estimate, and most of its directions carry no within-genome variation at all. | **Host-subspace projection.** L2-normalise, then project onto the leading principal directions *of the chromosome being analysed*, so the host defines its own frame of reference. (We also tried subtracting a long-wavelength positional trend to remove replichore structure. Measured over all 118 genomes it *hurts* by 6.9 F1 — see below.) |
| **Blunt boundaries.** Window-level scores can only say "somewhere in this 10 kb". SSG-LUGIA additionally credits each window to its *leading* 100 bp, which shifts every predicted island systematically leftward. | **Centred projection + CUSUM refinement.** Credit windows to their centres, then place each edge at the maximiser of a normalised mean-shift statistic on the tile-resolution anomaly track. |

We expected a fourth ingredient — fusing the embedding with the eleven compositional
statistics — to be the thing that made it work. **It isn't.** Measured over all 118 genomes at
three operating points, adding the compositional block to the embedding never helps, and at the
recall-oriented point it *hurts* significantly (−1.1 F1, *p* = 2×10⁻⁴). SAGE-GI is therefore the
embedding-only pipeline, and the fusion is reported as an ablation that earns nothing. See
[the negative results](#two-things-that-did-not-work).

## Results

Evaluation follows the IslandPick protocol exactly (Langille *et al.* 2008): accuracy is
measured on the curated island and backbone region sets only, and metrics are averaged
**over genomes**, not over pooled nucleotides, so large genomes cannot dominate.

### Detection, 118 chromosomes

| method | P (%) | R (%) | F1 (%) |
|---|---|---|---|
| **SAGE-GI** (DNABERT-S embedding) | **58.6** | **73.1** | **60.3** |
| SAGE-GI + composition (ablation) | 58.4 | 72.5 | 59.6 |
| hand-crafted composition, 11 features | 55.2 | 68.5 | 56.4 |
| SSG-LUGIA-F (published state of the art) | 55.1 | 65.2 | 55.2 |
| DNABERT-2 embedding alone | 47.2 | 62.1 | 48.5 |
| Nucleotide Transformer v2 alone | 43.4 | 54.0 | 43.2 |

SAGE-GI improves F1 by **+5.2 points** over the published method (paired Wilcoxon
*p* = 8.6×10⁻⁵, better on 76 of 118 genomes) — improving precision **and** recall
simultaneously, which the trade-off in this field usually forbids. Adding the compositional
block on top changes F1 by −0.7 (*p* = 0.18), so it is not part of the method.

### Which representation carries the signal?

Detector-free Mann–Whitney AUC (island tiles vs backbone tiles), all 118 chromosomes:

| representation | AUC | contrast |
|---|---|---|
| **DNABERT-S** (species-aware) | **0.759** | 1.17 |
| hand-crafted composition (11 features) | 0.745 | **1.87** |
| DNABERT-2 (general-purpose) | 0.655 | 1.09 |
| Nucleotide Transformer v2 (general-purpose) | 0.547 | 1.03 |

- **DNABERT-S vs DNABERT-2: +0.105 AUC, *p* = 8.2×10⁻¹⁴, better on 97/118.** The two models
  share architecture, size, tokenizer and corpus, differing *only* in the species-aware
  contrastive stage — so that stage is what makes a foundation model useful here.
- **DNABERT-S vs Nucleotide Transformer v2: +0.212, *p* = 7.6×10⁻¹⁹, better on 106/118.**
  A third model from a different family, trained on a different corpus, lands at 0.547 —
  barely above chance. Two independent general-purpose genomic foundation models both fall
  *below* hand-crafted composition; only the species-aware one beats it. The claim is
  therefore about the training objective, not about DNABERT-S.
- **DNABERT-S vs composition: +0.014, *p* = 0.49 — indistinguishable.** A frozen,
  task-agnostic model matches two decades of hand-tuning.

### The competitor comparison, and why its number needs reading carefully

[TreasureIsland](https://academic.oup.com/bioinformaticsadvances/article/4/1/vbae089/7695232)
is the nearest published method. Run over the same 118 chromosomes under the same protocol it
scores **F1 75.2** against our 60.3. That number should not be read as a like-for-like result,
for two reasons we measured rather than asserted.

**It is supervised.** Its representation is unsupervised, but the classifier is an SVM shipped
pre-trained, fitted to the Benbow compilation whose main component is the revised version of
this benchmark. Its own per-genome similarity check flags only 11 of our 118 chromosomes as
unlike its training data:

| | genomes | P | R | F1 |
|---|---|---|---|---|
| chromosomes it recognises | 107 | 79.5 | 92.4 | **82.9** |
| chromosomes it flags as unfamiliar | 11 | 0.0 | 0.0 | **0.0** |

On those 11 it returns **no predictions at all**. SAGE-GI scores **59.7** and **66.3** on the
same two subsets, having no training distribution to be inside or outside of. The unfamiliar
subset is exactly the case this project is about.

**The benchmark cannot see a third of its calls.** The IslandPick protocol grades only ~13% of
a chromosome. Recomputing precision with every base outside a curated island counted as a
negative — a lower bound, since the unscored majority contains uncurated islands — shows how
far each method depends on sequence the benchmark never looks at:

| method | called % | P official | P strict | gap |
|---|---|---|---|---|
| TreasureIsland | 29.8 | 79.5 | **8.5** | **71.0** |
| **SAGE-GI** | 11.0 | 59.1 | 13.9 | 45.2 |
| SSG-LUGIA-F | 10.2 | 55.8 | 13.8 | 42.0 |

Every method loses 42–46 points — that is the cost of a partially scored benchmark, and it
means absolute precision is optimistic across this whole literature. TreasureIsland loses 71.
Our own gap is slightly larger than the baseline's, because we call slightly more sequence.

### Two things that did not work

**Positional detrending** — one of the two ideas this project set out to test — **costs 6.9 F1**
(*p* = 5×10⁻¹²) at every bandwidth and dimensionality tried. Replichore structure spans
megabases while islands reach 200 kb, so no filter is simultaneously wide enough to model the
trend and narrow enough to treat an island as an outlier.

**Fusing learned and hand-crafted features** — the other idea — earns nothing. The two
representations do separate islands differently (composition by magnitude, the embedding by
rank), and we expected that to compound. Over 118 genomes at three operating points it never
does:

| operating point | fusion | embedding alone | Δ | *p* |
|---|---|---|---|---|
| default | 59.6 | **60.3** | −0.69 | 0.18 |
| precision-oriented | 57.9 | **58.0** | −0.08 | 0.56 |
| recall-oriented | 48.1 | **49.2** | −1.12 | **0.0002** |

Both are off by default and reported as negative results. We would rather publish them than
have someone else spend a month rediscovering them.

Raw per-genome tables are in [`results/`](results/); figures on the
[project site](https://saimoon-oman.github.io/SAGE-GI/).

## Install

```bash
git clone https://github.com/saimoon-oman/SAGE-GI.git
cd SAGE-GI
pip install -r requirements.txt
```

`transformers` is pinned to a 4.x release: DNABERT-2/S ship custom modelling code that
predates the meta-device loader introduced in transformers v5.

## Reproducing the study

Everything is resumable — re-running any step skips work already on disk.

```bash
# 1. Benchmark: fetch and parse the IslandPick additional files into tidy CSVs
python scripts/00_prepare_benchmarks.py

# 2. Genomes: 118 IslandPick chromosomes + 3 case-study genomes from NCBI (~500 MB)
python scripts/01_download_genomes.py

# 3. Compositional baseline: faithful SSG-LUGIA (P / F / R variants) on all 118 genomes
python scripts/02_run_ssg_lugia.py --workers 3

# 4. Tile embeddings.  This is the only GPU step.  On a laptop CPU it is ~70 h per model,
#    so run notebooks/01_extract_embeddings_gpu.ipynb on a free Colab T4 instead (~30 min
#    per model) and unzip the result into $SAGEGI_WORK/emb/.
python scripts/03_extract_embeddings.py --model DNABERT-S --tile 1000

# 5. SAGE-GI and every ablation
python scripts/04_run_sage_gi.py --configs configs/methods.json --workers 3

# 6. Hyper-parameter sweep with grouped 5-fold cross-validation
python scripts/05_tune_cv.py --channel sage --workers 3

# 7. Controlled synthetic-insertion benchmark and the three case-study genomes
python scripts/06_synthetic.py --workers 3
python scripts/07_case_studies.py

# 8. Figures and paper tables
python scripts/90_make_figures.py && python scripts/91_make_tables.py
```

Heavy intermediates (FASTA, embeddings, cached features) go to `$SAGEGI_WORK`
(default `C:/sagegi_work`, override with the environment variable). Only small tidy CSVs
land in `results/`.

## Using SAGE-GI on your own genome

```python
from sagegi.genome import read_fasta
from sagegi.embed import EmbedConfig, embed_genome
from sagegi.sage import SageConfig, run_sage_gi

seq  = read_fasta("my_genome.fna")
emb, _ = embed_genome(seq, EmbedConfig(model_path="zhihan1996/DNABERT-S", tile=1000))

res = run_sage_gi(len(seq), SageConfig(), emb=emb, tile=1000)
for start, end in res.islands:
    print(f"{start+1}\t{end+1}\t{end-start+1} bp")
```

No annotation, no reference genome, and no training on your organism is required.

## Repository layout

```
sagegi/            library
  genome.py          FASTA I/O, integer encoding, k-mer identifiers
  kmer.py            vectorised sliding-window k-mer counting
  features.py        faithful SSG-LUGIA compositional features
  embed.py           tile embedding extraction + closed-form window pooling
  store.py           compact embedding archives (PCA-reduced, lossless downstream)
  sage.py            SAGE-GI: host projection, anomaly scoring, CUSUM refinement
  anomaly.py         two-stage robust (MCD / Mahalanobis) anomaly detection
  postprocess.py     smoothing, window->nucleotide projection, island filtering
  evaluate.py        IslandPick protocol metrics + boundary metrics
  synthetic.py       chimeric-genome construction with exact ground truth
scripts/           numbered, runnable pipeline stages
notebooks/         GPU embedding extraction (Colab / Kaggle)
configs/           method and ablation definitions
data/benchmarks/   tidy IslandPick tables (positive, negative, baseline accuracies)
results/           per-genome metrics, intervals, sweeps
paper/             LaTeX sources (OUP template), Overleaf-ready
docs/              GitHub Pages site
wiki/              GitHub wiki pages (copy into the wiki repo)
```

## Data and provenance

| Asset | Source |
|---|---|
| IslandPick benchmark: 118 chromosomes, 771 curated islands (12.4 Mbp), 3,700 backbone regions (50.6 Mbp) | Langille MGI, Hsiao WWL, Brinkman FSL. *BMC Bioinformatics* 2008;**9**:329, Additional Files 2, 4, 6 |
| Baseline per-genome accuracies (SIGI-HMM, Centroid, IslandPath-DIMOB, PAI-IDA, IslandPath-DINUC, AlienHunter) | Same, Additional File 6 |
| Genome sequences | NCBI Entrez `nuccore` |
| DNABERT-S | [`zhihan1996/DNABERT-S`](https://huggingface.co/zhihan1996/DNABERT-S) — Zhou *et al.*, *Bioinformatics* 2025;**41**:i255–i264 |
| DNABERT-2 | [`zhihan1996/DNABERT-2-117M`](https://huggingface.co/zhihan1996/DNABERT-2-117M) — Zhou *et al.*, arXiv:2306.15006 |
| Reference implementation of the compositional baseline | [`nibtehaz/SSG-LUGIA`](https://github.com/nibtehaz/SSG-LUGIA) — Ibtehaz *et al.*, *Brief Bioinform* 2021;**22**:bbab116 |

## Citing

```bibtex
@article{oman2026sagegi,
  title   = {SAGE-GI: species-aware genomic language model embeddings improve
             unsupervised genomic island detection and boundary resolution},
  author  = {Oman, Saimoon Al Farshi and Bayzid, Md. Shamsuzzoha},
  journal = {Briefings in Bioinformatics},
  year    = {2026}
}
```

## Authors

- **Saimoon Al Farshi Oman** — Department of Computer Science and Engineering, Bangladesh
  University of Engineering and Technology (BUET), Dhaka, Bangladesh
- **Dr. Md. Shamsuzzoha Bayzid** — Department of Computer Science and Engineering, BUET
  (supervisor)

Developed as a CSE6406 *Bioinformatics Algorithms* project.

## License

MIT — see [LICENSE](LICENSE). The DNABERT models carry their own (Apache-2.0) licence.
