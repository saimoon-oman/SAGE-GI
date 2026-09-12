# SAGE-GI — Project Log

> **Purpose of this file.** A single place that records *what this project is*, *what has been
> done*, *what is running*, and *what remains*. Anyone (including a future me) should be able to
> read only this file and understand the full state of the work.

---

## 0. One-paragraph summary

**SAGE-GI** (*Species-Aware Genomic Embeddings for Genomic Island prediction*) is an unsupervised,
single-sequence, annotation-free pipeline that detects genomic islands (GIs) in bacterial genomes by
replacing the hand-crafted compositional features of **SSG-LUGIA** (Ibtehaz *et al.*, *Briefings in
Bioinformatics*, 2021) with frozen embeddings from **DNABERT-S** (Zhou *et al.*, *Bioinformatics* /
ISMB 2025), a genome foundation model explicitly trained to be *species-aware*. The core research
question: **do pretrained species-aware embeddings localise horizontally acquired DNA — and in
particular its boundaries — more accurately than hand-crafted compositional statistics?**

Target venue: **Briefings in Bioinformatics** (OUP), *Problem Solving Protocol* article type — the
same track as the SSG-LUGIA base paper.

Repository: <https://github.com/saimoon-oman/SAGE-GI> · Site: <https://saimoon-oman.github.io/SAGE-GI/>

---

## 1. Authors

| | |
|---|---|
| Saimoon Al Farshi Oman | Dept. of CSE, Bangladesh University of Engineering and Technology (BUET), Dhaka, Bangladesh — `1024052020@grad.cse.buet.ac.bd` |
| Dr. Md. Shamsuzzoha Bayzid | Dept. of CSE, BUET (supervisor) |

Origin: CSE6406 *Bioinformatics Algorithms* course project, Oct 2025 term.

---

## 2. Goal & contributions being claimed

1. **Tile-and-pool embedding extraction.** Naively embedding SSG-LUGIA's sliding windows
   (10 kb window / 100 bp step) costs ~100x the genome length in transformer forward passes
   (≈48,500 CPU-hours for IslandPick). SAGE-GI embeds *non-overlapping* 1 kb tiles once and
   mean-pools them into window embeddings — cost becomes **linear in genome length (1x)**, a ~100x
   reduction, with no loss of window resolution above 1 kb.
2. **Host-referenced whitening.** A species-aware embedding separates *different* genomes; used
   raw, within *one* genome its leading variance directions are dominated by nuisance structure
   (GC drift along the replichore, coding density, strand asymmetry). SAGE-GI centres and whitens
   tile embeddings against a robust estimate of the *host* distribution so that anomaly scoring
   operates on a residual "foreignness" axis. This is the step that makes a cross-genome embedding
   usable for intra-genome anomaly detection.
3. **Score-level fusion** of the compositional channel and the embedding channel, which are shown
   to be complementary rather than redundant.
4. **Change-point boundary refinement** replacing/augmenting SSG-LUGIA's median filter, evaluated
   with an explicit *boundary-distance* metric that prior GI work does not report.
5. **A controlled synthetic-insertion benchmark** with exact ground truth, giving detection
   sensitivity as a function of insert length, donor phylogenetic distance and compositional
   contrast.

---

## 3. Data

| Asset | Status | Location / provenance |
|---|---|---|
| IslandPick benchmark: 118 chromosomes, positive (GI) + negative (non-GI) regions | **acquired** | Langille *et al.* 2008, *BMC Bioinformatics* 9:329 — Additional Files 2 (positives), 4 (negatives), 6 (per-genome TP/FP/TN/FN for 6 baseline tools). Downloaded from Springer static content. |
| Per-genome accuracies of SIGI-HMM, Centroid, IslandPath-DIMOB, PAI-IDA, IslandPath-DINUC, AlienHunter | **acquired** | Additional File 6 (same source) — this is exactly the table SSG-LUGIA used for its Table 2. |
| 118 genome FASTA sequences | *in progress* | NCBI Entrez `nuccore` by accession |
| *S. typhi* CT18 (NC_003198.1) + 19 curated GIs | pending | Lu & Leong 2016; SSG-LUGIA sample data |
| *C. diphtheriae* NCTC13129, *P. aeruginosa* LESB58 | pending | NCBI |
| DNABERT-S / DNABERT-2 weights | **acquired** | HuggingFace `zhihan1996/DNABERT-S`, `zhihan1996/DNABERT-2-117M` (447 MB each) |

Total benchmark size: **484,996,900 bp** across 118 chromosomes (mean 4.11 Mbp, range 1.16–8.26 Mbp).

---

## 4. Compute situation (important)

Local machine: Intel i5-8250U, 4 physical cores, 12 GB RAM, **no GPU**, `torch 2.12.0+cpu`.

Measured DNABERT-S throughput on this machine: **~2,000 bp/s** (fp32, batch 16, 1 kb tiles).
- int8 dynamic quantization was tested and **rejected**: mean-pooled embedding cosine similarity to
  fp32 collapsed to 0.08–0.39, i.e. the embeddings are destroyed.
- Full IslandPick at 1x coverage ⇒ ~70 h **per model**, ~140 h for both. Not viable locally.

**Decision (user, 2026-08-31): run embedding extraction on a free Colab/Kaggle GPU.** On a T4 this
is ~3–5 GPU-hours for both models over all 118 genomes at full coverage. Notebooks are a project
deliverable so the run is reproducible.

In parallel a **local CPU pilot run** is used to develop and validate the downstream pipeline on
real embeddings, ordered by scientific priority (CT18 case study → stratified IslandPick subset →
remainder).

---

## 5. Progress

*Rewritten 2026-09-01. The original version of this section was written before any experiment
had run and had become misleading — it still listed the paper, the website and the whole
evaluation as "to do".*

### Complete

- [x] Framing, benchmark recovery, reference SSG-LUGIA implementation and its validation
      (reproduced the published F1 to within 0.09).
- [x] Genome acquisition: 118 benchmark chromosomes + case-study genomes.
- [x] Tile-and-pool embedding extractor, resumable, with Colab and Kaggle notebooks. Verified
      exact against direct windowing to 3e-8, and fp16-vs-fp32 to cosine 0.99999.
- [x] SAGE-GI: host-subspace projection, robust anomaly scoring, CUSUM boundary refinement,
      score fusion.
- [x] Full 118-genome evaluation, 17 configurations, with per-genome averaging.
- [x] Ablations: composition / DNABERT-2 / DNABERT-S / fusion; positional detrending (a
      negative result, reported as such); tile size; `fit_stride`.
- [x] Detector-free representation analysis (Mann-Whitney AUC), which is the cleanest result in
      the paper.
- [x] Grouped 5-fold cross-validation over genomes; paired Wilcoxon tests throughout.
- [x] Synthetic insertion benchmark.
- [x] TreasureIsland comparison, run directly rather than quoted.
- [x] Benchmark-wide functional enrichment — including the differential analysis that made us
      withdraw a claim.
- [x] Paper (LaTeX, Overleaf-ready), website (`docs/`), README, `wiki/`,
      `project_explanation_bangla.docx`, `ASSESSMENT.md`.
- [x] 20 unit tests.

### Running now

- [ ] `finalise.sh`: final detector run over all 118 genomes, then case studies, synthetic
      insertions, TreasureIsland, biology, tables, figures. Resumable; safe to interrupt.

### Waiting on the user

- [ ] Nucleotide Transformer v2 embeddings from Kaggle. Everything downstream is pre-wired —
      configurations, pipeline, tables — so this is a drop-in.

### Deliberately not done

- [ ] The revised 2019 IslandPick benchmark (104 genomes, 20 methods). This is the strongest
      remaining addition and is listed as such in `ASSESSMENT.md`; it is a scope decision, not
      an oversight.
- [ ] Fine-tuning DNABERT-S. Doing so would reintroduce the dependence on labelled data that
      the whole approach exists to avoid.

## 6. Changelog

- **2026-08-31** — Project started. Literature read, benchmark data recovered, models verified,
  compute strategy decided, repository initialised.

- **2026-08-31 (cont.)** — Core library complete and unit-tested:
  `genome`, `kmer`, `features` (faithful vectorised SSG-LUGIA), `anomaly`, `postprocess`,
  `evaluate`, `embed`, `store`, `sage`, `synthetic`. Verified `pool_windows` against a
  brute-force weighted mean (max error 3e-8); verified chimera construction, tile-provenance
  and coordinate mapping; verified positional detrending removes a genome-scale ramp
  (sigma 3.54 -> 0.03) while preserving a sharp inserted block (8.0 -> 8.0).
  Measured decimated MCD fitting: `fit_stride=10` is **6.6x faster** (136 s -> 21 s per
  genome) with **99.84 % identical calls** (Pearson 0.99999) — adopted for SAGE-GI, left off
  for the faithful SSG-LUGIA baseline. GPU extraction notebook written, validated and handed
  to the user. Local CPU embedding paused: with the baseline sweep occupying the cores the
  machine delivers ~300 bp/s, so embeddings will come from the GPU run.

- **2026-08-31 (afternoon)** — Deliverable scaffolding complete.
  * `README.md`, `LICENSE` (MIT), `requirements.txt`, `.gitignore`, `CITATION.cff`.
  * **20-case unit-test suite** (`tests/test_sagegi.py`) — all passing. It caught one real
    off-by-one in my own test expectation for `merge_close` (the library was right).
  * **Baseline faithfulness check.** On the first 18 IslandPick chromosomes our
    re-implementation gives SSG-LUGIA-F P 52.2 / R 73.3 / F1 56.4 against the published
    all-118 figures of 55.1 / 65.2 / 55.1, and SSG-LUGIA-R 38.9 / 95.8 / 50.1 against
    40.0 / 86.8 / 49.1. The re-implementation tracks the published algorithm closely.
  * **Figure system.** Palette taken from a CVD-validated categorical set; the skill's
    JS validator was ported to Python (`C:/sagegi_work/validate_palette.py`) and
    reproduces its documented figures exactly (adjacent CVD dE 9.1 protan, normal-vision
    19.6). Found that blue/orange/aqua/violet clears the stricter *all-pairs* gates, which
    is exactly the four representations compared here, so scatter plots use that set. An
    earlier attempt to distinguish SSG-LUGIA-P/F/R by invented orange tints FAILED the
    normal-vision floor and was replaced by family colour + axis labels.
  * **Paper.** Official OUP `oup-authoring-template` class from CTAN; compiles locally with
    MiKTeX, no errors, no undefined references. Introduction and Methods written. All
    result numbers come from `paper/numbers.tex`, regenerated by `scripts/91_make_tables.py`
    from `results/`, so the prose can never drift from the tables.

- **2026-08-31 (late afternoon)** — Deliverables 4, 5, 7 complete.
  * **Website** (`docs/`) — self-contained static site for GitHub Pages, theme-aware,
    responsive, no external requests. Verified in-browser: stylesheet loads, tokens apply,
    cards and figures render, responsive stacking works. Removed `backdrop-filter` from the
    sticky header after it broke compositing (and `color-mix()` alongside it) — a portability
    fix, not just a screenshot fix.
  * **Wiki** (`wiki/`) — 10 pages: Home, Installation, Quick Start, Method Overview,
    Benchmarks and Data, Pipeline Reference, API Reference, Reproducing the Paper,
    Troubleshooting, plus `_Sidebar`/`_Footer` and a README explaining how to publish to the
    separate wiki repository. All `[[links]]` verified to resolve.
  * **`project_explanation_bangla.docx`** — 13 pages, ~3,200 words, built with python-docx
    with the complex-script font slot (`w:cs`/`w:szCs`) set on every run so Bengali conjuncts
    render correctly in Nirmala UI. Verified by exporting to PDF through Word COM and reading
    the rendered pages. Covers the biology from scratch, both base papers, all three
    obstacles and fixes, the pipeline, the data, the contributions, the limitations, how to
    run it, and a 20-term glossary.
  * Armed a chained job: when the baseline finishes, the freed cores start a local CPU
    embedding pilot (9 genomes incl. CT18) so the SAGE-GI path can be validated on real
    embeddings regardless of when the GPU run lands.

---

- **2026-08-31 (evening)** — Baseline complete; three notebook bugs found and fixed.

  **SSG-LUGIA baseline finished on all 118 chromosomes, and it validates.** Averaged over
  genomes, our re-implementation lands within **0.09 F1 points** of the published Table 2 on
  every variant:

  | variant | ours (P / R / F1) | published (P / R / F1) | dF1 |
  |---|---|---|---|
  | SSG-LUGIA-P | 60.29 / 54.97 / 52.33 | 60.28 / 55.04 / 52.38 | −0.05 |
  | SSG-LUGIA-F | 55.28 / 65.18 / 55.17 | 55.10 / 65.22 / 55.08 | +0.09 |
  | SSG-LUGIA-R | 40.06 / 86.72 / 49.10 | 40.01 / 86.75 / 49.06 | +0.04 |

  Separately, the six published baselines recomputed from Additional File 6 match Langille's
  own table to within 0.12 on every metric, confirming the benchmark parsing and the metric
  code.

  **Boundary finding (needs care in the write-up).** Mean signed errors for SSG-LUGIA-F are
  start −19.4 kb, end +14.8 kb: predictions are not merely shifted left, they are roughly
  34 kb *wider* than the curated islands. Centred projection moves both edges right by
  ~4 kb (start −19.4 → −15.4) and improves median Jaccard 0.294 → 0.306, but the dominant
  error is over-extension, not translation. The paper must say this precisely rather than
  claim a pure "half-window leftward shift".

  **Bugs found by running the GPU notebook, reported by the user and reproduced locally:**
  1. *fp16 crash.* `RuntimeError: expected scalar type Float but found Half`. Two causes,
     both in the vendored model code: the extended attention mask is built with a hard-coded
     `.to(dtype=torch.float32)` (its "fp16 compatibility" comment is true only for the Triton
     kernel, which casts internally), and the ALiBi tensor is a plain attribute rather than a
     registered buffer, so `.half()` never converts it. Both feed `bias`, which promotes the
     scores back to fp32. Fixed with a forward pre-hook casting `bias` to the module dtype;
     verified fp16-vs-fp32 cosine **0.999990 min / 0.999996 mean**, and the ALiBi rebuild
     path (>512 tokens) works too.
  2. *DNABERT-2 will not load at all* on transformers 4.46.3 — independent of load order.
     `AutoModel.register` rejects it because its model classes inherit `config_class` from
     the stock `BertPreTrainedModel` while its own `configuration_bert.BertConfig` subclasses
     `PretrainedConfig` directly. Fixed by loading the remote class explicitly via
     `get_class_from_dynamic_module`, which skips registration; both checkpoints now load in
     one process.
  3. *`nvidia-smi` missing crashed cell 1* on any machine without the NVIDIA driver.

  Fixes applied to both `sagegi/embed.py` and the notebooks. Added an automatic fp32
  fallback if the precision check ever fails, so the run cannot be corrupted silently.
  Verified by executing every notebook cell locally, both models, end to end.

  Also: a **Kaggle** notebook variant, and the site moved `website/` → `docs/` because
  GitHub Pages branch deployment only offers root or `/docs`.

---

## 2026-08-31 (night) — GPU embeddings received; first honest look at the science

The user's Kaggle run produced both archives. Integrity verified: **120 genomes per model
(all 118 IslandPick chromosomes + 3 case studies), 494,025 tiles = 494 Mbp embedded, no
missing or malformed files, tile counts consistent with genome lengths.**

### Diagnostic analyses run before the main evaluation

These are detector-free: for each representation, score every tile by a fixed anomaly
statistic and ask how well it ranks curated island tiles above curated backbone tiles
(Mann-Whitney AUC). Differences are then attributable to the representation alone.

**On CT18 (19 curated islands):**

| representation | AUC | in-island | outside |
|---|---|---|---|
| compositional (11 features) | **0.817** | 4.94 | 2.88 |
| DNABERT-S, d=16, no detrend | 0.809 | 4.38 | 3.82 |
| DNABERT-S, d=24, no detrend | 0.805 | 5.42 | 4.72 |
| DNABERT-2, d=48, no detrend | 0.792 | 7.91 | 6.87 |
| DNABERT-S, d=24, detrend 250 kb | 0.743 | 5.31 | 4.76 |

Three things here contradict my working assumptions, and all of them have to go into the
paper as stated rather than as hoped:

1. **Positional detrending HURTS** (0.805 -> 0.743 on CT18). The 250 kb running median is
   removing signal, not just nuisance: CT18's islands are large (one is 134 kb) and
   island-dense regions occupy enough of the filter window to move even a median. The
   early-returning genomes of the main run agree — `SAGE-GI/no-detrend` is the best variant
   on most of them so far.
2. **DNABERT-2 is roughly as good as DNABERT-S.** The species-aware contrastive stage
   confers no clear advantage for this *intra*-genome task. That directly weakens the
   paper's original central hypothesis.
3. **Hand-crafted composition is at least as informative as either embedding**, and has a far
   better contrast ratio (1.72 vs ~1.15).

### A genuinely better anomaly statistic — which still does not help detection

The mean of the unit embeddings is pulled toward the minority of atypical tiles; the median
is not. Their difference estimates, with no labels, the direction along which foreign DNA
departs from the host bulk. Projecting each tile onto it gives a single well-conditioned
discriminant. Matched comparison on the 6 genomes whose full 768-d embeddings were archived:

| statistic | mean AUC | median AUC |
|---|---|---|
| **mean-minus-median projection (768-d)** | **0.842** | **0.899** |
| unscaled distance, PCA d=24 | 0.800 | 0.842 |
| cosine to centroid, raw 768-d | 0.789 | 0.834 |
| robust distance, PCA d=16 (what SAGE-GI currently uses) | 0.779 | 0.867 |
| robust distance, PCA d=24 | 0.762 | 0.838 |

Two caveats, both important:

* **It does not survive PCA compression** (0.842 -> 0.685). The discriminating direction
  lives in the *low-variance* directions, so the PCA-96 archive cannot reproduce it. Using it
  benchmark-wide would need the full 768-d embeddings re-exported.
* **More importantly, the AUC gain does not become an F1 gain.** Run through the full
  pipeline on those same 6 genomes (detrend off, refine off):

  | configuration | P | R | F1 | MABE (kb) |
  |---|---|---|---|---|
  | alien-direction + composition | 77.31 | 81.50 | **78.03** | 18.32 |
  | composition only | 77.09 | 81.26 | 77.88 | 18.23 |
  | PCA-24 + composition | 70.62 | 84.00 | 75.16 | 19.44 |
  | alien-direction alone (1-D) | 70.46 | 72.53 | 68.74 | 18.87 |

  +0.15 F1 over composition alone is noise. The embedding signal is largely **redundant**
  with hand-crafted composition.

One suggestive exception: on *P. aeruginosa* PAO1 (NC_002516.2), the genome where
composition does worst (F1 35.3), fusion helps materially (41.5). If that pattern holds
across the benchmark — embeddings helping precisely where composition fails — it is a real
and publishable finding, just a much narrower one than originally proposed. The full
118-genome run will settle it.

> **Status: the original hypothesis is not being supported by the data.** Full results
> pending; no claim will be written that the numbers do not support.

### Full 118-genome results — the earlier pessimism was wrong

The negative read above came from a single genome (CT18) with **positional detrending on**.
With all 118 chromosomes and detrending off, the picture reverses.

**Detector level** (IslandPick protocol, averaged over genomes; paired Wilcoxon vs SAGE-GI):

| method | P | R | F1 | dF1 vs SAGE-GI | p |
|---|---|---|---|---|---|
| **SAGE-GI (fusion, no detrend)** | **58.42** | **72.51** | **59.63** | — | — |
| Comp-only + centred | 55.17 | 68.54 | 56.36 | −3.27 | 1.2e-03 |
| SSG-LUGIA-F (published algorithm) | 55.10 | 65.22 | 55.17 | −4.47 | 1.0e-04 |
| Comp-only (= our SSG-LUGIA reproduction) | 54.99 | 65.27 | 55.03 | −4.60 | 4.5e-05 |
| SAGE-GI with 250 kb detrending | 50.23 | 68.02 | 52.72 | −6.91 | 5.0e-12 |
| DNABERT-S embedding alone | 48.99 | 63.65 | 49.81 | −9.82 | 5.8e-13 |
| SAGE-GI fused with DNABERT-2 instead | 43.85 | 57.98 | 44.93 | −14.70 | 1.5e-15 |
| DNABERT-2 embedding alone | 39.92 | 48.86 | 39.46 | −20.17 | 4.6e-18 |

SAGE-GI beats the published state of the art **on precision, recall and F1 simultaneously**
(+3.3 / +7.3 / +4.5), winning on 72 of 111 comparable genomes. Most methods in this field
buy one metric with another; this does not.

**Representation level** (detector-free Mann-Whitney AUC, all 118 genomes, no detrending):

| representation | best dims | AUC | contrast |
|---|---|---|---|
| DNABERT-S | 16 | **0.759** | 1.17 |
| compositional (11 hand-crafted features) | 11 | 0.745 | **1.87** |
| DNABERT-2 | 48 | 0.655 | 1.09 |

| contrast | dAUC | p | wins |
|---|---|---|---|
| DNABERT-S vs DNABERT-2 | **+0.105** | 8.2e-14 | 97/118 |
| compositional vs DNABERT-2 | +0.091 | 5.7e-07 | 82/118 |
| DNABERT-S vs compositional | +0.014 | 0.49 (n.s.) | 59/118 |

The precise statement the data supports:

1. **Species-aware pre-training is what makes a genome foundation model useful here.**
   DNABERT-S beats DNABERT-2 by 0.105 AUC and 10.4 F1 — the two models are identical in
   size, architecture and tokenizer, differing only in that contrastive stage.
2. **A frozen, task-agnostic DNABERT-S is as informative as decades of hand-tuned
   compositional statistics** (AUC difference not significant), which is the headline result.
3. **They are informative in different ways** — composition separates by magnitude
   (contrast 1.87 vs 1.17), the embedding by rank — which is why fusing them beats either.
4. **Positional detrending does not work.** Tested at 4 bandwidths x 4 dimensionalities x
   118 genomes, it is harmful everywhere. Reported as a negative result; default is now off.

Corrected the module docstring, the config registry and `SageConfig.detrend_bp` (now 0).
20/20 unit tests still pass.

---

## 2026-09-01 — tuning, case studies, biology; project essentially complete

### Tuning and cross-validation (`scripts/05_tune_cv.py`, 118 genomes x 108 configurations)

| objective | c1 | c2 | filter | min len | refine | in-sample P/R/F1 | 5-fold held-out P/R/F1 |
|---|---|---|---|---|---|---|---|
| F1 (SAGE-GI-F) | 0.15 | 0.05 | 400 | 10000 | yes | 58.42 / 72.51 / 59.63 | **57.69 / 72.14 / 59.10** |
| precision (SAGE-GI-P) | 0.10 | 0.05 | 400 | 10000 | yes | 60.50 / 64.89 / 57.94 | 59.84 / 64.15 / 57.42 |
| recall (SAGE-GI-R) | 0.20 | 0.25 | 200 | 8000 | no | 37.60 / 92.98 / 48.09 | 37.95 / 92.67 / 48.17 |

Two things matter here. The F1-optimal configuration is **exactly the defaults inherited from
SSG-LUGIA**, so the improvement is not a product of re-tuning. And the held-out F1 (59.10) sits
0.53 below in sample (59.63), with 4 of 5 folds choosing the same configuration — the gain
generalises to genomes the selection never saw.

### Case studies (whole-chromosome evaluation against curated island lists)

| genome | SSG-LUGIA-F | SAGE-GI | best published (quoted) |
|---|---|---|---|
| *S.* Typhi CT18 | 75.41 | **80.11** | GIHunter 74.4 |
| *C. diphtheriae* NCTC13129 | 19.28 (P variant) | 30.09 | — |
| *P. aeruginosa* LESB58 | 60.69 | 73.52 (DNABERT-S alone 76.35) | — |

On CT18 SAGE-GI would be the best method in the published comparison table. The two harder
genomes show the complementarity pattern directly: where composition collapses, the embedding
channel carries the result.

### Biological relevance (`scripts/09_biological_relevance.py`)

Using the 4,111 annotated genes of CT18 and a conservative keyword set, genes inside predicted
islands versus outside (odds ratio, one-sided Fisher):

| method | mobility | virulence | any marker |
|---|---|---|---|
| curated islands (reference) | 10.25 | 6.51 | 6.26 |
| **SAGE-GI** | **8.45** (p = 8e-21) | **4.54** (p = 3e-18) | **4.79** (p = 4e-31) |
| DNABERT-S alone | 8.36 | 4.36 | 4.69 |
| SSG-LUGIA-F | 6.63 | 4.24 | 4.14 |

SAGE-GI's calls are more functionally coherent than the baseline's, not merely more numerous —
the improvement in benchmark agreement is accompanied by an improvement in biology. No method
enriches for resistance genes, which is correct rather than a failure: CT18 carries its
multidrug resistance on the IncHI1 plasmid pHCM1, not on the chromosome.

### Deliverables brought in line with the corrected science

Positional detrending was presented as a contribution in the first drafts of every document.
It is now reported as a measured negative result in: `sagegi/sage.py` (docstring and default),
`configs/*.json`, the paper (Methods 3.2, Results, Discussion), `README.md`,
`wiki/Home.md`, `wiki/Method-Overview.md`, `wiki/API-Reference.md`, `wiki/Troubleshooting.md`,
`wiki/Reproducing-the-Paper.md`, `docs/index.html`, the pipeline schematic (Fig. 1) and the
Bangla document.

### Engineering notes from this session

* `\d2fused` in the generated `numbers.tex` silently broke the LaTeX build: `\d` is a real
  accent command, so `2fused` leaked into the preamble as text and produced
  "Missing \begin{document}" 200 lines later. `scripts/91_make_tables.py` now refuses any
  macro name that is not purely alphabetic.
* The generator pre-registers every expected macro as `TBD` and always writes every table file
  (with a visible placeholder), so an unfinished pipeline stage shows up in the PDF instead of
  breaking the compile.
* `sagegi.viz.save()` now mirrors figures into `paper/figures/` and `docs/assets/`
  automatically — that copy step had been silently forgotten after a re-run.
* Manuscript body is **4,871 words** against the BiB limit of 5,000
  (`scripts/92_wordcount.py`).

### Written this session

`ASSESSMENT.md` — an honest evaluation of publication prospects: measured evidence, three
simulated referee reports, calibrated acceptance estimates (~40% as-is, ~60-65% after the
three priority additions), and a prioritised list of what to add. The three priorities are a
TreasureIsland comparison, a third foundation model, and extending the biological analysis
beyond CT18.

---

## 2026-09-01 (afternoon) — incident, and the two comparisons the assessment called for

### Incident: the detector run was killed twice

Both times the background job died with no traceback — killed externally when the session was
torn down, not by any fault in the code. The first loss was the CPU embedding pilot; the
second was `04_run_sage_gi.py` at genome 46 of 118, costing about two hours.

Two fixes, both worth having anyway:

* **`04_run_sage_gi.py` is now resumable.** Per-genome results are cached under
  `$SAGEGI_WORK/sage_cache/`, keyed by accession *and* a SHA-1 digest of the configuration
  file, so editing a configuration invalidates the cache rather than silently reusing stale
  predictions. An interrupted run now resumes instead of restarting.
* **The pipeline is one sequential script** (`finalise.sh`) rather than parallel jobs racing
  for four cores. The earlier `after_chain.sh` used `pgrep -f` to wait for the detector run to
  finish; that never matched under this shell, so TreasureIsland started immediately and
  competed for memory (free RAM fell to 1.1 GB of 12.5 GB). Waiting on a *process name* was
  the wrong mechanism; sequencing the stages in one script is the right one.

### TreasureIsland: the comparison, and why the published numbers cannot be quoted

TreasureIsland (Banerjee *et al.* 2024) is the closest competitor — also unsupervised, also
single-sequence, also embedding-based, but with a Doc2Vec representation trained from scratch
rather than a pretrained foundation model.

**Its published F1 of 0.86 is not comparable to our 0.60.** It was evaluated on the Benbow
compilation, on 20 genomes, as a *classification of 566 pre-defined regions* (413 islands
against 153 non-islands). That is a balanced decision over given intervals. Ours is a
genome-wide scan in which the method must also decide where the intervals are. Putting the two
numbers side by side would be meaningless, and a referee would notice.

So `scripts/10_treasureisland.py` runs it on our 118 chromosomes and scores it with the same
code that scores everything else. It is installed in its own virtual environment (it pins its
own dependency stack) and invoked as a subprocess. It runs at roughly 17 kb/s, and reports its
own `out_of_distribution` flag, which matters: it was trained on Pseudomonadota and Firmicutes,
and IslandPick spans much more than that.

### A property of the benchmark that nobody reports

Measuring what fraction of each chromosome a method calls turned up something that belongs in
the paper:

| | % of genome called |
|---|---|
| curated islands (the positive set) | 2.4 |
| curated backbone (the negative set) | 10.7 |
| **unscored** | **87.0** |
| SSG-LUGIA-F | 10.2 |
| SAGE-GI | 12.9 |
| SSG-LUGIA-R | 29.6 |

The IslandPick protocol scores only **13.1%** of each chromosome. A method is not penalised
for false positives in the other 87%. This does not invalidate the benchmark — every method
here is scored identically — but it does mean the absolute numbers understate false-positive
rates, and a very liberal predictor is under-penalised. We now report the called fraction
alongside precision and recall.

### Benchmark-wide biology, including the question that matters most

`scripts/11_biology_benchmark.py` extends the CT18 enrichment analysis to the whole benchmark
using NCBI feature tables, and adds the differential test a referee will ask for: are the
islands SAGE-GI finds that SSG-LUGIA **misses** functionally coherent, or compositional noise?

First six genomes (full run in progress) — pooled odds ratios:

| region set | mobility | virulence | any marker |
|---|---|---|---|
| curated islands | 8.15 | 1.78 | 3.97 |
| both methods agree | 7.93 | 2.14 | 3.64 |
| **SAGE-GI only (the extra calls)** | **1.70** | **0.99** | **1.23** |

This is a genuinely uncomfortable result and it will go into the paper as stated. The extra
recall is not worthless — mobility enrichment of 1.70 is significant (p = 3e-5) — but it is far
weaker than the calls the two methods agree on. The F1 gain is therefore partly bought with
regions whose biological status is much less certain. Waiting for the full 118-genome run
before writing the claim.

### Prepared for the third foundation model

The user is running Nucleotide Transformer v2 (50M, multi-species) on Kaggle — a second
*general-purpose* genomic model, the control that would upgrade "DNABERT-S beats DNABERT-2"
into "species-aware pre-training is the ingredient". `notebooks/02_extract_third_model.ipynb`,
`configs/methods_nt.json`, `scripts/08_representation_analysis.py` and
`scripts/91_make_tables.py` are all wired for it; it slots in with no further changes.

---

## 2026-09-01 (evening) — the benchmark-wide biology, and a claim I had to withdraw

### CT18 was not representative

`scripts/09` (CT18, protein table) said SAGE-GI's islands were *more* functionally coherent
than the baseline's — mobility odds ratio 8.45 against 6.63. `scripts/11`, run over the whole
benchmark with NCBI feature tables, says the opposite:

| region set | mobility | virulence | any marker |
|---|---|---|---|
| curated islands (reference) | 6.14 | 1.30 | 3.28 |
| SSG-LUGIA-F | **5.25** | 2.36 | 3.28 |
| SAGE-GI | 4.55 | 1.86 | 2.74 |
| regions both methods call | 5.26 | 2.15 | 3.18 |
| **regions only SAGE-GI calls** | **1.63** | 1.18 | 1.36 |

The regions the two methods agree on are strongly enriched and essentially match the curated
islands. The regions SAGE-GI *adds* — the source of its higher recall — are enriched at only
1.63x. That is not noise (p = 1.2e-53 over 1,285 marker genes), but it is far thinner than the
agreed calls, and because SAGE-GI calls more sequence its pooled enrichment ends up *below*
the more conservative SSG-LUGIA-F.

**The CT18-based claim has been withdrawn from the paper** and replaced with the benchmark-wide
result plus a new subsection, "Where the extra recall comes from — and what it costs". The
honest framing: nucleotide F1 rewards recovering annotated island sequence and never asks
whether that sequence carries the functions that make islands interesting. High-recall users
should take SAGE-GI; users who need every call defensible should take the precision-oriented
operating point. Both are reported rather than presenting one number as if the choice did not
exist.

This is the single most important correction in the project so far, and it came from building
the analysis a sceptical referee would have demanded rather than waiting to be asked.

### The protocol scores only 13% of each chromosome

Now stated explicitly in Methods: curated islands cover 2.4% of an average chromosome and
curated backbone 10.7%, so 87% is unscored and predictions there are neither rewarded nor
penalised. Comparisons stay fair because every method is scored identically, but absolute
precision is optimistic and liberal predictors are under-penalised. The fraction of each
chromosome a method calls is now reported alongside (SSG-LUGIA-F 10.2%, SAGE-GI 12.9%).

### Nucleotide Transformer notebook: my bug, fixed

The user's run failed with
`ImportError: cannot import name 'find_pruneable_heads_and_indices' from
'transformers.pytorch_utils'`.

My fault: I wrote `pip install -U transformers` in that notebook on the assumption that NT-v2
uses the standard loader and therefore needs no pin. It does use `trust_remote_code`, and its
bundled `modeling_esm.py` imports two symbols that transformers 5 removed. Fixed by pinning
`transformers==4.46.3` (the same pin the DNABERT notebook uses), pinning the *model* revision
so an upstream edit cannot silently change the embeddings, adding an assert that fails fast if
the wrong version is imported, and adding the fp16 fidelity check with automatic fp32
fallback.

Verified locally end to end: loads under 4.46.3 (55.9 M parameters, 512-d embeddings), fp16
versus fp32 cosine 0.999999, full 1.16 Mbp genome embedded at 2,713 bp/s on CPU, and the
output `.npz` loads through `sagegi.store.load_tiles` as pipeline-compatible.

### Pipeline robustness

Two further problems, both mine, both fixed:

* Adding the NT-v2 configurations silently disabled the entire detector run.
  `04_run_sage_gi.py` selected genomes by intersecting *every* referenced embedding directory,
  so a configuration for a model with no embeddings yet made the intersection empty and the
  script reported "nothing to do". It now selects on cached compositional features and skips
  individual configurations per genome, printing a note for any model directory that is empty
  or partial.
* Two pipeline instances ended up running concurrently after a restart. Killed both, verified
  a clean process table, and relaunched exactly one.

## 2026-09-01 (late) — bringing the documents back into agreement with the evidence

The correction above was made in the paper first. Three other artefacts were still asserting
things the evidence no longer supports, which is exactly how a project ends up contradicting
itself in front of a referee. All three are now consistent.

### `project_explanation_bangla.docx` — rebuilt, and it had two errors, not one

The obvious one: finding 5 in §7.1 quoted the withdrawn CT18 figure ("8.45× against 6.63×") as
a clean win. It now gives the benchmark-wide result, including the 1.63× figure for the added
regions, and explains why the earlier claim was withdrawn.

The one I nearly missed: contribution #2 still described positional detrending as *"the step
without which the whole idea does not work."* That is the precise opposite of what the
118-genome ablation found, and it had survived because the earlier fix touched only the
methods narrative in §4.2, not the contributions list in §8. It now appears where it belongs —
as a clearly labelled negative result, with the reason it is off by default.

Two limitations were also added: the thinner enrichment of the regions SAGE-GI adds, and the
13% protocol coverage.

Verified after rebuilding: 129 paragraphs, 28 headings, 19,032 characters, all 134 Bengali runs
carrying a complex-script font slot (`w:cs`), and both stale strings confirmed absent.

### `ASSESSMENT.md` and the published page

The acceptance estimate moves from ≈40% to **≈55%**. The reviewer reports were rewritten
against the current manuscript rather than the one they were drafted for — leaving them as they
were would have made the document read as if TreasureIsland and the third model were still
missing.

Worth noting how Reviewer 2's objection behaves once the requested analysis exists: it does not
disappear, it sharpens. Having the benchmark-wide numbers turns "extend the biology" into
"reframe the paper around the representation result, which is unambiguous, rather than the F1
gain, which is qualified." That is a fair reading of our own evidence and deserves a decision
before submission rather than during review.

The estimate rose less than the completed work alone would suggest, and the reason should be
plain: three near-certain objections were removed, but the biology result weakened the
substantive claim. A careful referee will see a well-executed study with a more modest advance
than the abstract first suggested.

### In flight

`finalise.sh` running as a single instance: detector run at 28/118 genomes, then case studies,
synthetic insertions, TreasureIsland, biology, tables, figures. NT-v2 embeddings drop into
pre-wired configurations when they arrive.

### The abstract and Key Points were still selling the component that failed

Correcting the body of the paper was not enough. Two of the most-read parts of the manuscript
had not been touched since before the ablation, and both still advertised positional
detrending as a core contribution:

* The **abstract** named it as one of "three components [that] do the work", crediting it with
  removing "the replichore-scale variation that otherwise makes the detector fire on the origin
  of replication". The paper's own Methods, Results and Discussion all report that this does not
  happen and that the component costs 6.9 F1.
* **Key Point 1** said the embedding signal appears "only after the host's own replichore-scale
  structure is removed from them" — the same refuted claim, in the five bullets an editor reads
  first.

A referee who read the abstract and then Section 3.4 would have found the paper contradicting
itself on its own headline method. Both are rewritten. The abstract now leads with the
controlled DNABERT-S/DNABERT-2 comparison, which is the cleanest evidence in the study, and
states both qualifying findings — detrending, and the thinner enrichment of the added regions.
Key Point 1 now carries the controlled comparison; Key Point 5 carries the caveat.

This is worth recording as a process lesson: when a result is withdrawn, the abstract, the key
points, the contribution list and the plain-language summary are all separate copies of the
claim, and every one of them has to be chased down. Three of the four were stale here.

### Word count brought into compliance

Checked the journal's actual requirement rather than working from memory: Briefings in
Bioinformatics specifies **2,000–5,000 words** for a Problem Solving Protocol, and states no
abstract limit. The body had grown to 5,171 — over.

Cut to **4,932** without losing a single fact, by removing genuine duplication rather than
content. The detrending failure was being explained three times: in full in Methods, again in
Results, and a third time in Discussion, with one sentence appearing verbatim in two of them.
Discussion 4.1 and 4.2 were also restating numbers the Results tables already give instead of
interpreting them. Methods was left untouched, since cutting there would trade reproducibility
for word count.

The abstract went from 376 to 296 words. Also macro-ised a hard-coded "97 of 118" in the
introduction, which would otherwise have gone silently stale when the final run finishes.

The remaining 4,932 leaves less headroom than the 4,700 I would prefer, since revisions add
text. If a referee asks for material, Methods 2.1 and 2.11 are the passages to move to
supplementary.

## 2026-09-01 (evening) — Nucleotide Transformer arrives, and a cache bug worth the delay

### The embeddings validated cleanly

`NT-v2-50M_t1000.zip` delivered: 120 `.npz` files, 179 MB uncompressed. Validated before letting
the pipeline near it — accession set identical to DNABERT-S with nothing missing or extra, every
file loading through `sagegi.store.load_tiles`, all values finite, tile length 1,000 throughout,
96 stored principal components, and per-genome tile counts and genome lengths matching DNABERT-S
exactly on every spot check. No problems.

### A bug that had silently disabled cache invalidation

Wiring NT-v2 in exposed something worse than the integration itself. The detector cache is keyed
on accession plus a SHA-1 digest of the configuration file, and the header comment promises that
"editing a configuration invalidates it rather than silently reusing stale predictions."

It did not. The stem `NC_000913.2__4b5cd3ca35` was turned into a filename with
`Path.with_suffix(".rows.csv")` — and `with_suffix` replaces everything after the *last* dot, so
it read `.2__4b5cd3ca35` as the extension and threw it away. Every cache entry was landing at
`NC_000913.rows.csv`:

* the configuration digest was gone, so **editing `methods.json` invalidated nothing**;
* the accession version was gone, so two accessions differing only by version would have
  collided.

This is the kind of defect that does not announce itself. Nothing crashes; results are simply
computed under one configuration and reused under another.

**Whether it had actually corrupted anything was a question worth answering rather than
assuming.** I held out two genomes, deleted their cache, recomputed from scratch and compared
against the cached values across all seventeen shared configurations. Maximum absolute
difference in F1, precision and recall: **exactly zero**, on both genomes. The cached results
were computed under the current configuration after all — the only edit that had happened since
was the *addition* of the two NT entries, which leaves the other seventeen specifications
untouched. So no published number was ever wrong; the bug was a loaded gun that had not yet
gone off.

Fixed with an explicit `stem_path()` helper that concatenates rather than substitutes, and the
56 existing entries were renamed to the correct scheme instead of being discarded.

### Two further pipeline corrections

* **Cache completeness.** A genome cached before an embedding archive arrived contains rows only
  for the configurations that were runnable at the time. Reusing such an entry afterwards would
  have averaged NT-v2 over 56 fewer genomes than DNABERT-S and quietly biased the model
  comparison — the exact error the paper's central claim depends on not making. `one_genome` now
  compares the cache against what is runnable *now*, using a sidecar that records what was
  attempted so a configuration that failed at run time is not retried forever.
* **Incremental recompute.** The first version of that check re-ran all nineteen configurations
  to add two, turning a forty-minute job into four hours. Since configurations are independent
  given the same inputs, it now computes only the missing ones and appends.

* The representation analysis (`08`) was missing from `finalise.sh` and has been added; without
  it the NT-v2 AUC macros would have stayed at `TBD` and the third-model claim would have had no
  numbers behind it.

### Status

Pipeline running: detector run adding NT to the 56 cached genomes and computing the remaining 62
in full, then representation analysis over three models, case studies, synthetic insertions,
TreasureIsland, biology, tables and figures.

## 2026-09-01 (night) — the first fully correct run, and it breaks the fusion claim

### What the third model says

Nucleotide Transformer v2 answers the reviewer's question decisively, and in our favour.
All 118 chromosomes, 19 configurations, uniform coverage (2,242 rows, every configuration on
every genome — the cache-completeness fix doing exactly its job):

| configuration | F1 |
|---|---|
| DNABERT-S-only | **60.32** |
| SAGE-GI (fusion) | 59.63 |
| Comp-only+center | 56.36 |
| SSG-LUGIA-F (published baseline) | 55.17 |
| SAGE-GI(D2) | 51.54 |
| DNABERT-2-only | 48.45 |
| **NT-v2-only** | **43.18** |

Paired over genomes: DNABERT-S-only beats NT-v2-only by **+17.14 F1** (p = 1.1×10⁻¹⁵, ahead on
101/118) and DNABERT-2-only by +11.87 (p = 6.8×10⁻¹³, 97/118). Substituting either
general-purpose model into the full pipeline costs about 8 F1.

The claim is now much stronger than it was. **Two independent general-purpose genomic
foundation models both fall below hand-crafted composition; the species-aware one beats it.**
That is no longer a statement about DNABERT-S versus its sibling — it is a statement about what
kind of pre-training transfers.

### The number that breaks a claim

The same run says something the paper does not want to hear. SAGE-GI and DNABERT-S-only differ
in exactly one setting, `use_composition`; everything else — window geometry, dimensionality,
contamination, refinement — is identical. Adding the compositional channel to the embedding is
worth:

| metric | fusion | embedding only | delta | p | fusion wins |
|---|---|---|---|---|---|
| precision | 58.42 | 58.62 | −0.20 | 0.135 | 33/118 |
| recall | 72.51 | 73.13 | −0.62 | 0.492 | 21/118 |
| F1 | 59.63 | 60.32 | −0.69 | 0.181 | 34/118 |
| boundary error | 17.9 kb | 18.1 kb | −0.2 kb | 0.829 | 24/111 |

Nothing. On every metric the difference is insignificant and the point estimate slightly
favours the *simpler* method. **The frozen species-aware embedding, alone, beats the published
state of the art by +5.15 F1 (p = 8.6×10⁻⁵, ahead on 76/118) — and the compositional channel
adds no measurable value on top of it.**

This contradicts the paper in five places: the abstract ("the channels prove complementary"),
Key Point 3, the introduction contribution ("fusing them therefore beats either alone"),
Results §3.3 ("neither embedding is competitive on its own" — it is in fact the *best* single
channel), and the Discussion's complementarity argument.

### Why the old numbers hid this

The previous `numbers.tex` gave `\sgF = 52.72` and `\embSF = 49.81`. The current run gives
`SAGE-GI/detrend=250kb = 52.72` — exactly the old headline value. The old macros were generated
while detrending was still ON by default. When the 118-genome ablation showed detrending was
harmful and the default was flipped to 0, the *code* changed but the numbers were never
regenerated, so every quoted figure in the manuscript has been describing a configuration the
paper itself now disavows.

Two lessons, both mine to own. Regenerating derived numbers is part of changing a default, not a
follow-up task. And the cache-filename bug meant the config digest never invalidated anything,
so nothing forced the recomputation that would have surfaced this months earlier.

### Not yet decided

Whether to keep SAGE-GI as the fusion and report that composition adds nothing, or to redefine
the method as the embedding alone — simpler, and with the better point estimate. Holding that
until the case studies and the synthetic benchmark are in, because the fusion may still earn its
place where composition-based detection is hardest. That is the evidence that should decide it,
not a preference for the story already written.

## 2026-09-01 (night, later) — three process errors and a design correction

### A misdiagnosis I should record

I reported that `04_run_sage_gi.py` had "exited cleanly with exit code 0 having done nothing"
and that a silent `return` had swallowed the failure. That was wrong, and the way it was wrong
is worth writing down.

Launching with `nohup ... &` returns immediately. The completion notification I received was the
*launcher shell* exiting, not the job. I read it as the job finishing, looked for the output
file, found none — because the job was still running and writes only at the end — and concluded
the script had failed silently. I then started a second copy, and both ran against the same
cache for about fifteen minutes before I noticed two parent processes with the same command
line.

The defensive change stays: a run that produces no results now raises instead of exiting 0, and
that is genuinely better. But the premise was false, and "the harness said the launcher
finished" is not evidence about the job. Verify the worker processes, not the launcher.

### The synthetic sweep was mis-scoped by a default

Measured rather than estimated: **1 condition per minute, 2,025 remaining — 34 hours.**

The cause is a default rather than a decision. `06_synthetic.py` sets
`donors = [a for a in sorted(ready) if a not in hosts]` — *every* genome with embeddings. The
sweep was therefore 2 hosts x 118 donors x 5 insert lengths x 2 seeds = 2,360 conditions, or
236 replicates per insert length. What the experiment needs is coverage of the donor-host
*compositional distance* axis; the replication was incidental.

Computed the actual tetranucleotide distance from every candidate donor to both hosts and
sampled 24 donors evenly along it (`configs/synthetic_donors.txt`, distances in
`results/donor_distances.csv`):

* distance span retained in full: tnf_euclidean **0.0192-0.0811**, |ΔGC| **0.063-0.237**
* **480 conditions instead of 2,360** — about 3 hours instead of 34
* 48 replicates per insert length, which is ample

Endpoints are kept, so the "compositionally similar donor" regime the paper leans on is still
sampled at its extreme.

### The sweep is now resumable

It previously accumulated every result in memory and wrote one CSV at the end, so an
interruption at hour twelve would have lost twelve hours — and killing it did lose the 335
conditions already done. It now caches one CSV per condition under `$SAGEGI_WORK/synth_cache`,
keyed by host, donor, length, seed and a digest of the configuration file. A restart costs only
the unfinished conditions, and partial results can be inspected while the run is in progress.

### Current state

`finalise.sh` continued past the killed synthetic stage into TreasureIsland, and is running
through biology, tables and figures. The operating-point check is running alongside. A queued
script waits for the main chain to finish before starting the reduced synthetic sweep, so the
two do not compete for cores, and regenerates tables and figures afterwards.

## 2026-09-01 (late) — the fusion claim is withdrawn

### The operating-point check settles it

Adding the eleven compositional features to the species-aware embedding, measured over all 118
chromosomes at three operating points (configurations differing in exactly one setting,
`use_composition`):

| operating point | fusion | embedding only | delta | p | fusion wins |
|---|---|---|---|---|---|
| default | 59.63 | 60.32 | −0.69 | 0.181 | 34/118 |
| precision-oriented | 57.94 | 58.01 | −0.08 | 0.558 | 39/118 |
| recall-oriented | 48.09 | **49.21** | **−1.12** | **0.0002** | 35/118 |

It never helps, and at the recall-oriented point it is significantly *harmful*. Precision,
recall and boundary error agree at the default point (all p > 0.13, every point estimate
favouring the simpler method), and the three case studies agree too: the fusion wins on CT18
and loses on *C. diphtheriae* and LESB58.

The check also produced a clean determinism result. `SAGE-GI-P` and `SAGE-GI-R` were rerun in a
separate cache namespace from identical configuration text: maximum absolute difference from
the main run, over 118 genomes, **0.000000000000**.

### What changed in the paper

SAGE-GI is now the embedding-only pipeline. The compositional block is reported as an ablation
that earns nothing. Seven prose passages and four table-layer entries were rewritten:

* the abstract no longer claims the channels are complementary;
* Key Point 3 changes from "fusing beats either alone" to the learned representation
  *replacing* feature engineering;
* the introduction's third contribution says plainly that we expected complementarity and did
  not find it;
* Results 3.3 no longer says "neither embedding is competitive on its own" --- the species-aware
  one is the best single channel in the study;
* Results 3.4 gains "A component that earns nothing", carrying the three-operating-point result;
* Methods 2.6 is renamed from "Fusion and anomaly detection" and describes the concatenation as
  a variant rather than the method;
* the Discussion says LESB58 --- which we had read as our best evidence for complementarity ---
  actually scores *lower* fused than embedding-only.

`SageConfig.use_composition` now defaults to `False`, so the library default matches the
published method. Every configuration file sets the flag explicitly, so no cached result is
affected; the 21 tests still pass.

The framing throughout is **parsimony, not superiority**. At the default and precision points
the two are statistically tied, and claiming the simpler one is *better* would overreach in the
opposite direction. The defensible statement is that the compositional block does not earn its
place, so it goes.

### Word count

The rewrite pushed the body to 5,265 against BiB's 5,000. Rather than delete reproducibility
detail, three passages moved to a new `paper/supplementary.tex` --- decimated covariance
fitting, the hyper-parameter sweep bookkeeping, and the compositional feature definitions ---
each replaced by a one-line pointer. With further consolidation of the now-repetitive fusion
argument the body is **4,981**. That is compliant but tight; if a referee asks for additions,
Methods 2.1 is the next candidate to move.

### Still to do

The biology analysis was computed with the fused configuration as "our method" and must be
re-run against the embedding-only one (`scripts/11_biology_benchmark.py` now takes
`--main-method` for exactly this). Then tables, figures, README, website, wiki and the Bangla
document need the new numbers.

## 2026-09-01 (late) — propagating the withdrawal through every artefact

Withdrawing the fusion claim from the paper was the easy part. It appeared in eleven other
places, and two of them were figures rather than prose:

* **The pipeline figure** (Figure 1, and the same image on the project site) showed
  "4 · Fuse channels" as a coloured stage of the method. It is now a greyed box, "4 ·
  Composition (off by default)", saying that over 118 genomes it earns nothing.
* **The pipeline caption** in Methods carried *both* withdrawn claims at once — detrending as
  part of host referencing in step 2, and concatenation as part of the method in step 4. A
  reader comparing the caption to Section 3.4 would have found the paper contradicting itself
  in the first figure.
* **README, project site and wiki** — tables, narrative and the "why fusion works" card, which
  is now "and a second negative result".
* **The Bangla document** — the pipeline table, finding 3 in §7.1, contribution 3 and one
  limitation. Rebuilt and verified: 132 paragraphs, 28 headings, all 139 Bengali runs carrying
  a complex-script font, the new numbers present and both stale claims confirmed absent.
* **`ASSESSMENT.md` and the published page** — evidence table, a new section on the withdrawal,
  revised reviewer expectations, and the odds.

Also pointed the CT18 enrichment macros at whichever configuration `MAIN` names, so the
biology numbers cannot silently keep describing the fused method after the rename.

### Acceptance estimate: ≈55% → ≈60%

The third model is what moved it. Nucleotide Transformer v2 reaches AUC 0.547 — barely above
chance — against 0.759 for DNABERT-S (+0.212, p = 7.6×10⁻¹⁹, ahead on 106/118). Two
general-purpose genomic foundation models from different families both fall below hand-crafted
k-mer statistics, and only the species-aware one clears them. That converts the central claim
from a statement about one checkpoint into a statement about a training objective, which is the
sentence worth citing.

Set against it, the paper now reports three negative results (detrending, the fusion, the
thinner enrichment of added regions). Referees will respect that; an editor skimming for a
headline might read it as more cautionary than contributory. The abstract has been ordered to
lead with what works and place the qualifications after, and it should stay that way.

### Everything is chained to finish unattended

Three scripts, each waiting on the one before it, all verified running rather than assumed
running:

1. `finalise.sh` — TreasureIsland, then biology, tables, figures.
2. `after_final.sh` — waits, then the reduced synthetic sweep (480 conditions), then tables and
   figures.
3. `final_pass.sh` — waits, re-runs both enrichment analyses against the embedding-only method,
   regenerates tables and figures, checks the word count and builds the PDF.

The biology re-run matters: `finalise.sh` computed enrichment with the fused configuration as
"our method", which is no longer what SAGE-GI means.

## 2026-09-02 — TreasureIsland, and a benchmark property that explains it

### The comparison that looked like a defeat

TreasureIsland scores **F1 75.15** under our protocol, against 60.32 for SAGE-GI. Three
readings, two of them wrong, and the wrong ones are worth recording:

1. *"It beats us."* Incomplete.
2. *"It only over-calls."* **Wrong.** It calls 27% of each chromosome against a curated island
   content of 2.4%, which looked like the whole story — until the check: at a comparable
   called fraction, SSG-LUGIA-R reaches precision 40.1 and TreasureIsland 72.1. Over-calling
   does not explain a 32-point precision gap. This explanation was self-serving and it did not
   survive contact with the data.
3. The actual reason, in two independent parts.

### Part one: it is not an unsupervised method

The package ships `svm_embedding_*` files. Its *representation* is unsupervised, its
*classifier* is a supervised SVM shipped pre-trained. Our own runner's docstring described it
as "also unsupervised" — my error, and it would have gone into the paper.

It was trained on the Benbow compilation, whose principal component is *"104 genomes with
1,845 GEIs and 3,266 non-GEIs"* — the revised version of the benchmark we evaluate on. So a
large part of our test set is related to its training set.

The package exposes a per-genome similarity check, which turns that from an argument into a
measurement:

| | n | P | R | F1 |
|---|---|---|---|---|
| genomes it recognises | 107 | 79.5 | 92.4 | **82.9** |
| genomes it flags as unfamiliar | 11 | 0.0 | 0.0 | **0.0** |

On the 11 it returns *no predictions at all* — verified against `n_pred`, so this is the tool
declining rather than our runner failing. Its headline is the average of 82.9 on familiar
genomes and nothing on unfamiliar ones. SAGE-GI scores 59.7 and 66.3 on the same two subsets,
having no training distribution to be inside or outside of.

### Part two: the protocol cannot see a third of its predictions

`scripts/12_protocol_sensitivity.py` recomputes precision with every base outside a curated
island counted as a negative — a lower bound, since the unscored majority certainly contains
uncurated islands. The gap between that and the official figure is what the benchmark cannot
see:

| method | called % | P official | P strict | gap |
|---|---|---|---|---|
| TreasureIsland | 29.8 | 79.5 | **8.5** | **71.0** |
| SAGE-GI | 11.0 | 59.1 | 13.9 | 45.2 |
| SSG-LUGIA-F | 10.2 | 55.8 | 13.8 | 42.0 |
| Comp-only+center | 10.4 | 55.6 | 14.2 | 41.5 |

Every method loses ~42–46 points, which is the baseline cost of a protocol that grades 13% of
a chromosome. TreasureIsland loses 71. This is independent of the training-overlap argument,
so the two do not stand or fall together.

Worth stating against ourselves: SAGE-GI's gap (45.2) is slightly *larger* than the baseline's
(42.0), because we call slightly more sequence. The paper says so.

Both findings are now Results subsections, with every number generated by the pipeline.

### The biology reversed again, this time in our favour

Re-running both enrichment analyses against SAGE-GI *as now defined* (embedding-only) changes
the picture substantially:

| | fused (previous) | embedding-only (now) |
|---|---|---|
| SAGE-GI mobility OR | 4.55 — **below** the baseline | **5.95 — above it** |
| regions it adds | 1.63 | **2.96** |
| baseline (SSG-LUGIA-F) | 5.25 | 5.25 |
| curated islands | 6.14 | 6.14 |

The limitation written two rounds ago — "its pooled enrichment lands below the more
conservative SSG-LUGIA-F" — is no longer true. Dropping the compositional channel improved
biological coherence as well as F1. The qualification that survives is narrower: the regions
SAGE-GI *adds* are about half as enriched as the ones both methods agree on, and the baseline
remains ahead on **virulence** markers (2.36 against 1.85). Both are now stated.

This is the third time the biology conclusion has moved. The reason is that it was being
computed against whichever configuration was called "SAGE-GI" at the time, and the definition
changed. `scripts/11` now takes `--main-method` and the table macros follow `MAIN`, so this
particular failure cannot recur silently.

### Audit

`C:/sagegi_work/audit.py` checks what has actually gone wrong here before: macros resolving to
TBD, macros used in prose but never generated, withdrawn claims resurfacing in any of the six
documents, headline numbers drifting between paper and README, library defaults disagreeing
with the published method, non-uniform genome coverage, and the test suite. It reports **no
problems**: 21 tests pass, all 19 configurations cover 118 genomes uniformly, and the only
unresolved macros are the twelve synthetic ones while that sweep runs.

## 2026-09-02 — journal formatting, intervals, and an overclaim caught by checking my own claim

### A statistical claim I made, then had to test

Adding confidence intervals meant also saying something about multiple comparisons, so
Methods now states that p-values are Holm-corrected across the reported family and that every
significant comparison survives. That is an assertion about our own numbers, so
`C:/sagegi_work/check_holm.py` parses every p-value macro out of `numbers.tex`, applies Holm,
and reports whether any comparison changes verdict. **All 14 significant comparisons survive
across a family of 18.** The claim stands, and it is now checkable rather than asserted.

### The check found an overclaim

Sorting the p-values put `pstartAssign` at **0.75**, next to `startAssignDelta` at **−0.77**.
The paper was saying, in three places, that projecting each window to its centre rather than
its leading bases "recovers \\startAssignDelta{} F1" — which renders as "recovers −0.77 F1",
with the sign reading backwards, and quotes a p-value of 0.75 as though it supported the
claim.

Measured properly, over 502 matched islands:

| quantity | effect | p |
|---|---|---|
| F1 | +0.77 | 0.83 |
| median absolute boundary error | −162 bp | 0.60 |
| **signed start error** | **+1,508 bp toward zero** | **7.5×10⁻⁵** |

So the *mechanistic* claim is real and strongly supported — the leading-edge convention does
displace calls leftward, and centring removes about 1.5 kb of it — but the *accuracy* claim is
not. F1 and MABE do not move measurably.

Rewritten in all three places to say exactly that: it is a correctness fix, not an accuracy
gain, and "a metric that cannot see a systematic 1.5 kb displacement is a metric worth
supplementing" is now the point being made rather than a spurious F1 improvement. The Key
Point that quoted the F1 figure has been dropped.

Worth noting how this surfaced: not from reading the paper, but from sorting our own p-values
to check a *different* claim. The overclaim had survived several read-throughs because the
macro name looked plausible and the number was never displayed next to its p-value.

### Confidence intervals

Twelve p-values and no intervals was a standing invitation for a methods referee. Bootstrap
95% intervals (10,000 resamples over genomes) now accompany the comparisons the argument rests
on. Two of them do more work than any p-value:

| comparison | 95% CI |
|---|---|
| SAGE-GI − baseline (F1) | [2.3, 7.9] |
| **adding composition (F1)** | **[−2.3, 1.0]** — contains zero |
| DNABERT-S − DNABERT-2 (AUC) | [0.083, 0.127] |
| DNABERT-S − NT-v2 (AUC) | [0.182, 0.243] |
| **DNABERT-S − composition (AUC)** | **[−0.019, 0.047]** — contains zero |

An interval containing zero *and* excluding any effect of practical size is the defensible way
to assert equivalence. "We failed to reject the null" is not, and that is what the paper was
previously relying on for both of those claims.

### Journal formatting

Checked against the Briefings in Bioinformatics manuscript-preparation page rather than from
memory. Fixed:

* **Reference style** — the journal specifies PubMed/Vancouver ("Attwood T.K. ... Brief
  Bioinform 2000;1:45–59"). We were using `oup-plain`, which renders full author lists in a
  different order. Now `vancouver.bst`, which matches.
* **Key Points** — "should consist of 3-5 brief sentences". Ours were five long multi-clause
  paragraphs. Now five brief sentences.
* **Author description** — "about 30 words" per author, required, and absent. Added.
* **ORCID** — required for the submitting author. Fields added with placeholder iDs, flagged
  in a comment and by the audit.
* **Supplementary statement** — said S1–S3 after S4 was added.

Verified as already compliant: six keywords (limit six), figures at 600 dpi (requirement 600
for line art, 300 greyscale), Data availability, Funding and Conflict of interest sections all
present.

The audit now checks all of this, so a future edit that breaks a journal requirement fails
loudly instead of reaching a desk editor.

## 2026-09-02 — the synthetic benchmark was measuring a basis mismatch

The sweep finished (480/480, no failures) and the numbers made no sense:

| method | 5 kb | 10 kb | 20 kb | 40 kb | 80 kb |
|---|---|---|---|---|---|
| hand-crafted composition | 29.5 | 46.4 | 86.8 | 92.4 | 90.5 |
| **SAGE-GI (embedding)** | **9.2** | **11.8** | **15.5** | **21.4** | **33.7** |

On the real benchmark the embedding beats composition and the compositional block earns
nothing. On synthetic insertions the embedding recovered 18% of inserts against composition's
69%. Both cannot be true, and an 80 kb block from a distant genome should be trivially
detectable by *any* working method.

### The cause

Chimera tile embeddings are assembled by concatenating precomputed host and donor matrices.
The archive stores each chromosome **projected onto its own 96 principal directions**, and
does not retain the basis. So column *j* of the assembled matrix means "host PC *j*" for host
rows and "donor PC *j*" for donor rows --- different directions in the original 768-dimensional
space. The chimera's embedding is not an embedding of anything.

Diagnosed rather than assumed (`C:/sagegi_work/diagnose_synth.py`), using the most
compositionally distant donor in the sweep and the longest inserts:

| | AUC separating inserted/island tiles |
|---|---|
| 80 kb inserts in the assembled chimera | **0.273** |
| curated islands, single coherent basis | **0.784** |

An AUC of 0.27 is *below chance*: inserted tiles score as systematically **less** anomalous
than host tiles. PCA orders components by variance within each genome, so donor coordinates
have host-like magnitudes while meaning something else entirely --- the anomaly signal is not
weakened, it is destroyed.

### What this invalidates

Everything the synthetic benchmark said about the embedding channel. The compositional arm is
unaffected, because those features are recomputed on the real chimeric sequence, so the
construction and the boundary metric are still validated --- but the SAGE-GI-versus-composition
comparison is withdrawn, and with it the paper's claim that boundary conclusions rest
primarily on exact-coordinate data. They now rest on the curated benchmark, with the
uncertainty that implies. Figure 5 shows the compositional channel only, with a caption saying
why.

Repeating the experiment needs full-dimensional embeddings for the two hosts and 24 donors ---
about one GPU-hour --- and is the first thing to run next.

### Why it survived so long

Nothing in the code forbade it, and the result looked like a weak finding rather than an
invalid one. `assemble_tiles` now raises unless given full-dimensional embeddings, with the
measured AUC figures in the error path so the next person meets the evidence rather than a
rule, and `06_synthetic.py` passes the flag through. Test 22 covers it.

That is the general lesson from this project, arriving for the fourth time: a number that
merely looks disappointing gets explained; a number that looks impossible gets investigated.
The embedding arm had been quietly reported as a modest result for weeks.

## 2026-09-12 � full-dimensional re-embedding, synthetic repair, IEEE version

### Kaggle re-run and swap

The user ran the updated `01_extract_embeddings_kaggle.ipynb` (RAW_SUBSET = 31:
2 synthetic hosts + 24 donors + 5 previous). Both zips validated before touching
the pipeline: 120/120 `.npz` per model, all 26 needed genomes carrying finite
full-768-d `emb`, accession sets identical, and spot-checked `pcs` bitwise
identical to the previous archive (deterministic inference). Old archives moved to
`C:/sagegi_work/emb_backup_oldpc/`; new ones live in `emb/`.

### Two bugs the re-run exposed

1. **`load_tiles` preferred `pcs` even when `emb` was present.** The new archives
   carry both keys for RAW_SUBSET genomes, so the synthetic sweep silently kept
   getting PCA tiles and `assemble_tiles` kept raising. Fix: `prefer_full`
   parameter on `load_tiles` (`sagegi/store.py`), default False so every existing
   caller � and every published number � is untouched; `06_synthetic.py` is the
   one caller that passes True. Verified first on a single chimera before
   launching the sweep.
2. **Stale cache reuse trap.** `synth_cache` entries are keyed by config digest,
   which did not change, so the old invalid rows would have been reused
   silently. Old cache moved aside to `synth_cache_oldbasis/`; sweep recomputed
   fresh (480/480, zero failures, ~3h).

### Diagnosis, then and now

Replicating `diagnose_synth.py` exactly (AC_000091.1 + NC_009080.1, 80 kb):
chimera AUC **0.27 ? 0.991**. The basis mismatch is gone; the chimera embedding
is a valid embedding of the chimera.

### Sweep results (2 hosts x 24 donors x 5 lengths x 2 seeds)

Detection rate by insert length (Comp / DNABERT-2-only / DNABERT-S-only /
SAGE-GI-fused): 5 kb: 29.5 / 24.0 / 28.3 / 29.5; 10 kb: 46.4 / 40.1 / 42.4 /
48.4; 20 kb: 86.8 / 86.8 / 92.5 / 92.4; 40 kb: 92.4 / 91.8 / 93.8 / 93.9;
80 kb: 90.5 / 82.8 / 87.3 / 87.9. Overall: Comp 69.1%, DNABERT-S 68.9%,
DNABERT-2 65.1%. MABE: Comp 8.8 kb, DNABERT-S 6.5 kb, DNABERT-2 6.4 kb.
Comp-only reproduced bitwise (69.1/37.9), confirming determinism. The embedding
channel matches composition overall, exceeds it at 20-40 kb, and places edges
substantially better � the boundary claim now has exact-coordinate support.

### Table generator bug

`table4_ct18.tex` printed DNABERT-S-only twice instead of SAGE-GI: `MAIN` had
not been renamed with the method, so the table loop emitted the same method
under two names. Fixed in `91_make_tables.py` (`disp(m)` in the CT18 loop);
regenerated table now reads SAGE-GI 77.8 / SAGE-GI + composition 80.1,
consistent with the case-studies prose.

### Regeneration discipline

`91_make_tables.py` + `90_make_figures.py` re-ran. `numbers.tex` diff: only
syn* embedding macros changed (18.3?68.9 etc.) plus previously-unemitted edge*
macros; every other macro bitwise identical. Tables diff: only table4 (the fix
above). Fig5 now shows all four channels. Because the BiB paper is frozen at
its word limit, `paper/` keeps the Comp-only fig5 matching its caption; the
4-method fig5 lives in `paper-ieee/`. The figure script carries a comment
explaining this split. BiB still compiles clean (12 pp), word count 4,998,
audit NO PROBLEMS FOUND.

### IEEE two-column version (`paper-ieee/`)

Built for conference submission: IEEEtran conference class, same text/sections/
tables/figures/numbers/bib as the OUP version, adapted frontmatter (no ORCID,
journal metadata, Key Points), `IEEEtran.bst`, `\botrule`?`\bottomrule`,
`tablenotes`?footnotes. Compiles clean (11 pp, no undefined refs), Overleaf-ready.
After the synthetic repair it carries the real embedding results (rewritten �3.7,
fig5 caption, boundary pointer, limitations paragraph); the BiB version keeps the
withdrawal framing per the 5,000-word freeze.

## 2026-09-13 — full work check, BiB synthetic update, TCBB assessment

### Check: everything verified current

Audit clean (22 tests, 19 configs x 118 genomes uniform), both papers compile
(12 pp BiB, 11 pp IEEE), all IEEE macros resolve. Stale-claim grep over the repo:
wiki/README/docs already describe the repaired synthetic pipeline (a previous
session updated them, including `prefer_full=True`); no changes needed there. The
only stale text anywhere was the BiB paper's own withdrawal framing — updated
below, which is what "update everything" required.

### BiB paper now carries the real synthetic results too

The 5,000-word freeze is lifted by editing, not by exception: section 3.7
rewritten (withdrawal out, results in), fig5 caption replaced, boundary pointer
and limitations paragraph updated, 4-method fig5 regenerated into
`paper/figures`. Net effect: body **4,998 → 4,903** — still compliant with ~100
words of margin. Recompiled clean (12 pp), audit NO PROBLEMS FOUND. The
figure-script comment claiming paper/ keeps a Comp-only fig5 was removed as
outdated; both versions now show all four channels.

Deliberately unchanged: every non-synthetic number (regeneration diff proved only
syn* macros moved; Comp-only bitwise identical), all tables except the
already-fixed table4, the abstract, Key Points, and conclusions.

### ASSESSMENT_IEEE.md

New venue assessment for IEEE/ACM TCBB (the only top-tier IEEE journal in scope;
TBME/JBHI are clinical, TPAMI/TKDE out of field). Calibrated eventual acceptance
**approx 50-55% via major revision** — below BiB's 60-65% for three concrete
reasons: no venue-fit advantage, a transactions-level novelty bar ("engineering
around a pretrained model" attack), and the 2018 benchmark reading as expectation
rather than bonus. Three simulated reviews included (methods-favourable,
domain-sceptical, ML-positive); all three converge on the 2018 set +
fine-tuning headroom. Priority list: 2018 benchmark first (2-4 days, same GPU
sessions as fine-tuning), then detector/hyperparameter ablations (CPU-only), MAG
+ island-free smoke test, and transactions presentation items (bios, peak memory,
seed locations).
