# Pipeline Reference

Numbered scripts under `scripts/`. Each is independently runnable and **resumable** — every
stage caches to `$SAGEGI_WORK` and skips work already on disk.

---

## `00_prepare_benchmarks.py`

Downloads the IslandPick additional files from the publisher and parses the Excel into tidy
CSVs.

```bash
python scripts/00_prepare_benchmarks.py
```

**Writes** `data/benchmarks/`:

| File | Contents |
|---|---|
| `islandpick_positive.csv` | 771 curated islands: accession, organism, start, end, length |
| `islandpick_negative.csv` | 3,700 curated backbone regions, same columns |
| `islandpick_baselines.csv` | per-genome TP/FP/TN/FN for six published tools |
| `islandpick_accessions.json` | the 118 accessions |

## `01_download_genomes.py`

Fetches the 118 benchmark chromosomes plus three case-study genomes from NCBI Entrez, three
requests at a time under the 3 req/s limit, with retries.

```bash
python scripts/01_download_genomes.py
```

**Writes** `$SAGEGI_WORK/genomes/<accession>.fna` (~500 MB total).

## `02_run_ssg_lugia.py`

The compositional baseline: a faithful re-implementation of SSG-LUGIA in its P, F and R
variants, each also evaluated with centred window projection so the boundary bias can be
isolated.

```bash
python scripts/02_run_ssg_lugia.py --workers 3
```

| Flag | Default | Meaning |
|---|---|---|
| `--workers` | 3 | parallel processes (each pins BLAS to one thread) |
| `--limit` | 0 | process only the first *N* genomes |
| `--redo` | off | recompute cached features and scores |

**Writes** `$SAGEGI_WORK/comp/<acc>.npz` (feature matrix, window starts, per-variant scores)
and `results/ssg_lugia_{per_genome,boundaries,intervals}.csv`.

> This runs the exact published algorithm (`fit_stride = 1`), so it is slower than SAGE-GI by
> design. Budget roughly 2 h for all 118 chromosomes on four cores.

## `03_extract_embeddings.py`

Tile embedding extraction. Identical code path on CPU and GPU.

```bash
python scripts/03_extract_embeddings.py --model DNABERT-S --tile 1000
python scripts/03_extract_embeddings.py --model DNABERT-2 --tile 1000 --batch 64
```

| Flag | Default | Meaning |
|---|---|---|
| `--model` | `DNABERT-S` | `DNABERT-S` or `DNABERT-2` |
| `--tile` | 1000 | tile length in bp |
| `--batch` | 16 | tiles per forward pass (raise to 64 on a GPU) |
| `--pca-k` | 96 | principal components stored per genome |
| `--keep-raw` | – | accessions for which the full 768-d matrix is also kept |
| `--accessions` | all | restrict to specific accessions |

**Writes** `$SAGEGI_WORK/emb/<model>_t<tile>/<acc>.npz`.

Genomes are processed **smallest first** (after the case studies), so a partial run is still
useful. See **[[Reproducing the Paper]]** for the Colab route.

## `04_run_sage_gi.py`

Runs SAGE-GI and every ablation over all genomes that have both cached compositional features
and embeddings.

```bash
python scripts/04_run_sage_gi.py --configs configs/methods.json --workers 3
```

Configurations are named entries in a JSON file; each is a partial `SageConfig` plus an
`embedding` key naming the embedding directory (or `null` for a composition-only run). See
`configs/methods.json` for the 17 shipped configurations.

**Writes** `results/sage_gi_{per_genome,boundaries,intervals}.csv`.

## `05_tune_cv.py`

Hyper-parameter sweep with grouped 5-fold cross-validation.

```bash
python scripts/05_tune_cv.py --channel sage --workers 3
python scripts/05_tune_cv.py --channel comp
python scripts/05_tune_cv.py --channel emb
```

The expensive robust covariance fit depends only on the feature matrix and the two
contamination levels, so it is done once per `(genome, c1, c2)` and the post-processing
parameters are swept for free. The complete genome x configuration table is built once;
cross-validation then only *reads* it, so it costs nothing extra.

**Writes** `results/sweep_<tag>.csv` and `results/tuning_<tag>.json` (in-sample and held-out
performance, plus the configuration each fold chose).

## `06_synthetic.py`

The controlled synthetic-insertion benchmark.

```bash
python scripts/06_synthetic.py --workers 3
```

Splices tile-aligned donor blocks into host chromosomes at known coordinates, sweeping insert
length (5–80 kb) and donor–host compositional distance. Because inserts are tile-aligned, the
chimera embeddings are assembled by concatenating precomputed tiles — exactly, with no GPU
work per condition — but host and donor tiles must share one basis, so this step loads the
full-dimensional embeddings (`load_tiles(..., prefer_full=True)`); per-genome PCA archives
are refused by `assemble_tiles`. Compositional features *are* recomputed on the real
chimeric sequence.

**Writes** `results/synthetic_insertions.csv`.

## `90_make_figures.py` and `91_make_tables.py`

```bash
python scripts/90_make_figures.py
python scripts/90_make_figures.py --only fig3 fig4
python scripts/91_make_tables.py
```

Figures go to `figures/` as 600 dpi PDF plus PNG. Tables go to `paper/tables/*.tex`, and every
number quoted in the manuscript is written to `paper/numbers.tex` as a LaTeX macro — so the
prose can never drift from the data. Missing inputs are skipped with a note rather than
crashing, so both are safe to run while the pipeline is still filling in.

---

**Next:** **[[API Reference]]** · **[[Reproducing the Paper]]**
