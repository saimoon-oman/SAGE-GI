# Reproducing the Paper

Every number, figure and table in the manuscript is produced by the scripts below. Nothing is
hand-entered: `scripts/91_make_tables.py` writes `paper/numbers.tex`, and the LaTeX quotes
those macros, so the prose cannot drift from the data.

Total cost: about **2–3 h of CPU** plus **one Colab session** for the GPU step.

---

## Step 0 — install and verify

```bash
git clone https://github.com/saimoon-oman/SAGE-GI.git
cd SAGE-GI
pip install -r requirements.txt
python tests/test_sagegi.py        # expect 20/20
```

## Step 1 — benchmark tables (~1 min)

```bash
python scripts/00_prepare_benchmarks.py
```

Expected output:

```
positive :   771 islands   over 118 genomes  (12,376,375 bp)
negative :  3700 regions   over 118 genomes  (50,638,496 bp)
baselines:   708 rows      ['alien_hunter', 'centroid', 'islandpath_dimob',
                            'islandpath_dinuc', 'pai_ida', 'sigi_hmm']
```

If those three lines match, the benchmark is byte-for-byte what the paper used.

## Step 2 — genomes (~20 min, network-bound)

```bash
python scripts/01_download_genomes.py
```

121 FASTA files, about 500 MB, into `$SAGEGI_WORK/genomes/`. Safe to interrupt and re-run.

## Step 3 — compositional baseline (~2 h on four cores)

```bash
python scripts/02_run_ssg_lugia.py --workers 3
```

This runs the published SSG-LUGIA algorithm exactly, with no decimation of the covariance
fit, so it is deliberately the slowest stage. It caches per-genome features and scores, so a
re-run is nearly free.

**Sanity check.** Averaged over genomes, SSG-LUGIA-F should land near the published
P 55.1 / R 65.2 / F1 55.1 and SSG-LUGIA-R near 40.0 / 86.8 / 49.1. Exact agreement is not
expected — the published numbers come from a stochastic tuning procedure — but a large gap
means something is wrong.

## Step 4 — tile embeddings (the only GPU step)

### Recommended: free Colab

1. Open [`notebooks/01_extract_embeddings_gpu.ipynb`](https://github.com/saimoon-oman/SAGE-GI/blob/main/notebooks/01_extract_embeddings_gpu.ipynb)
   in Colab.
2. **Runtime → Change runtime type → T4 GPU.**
3. **Run all.**

The notebook downloads both DNABERT models and all 121 chromosomes itself, embeds every
chromosome, verifies that half precision does not perturb the embeddings (fp16-vs-fp32
cosine > 0.999), and zips the result to `Drive/MyDrive/SAGE-GI-embeddings/`. Roughly 20–40 min
per model. It is resumable, so a disconnect costs only the genome in flight.

4. Then open
   [`notebooks/02_extract_third_model.ipynb`](https://github.com/saimoon-oman/SAGE-GI/blob/main/notebooks/02_extract_third_model.ipynb)
   and run all. This produces the Nucleotide Transformer v2 embeddings, which are what turn the
   species-awareness result from a claim about one matched pair into a claim about training
   objectives.

   That notebook pins `transformers==4.46.3`, `tokenizers<0.21` and the model revision
   `81b29e5786726d891dbf929404ef20adca5b36f1`. The pins are not optional: NT-v2 loads through
   `trust_remote_code`, and its bundled `modeling_esm.py` imports two symbols that later
   `transformers` releases removed. Pinning the *model* revision as well means an upstream edit
   cannot silently change your embeddings.

Unzip all three archives into `$SAGEGI_WORK/emb/`:

```
$SAGEGI_WORK/emb/DNABERT-S_t1000/*.npz
$SAGEGI_WORK/emb/NT-v2-50M_t1000/*.npz
$SAGEGI_WORK/emb/DNABERT-2_t1000/*.npz
```

### Alternative: locally

```bash
python scripts/03_extract_embeddings.py --model DNABERT-S --tile 1000
python scripts/03_extract_embeddings.py --model DNABERT-2 --tile 1000
```

On a GPU this is the same ~30 min per model. On a laptop CPU expect **~70 h per model** at
about 2,000 bp/s — which is exactly why the notebook exists.

> **Do not try to save time with int8 quantisation.** We measured it: dynamic int8
> quantisation of the linear layers gives ~1.5x throughput but the mean-pooled embeddings
> collapse — cosine similarity to the fp32 embedding drops to 0.08–0.39. The embeddings are
> destroyed. fp16 on GPU is fine (cosine > 0.999); int8 is not.

## Step 5 — SAGE-GI and the ablations (~1 h)

```bash
python scripts/04_run_sage_gi.py --configs configs/methods.json --workers 3
```

17 configurations: the full model, the two embedding channels alone, the compositional
channel alone, and component ablations (with detrending, without refinement, leading-edge
projection, alternative subspace dimension).

```bash
python scripts/08_representation_analysis.py --workers 3
```

Step 5b measures the representations directly, without any detector: Mann-Whitney AUC of a
fixed anomaly statistic over island versus backbone tiles, swept across subspace dimension
and detrending bandwidth. This is what isolates the species-awareness effect from everything
the detector does.

## Step 6 — tuning and cross-validation

```bash
python scripts/05_tune_cv.py --channel sage --workers 3
python scripts/05_tune_cv.py --channel comp --workers 3
```

Produces both the in-sample numbers (the protocol SSG-LUGIA used, for comparability) and
grouped 5-fold cross-validated numbers, in which the configuration is chosen on four fifths
of the chromosomes and scored on the held-out fifth.

## Step 7 — synthetic insertions

```bash
python scripts/06_synthetic.py --workers 3
```

## Step 8 — figures, tables and the PDF

```bash
python scripts/90_make_figures.py
python scripts/91_make_tables.py

cd paper
pdflatex main && bibtex main && pdflatex main && pdflatex main
```

`91_make_tables.py` prints how many macros it resolved, e.g. `numbers.tex: 13/13 macros
resolved`. Any unresolved macro appears in the PDF as `TBD`, which makes a missing result
impossible to overlook.

### Overleaf

Upload the whole `paper/` folder and set `main.tex` as the main document. The OUP class
(`oup-authoring-template.cls`) and both `.bst` files are included, so nothing else is needed.

---

## Determinism

| Source of randomness | How it is controlled |
|---|---|
| `EllipticEnvelope` (MCD subsampling) | `random_state=3`, as in the reference implementation |
| PCA (`svd_solver="randomized"`) | `random_state=0` |
| Synthetic insert placement | explicit `numpy.random.default_rng(seed)`, seeds recorded in the output CSV |
| Cross-validation folds | `numpy.random.default_rng(0)` |

Re-running any stage reproduces the same numbers. Model inference is deterministic given a
fixed dtype; fp16 and fp32 differ in the last decimal, far below any reported precision (we
check this in the notebook).

---

**Next:** **[[Pipeline Reference]]** · **[[Troubleshooting]]**
