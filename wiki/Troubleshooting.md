# Troubleshooting

Things that actually went wrong while building this, and what fixed them.

---

## Model loading

### `RuntimeError: Tensor on device meta is not on the expected device cpu!`

You are on **transformers 5**. DNABERT-2/S ship custom modelling code written against the
4.x loader; transformers 5 initialises modules on the `meta` device first and the bundled
code creates real tensors during `__init__`.

```bash
pip install "transformers==4.46.3" "tokenizers<0.21"
```

### `ERROR: Failed building wheel for tokenizers`

There is no prebuilt `tokenizers` wheel for your Python version at the pinned version, so pip
tries to compile Rust. Either use a `tokenizers` release that has wheels for your Python
(3.13 needs ≥ 0.20.3), or use Python 3.11/3.12.

### `UserWarning: Unable to import Triton; defaulting MosaicBERT attention implementation to pytorch`

Harmless. At 1 kb tiles (~214 tokens) attention is a negligible share of the compute.

On a GPU where Triton *is* importable, the bundled `flash_attn_triton.py` may fail with newer
Triton releases. The provided notebook sidesteps this by forcing the PyTorch path:

```python
import sys
for name, mod in list(sys.modules.items()):
    if name.endswith("bert_layers") and hasattr(mod, "flash_attn_qkvpacked_func"):
        mod.flash_attn_qkvpacked_func = None
```

### `OMP: Error #15: Initializing libiomp5md.dll, but found libiomp5md.dll already initialized`

Anaconda's MKL and PyTorch each ship an OpenMP runtime. The standard workaround:

```bash
export KMP_DUPLICATE_LIB_OK=TRUE
```

---

## Performance

### Embedding extraction is impossibly slow

Check you are actually on a GPU:

```python
import torch; print(torch.cuda.is_available(), torch.__version__)
```

`2.x+cpu` means you installed the CPU build. On CPU, expect ~2,000 bp/s — about 70 h for the
whole benchmark per model. Use the Colab notebook instead.

### Can I quantise the model to speed up CPU inference?

**No.** We measured dynamic int8 quantisation of the linear layers: ~1.5x throughput, but
cosine similarity between the int8 and fp32 mean-pooled embeddings collapses to **0.08–0.39**.
The embeddings are destroyed and every downstream number becomes meaningless. fp16 on GPU is
safe (cosine > 0.999).

### The baseline sweep is taking hours

It is meant to. `02_run_ssg_lugia.py` runs the published algorithm exactly, fitting the
Minimum Covariance Determinant on every window. SAGE-GI itself uses `fit_stride=10`, which is
6.6x faster and agrees on 99.84% of window calls — but the baseline deliberately does not, so
that it *is* the published method.

If you only need the baseline approximately, `--workers 4` and editing `fit_stride` into the
model dicts will cut it to ~20 min.

### It is using too much memory

Each worker holds a `(n_windows, 256)` tetranucleotide matrix — about 170 MB for an 8 Mbp
chromosome. With `--workers 3` peak usage is roughly 3 GB. Reduce `--workers` on a small
machine.

---

## Results that look wrong

### F1 is 0 for some genomes

Expected, and not necessarily a bug. Some IslandPick chromosomes have only one or two curated
islands totalling a few kb, so a single missed call takes recall to zero. Look at the
per-genome table (`results/*_per_genome.csv`) rather than the mean before concluding anything.

### Predictions are systematically shifted left

That is the `assign="start"` projection, which reproduces the reference implementation's
behaviour: it credits each window to its *leading* `dw` bases while the window extends `w`
bases to the right, so islands drift left by about `w/2`. Use `assign="center"` (the SAGE-GI
default). The effect is quantified by `median_signed_start` in the boundary tables.

### Everything is called an island near one place in the genome

Check `detrend_bp` is 0 (the default). Positional detrending was measured to *hurt* across
the benchmark, and switching it on can concentrate calls in the wrong places. If calls still
cluster in one region, plot `res.scores` against `res.window_starts` — a single broad peak
usually means the chromosome has strong compositional structure of its own (a large prophage
or a recent duplication) rather than many islands.

### My genome has no islands at all

Check `res.scores`. If the whole track sits above zero, the contamination parameters are too
strict for this genome; try `SAGE-GI-R` settings (`contamination_model1=0.20`,
`contamination_model2=0.25`). Also check the chromosome is long enough — with
`min_island_len=10000` and a 10 kb window, anything under ~500 kb has very little to work
with.

---

## Data

### NCBI downloads time out

`01_download_genomes.py` retries five times with backoff and throttles to three concurrent
requests. If it still fails, you are likely rate-limited; wait and re-run — completed genomes
are skipped.

### The Springer additional files will not download

`00_prepare_benchmarks.py` sets a browser user-agent because the publisher's CDN rejects the
default. If it still fails, download Additional Files 2, 4 and 6 from the
[article page](https://doi.org/10.1186/1471-2105-9-329) by hand into
`$SAGEGI_WORK/raw/langille_AF{2,4,6}.xls` and re-run.

### `xlrd` refuses to open the file

`xlrd` ≥ 2.0 reads only legacy `.xls`, which is exactly what these files are. If you get a
format error, the download was truncated — delete it and re-run.

---

## LaTeX

### `File 'flushend.sty' not found`

The OUP class needs `sttools`. On MiKTeX:

```bash
miktex packages install sttools
initexmf --set-config-value "[MPM]AutoInstall=1"
```

On TeX Live: `tlmgr install sttools`.

### Numbers in the PDF read `TBD`

`paper/numbers.tex` has unresolved macros because the corresponding result table does not
exist yet. Run the missing pipeline stage, then `python scripts/91_make_tables.py`. The
script prints how many macros it resolved.

---

Still stuck? Open an issue at
<https://github.com/saimoon-oman/SAGE-GI/issues> with the command you ran, the full
traceback, and the output of `python -c "import torch, transformers, sklearn; print(torch.__version__, transformers.__version__, sklearn.__version__)"`.
