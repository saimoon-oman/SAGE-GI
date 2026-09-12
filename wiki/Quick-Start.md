# Quick Start

## Call it on one genome

```python
from sagegi.genome import read_fasta
from sagegi.embed import EmbedConfig, embed_genome
from sagegi.sage import SageConfig, run_sage_gi

seq = read_fasta("my_genome.fna")

# 1. embed 1 kb tiles once (the only slow step; uses a GPU if one is visible)
emb, _ = embed_genome(seq, EmbedConfig(model_path="zhihan1996/DNABERT-S", tile=1000))

# 2. detect
res = run_sage_gi(len(seq), SageConfig(), emb=emb, tile=1000)

for start, end in res.islands:
    print(f"{start+1}\t{end+1}\t{end-start+1} bp")
```

Coordinates returned by the library are **0-based inclusive**; add 1 for the 1-based
convention used by GenBank and by the printed output above.

## Choosing an operating point

Genomic island calling has an unavoidable precision/recall trade-off, so — following
SSG-LUGIA — three tunings are provided rather than one:

| Variant | Use when |
|---|---|
| `SAGE-GI-P` | you will follow up every call experimentally and false positives are costly |
| `SAGE-GI-F` | balanced; the sensible default |
| `SAGE-GI-R` | you are screening and would rather not miss anything |

```python
cfg = SageConfig(contamination_model1=0.075, contamination_model2=0.075)  # precision-oriented
cfg = SageConfig(contamination_model1=0.20,  contamination_model2=0.25)   # recall-oriented
```

## Composition only (no GPU, no model download)

```python
res = run_sage_gi(len(seq), SageConfig(use_embedding=False), 
                  comp_X=comp.X, comp_starts=comp.starts)
```

This is the SSG-LUGIA feature set inside SAGE-GI's post-processing. Useful as a fast
first pass, or when you have no GPU budget at all.

## Inspecting the anomaly track

`run_sage_gi` returns more than coordinates:

```python
res.scores          # per-window decision score; < 0 means alien
res.window_starts   # window start coordinate for each score
res.tile_signal     # tile-resolution anomaly track used for boundary refinement
res.coarse_islands  # intervals before CUSUM refinement
```

Plotting `res.scores` against `res.window_starts` is the fastest way to see whether a genome
has a clean signal or is compositionally noisy.

## Speed

| Setting | 5 Mbp chromosome |
|---|---|
| T4 GPU, fp16, batch 64 | well under a minute |
| Laptop CPU (4 cores), fp32 | ~45 minutes |
| Composition only | ~10 seconds |

If you have many genomes and no local GPU, run
[`notebooks/01_extract_embeddings_gpu.ipynb`](https://github.com/saimoon-oman/SAGE-GI/blob/main/notebooks/01_extract_embeddings_gpu.ipynb)
on free Colab — it is resumable and writes one compact `.npz` per genome.

**Next:** **[[Method Overview]]** · **[[Pipeline Reference]]**
