# API Reference

The library is `sagegi/`. Everything is plain NumPy / scikit-learn; nothing is stateful.

---

## `sagegi.genome`

| Function | Purpose |
|---|---|
| `read_fasta(path) -> str` | concatenated, upper-cased sequence of a single-record FASTA |
| `encode(seq) -> int8[]` | A,T,C,G to 0,1,2,3; every other symbol to −1 |
| `kmer_ids(code, k) -> int32[]` | base-4 id of the k-mer at each position, −1 if it spans an ambiguity |
| `window_starts(L, w, dw)` | window offsets, matching the reference implementation's loop exactly |

`BASES = "ATCG"` — the index order is inherited from the reference SSG-LUGIA code and must
not be changed, or the amino-acid mapping shifts.

## `sagegi.kmer`

| Function | Purpose |
|---|---|
| `window_kmer_counts(code, k, starts, w)` | `(n_windows, 4**k)` counts |
| `window_base_counts(code, starts, w)` | `(n_windows, 4)` A/T/C/G counts |

Counts are obtained by binary search over per-k-mer sorted position lists rather than a
Python deque: identical results, roughly two orders of magnitude faster.

## `sagegi.features`

```python
wf = extract_features(code, params, keep_raw=False)
wf.X        # (n_windows, 11) feature matrix
wf.starts   # window start coordinates
wf.gc       # raw GC content per window
```

The 11 columns, in order:

```
gc_content, gc_skew,
H(dinucleotide bias), H(amino-acid usage), H(4-mer frequency),
PCA_dinuc_1..2, PCA_aa_1..2, PCA_4mer_1..2
```

Two definitions are reproduced **exactly as published**, including their non-standard form,
so that this block *is* SSG-LUGIA rather than a lookalike:

- `karlin_dinucleotide(..., "normalized")` — `n(XY) / (sqrt(n(X)) * sqrt(n(Y)) + 1e-6)`,
  computed on raw counts inside the window.
- `offset_entropy(f)` — `sum_i p_i log2 p_i` evaluated on `p = f / sum(f) + 1`; note the `+1`
  offset and the absent negation.

## `sagegi.embed`

```python
cfg = EmbedConfig(model_path=..., tile=1000, batch=16, device="auto")
emb, starts = embed_genome(seq, cfg)          # (n_tiles, 768) float16
pooled      = pool_windows(emb, tile, tile, window_starts, w)
```

`pool_windows` is the closed-form overlap-weighted mean (see **[[Method Overview]]**, §3.1).
It requires non-overlapping tiles and is verified against brute force in the test suite.

## `sagegi.store`

Compact archives. The full 768-d float16 embeddings for all 118 chromosomes come to ~750 MB
per model, which is awkward to move off a Colab session. Because SAGE-GI only uses the
leading principal directions of a genome's own tiles, storing L2-normalised **PCA-96 scores**
loses nothing: PCA components are nested, so slicing the stored matrix to 24 columns is
*identical* to having run PCA with `n_components=24`.

```python
Z, tile, genome_len, is_pcs = load_tiles(path)
```

`is_pcs=True` tells the caller that PCA has already been applied and it should slice columns
rather than re-fit.

## `sagegi.sage`

```python
# the published method: embedding only, no compositional features needed
cfg = SageConfig(n_components=24, assign="center", refine=True, fit_stride=10)
res = run_sage_gi(genome_len, cfg, emb=Z, tile=1000, emb_is_pcs=is_pcs)

# the ablation, if you want to reproduce it
cfg = SageConfig(use_composition=True)
res = run_sage_gi(genome_len, cfg, emb=Z, tile=1000,
                  comp_X=comp.X, comp_starts=comp.starts, emb_is_pcs=is_pcs)
```

`SageConfig` fields worth knowing:

| Field | Default | Effect |
|---|---|---|
| `n_components` | 24 | host subspace dimension |
| `detrend_bp` | 0 | running-median bandwidth for positional detrending. **Leave at 0** — detrending costs 6.9 F1 across the benchmark; the parameter exists so the ablation stays reproducible |
| `use_embedding` | `True` | the species-aware embedding channel |
| `use_composition` | `False` | **Leave at `False`** — adding the eleven compositional features earns nothing over 118 genomes at three operating points, and costs 1.1 F1 at the recall-oriented one. The flag exists so the ablation stays reproducible |
| `contamination_model1` / `2` | 0.15 / 0.05 | the two anomaly-detection stages |
| `assign` | `"center"` | `"start"` reproduces the reference's leftward bias |
| `refine` | `True` | CUSUM boundary refinement |
| `fit_stride` | 10 | fit the robust covariance on every k-th window |

Returns a `SageResult` with `.islands`, `.scores`, `.window_starts`, `.tile_signal` and
`.coarse_islands` (before refinement).

Helpers: `detrend`, `tile_matrix`, `refine_boundaries`, `_cusum_change_point`.

## `sagegi.anomaly`

```python
labels, scores = two_stage_envelope(X, params)     # label −1 / score < 0 means alien
```

`params["fit_stride"]` decimates the *fitting* sample only; every window is always scored. At
the standard geometry consecutive windows overlap by 99%, so `fit_stride=10` is 6.6x faster
and agrees on 99.84% of calls.

## `sagegi.postprocess`

`median_smooth`, `project_to_nucleotides` (`assign="start"` or `"center"`), `runs`,
`drop_short`, `merge_close`, `intervals_to_track`.

Note the `merge_close` convention: `gap` counts the **zero bases between** two runs, so
`gap=0` merges nothing and `gap=2` merges runs separated by up to two bases.

## `sagegi.evaluate`

```python
m  = evaluate_genome(pred, pos, neg, genome_len)   # IslandPick protocol
be = boundary_errors(pred, pos)                    # one row per curated island
summarise_boundaries(be)                           # MABE, median Jaccard, within-tolerance
macro_average(per_genome_df)                       # average over genomes, not nucleotides
```

## `sagegi.synthetic`

`make_chimera`, `assemble_tiles`, `negative_regions`, `shift_host_intervals`,
`composition_distance`. See **[[Method Overview]]** and `scripts/06_synthetic.py`.

## `sagegi.viz`

Figure style and the colour-vision-validated palette. `CAT8` is adjacent-pair safe (bars,
lines); `CAT4` clears the stricter all-pairs gates and is used for scatter plots. Colour
encodes the *feature representation family*, never the rank of a bar.

---

**Next:** **[[Pipeline Reference]]** · **[[Troubleshooting]]**
