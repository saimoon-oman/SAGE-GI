# SAGE-GI

**S**pecies-**A**ware **G**enomic **E**mbeddings for **G**enomic **I**sland prediction.

SAGE-GI finds genomic islands — blocks of horizontally acquired DNA — in a **single,
unannotated bacterial chromosome**, using frozen embeddings from
[DNABERT-S](https://huggingface.co/zhihan1996/DNABERT-S), a genome foundation model trained
so that sequences from different species separate in embedding space.

No annotation file. No reference genome. No training on your organism.

---

## Start here

| If you want to… | Go to |
|---|---|
| Install it and call it on a FASTA file | **[[Installation]]** → **[[Quick Start]]** |
| Understand *why* it is built this way | **[[Method Overview]]** |
| Re-run every experiment in the paper | **[[Reproducing the Paper]]** |
| Look up a function or a script flag | **[[API Reference]]** · **[[Pipeline Reference]]** |
| Know exactly what data was used | **[[Benchmarks and Data]]** |
| Fix something that broke | **[[Troubleshooting]]** |

---

## The one-paragraph version

Composition-based genomic island detectors describe each window of a genome with statistics
chosen by hand — GC skew, dinucleotide bias, codon usage, tetranucleotide frequency. A
genome foundation model could supply a *learned* description instead, and DNABERT-S is
trained on precisely the right signal: its contrastive objective exists to make sequences
from different species separate. But three things stand in the way, and fixing them is what
SAGE-GI is:

1. **Cost.** A 10 kb window stepped every 100 bp means every base is embedded ~100 times.
   SAGE-GI embeds non-overlapping 1 kb tiles *once* and reconstructs each window as the
   overlap-weighted mean of the tiles it covers — a closed form, exact to 3×10⁻⁸, that makes
   cost linear in genome length.
2. **Conditioning.** A 768-dimensional embedding cannot support a robust covariance
   estimate. SAGE-GI L2-normalises and projects onto the leading principal directions *of the
   chromosome being analysed*, so the host defines its own frame of reference. (We also tried
   removing a long-wavelength positional trend, expecting replichore structure to dominate.
   Measured across 118 genomes it *costs* 6.9 F1 — a negative result reported in full.)
3. **Boundaries.** Nucleotide F1 rewards *covering* an island and says nothing about *where*
   its edges are. SAGE-GI reports boundary error explicitly and places each edge with a
   CUSUM change-point search on the tile-resolution anomaly track.

We expected the learned and hand-crafted channels to be complementary, and built SAGE-GI to
fuse them. Measured over 118 genomes at three operating points, the compositional block earns
nothing on top of the embedding — and at the recall-oriented point it *hurts* (−1.1 F1,
*p* = 2×10⁻⁴). **SAGE-GI is therefore the embedding-only pipeline**, and the fusion is reported
as a second negative result.

---

## At a glance

| | |
|---|---|
| **Input** | one FASTA chromosome |
| **Output** | island start/end coordinates (and optionally the anomaly track) |
| **Supervision** | none — fully unsupervised |
| **External requirements** | none at inference; model weights are downloaded once |
| **Benchmark** | IslandPick — 118 chromosomes, 485 Mbp, 771 curated islands |
| **Runtime** | ≈1 min per 5 Mbp chromosome on a T4 GPU; ≈45 min on a laptop CPU |
| **Licence** | MIT (models are Apache-2.0) |

## Links

- Repository — <https://github.com/saimoon-oman/SAGE-GI>
- Project site — <https://saimoon-oman.github.io/SAGE-GI/>
- GPU extraction notebook — [`notebooks/01_extract_embeddings_gpu.ipynb`](https://github.com/saimoon-oman/SAGE-GI/blob/main/notebooks/01_extract_embeddings_gpu.ipynb)
- Progress log — [`log.md`](https://github.com/saimoon-oman/SAGE-GI/blob/main/log.md)

## Citing

```bibtex
@article{oman2026sagegi,
  title   = {SAGE-GI: species-aware genomic language model embeddings sharpen
             unsupervised genomic island detection and boundary resolution},
  author  = {Oman, Saimoon Al Farshi and Bayzid, Md. Shamsuzzoha},
  journal = {Briefings in Bioinformatics},
  year    = {2026}
}
```
