"""Compositional feature extraction (faithful, vectorised SSG-LUGIA re-implementation).

Every quantity here reproduces ``feature_extraction.py`` of the reference SSG-LUGIA
release (github.com/nibtehaz/SSG-LUGIA) bit-for-bit, including its two non-standard
choices, which we keep deliberately so that our baseline *is* SSG-LUGIA and not a
lookalike:

* the *normalised* Karlin dinucleotide bias  ``n(XY) / (sqrt(n(X)) sqrt(n(Y)) + 1e-6)``
  where ``n(.)`` are raw counts inside the window;
* the "entropy" statistic ``sum_i p_i log2 p_i`` evaluated on ``p = f / sum(f) + 1``
  (offset by one and unnegated).

The resulting window feature vector has 11 dimensions, in this order::

    [GC content, GC skew,
     H(dinucleotide bias), H(amino-acid usage), H(4-mer frequency),
     PC1..PC2 of dinucleotide bias, PC1..PC2 of amino-acid usage, PC1..PC2 of 4-mer freq]
"""
from __future__ import annotations
from dataclasses import dataclass

import numpy as np
from sklearn.decomposition import PCA

from .genome import BASES, window_starts
from .kmer import window_base_counts, window_kmer_counts

FEATURE_NAMES = ["gc_content", "gc_skew", "H_dinuc", "H_aa", "H_4mer",
                 "PCA_dinuc_1", "PCA_dinuc_2", "PCA_aa_1", "PCA_aa_2",
                 "PCA_4mer_1", "PCA_4mer_2"]

_CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L",
    "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V",
    "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P",
    "CCC": "P", "CCA": "P", "CCG": "P", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*",
    "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "TGT": "C",
    "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G",
    "GGG": "G",
}


def _codon_to_aa_index() -> np.ndarray:
    """Map each of the 64 base-4 3-mer ids to one of 21 amino-acid / stop classes."""
    aas = sorted(set(_CODON_TABLE.values()))
    aa_ix = {a: i for i, a in enumerate(aas)}
    out = np.zeros(64, dtype=np.int64)
    for cid in range(64):
        codon = "".join(BASES[(cid >> (2 * (2 - j))) & 3] for j in range(3))
        out[cid] = aa_ix[_CODON_TABLE[codon]]
    return out


_AA_IX = _codon_to_aa_index()
N_AA = int(_AA_IX.max()) + 1


def offset_entropy(freq: np.ndarray) -> np.ndarray:
    """SSG-LUGIA's entropy statistic, row-wise. ``freq`` is (n_windows, n_bins)."""
    s = freq.sum(axis=1, keepdims=True)
    p = np.divide(freq, s, out=np.zeros_like(freq, dtype=np.float64),
                  where=np.abs(s) > 1e-6) + 1.0
    h = (p * np.log2(p)).sum(axis=1)
    h[np.abs(s[:, 0]) <= 1e-6] = 0.0
    return h


def karlin_dinucleotide(base_counts: np.ndarray, di_counts: np.ndarray,
                        mode: str = "normalized") -> np.ndarray:
    """Karlin dinucleotide bias for every window. ``di_counts`` is (n_windows, 16)."""
    di = di_counts.astype(np.float64)
    if mode == "raw":
        return di
    x = base_counts[:, :, None].astype(np.float64)     # first base of the dimer
    y = base_counts[:, None, :].astype(np.float64)     # second base
    denom = (x * y) if mode == "original" else np.sqrt(x) * np.sqrt(y)
    return di / (denom.reshape(len(di), 16) + 1e-6)


@dataclass
class WindowFeatures:
    X: np.ndarray                 # (n_windows, 11) feature matrix
    starts: np.ndarray            # window start coordinates (0-based)
    w: int
    dw: int
    gc: np.ndarray                # convenience: raw GC content per window
    dinuc: np.ndarray             # (n_windows, 16) Karlin bias, pre-PCA
    aa: np.ndarray                # (n_windows, 21) amino-acid usage, pre-PCA
    kmer4: np.ndarray             # (n_windows, 256) 4-mer counts, pre-PCA


def extract_features(code: np.ndarray, params: dict, keep_raw: bool = False) -> WindowFeatures:
    """Compute the 11-dimensional SSG-LUGIA feature matrix for one genome."""
    w, dw = params["w"], params["dw"]
    starts = window_starts(len(code), w, dw)

    nt = window_base_counts(code, starts, w)                       # (n, 4) A,T,C,G
    di = window_kmer_counts(code, 2, starts, w)                    # (n, 16)
    tri = window_kmer_counts(code, 3, starts, w)                   # (n, 64)
    tet = window_kmer_counts(code, 4, starts, w)                   # (n, 256)

    g, c = nt[:, 3].astype(np.float64), nt[:, 2].astype(np.float64)
    gc_content = (g + c) / w
    with np.errstate(divide="ignore", invalid="ignore"):
        gc_skew = np.where((g + c) > 0, (g - c) / (g + c), 0.0)

    dinuc = karlin_dinucleotide(nt, di, params.get("karlin_mode", "normalized"))
    aa = np.zeros((len(starts), N_AA), dtype=np.float64)
    np.add.at(aa.T, _AA_IX, tri.T.astype(np.float64))

    cols = [gc_content, gc_skew]
    if params.get("entropy_features", True):
        cols += [offset_entropy(dinuc), offset_entropy(aa), offset_entropy(tet.astype(np.float64))]

    blocks = []
    if params.get("pca_dn", 2) > 0:
        blocks.append(PCA(n_components=params["pca_dn"], svd_solver="full").fit_transform(dinuc))
    blocks.append(PCA(n_components=params["pca_amino_acid"], svd_solver="full").fit_transform(aa))
    blocks.append(PCA(n_components=params["pca_kmer4"], svd_solver="full")
                  .fit_transform(tet.astype(np.float64)))

    X = np.column_stack(cols + blocks)
    return WindowFeatures(X=X, starts=starts, w=w, dw=dw, gc=gc_content,
                          dinuc=dinuc if keep_raw else np.empty(0),
                          aa=aa if keep_raw else np.empty(0),
                          kmer4=tet if keep_raw else np.empty(0))
