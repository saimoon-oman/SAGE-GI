"""Genome I/O and integer encoding."""
from __future__ import annotations
from pathlib import Path

import numpy as np

BASES = "ATCG"                       # index order matches the reference SSG-LUGIA implementation
_LUT = np.full(256, -1, dtype=np.int8)
for _i, _b in enumerate(BASES):
    _LUT[ord(_b)] = _i
    _LUT[ord(_b.lower())] = _i


def read_fasta(path: str | Path) -> str:
    """Return the concatenated sequence of a single-record FASTA file, upper-cased."""
    seq = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq.append(line.strip())
    return "".join(seq).upper()


def encode(seq: str) -> np.ndarray:
    """Map a DNA string to int8 codes 0..3 (A,T,C,G); every other symbol becomes -1."""
    return _LUT[np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)]


def kmer_ids(code: np.ndarray, k: int) -> np.ndarray:
    """Base-4 identifier of the k-mer starting at each position; -1 if it spans an ambiguity.

    Returns an array of length ``len(code) - k + 1``.
    """
    n = len(code) - k + 1
    if n <= 0:
        return np.empty(0, dtype=np.int32)
    ids = np.zeros(n, dtype=np.int32)
    bad = np.zeros(n, dtype=bool)
    for j in range(k):
        c = code[j:j + n]
        ids = ids * 4 + np.maximum(c, 0).astype(np.int32)
        bad |= (c < 0)
    ids[bad] = -1
    return ids


def window_starts(genome_len: int, w: int, dw: int) -> np.ndarray:
    """Window start offsets exactly as produced by the reference implementation.

    The reference loops ``for st in range(0, len(seq), dw)`` and breaks as soon as
    ``st + w >= len(seq)``, so the last accepted start satisfies ``st + w < len(seq)``.
    """
    n = max(0, -(-(genome_len - w) // dw))            # ceil((L-w)/dw)
    return np.arange(n, dtype=np.int64) * dw
