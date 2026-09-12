"""Fast sliding-window k-mer counting.

The reference SSG-LUGIA implementation keeps running counts with a deque, which is O(L)
but pure Python and takes minutes per genome.  We obtain identical counts in a fully
vectorised way: for every k-mer we store the sorted list of positions at which it starts
and answer "how many occurrences begin inside this window" with two `searchsorted` calls.
Cost is O(4^k * n_windows * log L) with numpy-level constants.
"""
from __future__ import annotations

import numpy as np

from .genome import kmer_ids


def window_kmer_counts(code: np.ndarray, k: int, starts: np.ndarray, w: int) -> np.ndarray:
    """Counts of every k-mer whose first base lies in ``[s, s + w - k + 1)`` for each start.

    Matches the reference implementation, which counts k-mers at offsets
    ``0 .. len(window) - k`` of each window.
    """
    ids = kmer_ids(code, k)
    n_k = 4 ** k
    out = np.zeros((len(starts), n_k), dtype=np.int32)
    order = np.argsort(ids, kind="stable")
    ids_sorted = ids[order]
    # positions of each k-mer id inside the sorted array
    bounds = np.searchsorted(ids_sorted, np.arange(-1, n_k + 1))
    ends = starts + (w - k + 1)
    for kid in range(n_k):
        lo, hi = bounds[kid + 1], bounds[kid + 2]
        if hi <= lo:
            continue
        pos = np.sort(order[lo:hi])
        out[:, kid] = (np.searchsorted(pos, ends, side="left")
                       - np.searchsorted(pos, starts, side="left"))
    return out


def window_base_counts(code: np.ndarray, starts: np.ndarray, w: int) -> np.ndarray:
    """(n_windows, 4) counts of A,T,C,G in each window."""
    out = np.zeros((len(starts), 4), dtype=np.int32)
    ends = starts + w
    for b in range(4):
        pos = np.flatnonzero(code == b)
        out[:, b] = (np.searchsorted(pos, ends, side="left")
                     - np.searchsorted(pos, starts, side="left"))
    return out
