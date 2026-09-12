"""Post-processing: median smoothing, window→nucleotide projection, small-island removal."""
from __future__ import annotations

import numpy as np
from scipy.ndimage import median_filter


def median_smooth(ys: np.ndarray, wlen: int) -> np.ndarray:
    """Median filter matching the reference implementation, including its untouched edges."""
    if wlen <= 1 or len(ys) <= wlen:
        return np.asarray(ys, dtype=float).copy()
    half = wlen // 2
    out = median_filter(np.asarray(ys, dtype=float), size=wlen, mode="nearest")
    out[:half] = ys[:half]
    out[len(ys) - half:] = ys[len(ys) - half:]
    return out


def project_to_nucleotides(alien: np.ndarray, genome_len: int, w: int, dw: int,
                           assign: str = "start") -> np.ndarray:
    """Turn a per-window alien/native call into a per-nucleotide 0/1 track.

    ``assign="start"`` reproduces SSG-LUGIA, which credits window *i* to nucleotides
    ``[i*dw, (i+1)*dw)`` — the leading ``dw`` bases of the window.  Because the window
    extends ``w`` bases to the right of that anchor, predictions inherit a systematic
    left shift of roughly ``w/2``.  ``assign="center"`` instead credits the window to
    ``[i*dw + w/2 - dw/2, i*dw + w/2 + dw/2)``, removing that bias.
    """
    y = np.zeros(genome_len, dtype=np.uint8)
    idx = np.flatnonzero(alien)
    if len(idx) == 0:
        return y
    if assign == "start":
        lo, hi = idx * dw, (idx + 1) * dw
    elif assign == "center":
        shift = w // 2 - dw // 2
        lo, hi = idx * dw + shift, (idx + 1) * dw + shift
    else:
        raise ValueError(f"unknown assign mode {assign!r}")
    lo = np.clip(lo, 0, genome_len)
    hi = np.clip(hi, 0, genome_len)
    for a, b in zip(lo, hi):
        y[a:b] = 1
    return y


def runs(y: np.ndarray) -> list[tuple[int, int]]:
    """Maximal runs of 1s as inclusive 0-based ``(start, end)`` pairs."""
    if len(y) == 0:
        return []
    d = np.diff(np.concatenate(([0], (y > 0).view(np.uint8), [0])).astype(np.int8))
    return list(zip(np.flatnonzero(d == 1).tolist(),
                    (np.flatnonzero(d == -1) - 1).tolist()))


def drop_short(intervals: list[tuple[int, int]], min_len: int) -> list[tuple[int, int]]:
    return [(s, e) for s, e in intervals if (e - s + 1) >= min_len]


def merge_close(intervals: list[tuple[int, int]], gap: int) -> list[tuple[int, int]]:
    """Merge intervals separated by at most ``gap`` bases."""
    if not intervals:
        return []
    out = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s - out[-1][1] - 1 <= gap:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def intervals_to_track(intervals, genome_len: int) -> np.ndarray:
    y = np.zeros(genome_len, dtype=np.uint8)
    for s, e in intervals:
        y[max(0, s):min(genome_len, e + 1)] = 1
    return y
