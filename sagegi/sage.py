r"""SAGE-GI - Species-Aware Genomic Embeddings for Genomic Island prediction.

Pipeline
--------
1. **L2 normalisation** of tile embeddings.  DNABERT-S is trained with a cosine-similarity
   contrastive objective, so the unit sphere is the geometry its species signal lives in.
2. **Host subspace reduction.**  A 768-dimensional embedding cannot support a robust
   covariance estimate (MCD is O(d^3) and needs n >> d), and most directions carry no
   within-genome variation at all.  We project onto the leading ``n_components`` principal
   directions *of the genome under analysis* - the host defines its own subspace.
3. **Positional detrending - available but OFF by default, because it does not work.**
   The idea was that a species-aware embedding, trained to separate *different* genomes,
   would inside a single chromosome be dominated by host structure irrelevant to horizontal
   transfer (the GC-skew inversion at the origin and terminus of replication, replichore
   asymmetry, slow compositional drift), and that subtracting a long-wavelength running
   median would remove it while leaving 10-200 kb islands intact.

   Measured over all 118 IslandPick chromosomes, it costs **6.9 F1** (52.7 with a 250 kb
   filter versus 59.6 without; paired Wilcoxon p = 5e-12). Genomic islands are simply not
   short enough relative to any filter wide enough to capture replichore structure - a
   134 kb island moves even a median inside a 250 kb window - so the filter removes signal
   along with the trend. ``detrend_bp`` is retained so the ablation is reproducible; leave
   it at 0.

   What *is* needed is the host-subspace projection in step 2, which conditions the problem
   without touching the genome axis.
4. **Window pooling** onto SSG-LUGIA's (10 kb, 100 bp) window grid, with fractional
   overlap weights (see :func:`sagegi.embed.pool_windows`).
5. **Fusion** with the compositional feature block.  Because the detector scores windows by
   Mahalanobis distance, which is affine-invariant, the two blocks can simply be
   concatenated: no scaling constant has to be chosen.
6. **Two-stage robust anomaly detection**, identical in form to SSG-LUGIA.
7. **Median smoothing, centred window-to-nucleotide projection, small-island removal.**
8. **CUSUM boundary refinement** at tile resolution, which is what converts a coarse
   10 kb-window call into a sharp coordinate.
"""
from __future__ import annotations
from dataclasses import asdict, dataclass, field

import numpy as np
from scipy.ndimage import median_filter
from sklearn.decomposition import PCA

from .anomaly import two_stage_envelope
from .embed import pool_windows
from .postprocess import (drop_short, median_smooth, merge_close,
                          project_to_nucleotides, runs)


# --------------------------------------------------------------------------------- config
@dataclass
class SageConfig:
    # window geometry (kept identical to SSG-LUGIA for commensurability)
    w: int = 10000
    dw: int = 100
    # embedding channel
    n_components: int = 24
    detrend_bp: int = 0          # positional detrending is OFF by default: the
                                 # 118-genome ablation shows it costs 6.9 F1
    l2_normalise: bool = True
    # feature composition
    use_embedding: bool = True
    use_composition: bool = False   # the compositional channel is OFF by default: over 118
                                    # genomes and three operating points it never helps, and
                                    # at the recall-oriented point it costs 1.1 F1 (p=2e-4)
    # anomaly detection
    contamination_model1: float = 0.15
    contamination_model2: float = 0.05
    support_fraction_model1: float = 0.75
    support_fraction_model2: float = 0.9
    fit_stride: int = 10          # fit the robust covariance on every k-th window (see anomaly.py)
    # post-processing
    median_filter_window_len: int = 400
    min_island_len: int = 10000
    merge_gap: int = 5000
    assign: str = "center"
    # boundary refinement
    refine: bool = True
    refine_radius: int = 15000
    refine_context: int = 20000

    def to_dict(self) -> dict:
        return asdict(self)


# ------------------------------------------------------------------- embedding preparation
def detrend(Z: np.ndarray, width_tiles: int) -> np.ndarray:
    """Subtract a long-wavelength running median along the genome axis."""
    if width_tiles < 3 or width_tiles >= len(Z):
        return Z - np.median(Z, axis=0, keepdims=True)
    trend = median_filter(Z, size=(width_tiles, 1), mode="nearest")
    return Z - trend


def tile_matrix(emb: np.ndarray, cfg: SageConfig, tile: int,
                is_pcs: bool = False) -> np.ndarray:
    """L2-normalise, then PCA, then positional detrending. Returns (n_tiles, n_components).

    ``is_pcs=True`` means ``emb`` already holds L2-normalised PCA scores (see
    :mod:`sagegi.store`); because PCA components are nested, slicing the leading columns is
    exactly equivalent to having fitted PCA with ``n_components`` in the first place.
    """
    if is_pcs:
        Z = np.asarray(emb, dtype=np.float32)[:, :cfg.n_components]
    else:
        E = emb.astype(np.float32)
        if cfg.l2_normalise:
            E = E / np.maximum(np.linalg.norm(E, axis=1, keepdims=True), 1e-8)
        d = min(cfg.n_components, E.shape[1], max(2, E.shape[0] - 1))
        Z = PCA(n_components=d, svd_solver="randomized", random_state=0).fit_transform(E)
    if cfg.detrend_bp > 0:
        Z = detrend(Z, max(3, cfg.detrend_bp // tile))
    return Z


# ---------------------------------------------------------------------- boundary refinement
def _cusum_change_point(a: np.ndarray) -> int:
    """Index k maximising the normalised two-sample mean-difference statistic of ``a``."""
    n = len(a)
    if n < 4:
        return n // 2
    c = np.cumsum(a, dtype=np.float64)
    k = np.arange(1, n, dtype=np.float64)
    left = c[:-1] / k
    right = (c[-1] - c[:-1]) / (n - k)
    stat = np.sqrt(k * (n - k) / n) * np.abs(left - right)
    return int(np.argmax(stat) + 1)


def refine_boundaries(intervals, signal: np.ndarray, tile: int, cfg: SageConfig,
                      genome_len: int):
    r"""Sharpen each island's coordinates with a local CUSUM change-point search.

    ``signal`` is a tile-resolution anomaly track (larger = more alien).  Around each coarse
    boundary we take a window spanning ``refine_context`` bases of putative host on the
    outside and ``refine_radius`` on the inside, and place the boundary at the maximiser of
    the normalised mean-shift statistic - i.e. where the anomaly track actually steps.
    """
    out = []
    n = len(signal)
    for s, e in intervals:
        lo = max(0, (s - cfg.refine_context) // tile)
        hi = min(n, (s + cfg.refine_radius) // tile + 1)
        new_s = s
        if hi - lo >= 4:
            new_s = int((lo + _cusum_change_point(signal[lo:hi])) * tile)

        lo2 = max(0, (e - cfg.refine_radius) // tile)
        hi2 = min(n, (e + cfg.refine_context) // tile + 1)
        new_e = e
        if hi2 - lo2 >= 4:
            new_e = int((lo2 + _cusum_change_point(signal[lo2:hi2])) * tile) - 1

        if new_e <= new_s:                       # refinement collapsed: keep the original
            new_s, new_e = s, e
        out.append((max(0, new_s), min(genome_len - 1, new_e)))
    return out


# --------------------------------------------------------------------------------- pipeline
@dataclass
class SageResult:
    islands: list
    scores: np.ndarray
    window_starts: np.ndarray
    tile_signal: np.ndarray = field(default_factory=lambda: np.empty(0))
    coarse_islands: list = field(default_factory=list)


def run_sage_gi(genome_len: int, cfg: SageConfig, emb: np.ndarray | None = None,
                tile: int = 1000, comp_X: np.ndarray | None = None,
                comp_starts: np.ndarray | None = None,
                emb_is_pcs: bool = False) -> SageResult:
    """Predict genomic islands. Supply ``emb`` and/or ``comp_X``; at least one is required."""
    blocks, starts = [], comp_starts
    Zt = None

    if cfg.use_embedding:
        if emb is None:
            raise ValueError("use_embedding=True but no embeddings supplied")
        Zt = tile_matrix(emb, cfg, tile, is_pcs=emb_is_pcs)
        if starts is None:
            n_win = max(0, -(-(genome_len - cfg.w) // cfg.dw))
            starts = np.arange(n_win, dtype=np.int64) * cfg.dw
        blocks.append(pool_windows(Zt, tile, tile, starts, cfg.w))

    if cfg.use_composition:
        if comp_X is None:
            raise ValueError("use_composition=True but no compositional features supplied")
        blocks.append(np.asarray(comp_X, dtype=np.float32))

    if not blocks:
        raise ValueError("no feature channel enabled")

    n = min(len(b) for b in blocks)
    X = np.column_stack([b[:n] for b in blocks]).astype(np.float64)
    starts = np.asarray(starts)[:n]

    _, scores = two_stage_envelope(X, cfg.to_dict())

    smoothed = median_smooth(scores, cfg.median_filter_window_len)
    track = project_to_nucleotides(smoothed < 0, genome_len, cfg.w, cfg.dw, cfg.assign)
    coarse = drop_short(merge_close(runs(track), cfg.merge_gap), cfg.min_island_len)

    islands, tile_sig = coarse, np.empty(0)
    if cfg.refine and Zt is not None and len(coarse):
        # tile-resolution anomaly track: robust distance from the host centre in the
        # detrended embedding subspace
        centre = np.median(Zt, axis=0)
        mad = np.median(np.abs(Zt - centre), axis=0) * 1.4826 + 1e-8
        tile_sig = np.linalg.norm((Zt - centre) / mad, axis=1)
        tile_sig = median_filter(tile_sig, size=5, mode="nearest")
        islands = refine_boundaries(coarse, tile_sig, tile, cfg, genome_len)
        islands = drop_short(merge_close(sorted(islands), cfg.merge_gap), cfg.min_island_len)

    return SageResult(islands=islands, scores=scores, window_starts=starts,
                      tile_signal=tile_sig, coarse_islands=coarse)
