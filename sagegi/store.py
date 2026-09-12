r"""Storage format for tile embeddings.

Full 768-dimensional float16 embeddings for all 118 IslandPick chromosomes come to
~750 MB per model, which is awkward to move off a Colab session.  Because SAGE-GI only
ever uses the leading principal directions of a genome's own tiles, we can store the
L2-normalised embeddings already projected onto their first ``k`` principal components
(``k = 96`` by default) and lose nothing: PCA components are nested, so slicing the stored
matrix to 24 columns is *identical* to having run PCA with ``n_components=24``.

An archive therefore contains either

``emb``  : (n_tiles, 768) float16   raw mean-pooled hidden states, or
``pcs``  : (n_tiles, k)  float32    L2-normalised then PCA-projected tile scores,

plus ``starts``, ``tile``, ``step``, ``genome_len`` and, for ``pcs``, the
``explained_variance_ratio``.  :func:`load_tiles` returns a usable matrix either way.
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
from sklearn.decomposition import PCA


def l2(E: np.ndarray) -> np.ndarray:
    E = E.astype(np.float32)
    return E / np.maximum(np.linalg.norm(E, axis=1, keepdims=True), 1e-8)


def reduce_embeddings(emb: np.ndarray, k: int = 96, seed: int = 0):
    """L2-normalise then project onto the leading ``k`` principal components."""
    E = l2(emb)
    k = int(min(k, E.shape[1], max(2, E.shape[0] - 1)))
    p = PCA(n_components=k, svd_solver="randomized", random_state=seed)
    return p.fit_transform(E).astype(np.float32), p.explained_variance_ratio_.astype(np.float32)


def save_reduced(path, emb, starts, tile, step, genome_len, model, k: int = 96) -> None:
    pcs, evr = reduce_embeddings(emb, k)
    np.savez_compressed(path, pcs=pcs, explained_variance_ratio=evr, starts=starts,
                        tile=tile, step=step, genome_len=genome_len, model=model)


def load_tiles(path: str | Path, l2_normalise: bool = True, prefer_full: bool = False):
    """Return ``(Z, tile, genome_len, is_pcs)`` where ``Z`` is (n_tiles, d) float32.

    ``is_pcs`` tells the caller that PCA has already been applied, so it should slice
    columns rather than re-fit.

    Archives in ``RAW_SUBSET`` carry both ``pcs`` and the full-dimensional ``emb``.
    By default the compact ``pcs`` is returned (every existing caller depends on
    this); pass ``prefer_full=True`` to get the full-dimensional matrix instead,
    which chimera assembly requires so that host and donor tiles share one basis.
    """
    z = np.load(path, allow_pickle=False)
    tile = int(z["tile"])
    genome_len = int(z["genome_len"])
    if prefer_full and "emb" in z.files:
        E = z["emb"]
        return (l2(E) if l2_normalise else E.astype(np.float32)), tile, genome_len, False
    if "pcs" in z.files:
        return z["pcs"].astype(np.float32), tile, genome_len, True
    E = z["emb"]
    return (l2(E) if l2_normalise else E.astype(np.float32)), tile, genome_len, False
