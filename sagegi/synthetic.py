r"""Controlled synthetic-insertion benchmark.

Curated benchmarks such as IslandPick give *approximate* island boundaries inferred by
comparative genomics, so they can never settle how precisely a method localises an edge.
We therefore also build chimeric chromosomes in which foreign DNA is spliced into a native
host at coordinates we choose, and where the ground truth is exact by construction.

Design
------
Blocks of a **donor** chromosome are inserted into a **host** chromosome at random,
well-separated positions that avoid the host's own annotated islands.  Sweeping the insert
length and the donor-host compositional/phylogenetic distance yields detection-sensitivity
surfaces for each feature representation.

Why this is cheap *and* exact
-----------------------------
Insert offsets and lengths are constrained to be multiples of the embedding tile length.
The tiling of the chimeric chromosome is therefore an exact interleaving of host tiles and
donor tiles - every tile of the chimera is a tile that already exists in the host's or the
donor's precomputed embedding matrix.  So the chimera's embeddings can be assembled by
concatenation instead of being recomputed on a GPU, with no approximation whatsoever.
Compositional features are recomputed on the real chimeric sequence, so the 10 kb windows
correctly see the junctions.
"""
from __future__ import annotations
from dataclasses import dataclass

import numpy as np

from .genome import encode, kmer_ids


# ------------------------------------------------------------------ compositional distance
def tetranucleotide_profile(code: np.ndarray) -> np.ndarray:
    ids = kmer_ids(code, 4)
    ids = ids[ids >= 0]
    v = np.bincount(ids, minlength=256).astype(np.float64)
    return v / max(v.sum(), 1.0)


def gc_fraction(code: np.ndarray) -> float:
    valid = code >= 0
    return float(np.count_nonzero((code == 2) | (code == 3)) / max(np.count_nonzero(valid), 1))


def composition_distance(host_code: np.ndarray, donor_code: np.ndarray) -> dict:
    """Distance summaries between a host and a donor chromosome."""
    ph, pd_ = tetranucleotide_profile(host_code), tetranucleotide_profile(donor_code)
    cos = float(ph @ pd_ / (np.linalg.norm(ph) * np.linalg.norm(pd_) + 1e-12))
    return dict(delta_gc=abs(gc_fraction(host_code) - gc_fraction(donor_code)),
                tnf_euclidean=float(np.linalg.norm(ph - pd_)),
                tnf_manhattan=float(np.abs(ph - pd_).sum()),
                tnf_cosine_dist=1.0 - cos)


# --------------------------------------------------------------------------- construction
@dataclass
class Chimera:
    seq: str
    inserts: list                  # [(start, end)] inclusive, 0-based, in chimera coordinates
    host_len: int
    tile: int
    donor_tile_index: list         # donor tile indices used, per insert
    host_blocks: list              # (host_start_tile, n_tiles) segments retained, in order
    host_offsets: list             # host coordinates at which each block was spliced in
    insert_len: int


def _forbidden_mask(host_len: int, avoid, margin: int) -> np.ndarray:
    m = np.zeros(host_len, dtype=bool)
    for s, e in avoid:
        m[max(0, int(s) - margin):min(host_len, int(e) + 1 + margin)] = True
    return m


def make_chimera(host_seq: str, donor_seq: str, insert_len: int, n_inserts: int,
                 tile: int = 1000, avoid=(), min_separation: int = 100_000,
                 edge_margin: int = 50_000, rng: np.random.Generator | None = None) -> Chimera:
    """Splice ``n_inserts`` blocks of ``donor_seq`` into ``host_seq`` at tile-aligned sites."""
    rng = rng or np.random.default_rng(0)
    if insert_len % tile:
        raise ValueError("insert_len must be a multiple of the tile length")
    host_len = len(host_seq)
    banned = _forbidden_mask(host_len, avoid, margin=25_000)
    banned[:edge_margin] = True
    banned[host_len - edge_margin:] = True

    # candidate tile-aligned host offsets
    cands = np.arange(edge_margin // tile, (host_len - edge_margin) // tile) * tile
    cands = np.array([c for c in cands if not banned[c]], dtype=np.int64)
    rng.shuffle(cands)

    chosen: list[int] = []
    for c in cands:
        if len(chosen) == n_inserts:
            break
        if all(abs(int(c) - p) >= min_separation for p in chosen):
            chosen.append(int(c))
    chosen.sort()
    if len(chosen) < n_inserts:
        raise ValueError(f"host too small/constrained: placed {len(chosen)}/{n_inserts}")

    # donor blocks, tile aligned, non-overlapping
    n_donor_tiles = len(donor_seq) // tile
    need = insert_len // tile
    if n_donor_tiles < need * n_inserts:
        raise ValueError("donor genome too short for the requested inserts")
    slots = rng.permutation(n_donor_tiles - need)[:n_inserts]

    pieces, inserts, donor_idx, host_blocks = [], [], [], []
    cursor, out_pos = 0, 0
    for host_off, slot in zip(chosen, slots):
        pieces.append(host_seq[cursor:host_off])
        host_blocks.append((cursor // tile, (host_off - cursor) // tile))
        out_pos += host_off - cursor
        d0 = int(slot) * tile
        pieces.append(donor_seq[d0:d0 + insert_len])
        inserts.append((out_pos, out_pos + insert_len - 1))
        donor_idx.append(list(range(int(slot), int(slot) + need)))
        out_pos += insert_len
        cursor = host_off
    pieces.append(host_seq[cursor:])
    host_blocks.append((cursor // tile, (len(host_seq) - cursor) // tile))

    return Chimera(seq="".join(pieces), inserts=inserts, host_len=host_len, tile=tile,
                   donor_tile_index=donor_idx, host_blocks=host_blocks,
                   host_offsets=chosen, insert_len=insert_len)


def assemble_tiles(ch: Chimera, host_tiles: np.ndarray, donor_tiles: np.ndarray,
                   host_is_pcs: bool = False, allow_mixed_basis: bool = False) -> np.ndarray:
    """Exact tile matrix of the chimera, assembled from precomputed host and donor tiles.

    The tile alignment is exact: inserts are tile-aligned with tile-multiple lengths, so no
    tile of the chimera ever straddles a junction.

    The *representation* is only exact if both matrices live in the same coordinate system.
    Tile archives store each chromosome projected onto its own principal directions
    (``is_pcs=True``) and do not retain the basis, so concatenating two such matrices gives
    a chimera in which column ``j`` means a different direction for host and donor rows.
    The effect is not subtle: 80 kb insertions from the most compositionally distant donor
    available separate from host tiles at AUC 0.27, below chance, where real islands in a
    coherent basis reach 0.78.

    Pass full-dimensional embeddings (``host_is_pcs=False``), or set
    ``allow_mixed_basis=True`` if you are deliberately reproducing the earlier, invalid
    construction.
    """
    if host_is_pcs and not allow_mixed_basis:
        raise ValueError(
            "assemble_tiles received per-genome PCA-reduced tiles. Host and donor "
            "coordinates are then expressed in different bases and the assembled matrix is "
            "not a valid embedding of the chimera. Re-extract full-dimensional embeddings "
            "for the host and donor genomes, or pass allow_mixed_basis=True to reproduce "
            "the earlier invalid behaviour deliberately.")
    out = []
    for i, (h0, hn) in enumerate(ch.host_blocks):
        if hn > 0:
            out.append(host_tiles[h0:h0 + hn])
        if i < len(ch.donor_tile_index):
            out.append(donor_tiles[ch.donor_tile_index[i]])
    Z = np.concatenate(out, axis=0)
    n_expected = len(ch.seq) // ch.tile
    return Z[:n_expected]


def negative_regions(ch: Chimera, avoid_in_chimera, flank: int = 20_000):
    """Host-derived stretches that should *not* be called: everything except the inserts,
    the host's own annotated islands, and a flank around every insert."""
    n = len(ch.seq)
    mask = np.ones(n, dtype=bool)
    for s, e in ch.inserts:
        mask[max(0, s - flank):min(n, e + 1 + flank)] = False
    for s, e in avoid_in_chimera:
        mask[max(0, int(s)):min(n, int(e) + 1)] = False
    d = np.diff(np.concatenate(([0], mask.view(np.uint8), [0])).astype(np.int8))
    return list(zip(np.flatnonzero(d == 1).tolist(), (np.flatnonzero(d == -1) - 1).tolist()))


def to_chimera_coords(ch: Chimera, x: int) -> int:
    """Map a host coordinate to its position in the chimera."""
    k = int(np.searchsorted(np.asarray(ch.host_offsets, dtype=np.int64), int(x), side="right"))
    return int(x) + k * ch.insert_len


def shift_host_intervals(ch: Chimera, intervals):
    """Map host intervals into chimera coordinates (inserts push everything rightwards)."""
    return [(to_chimera_coords(ch, int(s)), to_chimera_coords(ch, int(e))) for s, e in intervals]
