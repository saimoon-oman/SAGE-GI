"""Unit tests for the SAGE-GI library.

Run with ``python -m pytest tests -q`` (or ``python tests/test_sagegi.py`` for a plain run).
These cover the parts where a silent error would corrupt every downstream number:
k-mer counting, window pooling, the SSG-LUGIA feature definitions, detrending, the
change-point search, chimera construction, and the evaluation metrics.
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from sagegi.embed import pool_windows, tile_bounds                  # noqa: E402
from sagegi.evaluate import boundary_errors, confusion, evaluate_genome, rates  # noqa: E402
from sagegi.features import extract_features, karlin_dinucleotide, offset_entropy  # noqa: E402
from sagegi.genome import encode, kmer_ids, window_starts           # noqa: E402
from sagegi.kmer import window_base_counts, window_kmer_counts      # noqa: E402
from sagegi.postprocess import drop_short, merge_close, project_to_nucleotides, runs  # noqa: E402
from sagegi.sage import SageConfig, _cusum_change_point, detrend    # noqa: E402
from sagegi.synthetic import (assemble_tiles, make_chimera,         # noqa: E402
                              negative_regions, shift_host_intervals)

RNG = np.random.default_rng(20260831)


def random_seq(n: int, p=(0.25, 0.25, 0.25, 0.25)) -> str:
    return "".join(RNG.choice(list("ATCG"), n, p=list(p)))


# --------------------------------------------------------------------------- encoding
def test_encode_and_kmer_ids():
    assert list(encode("ATCG")) == [0, 1, 2, 3]
    assert encode("ATNG")[2] == -1
    ids = kmer_ids(encode("ATCG"), 2)
    assert list(ids) == [0 * 4 + 1, 1 * 4 + 2, 2 * 4 + 3]
    assert kmer_ids(encode("ANCG"), 2)[0] == -1      # a k-mer spanning an ambiguity is invalid


def test_window_starts_matches_reference_loop():
    """The reference implementation breaks as soon as ``st + w >= len(seq)``."""
    for L, w, dw in [(1000, 100, 10), (1005, 100, 10), (999, 100, 100), (10000, 10000, 100)]:
        expected = []
        for st in range(0, L, dw):
            if st + w >= L:
                break
            expected.append(st)
        assert list(window_starts(L, w, dw)) == expected, (L, w, dw)


# ------------------------------------------------------------------------ k-mer counting
def test_window_kmer_counts_matches_bruteforce():
    seq = random_seq(5000)
    code = encode(seq)
    w, dw, k = 500, 137, 3
    starts = window_starts(len(code), w, dw)
    fast = window_kmer_counts(code, k, starts, w)
    for i in RNG.choice(len(starts), 12, replace=False):
        s = int(starts[i])
        slow = np.zeros(4 ** k, dtype=int)
        for j in range(s, s + w - k + 1):
            kid = kmer_ids(code[j:j + k], k)[0]
            if kid >= 0:
                slow[kid] += 1
        assert np.array_equal(fast[i], slow)


def test_window_base_counts_matches_python_count():
    seq = random_seq(4000)
    code = encode(seq)
    starts = window_starts(len(code), 700, 91)
    fast = window_base_counts(code, starts, 700)
    for i in RNG.choice(len(starts), 10, replace=False):
        s = int(starts[i])
        sub = seq[s:s + 700]
        assert list(fast[i]) == [sub.count(b) for b in "ATCG"]


# ------------------------------------------------------------------------- feature defs
def test_offset_entropy_matches_reference_formula():
    f = RNG.random((5, 16)) * 100
    out = offset_entropy(f)
    for i in range(5):
        p = f[i] / f[i].sum() + 1.0
        assert np.isclose(out[i], (p * np.log2(p)).sum())
    assert offset_entropy(np.zeros((1, 16)))[0] == 0.0     # guard against divide-by-zero


def test_karlin_dinucleotide_normalised_formula():
    nt = np.array([[100, 200, 300, 400]])
    di = np.arange(16, dtype=np.int32)[None, :]
    out = karlin_dinucleotide(nt, di, "normalized")
    for a in range(4):
        for b in range(4):
            expect = di[0, a * 4 + b] / (np.sqrt(nt[0, a]) * np.sqrt(nt[0, b]) + 1e-6)
            assert np.isclose(out[0, a * 4 + b], expect)


def test_extract_features_shape_and_finiteness():
    code = encode(random_seq(60_000))
    params = dict(w=10000, dw=1000, karlin_mode="normalized", pca_dn=2,
                  pca_amino_acid=2, pca_kmer4=2, entropy_features=True)
    wf = extract_features(code, params)
    assert wf.X.shape == (len(wf.starts), 11)
    assert np.isfinite(wf.X).all()
    assert (wf.gc >= 0).all() and (wf.gc <= 1).all()


# ------------------------------------------------------------------------ window pooling
def test_pool_windows_equals_bruteforce_weighted_mean():
    for tile, w, step in [(1000, 10000, 100), (500, 10000, 137), (2000, 7000, 250)]:
        n, H = 90, 5
        e = RNG.normal(size=(n, H)).astype(np.float16)
        ef = e.astype(np.float64)
        starts = np.arange(0, n * tile - w + 1, step)
        fast = pool_windows(e, tile, tile, starts, w)
        for i in RNG.choice(len(starts), 15, replace=False):
            s, t = int(starts[i]), int(starts[i]) + w
            num, den = np.zeros(H), 0.0
            for j in range(n):
                ov = min(t, (j + 1) * tile) - max(s, j * tile)
                if ov > 0:
                    num += ef[j] * ov
                    den += ov
            assert np.allclose(fast[i], num / den, atol=1e-4), (tile, w, step)


def test_tile_bounds_drops_trailing_partial_tile():
    assert list(tile_bounds(2500, 1000, 1000)) == [0, 1000]
    assert list(tile_bounds(3000, 1000, 1000)) == [0, 1000, 2000]


# ---------------------------------------------------------------------------- detrending
def test_detrend_removes_slow_trend_but_keeps_a_block():
    t = np.arange(4000)
    Z = np.column_stack([np.sin(2 * np.pi * t / 4000) * 5.0, np.zeros(4000)])
    Z[1800:1900, 1] += 8.0
    D = detrend(Z, 250)
    assert D[:, 0].std() < 0.15 * Z[:, 0].std()            # genome-scale ramp is gone
    assert np.isclose(D[1800:1900, 1].mean(), 8.0, atol=0.1)   # sharp block survives intact


# ------------------------------------------------------------------------- change points
def test_cusum_finds_the_step():
    for k in (17, 50, 83):
        a = np.concatenate([np.zeros(k), np.ones(100 - k) * 4.0])
        assert abs(_cusum_change_point(a) - k) <= 1
    assert _cusum_change_point(np.zeros(3)) == 1           # degenerate input must not crash


# ------------------------------------------------------------------------ post-processing
def test_runs_merge_and_drop():
    y = np.array([0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0], dtype=np.uint8)
    assert runs(y) == [(1, 2), (5, 5), (7, 9)]
    # ``gap`` counts the zero bases between two runs: 3,4 separate the first pair (gap 2)
    # and 6 separates the second pair (gap 1).
    assert merge_close(runs(y), gap=0) == [(1, 2), (5, 5), (7, 9)]
    assert merge_close(runs(y), gap=1) == [(1, 2), (5, 9)]
    assert merge_close(runs(y), gap=2) == [(1, 9)]
    assert drop_short([(1, 2), (7, 9)], 3) == [(7, 9)]


def test_projection_modes_differ_by_half_a_window():
    alien = np.zeros(200, dtype=bool)
    alien[100:110] = True
    a = runs(project_to_nucleotides(alien, 100_000, 10000, 100, "start"))[0]
    b = runs(project_to_nucleotides(alien, 100_000, 10000, 100, "center"))[0]
    assert b[0] - a[0] == 10000 // 2 - 100 // 2


# ---------------------------------------------------------------------------- evaluation
def test_confusion_and_rates_follow_islandpick_protocol():
    pred = np.zeros(1000, dtype=np.uint8); pred[100:300] = 1
    pos = np.zeros(1000, dtype=np.uint8); pos[150:350] = 1
    neg = np.zeros(1000, dtype=np.uint8); neg[500:700] = 1
    cm = confusion(pred, pos, neg)
    assert cm == dict(TP=150, FP=0, TN=200, FN=50)
    r = rates(cm)
    assert np.isclose(r["precision"], 100.0)
    assert np.isclose(r["recall"], 75.0)
    # sequence in neither curated set is ignored entirely
    assert cm["TP"] + cm["FP"] + cm["TN"] + cm["FN"] == 400


def test_boundary_errors_pick_best_overlapping_prediction():
    df = boundary_errors([(90, 220), (400, 500)], [(100, 200), (900, 1000)])
    assert bool(df.matched.iloc[0]) and not bool(df.matched.iloc[1])
    assert df.start_err.iloc[0] == 10 and df.end_err.iloc[0] == 20
    assert df.signed_start_err.iloc[0] == -10               # predicted edge sits to the left


def test_evaluate_genome_end_to_end():
    m = evaluate_genome([(0, 99)], [(50, 149)], [(200, 299)], 400)
    assert m["TP"] == 50 and m["FN"] == 50 and m["FP"] == 0 and m["TN"] == 100


# ----------------------------------------------------------------------------- synthetic
def test_chimera_is_exact_and_tile_aligned():
    host, donor = random_seq(300_000), random_seq(200_000)
    ch = make_chimera(host, donor, 20000, 3, tile=1000, min_separation=60_000,
                      edge_margin=30_000, rng=np.random.default_rng(2))
    assert len(ch.seq) == len(host) + 3 * 20000
    for (s, e), idx in zip(ch.inserts, ch.donor_tile_index):
        assert ch.seq[s:e + 1] == donor[idx[0] * 1000:(idx[-1] + 1) * 1000]
        assert s % 1000 == 0 and (e + 1) % 1000 == 0

    # tiles assembled from precomputed host/donor matrices must match insert provenance
    host_tiles = np.arange(len(host) // 1000, dtype=float)[:, None]
    donor_tiles = -(np.arange(len(donor) // 1000, dtype=float)[:, None] + 1)
    Z = assemble_tiles(ch, host_tiles, donor_tiles)
    assert len(Z) == len(ch.seq) // 1000
    from_donor = Z[:, 0] < 0
    truth = np.zeros(len(Z), bool)
    for s, e in ch.inserts:
        truth[s // 1000:(e + 1) // 1000] = True
    assert np.array_equal(from_donor, truth)


def test_host_coordinates_shift_by_preceding_inserts():
    host, donor = random_seq(300_000), random_seq(200_000)
    ch = make_chimera(host, donor, 20000, 3, tile=1000, min_separation=60_000,
                      edge_margin=30_000, rng=np.random.default_rng(2))
    for x in (0, 1000, 150_000, 299_999):
        k = sum(1 for o in ch.host_offsets if o <= x)
        assert shift_host_intervals(ch, [(x, x)])[0][0] == x + k * 20000


def test_negative_regions_exclude_inserts_and_flanks():
    host, donor = random_seq(300_000), random_seq(200_000)
    ch = make_chimera(host, donor, 20000, 3, tile=1000, min_separation=60_000,
                      edge_margin=30_000, rng=np.random.default_rng(2))
    neg = negative_regions(ch, [], flank=20000)
    for ns, ne in neg:
        for s, e in ch.inserts:
            assert ne < s - 20000 or ns > e + 20000


def test_sage_config_roundtrips():
    d = SageConfig(n_components=8).to_dict()
    assert d["n_components"] == 8 and d["fit_stride"] >= 1
    assert SageConfig(**{k: v for k, v in d.items()}).n_components == 8



def test_cache_stem_preserves_version_and_config_key():
    """Cache filenames must keep both the accession version and the config digest.

    `Path.with_suffix` reads `.2__4b5cd3ca` in `NC_000913.2__4b5cd3ca` as the extension and
    replaces it, which collapsed every entry to `NC_000913.rows.csv`: editing a
    configuration then invalidated nothing, and two accessions differing only by version
    would have overwritten each other.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "run_sage_gi",
        Path(__file__).resolve().parent.parent / "scripts" / "04_run_sage_gi.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    stem = Path("/cache") / "NC_000913.2__4b5cd3ca35"
    got = mod.stem_path(stem, ".rows.csv").name
    assert got == "NC_000913.2__4b5cd3ca35.rows.csv", got

    # the two properties the bug destroyed, stated directly
    assert "4b5cd3ca35" in got, "configuration digest must survive"
    assert ".2__" in got, "accession version must survive"

    # different configurations and different versions must not collide
    a = mod.stem_path(Path("/c") / "NC_000913.2__aaaaaaaaaa", ".rows.csv").name
    b = mod.stem_path(Path("/c") / "NC_000913.2__bbbbbbbbbb", ".rows.csv").name
    c = mod.stem_path(Path("/c") / "NC_000913.3__aaaaaaaaaa", ".rows.csv").name
    assert len({a, b, c}) == 3, (a, b, c)



def test_assemble_tiles_rejects_mixed_pca_bases():
    """Chimera tiles must not be assembled from per-genome PCA coordinates.

    Archived embeddings are projected onto each chromosome's own principal directions and
    the basis is not stored, so concatenating host and donor tiles yields a matrix whose
    column j means different directions for different rows. Measured effect: 80 kb inserts
    from the most distant available donor separated at AUC 0.27 -- below chance -- while
    real islands in a coherent basis reach 0.78. The construction has to refuse this.
    """
    from sagegi.synthetic import assemble_tiles, make_chimera

    rng = np.random.default_rng(0)
    # large enough to clear the placement margins make_chimera enforces
    host = "".join(rng.choice(list("ACGT"), 400_000))
    donor = "".join(rng.choice(list("ACGT"), 400_000))
    ch = make_chimera(host, donor, 5_000, 1, tile=1000, avoid=[], rng=np.random.default_rng(1))

    hz = rng.normal(size=(400, 8)).astype(np.float32)
    dz = rng.normal(size=(400, 8)).astype(np.float32)

    # full-dimensional embeddings are fine
    out = assemble_tiles(ch, hz, dz, host_is_pcs=False)
    assert out.shape[1] == 8

    # per-genome PCA coordinates must be refused
    try:
        assemble_tiles(ch, hz, dz, host_is_pcs=True)
    except ValueError as e:
        assert "different bases" in str(e) or "PCA-reduced" in str(e)
    else:
        raise AssertionError("assemble_tiles accepted mixed PCA bases")

    # ...unless the caller says so explicitly
    out = assemble_tiles(ch, hz, dz, host_is_pcs=True, allow_mixed_basis=True)
    assert out.shape[1] == 8


if __name__ == "__main__":
    import traceback
    fns = [(n, f) for n, f in sorted(globals().items())
           if n.startswith("test_") and callable(f)]
    failed = 0
    for name, fn in fns:
        try:
            fn()
            print(f"  PASS  {name}")
        except Exception:                                   # noqa: BLE001
            failed += 1
            print(f"  FAIL  {name}")
            traceback.print_exc()
    print(f"\n{len(fns) - failed}/{len(fns)} tests passed")
    sys.exit(1 if failed else 0)
