#!/usr/bin/env python
"""Detector-free analysis of how well each representation separates islands from backbone.

Every published comparison of genomic-island methods entangles two very different things:
how informative the *features* are, and how well the *detector and its post-processing* are
tuned. This script measures only the first. For each chromosome and each representation we
build a single anomaly statistic per tile -- the robust (median/MAD-standardised) distance
from the chromosome's own centre -- and ask how well it ranks curated island tiles above
curated backbone tiles, as a Mann-Whitney AUC.

An AUC of 0.5 means the representation carries no island signal at all; 1.0 means perfect
separation. Because the statistic is fixed and identical across representations, differences
in AUC are attributable to the representation alone.

We also report the *contrast ratio* (mean statistic inside islands / mean outside), which
matters separately: a representation can rank well while compressing islands and backbone
into a narrow band, which makes thresholding hard for any downstream detector.

    python scripts/08_representation_analysis.py --workers 3
"""
from __future__ import annotations
import argparse, json, os, sys, traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import median_filter
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _common import intervals_for, load_benchmark                    # noqa: E402
from sagegi.features import extract_features                         # noqa: E402
from sagegi.genome import encode, read_fasta                         # noqa: E402
from sagegi.paths import BENCH, EMB, GENOMES, RESULTS                # noqa: E402
from sagegi.sage import SageConfig, tile_matrix                      # noqa: E402
from sagegi.store import load_tiles                                  # noqa: E402

TILE = 1000
DIMS = [8, 16, 24, 48]
DETRENDS = [0, 100_000, 250_000, 1_000_000]
#: any model directory present under $SAGEGI_WORK/emb is analysed; the third
#: model is a general-purpose control for the species-awareness claim.
MODELS = ["DNABERT-S", "DNABERT-2", "NT-v2-50M"]

GLEN = json.loads((BENCH / "genome_lengths.json").read_text())
POS: dict = {}
NEG: dict = {}


def _init(pos_d, neg_d):
    global POS, NEG
    POS, NEG = pos_d, neg_d
    for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ[v] = "1"


def robust_stat(Z: np.ndarray) -> np.ndarray:
    c = np.median(Z, axis=0)
    mad = np.median(np.abs(Z - c), axis=0) * 1.4826 + 1e-8
    return median_filter(np.linalg.norm((Z - c) / mad, axis=1), size=5, mode="nearest")


def auc(sig: np.ndarray, pos: np.ndarray, neg: np.ndarray):
    """Mann-Whitney AUC over curated tiles only; backbone tiles are the negatives."""
    keep = pos | neg
    if pos.sum() < 3 or neg.sum() < 3:
        return np.nan, np.nan, np.nan
    s = sig[keep]
    p = pos[keep]
    r = rankdata(s)
    n1, n0 = int(p.sum()), int((~p).sum())
    a = (r[p].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
    return a, float(sig[pos].mean()), float(sig[neg].mean())


def tiles_mask(intervals, n_tiles):
    m = np.zeros(n_tiles, dtype=bool)
    for s, e in intervals:
        m[max(0, s // TILE):min(n_tiles, e // TILE + 1)] = True
    return m


def one_genome(acc: str):
    rows = []
    n_ref = None
    pos_iv, neg_iv = POS[acc], NEG[acc]

    for model in MODELS:
        f = EMB / f"{model}_t{TILE}" / f"{acc}.npz"
        if not f.exists():
            continue
        Z0, tl, _, is_pcs = load_tiles(f)
        n_ref = len(Z0)
        pm, nm = tiles_mask(pos_iv, n_ref), tiles_mask(neg_iv, n_ref)
        for d in DIMS:
            for dt in DETRENDS:
                Z = tile_matrix(Z0, SageConfig(n_components=d, detrend_bp=dt), tl,
                                is_pcs=is_pcs)
                a, mi, mo = auc(robust_stat(Z), pm, nm)
                rows.append(dict(accession=acc, representation=model, n_components=d,
                                 detrend_bp=dt, auc=a, mean_in=mi, mean_out=mo,
                                 contrast=mi / mo if mo else np.nan))

    # compositional reference, evaluated on the same tile grid
    code = encode(read_fasta(GENOMES / f"{acc}.fna"))
    wf = extract_features(code, dict(w=10000, dw=TILE, karlin_mode="normalized",
                                     pca_dn=2, pca_amino_acid=2, pca_kmer4=2,
                                     entropy_features=True))
    n = len(wf.X) if n_ref is None else min(n_ref, len(wf.X))
    pm, nm = tiles_mask(pos_iv, n), tiles_mask(neg_iv, n)
    a, mi, mo = auc(robust_stat(wf.X[:n]), pm, nm)
    rows.append(dict(accession=acc, representation="compositional", n_components=11,
                     detrend_bp=0, auc=a, mean_in=mi, mean_out=mo,
                     contrast=mi / mo if mo else np.nan))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    pos, neg = load_benchmark()
    have = {p.stem for p in (EMB / f"DNABERT-S_t{TILE}").glob("*.npz")}
    accs = sorted(have & set(pos.accession))
    if args.limit:
        accs = accs[:args.limit]
    pos_d = {a: intervals_for(pos, a) for a in accs}
    neg_d = {a: intervals_for(neg, a) for a in accs}
    print(f"{len(accs)} genomes", flush=True)

    frames = []
    with ProcessPoolExecutor(args.workers, initializer=_init,
                             initargs=(pos_d, neg_d)) as ex:
        futs = {ex.submit(one_genome, a): a for a in accs}
        for i, f in enumerate(as_completed(futs), 1):
            try:
                frames.append(f.result())
            except Exception:                                  # noqa: BLE001
                print(f"{futs[f]} FAILED"); traceback.print_exc()
            if i % 20 == 0 or i == len(accs):
                print(f"  {i}/{len(accs)}", flush=True)

    df = pd.concat(frames, ignore_index=True)
    df.to_csv(RESULTS / "representation_auc.csv", index=False)

    g = (df.groupby(["representation", "n_components", "detrend_bp"])
           [["auc", "contrast"]].mean().reset_index()
           .sort_values("auc", ascending=False))
    print("\n=== island-vs-backbone separation, averaged over genomes ===")
    print(f"{'representation':16s} {'dims':>5s} {'detrend':>9s} {'AUC':>7s} {'contrast':>9s}")
    for r in g.itertuples():
        dt = "off" if r.detrend_bp == 0 else f"{r.detrend_bp//1000}kb"
        print(f"{r.representation:16s} {r.n_components:5d} {dt:>9s} "
              f"{r.auc:7.3f} {r.contrast:9.3f}")


if __name__ == "__main__":
    main()
