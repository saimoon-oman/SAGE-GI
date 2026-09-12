#!/usr/bin/env python
"""Run the compositional (SSG-LUGIA) baseline over every available IslandPick genome.

For each chromosome we extract the 11-dimensional window feature matrix once and reuse it
for all model variants.  Per-genome outputs:

* ``$SAGEGI_WORK/comp/<acc>.npz``  – feature matrix, window starts, anomaly scores per variant
* ``results/ssg_lugia_per_genome.csv`` – nucleotide-level metrics per genome and variant
* ``results/ssg_lugia_boundaries.csv`` – per-island boundary errors
* ``results/ssg_lugia_intervals.csv``  – the predicted island coordinates themselves
"""
from __future__ import annotations
import argparse, os, sys, time, traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _common import (available_accessions, genome_codes, intervals_for,  # noqa: E402
                     load_benchmark)
from sagegi import models                                                # noqa: E402
from sagegi.paths import RESULTS, WORK                                   # noqa: E402
from sagegi.features import extract_features                             # noqa: E402
from sagegi.anomaly import two_stage_envelope                            # noqa: E402
from sagegi.pipeline import predict_from_scores                          # noqa: E402
from sagegi.evaluate import boundary_errors, evaluate_genome, summarise_boundaries  # noqa: E402

COMP = WORK / "comp"
COMP.mkdir(parents=True, exist_ok=True)
VARIANTS = ["SSG-LUGIA-P", "SSG-LUGIA-F", "SSG-LUGIA-R"]


def one_genome(acc: str, redo: bool = False):
    t0 = time.time()
    code, L = genome_codes(acc)
    cache = COMP / f"{acc}.npz"
    params0 = models.get(VARIANTS[0])

    if cache.exists() and not redo:
        z = np.load(cache)
        X, starts = z["X"].astype(np.float64), z["starts"]
        scores = {v: z[f"score_{v}"] for v in VARIANTS if f"score_{v}" in z}
    else:
        wf = extract_features(code, params0)
        X, starts, scores = wf.X, wf.starts, {}

    rows, brows, irows = [], [], []
    for v in VARIANTS:
        p = models.get(v)
        if v not in scores:
            scores[v] = two_stage_envelope(X, p)[1]
        for assign in ("start", "center"):
            pp = {**p, "assign": assign}
            isl = predict_from_scores(scores[v], L, pp)
            name = v if assign == "start" else f"{v}+center"
            m = evaluate_genome(isl, POS[acc], NEG[acc], L)
            be = boundary_errors(isl, POS[acc])
            rows.append(dict(accession=acc, method=name, genome_len=L,
                             n_pred=len(isl),
                             pred_bp=int(sum(e - s + 1 for s, e in isl)), **m,
                             **summarise_boundaries(be)))
            be.insert(0, "method", name); be.insert(0, "accession", acc)
            brows.append(be)
            for s, e in isl:
                irows.append(dict(accession=acc, method=name, start=s + 1, end=e + 1))

    np.savez_compressed(cache, X=X.astype(np.float32), starts=starts,
                        **{f"score_{v}": scores[v].astype(np.float32) for v in VARIANTS})
    return (pd.DataFrame(rows), pd.concat(brows, ignore_index=True),
            pd.DataFrame(irows), time.time() - t0)


POS: dict = {}
NEG: dict = {}


def _init(pos_d, neg_d):
    global POS, NEG
    POS, NEG = pos_d, neg_d
    os.environ["OMP_NUM_THREADS"] = "1"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--redo", action="store_true")
    args = ap.parse_args()

    pos, neg = load_benchmark()
    accs = [a for a in available_accessions() if a in set(pos.accession)]
    if args.limit:
        accs = accs[:args.limit]
    pos_d = {a: intervals_for(pos, a) for a in accs}
    neg_d = {a: intervals_for(neg, a) for a in accs}
    print(f"{len(accs)} genomes, {args.workers} workers", flush=True)

    R, B, I = [], [], []
    with ProcessPoolExecutor(args.workers, initializer=_init,
                             initargs=(pos_d, neg_d)) as ex:
        futs = {ex.submit(one_genome, a, args.redo): a for a in accs}
        for i, f in enumerate(as_completed(futs), 1):
            acc = futs[f]
            try:
                r, b, iv, dt = f.result()
                R.append(r); B.append(b); I.append(iv)
                f1 = r[r.method == "SSG-LUGIA-F"].f1.iloc[0]
                print(f"[{i}/{len(accs)}] {acc}  F1={f1:6.2f}  {dt:5.1f}s", flush=True)
            except Exception:                                   # noqa: BLE001
                print(f"[{i}/{len(accs)}] {acc}  FAILED", flush=True)
                traceback.print_exc()

    if not R:
        return
    pd.concat(R, ignore_index=True).to_csv(RESULTS / "ssg_lugia_per_genome.csv", index=False)
    pd.concat(B, ignore_index=True).to_csv(RESULTS / "ssg_lugia_boundaries.csv", index=False)
    pd.concat(I, ignore_index=True).to_csv(RESULTS / "ssg_lugia_intervals.csv", index=False)
    df = pd.concat(R, ignore_index=True)
    print(df.groupby("method")[["precision", "recall", "f1", "mabe"]].mean().round(2))


if __name__ == "__main__":
    main()
