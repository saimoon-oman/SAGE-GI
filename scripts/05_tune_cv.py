#!/usr/bin/env python
"""Hyper-parameter sweep with honest, grouped cross-validation.

Two things make this affordable.

1. **Refit only when you must.**  Within one genome the expensive step is the
   Minimum-Covariance-Determinant fit, which depends only on the feature matrix and the two
   contamination levels.  Everything after it - median smoothing, thresholding, projection,
   minimum-island filtering - is cheap.  So we fit once per ``(genome, c1, c2)`` and then
   sweep the post-processing parameters for free.

2. **Score the whole grid once, then cross-validate on paper.**  We build the complete
   ``genome x configuration`` metric table a single time.  Because model selection only ever
   *reads* that table, k-fold cross-validation costs nothing extra: for each fold we choose
   the configuration that maximises the objective on the training genomes and read off its
   value on the held-out genomes.

SSG-LUGIA tuned on repeated random samples of 20 of the 118 genomes and reported the
resulting numbers on the same benchmark.  We report both: the same in-sample protocol, for
comparability, *and* grouped 5-fold cross-validated numbers, which is what actually
estimates generalisation to a new genome.

    python scripts/05_tune_cv.py --channel sage --workers 3
"""
from __future__ import annotations
import argparse, itertools, json, os, sys, time, traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _common import intervals_for, load_benchmark                        # noqa: E402
from sagegi.anomaly import two_stage_envelope                            # noqa: E402
from sagegi.embed import pool_windows                                    # noqa: E402
from sagegi.evaluate import (boundary_errors, evaluate_genome,           # noqa: E402
                             summarise_boundaries)
from sagegi.paths import BENCH, EMB, RESULTS, WORK                       # noqa: E402
from sagegi.sage import SageConfig, refine_boundaries, tile_matrix       # noqa: E402
from sagegi.postprocess import (drop_short, median_smooth, merge_close,  # noqa: E402
                                project_to_nucleotides, runs)
from sagegi.store import load_tiles                                      # noqa: E402
from scipy.ndimage import median_filter                                  # noqa: E402

COMP = WORK / "comp"
SWEEP = WORK / "sweep"
SWEEP.mkdir(parents=True, exist_ok=True)
GLEN = json.loads((BENCH / "genome_lengths.json").read_text())

# ------------------------------------------------------------------------------- the grid
CONTAM = [(0.05, 0.05), (0.075, 0.075), (0.10, 0.05), (0.10, 0.10),
          (0.15, 0.05), (0.15, 0.10), (0.15, 0.20), (0.20, 0.10), (0.20, 0.25)]
FILTER_LEN = [200, 400, 600]
MIN_LEN = [8000, 10000]

POS: dict = {}
NEG: dict = {}
ARGS: dict = {}


def _init(pos_d, neg_d, args_d):
    global POS, NEG, ARGS
    POS, NEG, ARGS = pos_d, neg_d, args_d
    for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ[v] = "1"


def build_features(acc: str, channel: str, emb_dir: str, n_components: int,
                   detrend_bp: int):
    """Feature matrix, window starts and (for refinement) the tile-level signal."""
    z = np.load(COMP / f"{acc}.npz")
    comp_X, starts = z["X"].astype(np.float32), z["starts"]
    glen = int(GLEN[acc])
    if channel == "comp":
        return comp_X.astype(np.float64), starts, glen, None, 1000

    Z, tile, _, is_pcs = load_tiles(EMB / emb_dir / f"{acc}.npz")
    cfg = SageConfig(n_components=n_components, detrend_bp=detrend_bp)
    Zt = tile_matrix(Z, cfg, tile, is_pcs=is_pcs)
    pooled = pool_windows(Zt, tile, tile, starts, cfg.w)
    if channel == "emb":
        X = pooled.astype(np.float64)
    elif channel == "sage":
        n = min(len(pooled), len(comp_X))
        X = np.column_stack([pooled[:n], comp_X[:n]]).astype(np.float64)
        starts = starts[:n]
    else:
        raise ValueError(channel)
    return X, starts, glen, Zt, tile


def tile_signal(Zt: np.ndarray) -> np.ndarray:
    centre = np.median(Zt, axis=0)
    mad = np.median(np.abs(Zt - centre), axis=0) * 1.4826 + 1e-8
    return median_filter(np.linalg.norm((Zt - centre) / mad, axis=1), size=5, mode="nearest")


def sweep_genome(acc: str):
    a = ARGS
    cache = SWEEP / f"{a['tag']}__{acc}.csv"
    if cache.exists() and not a["redo"]:
        return pd.read_csv(cache)

    X, starts, glen, Zt, tile = build_features(acc, a["channel"], a["emb_dir"],
                                               a["n_components"], a["detrend_bp"])
    sig = tile_signal(Zt) if Zt is not None else None
    rows = []
    for c1, c2 in CONTAM:
        params = dict(contamination_model1=c1, contamination_model2=c2,
                      support_fraction_model1=0.75, support_fraction_model2=0.9,
                      fit_stride=a["fit_stride"])
        try:
            _, scores = two_stage_envelope(X, params)
        except Exception:                                     # noqa: BLE001
            traceback.print_exc()
            continue
        for flen in FILTER_LEN:
            sm = median_smooth(scores, flen)
            track = project_to_nucleotides(sm < 0, glen, 10000, 100, a["assign"])
            base = merge_close(runs(track), 5000)
            for mlen in MIN_LEN:
                isl = drop_short(base, mlen)
                for refine in ([False, True] if sig is not None else [False]):
                    if refine and isl:
                        cfg = SageConfig()
                        ri = refine_boundaries(isl, sig, tile, cfg, glen)
                        ri = drop_short(merge_close(sorted(ri), 5000), mlen)
                    else:
                        ri = isl
                    m = evaluate_genome(ri, POS[acc], NEG[acc], glen)
                    b = summarise_boundaries(boundary_errors(ri, POS[acc]))
                    rows.append(dict(accession=acc, c1=c1, c2=c2, filter_len=flen,
                                     min_len=mlen, refine=refine, n_pred=len(ri), **m,
                                     mabe=b["mabe"], median_jaccard=b["median_jaccard"]))
    df = pd.DataFrame(rows)
    df.to_csv(cache, index=False)
    return df


# ----------------------------------------------------------------------------- selection
CONFIG_KEYS = ["c1", "c2", "filter_len", "min_len", "refine"]


def select(train: pd.DataFrame, objective: str) -> tuple:
    """Pick the configuration maximising an objective, averaged over training genomes."""
    g = train.groupby(CONFIG_KEYS)[["precision", "recall", "f1"]].mean().reset_index()
    if objective == "f1":
        score = g.f1
    elif objective == "precision":                 # best precision at recall >= 50 %
        ok = g.recall >= 50
        score = np.where(ok, g.precision, g.precision - 1000)
    elif objective == "recall":                    # best recall at precision >= 35 %
        ok = g.precision >= 35
        score = np.where(ok, g.recall, g.recall - 1000)
    else:
        raise ValueError(objective)
    best = g.iloc[int(np.argmax(score))]
    return tuple(best[k] for k in CONFIG_KEYS)


def cross_validate(df: pd.DataFrame, objective: str, k: int = 5, seed: int = 0):
    accs = np.array(sorted(df.accession.unique()))
    rng = np.random.default_rng(seed)
    fold = rng.permutation(len(accs)) % k
    out, picks = [], []
    for f in range(k):
        te = set(accs[fold == f]); tr = set(accs[fold != f])
        conf = select(df[df.accession.isin(tr)], objective)
        sel = df[df.accession.isin(te)]
        for key, val in zip(CONFIG_KEYS, conf):
            sel = sel[sel[key] == val]
        out.append(sel.assign(fold=f))
        picks.append(dict(fold=f, **dict(zip(CONFIG_KEYS, conf))))
    return pd.concat(out, ignore_index=True), pd.DataFrame(picks)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--channel", default="sage", choices=["comp", "emb", "sage"])
    ap.add_argument("--emb-dir", default="DNABERT-S_t1000")
    ap.add_argument("--n-components", type=int, default=24)
    ap.add_argument("--detrend-bp", type=int, default=250000)
    ap.add_argument("--assign", default="center")
    ap.add_argument("--fit-stride", type=int, default=10)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--tag", default=None)
    ap.add_argument("--redo", action="store_true")
    ap.add_argument("--accessions", nargs="*", default=None)
    args = ap.parse_args()
    tag = args.tag or f"{args.channel}_{args.emb_dir if args.channel != 'comp' else 'x'}"

    pos, neg = load_benchmark()
    have_comp = {p.stem for p in COMP.glob("*.npz")}
    pool = have_comp if args.channel == "comp" else \
        have_comp & {p.stem for p in (EMB / args.emb_dir).glob("*.npz")}
    accs = args.accessions or sorted(pool & set(pos.accession))
    if not accs:
        print("no genomes ready for this channel"); return

    pos_d = {a: intervals_for(pos, a) for a in accs}
    neg_d = {a: intervals_for(neg, a) for a in accs}
    grid_n = len(CONTAM) * len(FILTER_LEN) * len(MIN_LEN) * (2 if args.channel != "comp" else 1)
    print(f"tag={tag}  {len(accs)} genomes x {grid_n} configurations", flush=True)

    a_d = dict(channel=args.channel, emb_dir=args.emb_dir, n_components=args.n_components,
               detrend_bp=args.detrend_bp, assign=args.assign, fit_stride=args.fit_stride,
               tag=tag, redo=args.redo)

    frames, t0 = [], time.time()
    with ProcessPoolExecutor(args.workers, initializer=_init,
                             initargs=(pos_d, neg_d, a_d)) as ex:
        futs = {ex.submit(sweep_genome, x): x for x in accs}
        for i, f in enumerate(as_completed(futs), 1):
            try:
                frames.append(f.result())
                print(f"[{i}/{len(accs)}] {futs[f]}  ({(time.time()-t0)/60:.1f} min)", flush=True)
            except Exception:                                 # noqa: BLE001
                print(f"[{i}/{len(accs)}] {futs[f]} FAILED", flush=True)
                traceback.print_exc()

    df = pd.concat(frames, ignore_index=True)
    df.to_csv(RESULTS / f"sweep_{tag}.csv", index=False)

    report = []
    for obj in ("f1", "precision", "recall"):
        insample_cfg = select(df, obj)
        sel = df.copy()
        for key, val in zip(CONFIG_KEYS, insample_cfg):
            sel = sel[sel[key] == val]
        cv, picks = cross_validate(df, obj)
        report.append(dict(
            objective=obj, config=dict(zip(CONFIG_KEYS, [float(x) for x in insample_cfg])),
            in_sample={k: round(float(sel[k].mean()), 3)
                       for k in ("precision", "recall", "f1", "mabe")},
            cv_heldout={k: round(float(cv[k].mean()), 3)
                        for k in ("precision", "recall", "f1", "mabe")},
            fold_choices=picks.to_dict("records")))
        print(f"\n[{obj}] chosen {dict(zip(CONFIG_KEYS, insample_cfg))}")
        print(f"   in-sample : {report[-1]['in_sample']}")
        print(f"   5-fold CV : {report[-1]['cv_heldout']}")

    (RESULTS / f"tuning_{tag}.json").write_text(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
