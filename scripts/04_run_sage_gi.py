#!/usr/bin/env python
"""Run SAGE-GI (and its ablations) over every genome that has tile embeddings.

Configurations are named entries in a JSON file (see ``configs/methods.json``); each is a
partial :class:`sagegi.sage.SageConfig` plus an ``embedding`` key naming the embedding
directory to use (or ``null`` for a composition-only run).

    python scripts/04_run_sage_gi.py --configs configs/methods.json --workers 3
"""
from __future__ import annotations
import argparse, json, os, sys, time, traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _common import intervals_for, load_benchmark                         # noqa: E402
from sagegi.paths import BENCH, EMB, RESULTS, WORK                        # noqa: E402
from sagegi.sage import SageConfig, run_sage_gi                           # noqa: E402
from sagegi.store import load_tiles                                       # noqa: E402
from sagegi.evaluate import (boundary_errors, evaluate_genome,            # noqa: E402
                             summarise_boundaries)

COMP = WORK / "comp"
#: Per-genome results are cached so an interrupted run resumes instead of starting over.
#: The cache key includes a digest of the configuration file, so editing a configuration
#: invalidates it rather than silently reusing stale predictions.
CACHE = WORK / "sage_cache"
CACHE.mkdir(parents=True, exist_ok=True)
GLEN = json.loads((BENCH / "genome_lengths.json").read_text())
CFG_KEY = ""
POS: dict = {}
NEG: dict = {}
CFGS: dict = {}


def _init(pos_d, neg_d, cfgs, cfg_key):
    global POS, NEG, CFGS, CFG_KEY
    POS, NEG, CFGS, CFG_KEY = pos_d, neg_d, cfgs, cfg_key
    for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ[v] = "1"


def load_inputs(acc: str, emb_dirs: set[str]):
    """Compositional features plus each requested embedding matrix for one genome."""
    z = np.load(COMP / f"{acc}.npz")
    comp_X, comp_starts = z["X"].astype(np.float32), z["starts"]
    tiles = {}
    for d in emb_dirs:
        p = EMB / d / f"{acc}.npz"
        if p.exists():
            tiles[d] = load_tiles(p)
    return comp_X, comp_starts, tiles


def stem_path(stem, ext: str) -> Path:
    """``Path.with_suffix`` cannot be used here: an accession such as ``NC_000913.2``
    makes ``.2__<key>`` look like the extension, so the version *and* the configuration
    digest were being stripped from the cache filename -- which silently defeated the
    invalidation this cache documents."""
    return stem.parent / (stem.name + ext)


def runnable_configs(acc: str) -> set[str]:
    """Configuration names whose inputs exist for this genome *right now*.

    A configuration is skipped when the embeddings it needs have not been extracted
    yet, so what a cache entry contains depends on when it was written. Comparing
    against this set is what stops a genome cached before an embedding archive arrived
    from being reused afterwards, which would otherwise average one model over fewer
    genomes than another and silently bias the comparison.
    """
    out = set()
    for name, spec in CFGS.items():
        d = spec.get("embedding")
        if d is None or (EMB / d / f"{acc}.npz").exists():
            out.add(name)
    return out


def one_genome(acc: str):
    stem = CACHE / f"{acc}__{CFG_KEY}"
    want = runnable_configs(acc)
    todo, old = want, None
    if (stem_path(stem, ".rows.csv")).exists():
        cached = pd.read_csv(stem_path(stem, ".rows.csv"))
        marker = stem_path(stem, ".attempted.json")
        # A configuration that raised at run time must not force a recompute on every
        # later invocation, so we compare against what was *attempted*, not what
        # succeeded. Caches written before this marker existed fall back to their rows.
        attempted = (set(json.loads(marker.read_text())) if marker.exists()
                     else set(cached["method"]) if len(cached) else set())
        todo = want - attempted
        cached_bnd = (pd.read_csv(stem_path(stem, ".bnd.csv"))
                      if stem_path(stem, ".bnd.csv").exists() else pd.DataFrame())
        cached_iv = pd.read_csv(stem_path(stem, ".iv.csv"))
        if not todo:
            return cached, cached_bnd, cached_iv
        # Configurations are independent given the same inputs, so running only the new
        # ones and appending is equivalent to re-running all of them -- and avoids
        # redoing seventeen configurations to add two.
        old = (cached, cached_bnd, cached_iv, attempted)
        print(f"  {acc}: adding {len(todo)} configuration(s) "
              f"({', '.join(sorted(todo))}) to {len(attempted)} cached", flush=True)

    emb_dirs = {c["embedding"] for n, c in CFGS.items()
                if n in todo and c.get("embedding")}
    comp_X, comp_starts, tiles = load_inputs(acc, emb_dirs)
    genome_len = int(GLEN[acc])

    rows, brows, irows = [], [], []
    for name, spec in CFGS.items():
        if name not in todo:
            continue
        spec = dict(spec)
        emb_dir = spec.pop("embedding", None)
        cfg = SageConfig(**spec)
        if cfg.use_embedding:
            if emb_dir not in tiles:
                continue
            Z, tile, glen, is_pcs = tiles[emb_dir]
        else:
            Z, tile, is_pcs = None, 1000, False
            glen = genome_len
        t0 = time.time()
        try:
            res = run_sage_gi(glen, cfg, emb=Z, tile=tile, comp_X=comp_X,
                              comp_starts=comp_starts, emb_is_pcs=is_pcs)
        except Exception:                                       # noqa: BLE001
            traceback.print_exc()
            continue
        m = evaluate_genome(res.islands, POS[acc], NEG[acc], glen)
        be = boundary_errors(res.islands, POS[acc])
        rows.append(dict(accession=acc, method=name, genome_len=glen,
                         n_pred=len(res.islands),
                         pred_bp=int(sum(e - s + 1 for s, e in res.islands)),
                         seconds=round(time.time() - t0, 1), **m,
                         **summarise_boundaries(be)))
        be.insert(0, "method", name); be.insert(0, "accession", acc)
        brows.append(be)
        for s, e in res.islands:
            irows.append(dict(accession=acc, method=name, start=s + 1, end=e + 1))
    R, B, I = (pd.DataFrame(rows),
               pd.concat(brows, ignore_index=True) if brows else pd.DataFrame(),
               pd.DataFrame(irows))
    if old is not None:
        cached, cached_bnd, cached_iv, attempted = old
        R = pd.concat([cached, R], ignore_index=True) if len(R) else cached
        B = (pd.concat([cached_bnd, B], ignore_index=True)
             if len(B) and len(cached_bnd) else (B if len(B) else cached_bnd))
        I = pd.concat([cached_iv, I], ignore_index=True) if len(I) else cached_iv
        want = want | attempted
    R.to_csv(stem_path(stem, ".rows.csv"), index=False)
    I.to_csv(stem_path(stem, ".iv.csv"), index=False)
    if len(B):
        B.to_csv(stem_path(stem, ".bnd.csv"), index=False)
    # Record what could have run, so a later invocation can tell a configuration that
    # was absent from one that was tried and failed.
    stem_path(stem, ".attempted.json").write_text(json.dumps(sorted(want)))
    return R, B, I


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--configs", default="configs/methods.json")
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--tag", default="sage_gi")
    ap.add_argument("--accessions", nargs="*", default=None)
    args = ap.parse_args()

    cfg_text = Path(args.configs).read_text()
    cfgs = json.loads(cfg_text)
    cfg_key = __import__("hashlib").sha1(cfg_text.encode()).hexdigest()[:10]
    pos, neg = load_benchmark()

    # A genome is runnable if its compositional features are cached. Individual
    # configurations are skipped per genome when the embeddings they need are absent
    # (see one_genome), so adding a configuration for a model that has not been extracted
    # yet must not disable the whole run -- which is exactly what intersecting every
    # embedding directory used to do.
    emb_dirs = {c["embedding"] for c in cfgs.values() if c.get("embedding")}
    have_comp = {p.stem for p in COMP.glob("*.npz")}
    accs = args.accessions or sorted(have_comp & set(pos.accession))
    if not accs:
        print("nothing to do: no cached compositional features "
              "(run scripts/02_run_ssg_lugia.py first)")
        return
    for d in sorted(emb_dirs):
        n = len(list((EMB / d).glob("*.npz")))
        if n == 0:
            print(f"  note: no embeddings under {d} -- configurations using it are skipped")
        elif n < len(accs):
            print(f"  note: {d} covers {n}/{len(accs)} genomes")

    pos_d = {a: intervals_for(pos, a) for a in accs}
    neg_d = {a: intervals_for(neg, a) for a in accs}
    print(f"{len(accs)} genomes x {len(cfgs)} configurations", flush=True)

    R, B, I = [], [], []
    done = len(list(CACHE.glob(f"*__{cfg_key}.rows.csv")))
    if done:
        print(f"resuming: {done} genomes already cached for this configuration", flush=True)
    with ProcessPoolExecutor(args.workers, initializer=_init,
                             initargs=(pos_d, neg_d, cfgs, cfg_key)) as ex:
        futs = {ex.submit(one_genome, a): a for a in accs}
        for i, f in enumerate(as_completed(futs), 1):
            acc = futs[f]
            try:
                r, b, iv = f.result()
                R.append(r); I.append(iv)
                if len(b):
                    B.append(b)
                best = r.sort_values("f1", ascending=False).head(1)
                print(f"[{i}/{len(accs)}] {acc}  best={best.method.iloc[0]} "
                      f"F1={best.f1.iloc[0]:.2f}", flush=True)
            except Exception:                                   # noqa: BLE001
                print(f"[{i}/{len(accs)}] {acc} FAILED", flush=True)
                traceback.print_exc()

    # Exiting 0 having written nothing looks identical to success in a log, and a caller
    # that only checks the exit status will happily report results that were never
    # computed. Fail loudly instead.
    if not R:
        raise SystemExit(
            f"no results produced for any of {len(accs)} genomes -- every task failed or "
            f"the worker pool never started. Nothing was written to "
            f"{RESULTS / f'{args.tag}_per_genome.csv'}.")
    df = pd.concat(R, ignore_index=True)
    df.to_csv(RESULTS / f"{args.tag}_per_genome.csv", index=False)
    pd.concat(I, ignore_index=True).to_csv(RESULTS / f"{args.tag}_intervals.csv", index=False)
    if B:
        pd.concat(B, ignore_index=True).to_csv(RESULTS / f"{args.tag}_boundaries.csv",
                                               index=False)
    print(df.groupby("method")[["precision", "recall", "f1", "mabe", "median_jaccard"]]
          .mean().round(2).to_string())


if __name__ == "__main__":
    main()
