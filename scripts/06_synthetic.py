#!/usr/bin/env python
"""Controlled synthetic-insertion experiment.

Splices tile-aligned blocks of a donor chromosome into a host chromosome at known
coordinates and measures, for each feature representation, how detection sensitivity and
boundary accuracy vary with

* insert length      (5, 10, 20, 40, 80 kb),
* donor-host compositional distance (tetranucleotide profile distance, |dGC|),
* replicate (different random placements).

Ground truth here is exact, so unlike IslandPick it can settle boundary questions.
Because inserts are tile-aligned, the chimera's tile embeddings are assembled by
concatenating precomputed host and donor tiles - no GPU work is needed per condition.
Compositional features *are* recomputed on the real chimeric sequence, so its 10 kb windows
see the junctions exactly as they would in a real analysis.

    python scripts/06_synthetic.py --workers 3
"""
from __future__ import annotations
import argparse, json, os, sys, traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _common import intervals_for, load_benchmark                        # noqa: E402
from sagegi.evaluate import boundary_errors, evaluate_genome, summarise_boundaries  # noqa: E402
from sagegi.features import extract_features                             # noqa: E402
from sagegi.genome import encode, read_fasta                             # noqa: E402
from sagegi.paths import EMB, GENOMES, RESULTS, WORK                     # noqa: E402
from sagegi.sage import SageConfig, run_sage_gi                          # noqa: E402
from sagegi.store import load_tiles                                      # noqa: E402
from sagegi.synthetic import (assemble_tiles, composition_distance,      # noqa: E402
                              make_chimera, negative_regions, shift_host_intervals)

INSERT_LENS = [5000, 10000, 20000, 40000, 80000]
N_INSERTS = 6
SEEDS = [0, 1]
CFG: dict = {}
#: One CSV per finished condition. The sweep is long enough that losing it to an
#: interruption is a real cost, and per-condition files also make the partial result
#: readable while the run is still in progress.
CACHE = WORK / "synth_cache"


def cache_path(host: str, donor: str, insert_len: int, seed: int, cfg_key: str) -> Path:
    return CACHE / f"{host}__{donor}__{insert_len}__{seed}__{cfg_key}.csv"


def _init(cfgs):
    global CFG
    CFG = cfgs
    for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ[v] = "1"


def one_condition(host: str, donor: str, insert_len: int, seed: int, emb_dirs: dict,
                  host_islands, cfg_key: str = ""):
    cp = cache_path(host, donor, insert_len, seed, cfg_key)
    if cp.exists():
        return pd.read_csv(cp)
    host_seq = read_fasta(GENOMES / f"{host}.fna")
    donor_seq = read_fasta(GENOMES / f"{donor}.fna")
    ch = make_chimera(host_seq, donor_seq, insert_len, N_INSERTS, tile=1000,
                      avoid=host_islands, rng=np.random.default_rng(seed))

    code = encode(ch.seq)
    comp = extract_features(code, dict(w=10000, dw=100, karlin_mode="normalized",
                                       pca_dn=2, pca_amino_acid=2, pca_kmer4=2,
                                       entropy_features=True))
    dist = composition_distance(encode(host_seq), encode(donor_seq))

    pos = ch.inserts
    neg = negative_regions(ch, shift_host_intervals(ch, host_islands), flank=20000)
    glen = len(ch.seq)

    tiles = {}
    for tag, d in emb_dirs.items():
        # Chimera assembly needs host and donor tiles in one coherent basis, so
        # this is the one place that asks for full-dimensional embeddings; every
        # other caller keeps the default compact form and is unaffected.
        hz, tile, _, is_pcs = load_tiles(EMB / d / f"{host}.npz", prefer_full=True)
        dz, _, _, _ = load_tiles(EMB / d / f"{donor}.npz", prefer_full=True)
        tiles[tag] = (assemble_tiles(ch, hz, dz, host_is_pcs=is_pcs), tile, is_pcs)

    rows = []
    for name, spec in CFG.items():
        spec = dict(spec)
        tag = spec.pop("embedding", None)
        cfg = SageConfig(**spec)
        if cfg.use_embedding and tag not in tiles:
            continue
        Z, tile, is_pcs = tiles[tag] if cfg.use_embedding else (None, 1000, False)
        try:
            res = run_sage_gi(glen, cfg, emb=Z, tile=tile, comp_X=comp.X,
                              comp_starts=comp.starts, emb_is_pcs=is_pcs)
        except Exception:                                    # noqa: BLE001
            traceback.print_exc()
            continue
        m = evaluate_genome(res.islands, pos, neg, glen)
        be = boundary_errors(res.islands, pos)
        b = summarise_boundaries(be)
        rows.append(dict(host=host, donor=donor, insert_len=insert_len, seed=seed,
                         method=name, n_inserts=len(pos), n_pred=len(res.islands),
                         detected=int(be.matched.sum()),
                         detection_rate=100 * float(be.matched.mean()), **m, **b, **dist))
    out = pd.DataFrame(rows)
    CACHE.mkdir(parents=True, exist_ok=True)
    out.to_csv(cp, index=False)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--configs", default="configs/methods_synthetic.json")
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--hosts", nargs="*", default=None)
    ap.add_argument("--donors", nargs="*", default=None)
    args = ap.parse_args()

    cfgs = json.loads(Path(args.configs).read_text())
    emb_dirs = {t: t for t in {c["embedding"] for c in cfgs.values() if c.get("embedding")}}
    ready = set.intersection(*[{p.stem for p in (EMB / d).glob("*.npz")}
                              for d in emb_dirs.values()]) if emb_dirs else set()
    pos_df, _ = load_benchmark()

    hosts = args.hosts or sorted(ready & set(pos_df.accession))[:2]
    donors = args.donors or [a for a in sorted(ready) if a not in hosts]
    if not hosts or not donors:
        print(f"need embeddings for at least one host and one donor "
              f"(ready: {sorted(ready)})")
        return
    print(f"hosts={hosts}\ndonors={donors}", flush=True)

    isl = {h: intervals_for(pos_df, h) for h in hosts}
    cfg_key = __import__("hashlib").sha1(
        Path(args.configs).read_text().encode()).hexdigest()[:10]
    jobs = list(product(hosts, donors, INSERT_LENS, SEEDS))
    CACHE.mkdir(parents=True, exist_ok=True)
    done = sum(cache_path(h, d, L, sd, cfg_key).exists() for h, d, L, sd in jobs)
    print(f"{len(jobs)} conditions x {len(cfgs)} methods"
          + (f" ({done} already cached)" if done else ""), flush=True)

    frames = []
    with ProcessPoolExecutor(args.workers, initializer=_init, initargs=(cfgs,)) as ex:
        futs = {ex.submit(one_condition, h, d, L, s, emb_dirs, isl[h], cfg_key): (h, d, L, s)
                for h, d, L, s in jobs}
        for i, f in enumerate(as_completed(futs), 1):
            key = futs[f]
            try:
                frames.append(f.result())
                print(f"[{i}/{len(jobs)}] {key}", flush=True)
            except Exception:                                # noqa: BLE001
                print(f"[{i}/{len(jobs)}] {key} FAILED", flush=True)
                traceback.print_exc()

    if not frames:
        raise SystemExit(
            f"no conditions produced results out of {len(jobs)} -- every task failed or the "
            f"worker pool never started; nothing was written.")
    df = pd.concat(frames, ignore_index=True)
    df.to_csv(RESULTS / "synthetic_insertions.csv", index=False)
    print(df.groupby(["method", "insert_len"])[["detection_rate", "precision", "mabe"]]
          .mean().round(2).to_string())


if __name__ == "__main__":
    main()
