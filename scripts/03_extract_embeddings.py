#!/usr/bin/env python
"""Extract frozen tile embeddings for a set of genomes with DNABERT-S or DNABERT-2.

Resumable: one ``.npz`` per (model, tile length, genome); existing files are skipped.
Runs unchanged on CPU and on a CUDA GPU (fp16 is used automatically on CUDA).

    python scripts/03_extract_embeddings.py --model DNABERT-S --tile 1000
    python scripts/03_extract_embeddings.py --model DNABERT-2 --tile 1000 --accessions NC_003198.1
"""
from __future__ import annotations
import argparse, json, sys, time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _common import available_accessions                       # noqa: E402
from sagegi.embed import EmbedConfig, embed_genome, load_model      # noqa: E402
from sagegi.store import save_reduced                              # noqa: E402
from sagegi.genome import read_fasta                           # noqa: E402
from sagegi.paths import EMB, GENOMES, MODELS, BENCH           # noqa: E402

PRIORITY = ["NC_003198.1", "NC_002935.2", "NC_011770.1"]        # case studies go first


def order_accessions(accs: list[str]) -> list[str]:
    """Case studies first, then the rest smallest-genome-first so partial runs are useful."""
    rest = [a for a in accs if a not in PRIORITY]
    sizes = {a: (GENOMES / f"{a}.fna").stat().st_size for a in rest}
    return [a for a in PRIORITY if a in accs] + sorted(rest, key=sizes.get)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="DNABERT-S", choices=["DNABERT-S", "DNABERT-2"])
    ap.add_argument("--model-path", default=None)
    ap.add_argument("--tile", type=int, default=1000)
    ap.add_argument("--batch", type=int, default=16)
    ap.add_argument("--device", default="auto")
    ap.add_argument("--accessions", nargs="*", default=None)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", default=None)
    ap.add_argument("--pca-k", type=int, default=96)
    ap.add_argument("--keep-raw", nargs="*", default=[])
    args = ap.parse_args()

    cfg = EmbedConfig(model_path=args.model_path or str(MODELS / args.model),
                      tile=args.tile, batch=args.batch, device=args.device)
    outdir = Path(args.out) if args.out else EMB / f"{args.model}_t{args.tile}"
    outdir.mkdir(parents=True, exist_ok=True)

    accs = args.accessions or order_accessions(available_accessions())
    todo = [a for a in accs if not (outdir / f"{a}.npz").exists()]
    if args.limit:
        todo = todo[:args.limit]
    print(f"{args.model} tile={args.tile} -> {outdir}", flush=True)
    print(f"{len(accs)} genomes known, {len(todo)} still to do", flush=True)
    if not todo:
        return

    tok, model, dev = load_model(cfg)
    print(f"device={dev} dtype={next(model.parameters()).dtype}", flush=True)

    total_bp = 0
    t_start = time.time()
    for i, acc in enumerate(todo, 1):
        seq = read_fasta(GENOMES / f"{acc}.fna")
        t0 = time.time()
        emb, starts = embed_genome(seq, cfg, tok, model, dev)
        dt = time.time() - t0
        if acc in set(args.keep_raw):
            np.savez_compressed(outdir / f"{acc}.npz", emb=emb, starts=starts, tile=cfg.tile,
                                step=cfg.step, genome_len=len(seq), model=args.model)
        else:
            save_reduced(outdir / f"{acc}.npz", emb, starts, cfg.tile, cfg.step,
                         len(seq), args.model, k=args.pca_k)
        total_bp += len(seq)
        rate = total_bp / max(1e-9, time.time() - t_start)
        eta = sum((GENOMES / f"{a}.fna").stat().st_size for a in todo[i:]) / max(rate, 1)
        print(f"[{i}/{len(todo)}] {acc} {len(seq):>9,} bp  {len(starts):>6,} tiles  "
              f"{dt:7.1f}s  {len(seq)/dt:7.0f} bp/s  ETA {eta/3600:5.2f} h", flush=True)


if __name__ == "__main__":
    main()
