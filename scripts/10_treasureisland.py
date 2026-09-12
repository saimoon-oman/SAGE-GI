#!/usr/bin/env python
"""Run TreasureIsland on the IslandPick benchmark and score it under our protocol.

TreasureIsland (Banerjee, Eulenstein & Friedberg, *Bioinformatics Advances* 2024) is the
closest published competitor: single-sequence, and replacing hand-crafted composition with a
learned Doc2Vec embedding rather than a pretrained foundation model.

It is **not** an unsupervised method, and the distinction matters for how its numbers are
read. Its *representation* is learned without labels, but the classifier on top is a
supervised SVM, shipped pre-trained (`svm_embedding_*` in the package). It was trained on the
Benbow compilation -- 1,329 islands from 145 genomes and 1,240 non-islands from 72 -- whose
main component is 104 genomes with 1,845 islands and 3,266 non-islands, that is, the revised
version of the IslandPick benchmark we evaluate on. Running the shipped model over these 118
chromosomes therefore scores a supervised classifier on data related to its own training set.

Its published numbers are separately incomparable with ours: it was evaluated on 20 genomes
as a classification of 566 pre-defined regions (413 islands against 153 non-islands), a
balanced decision over given intervals, where ours is a genome-wide scan that must also
decide *where* the intervals are. So we run it ourselves, on the same 118 chromosomes, score
it with the same code as every other method, and record the tool's own
`out_of_distribution` flag per genome so the training overlap can be measured rather than
argued about.

TreasureIsland pulls in its own pinned dependency stack, so it is installed in a separate
virtual environment and invoked as a subprocess.

Setup (once)::

    python -m venv $SAGEGI_WORK/ti_env
    $SAGEGI_WORK/ti_env/Scripts/python -m pip install treasureisland

Then::

    python scripts/10_treasureisland.py --workers 3
"""
from __future__ import annotations
import argparse, json, os, subprocess, sys, textwrap, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _common import intervals_for, load_benchmark                    # noqa: E402
from sagegi.evaluate import (boundary_errors, evaluate_genome,       # noqa: E402
                             summarise_boundaries)
from sagegi.paths import BENCH, GENOMES, RESULTS, WORK               # noqa: E402

TI_PY = WORK / "ti_env" / "Scripts" / "python.exe"
if not TI_PY.exists():                                   # POSIX layout
    TI_PY = WORK / "ti_env" / "bin" / "python"
CACHE = WORK / "ti_pred"
CACHE.mkdir(parents=True, exist_ok=True)
GLEN = json.loads((BENCH / "genome_lengths.json").read_text())

RUNNER = textwrap.dedent("""
    import json, sys, warnings
    warnings.filterwarnings('ignore')
    from treasureisland.Predictor import Predictor
    fasta, out_json, workdir = sys.argv[1], sys.argv[2], sys.argv[3]
    p = Predictor(fasta, workdir)
    pred = p.predict()
    rows, ood = [], bool(any(p.out_of_distribution))
    for key, regions in pred.items():
        for r in regions:
            rows.append([str(r[0]), int(float(r[1])), int(float(r[2])), float(r[3])])
    json.dump({'regions': rows, 'out_of_distribution': ood}, open(out_json, 'w'))
""")


def run_one(acc: str, redo: bool = False):
    out = CACHE / f"{acc}.json"
    if out.exists() and not redo:
        return json.loads(out.read_text())
    fasta = GENOMES / f"{acc}.fna"
    script = CACHE / "_runner.py"
    script.write_text(RUNNER)
    t = time.time()
    env = dict(os.environ, PYTHONIOENCODING="utf-8", KMP_DUPLICATE_LIB_OK="TRUE",
               OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    proc = subprocess.run([str(TI_PY), str(script), str(fasta), str(out),
                           str(CACHE / f"_work_{acc}")],
                          capture_output=True, text=True, encoding="utf-8",
                          errors="replace", env=env, timeout=7200)
    if not out.exists():
        raise RuntimeError(f"{acc}: TreasureIsland failed\n{proc.stderr[-1500:]}")
    d = json.loads(out.read_text())
    d["seconds"] = round(time.time() - t, 1)
    out.write_text(json.dumps(d))
    return d


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--redo", action="store_true")
    args = ap.parse_args()

    if not TI_PY.exists():
        print(f"TreasureIsland environment not found at {TI_PY}\n"
              "create it with:  python -m venv $SAGEGI_WORK/ti_env && "
              "$SAGEGI_WORK/ti_env/Scripts/python -m pip install treasureisland")
        return

    pos, neg = load_benchmark()
    accs = sorted(set(pos.accession) & {p.stem for p in GENOMES.glob("*.fna")})
    if args.limit:
        accs = accs[:args.limit]
    todo = [a for a in accs if args.redo or not (CACHE / f"{a}.json").exists()]
    print(f"{len(accs)} genomes, {len(todo)} to predict", flush=True)

    if todo:
        with ThreadPoolExecutor(args.workers) as ex:
            futs = {ex.submit(run_one, a, args.redo): a for a in todo}
            for i, f in enumerate(as_completed(futs), 1):
                acc = futs[f]
                try:
                    d = f.result()
                    print(f"[{i}/{len(todo)}] {acc}  {len(d['regions']):3d} regions  "
                          f"{d.get('seconds', 0):6.1f}s"
                          f"{'  [out of distribution]' if d['out_of_distribution'] else ''}",
                          flush=True)
                except Exception as exc:                      # noqa: BLE001
                    print(f"[{i}/{len(todo)}] {acc} FAILED: {exc}", flush=True)

    rows, ivrows = [], []
    for acc in accs:
        f = CACHE / f"{acc}.json"
        if not f.exists():
            continue
        d = json.loads(f.read_text())
        L = int(GLEN[acc])
        iv = [(max(0, s - 1), min(L - 1, e - 1)) for _, s, e, _ in d["regions"] if e > s]
        m = evaluate_genome(iv, intervals_for(pos, acc), intervals_for(neg, acc), L)
        b = summarise_boundaries(boundary_errors(iv, intervals_for(pos, acc)))
        rows.append(dict(accession=acc, method="TreasureIsland", genome_len=L,
                         n_pred=len(iv), pred_bp=int(sum(e - s + 1 for s, e in iv)),
                         out_of_distribution=d["out_of_distribution"], **m, **b))
        for s, e in iv:
            ivrows.append(dict(accession=acc, method="TreasureIsland",
                               start=s + 1, end=e + 1))

    if not rows:
        print("no predictions to score"); return
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "treasureisland_per_genome.csv", index=False)
    pd.DataFrame(ivrows).to_csv(RESULTS / "treasureisland_intervals.csv", index=False)

    print(f"\nTreasureIsland over {len(df)} genomes "
          f"({int(df.out_of_distribution.sum())} flagged out of distribution by the tool):")
    print(f"  precision {df.precision.mean():6.2f}   recall {df.recall.mean():6.2f}   "
          f"F1 {df.f1.mean():6.2f}   genome fraction called "
          f"{100 * (df.pred_bp / df.genome_len).mean():5.1f}%")
    ind = df[~df.out_of_distribution]
    if len(ind) and len(ind) < len(df):
        print(f"  in-distribution only ({len(ind)} genomes): "
              f"P {ind.precision.mean():.2f}  R {ind.recall.mean():.2f}  "
              f"F1 {ind.f1.mean():.2f}")


if __name__ == "__main__":
    main()
