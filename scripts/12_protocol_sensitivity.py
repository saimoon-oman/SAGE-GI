#!/usr/bin/env python
"""How much does the IslandPick protocol's partial coverage flatter each method?

The protocol scores only the curated island and backbone regions -- about 13% of an average
chromosome. A prediction falling in the other 87% is neither rewarded nor penalised, so a
method that calls a lot of sequence is under-penalised. That is a property of the benchmark
rather than of any method, but it bounds how the numbers should be read, and it is the
reason a method calling 27% of a chromosome can post a high F1.

This script recomputes precision under two scorings:

* **official** -- the published protocol: TP = |Y & P|, FP = |Y & N|, ignoring everything
  outside the curated sets;
* **strict**   -- every base outside a curated island counts as a negative, so predictions
  landing in unscored sequence are penalised.

Neither is "correct". The official scoring is what the literature reports and what makes our
numbers comparable with published ones; the strict scoring is a lower bound that answers
"what if the curated islands really are all the islands there are?". The truth lies between,
because the unscored majority certainly contains real islands nobody has curated. Reporting
the gap is what stops a reader over-reading either number.

    python scripts/12_protocol_sensitivity.py
"""
from __future__ import annotations
import json, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _common import intervals_for, load_benchmark                        # noqa: E402
from sagegi.paths import BENCH, RESULTS                                  # noqa: E402

GLEN = json.loads((BENCH / "genome_lengths.json").read_text())


def mask_from(intervals, n: int) -> np.ndarray:
    m = np.zeros(n, dtype=bool)
    for s, e in intervals:
        m[max(0, s):min(n, e + 1)] = True
    return m


def load_intervals(path: Path, method: str | None = None):
    if not path.exists():
        return None
    d = pd.read_csv(path)
    if method is not None:
        d = d[d.method == method]
    return d


def score(pred_iv, pos_iv, neg_iv, L):
    """Precision under the official partial-coverage protocol and under strict scoring."""
    P = mask_from(pos_iv, L)
    N = mask_from(neg_iv, L)
    Y = mask_from(pred_iv, L)
    tp = int((Y & P).sum())
    fp_official = int((Y & N).sum())
    fp_strict = int((Y & ~P).sum())          # everything outside a curated island
    prec_off = 100 * tp / (tp + fp_official) if tp + fp_official else np.nan
    prec_str = 100 * tp / (tp + fp_strict) if tp + fp_strict else np.nan
    return prec_off, prec_str, int(Y.sum()) / L


def main() -> None:
    pos, neg = load_benchmark()
    ours = load_intervals(RESULTS / "sage_gi_intervals.csv")
    ssg = load_intervals(RESULTS / "ssg_lugia_intervals.csv")
    ti = load_intervals(RESULTS / "treasureisland_intervals.csv")

    sources = []
    if ours is not None:
        for m in ["DNABERT-S-only", "SAGE-GI", "SAGE-GI-P", "Comp-only+center"]:
            if m in set(ours.method):
                sources.append((m, ours[ours.method == m]))
    if ssg is not None:
        for m in ["SSG-LUGIA-F", "SSG-LUGIA-P"]:
            if m in set(ssg.method):
                sources.append((m, ssg[ssg.method == m]))
    if ti is not None and len(ti):
        sources.append(("TreasureIsland", ti))

    if not sources:
        raise SystemExit("no interval files found -- run the detector first")

    rows = []
    for name, d in sources:
        for acc, g in d.groupby("accession"):
            if acc not in GLEN:
                continue
            L = int(GLEN[acc])
            pred = [(int(s) - 1, int(e) - 1) for s, e in zip(g.start, g.end)]
            po, ps, frac = score(pred, intervals_for(pos, acc), intervals_for(neg, acc), L)
            rows.append(dict(method=name, accession=acc, precision_official=po,
                             precision_strict=ps, called_frac=100 * frac))

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "protocol_sensitivity.csv", index=False)

    agg = (df.groupby("method")[["precision_official", "precision_strict", "called_frac"]]
             .mean().round(2))
    agg["drop"] = (agg.precision_official - agg.precision_strict).round(2)
    agg = agg.sort_values("precision_official", ascending=False)
    print(agg.to_string())
    print()
    print("`drop` is how much of each method's reported precision depends on the 87% of")
    print("sequence the protocol never scores. A larger drop means a method places more of")
    print("its predictions where the benchmark cannot see them.")


if __name__ == "__main__":
    main()
