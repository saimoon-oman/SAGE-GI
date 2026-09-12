r"""Evaluation under the IslandPick protocol, plus boundary-level metrics.

IslandPick protocol (Langille et al. 2008; also used by SSG-LUGIA)
-----------------------------------------------------------------
Accuracy is measured only on the two *curated* parts of each chromosome, never on the
whole genome:

* **positive set** – regions IslandPick's comparative-genomics procedure calls islands;
* **negative set** – backbone regions conserved across all comparison genomes.

Hence ``TP = |pred ∩ pos|``, ``FN = |pos \ pred|``, ``FP = |pred ∩ neg|``,
``TN = |neg \ pred|``; sequence belonging to neither set is ignored.  Using the same
convention keeps our numbers directly comparable with the published tables.

Boundary metrics
----------------
Nucleotide-level P/R/F1 rewards covering an island but says nothing about *where* a method
thinks the island begins and ends.  We therefore also report, for every curated island that
a method overlaps at all, the absolute distance between the true and predicted start and end
coordinates, and summarise these as a median absolute boundary error (MABE) and as the
fraction of boundaries recovered within a tolerance.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


# ----------------------------------------------------------------------------- set metrics
def track_from_intervals(intervals, genome_len: int) -> np.ndarray:
    y = np.zeros(genome_len, dtype=np.uint8)
    for s, e in intervals:
        y[max(0, int(s)):min(genome_len, int(e) + 1)] = 1
    return y


def confusion(pred: np.ndarray, pos: np.ndarray, neg: np.ndarray) -> dict:
    p, n = pos.astype(bool), neg.astype(bool)
    q = pred.astype(bool)
    tp = int(np.count_nonzero(q & p))
    fn = int(np.count_nonzero(~q & p))
    fp = int(np.count_nonzero(q & n))
    tn = int(np.count_nonzero(~q & n))
    return dict(TP=tp, FP=fp, TN=tn, FN=fn)


def rates(cm: dict) -> dict:
    tp, fp, tn, fn = cm["TP"], cm["FP"], cm["TN"], cm["FN"]
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    acc = (tp + tn) / (tp + tn + fp + fn) if (tp + tn + fp + fn) else 0.0
    mcc_den = np.sqrt(float(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = ((tp * tn - fp * fn) / mcc_den) if mcc_den > 0 else 0.0
    return dict(precision=100 * prec, recall=100 * rec, f1=100 * f1,
                accuracy=100 * acc, mcc=mcc, **cm)


def evaluate_genome(pred_intervals, pos_intervals, neg_intervals, genome_len: int) -> dict:
    pred = track_from_intervals(pred_intervals, genome_len)
    pos = track_from_intervals(pos_intervals, genome_len)
    neg = track_from_intervals(neg_intervals, genome_len)
    return rates(confusion(pred, pos, neg))


# -------------------------------------------------------------------------- boundary metrics
def boundary_errors(pred_intervals, pos_intervals) -> pd.DataFrame:
    """For each curated island, boundary errors against its best-overlapping prediction."""
    rows = []
    preds = [(int(s), int(e)) for s, e in pred_intervals]
    for ts, te in ((int(s), int(e)) for s, e in pos_intervals):
        best, best_ov = None, 0
        for ps, pe in preds:
            ov = min(te, pe) - max(ts, ps) + 1
            if ov > best_ov:
                best, best_ov = (ps, pe), ov
        if best is None:
            rows.append(dict(true_start=ts, true_end=te, matched=False,
                             start_err=np.nan, end_err=np.nan, mean_err=np.nan,
                             signed_start_err=np.nan, signed_end_err=np.nan, jaccard=0.0))
            continue
        ps, pe = best
        union = max(te, pe) - min(ts, ps) + 1
        rows.append(dict(true_start=ts, true_end=te, matched=True,
                         start_err=abs(ps - ts), end_err=abs(pe - te),
                         mean_err=(abs(ps - ts) + abs(pe - te)) / 2,
                         signed_start_err=ps - ts, signed_end_err=pe - te,
                         jaccard=best_ov / union))
    return pd.DataFrame(rows)


def summarise_boundaries(df: pd.DataFrame, tolerances=(1000, 2500, 5000)) -> dict:
    m = df[df.matched]
    if len(m) == 0:
        out = dict(n_matched=0, mabe=np.nan, median_jaccard=np.nan,
                   median_signed_start=np.nan, median_signed_end=np.nan)
        out.update({f"within_{t}": np.nan for t in tolerances})
        return out
    errs = np.concatenate([m.start_err.values, m.end_err.values])
    out = dict(n_matched=int(len(m)),
               mabe=float(np.median(errs)),
               median_jaccard=float(m.jaccard.median()),
               median_signed_start=float(m.signed_start_err.median()),
               median_signed_end=float(m.signed_end_err.median()))
    out.update({f"within_{t}": float(100 * np.mean(errs <= t)) for t in tolerances})
    return out


# ------------------------------------------------------------------------------ aggregation
def macro_average(per_genome: pd.DataFrame, cols=("precision", "recall", "f1", "accuracy")) -> dict:
    """Average metrics over genomes, not over pooled nucleotides.

    Ibtehaz et al. stress this point: pooling TP/FP/TN/FN across the whole benchmark lets
    the largest genomes and the largest islands dominate.  We follow their convention.
    """
    out = {c: float(per_genome[c].mean()) for c in cols if c in per_genome}
    if {"precision", "recall"} <= set(out):
        p, r = out["precision"], out["recall"]
        out["f1_of_means"] = 2 * p * r / (p + r) if (p + r) else 0.0
    return out
