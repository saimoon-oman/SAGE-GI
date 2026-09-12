#!/usr/bin/env python
"""Case studies on three well-characterised genomes.

Unlike IslandPick, these genomes have a curated list of islands covering the *whole*
chromosome, so evaluation is genome-wide: everything outside a reference island counts as
negative. This is the protocol used by Lu & Leong (2016) and by SSG-LUGIA's Table 3, which
is why we can quote their numbers alongside ours.

Genomes
-------
NC_003198.1  Salmonella enterica serovar Typhi CT18   -- 19 curated islands, 10 of them
             pathogenicity islands; the standard qualitative case study
NC_002935.2  Corynebacterium diphtheriae NCTC13129    -- 23 curated islands
NC_011770.1  Pseudomonas aeruginosa LESB58            -- 11 curated islands

Reference coordinates come from the GI-SVM / GI-review repository of Lu & Leong,
github.com/icelu/GI_Prediction, which is the source SSG-LUGIA cites.

    python scripts/07_case_studies.py
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from sagegi import models                                            # noqa: E402
from sagegi.anomaly import two_stage_envelope                        # noqa: E402
from sagegi.evaluate import (boundary_errors, rates, summarise_boundaries,  # noqa: E402
                             track_from_intervals)
from sagegi.features import extract_features                         # noqa: E402
from sagegi.genome import encode, read_fasta                         # noqa: E402
from sagegi.paths import BENCH, EMB, GENOMES, RESULTS, WORK          # noqa: E402
from sagegi.pipeline import predict_from_scores                      # noqa: E402
from sagegi.sage import SageConfig, run_sage_gi                      # noqa: E402
from sagegi.store import load_tiles                                  # noqa: E402

COMP = WORK / "comp"

CASES = {
    "NC_003198.1": "Salmonella enterica serovar Typhi CT18",
    "NC_002935.2": "Corynebacterium diphtheriae NCTC13129",
    "NC_011770.1": "Pseudomonas aeruginosa LESB58",
}

#: Published performance on CT18, quoted from Table 3 of Ibtehaz et al. (2021), which in turn
#: collected them from Lu & Leong (2016). Values are (recall, precision, F1) in per cent.
CT18_PUBLISHED = {
    "GI-SVM": (89.5, 44.6, 59.6),
    "EGID": (77.9, 53.5, 63.4),
    "SIGI-HMM": (24.1, 55.6, 33.7),
    "IslandViewer": (65.4, 67.0, 66.2),
    "GIHunter": (82.7, 67.6, 74.4),
    "IslandPath-DIMOB": (55.3, 78.8, 65.0),
    "tRNAcc": (28.6, 99.3, 44.4),
    "IslandPick": (6.0, 100.0, 11.4),
}


def genome_wide_metrics(pred, ref, genome_len):
    """Whole-chromosome confusion: everything outside a reference island is negative."""
    p = track_from_intervals(pred, genome_len).astype(bool)
    r = track_from_intervals(ref, genome_len).astype(bool)
    return rates(dict(TP=int(np.count_nonzero(p & r)), FP=int(np.count_nonzero(p & ~r)),
                      TN=int(np.count_nonzero(~p & ~r)), FN=int(np.count_nonzero(~p & r))))


def load_reference(acc):
    df = pd.read_csv(BENCH / "case_study_islands.csv")
    sub = df[df.accession == acc]
    return [(int(s) - 1, int(e) - 1) for s, e in zip(sub.start, sub.end)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--configs", default="configs/methods.json")
    ap.add_argument("--emb-dir", default="DNABERT-S_t1000")
    args = ap.parse_args()

    cfgs = json.loads(Path(args.configs).read_text())
    rows, brows, ivrows = [], [], []

    for acc, organism in CASES.items():
        fasta = GENOMES / f"{acc}.fna"
        if not fasta.exists():
            print(f"skip {acc}: no sequence"); continue
        ref = load_reference(acc)
        if not ref:
            print(f"skip {acc}: no reference islands"); continue
        seq = read_fasta(fasta)
        code = encode(seq)
        L = len(seq)
        print(f"\n=== {acc}  {organism}  {L:,} bp  {len(ref)} reference islands ===",
              flush=True)

        cache = COMP / f"{acc}.npz"
        if cache.exists():
            z = np.load(cache)
            X, starts = z["X"].astype(np.float64), z["starts"]
        else:
            wf = extract_features(code, models.get("SSG-LUGIA-F"))
            X, starts = wf.X, wf.starts

        # --- SSG-LUGIA baseline, exactly as published
        for v in ("SSG-LUGIA-P", "SSG-LUGIA-F", "SSG-LUGIA-R"):
            p = models.get(v)
            key = f"score_{v}"
            scores = z[key].astype(np.float64) if (cache.exists() and key in z) \
                else two_stage_envelope(X, p)[1]
            isl = predict_from_scores(scores, L, p)
            m = genome_wide_metrics(isl, ref, L)
            be = boundary_errors(isl, ref)
            rows.append(dict(accession=acc, organism=organism, method=v, n_pred=len(isl),
                             **m, **summarise_boundaries(be)))
            print(f"  {v:22s} P={m['precision']:6.2f} R={m['recall']:6.2f} "
                  f"F1={m['f1']:6.2f}", flush=True)
            for s, e in isl:
                ivrows.append(dict(accession=acc, method=v, start=s + 1, end=e + 1))

        # --- SAGE-GI and ablations, where embeddings exist
        tiles = {}
        for d in {c["embedding"] for c in cfgs.values() if c.get("embedding")}:
            f = EMB / d / f"{acc}.npz"
            if f.exists():
                tiles[d] = load_tiles(f)
        for name, spec in cfgs.items():
            spec = dict(spec)
            d = spec.pop("embedding", None)
            cfg = SageConfig(**spec)
            if cfg.use_embedding and d not in tiles:
                continue
            if cfg.use_embedding:
                Z, tile, _, is_pcs = tiles[d]
            else:
                Z, tile, is_pcs = None, 1000, False
            res = run_sage_gi(L, cfg, emb=Z, tile=tile, comp_X=X.astype(np.float32),
                              comp_starts=starts, emb_is_pcs=is_pcs)
            m = genome_wide_metrics(res.islands, ref, L)
            be = boundary_errors(res.islands, ref)
            rows.append(dict(accession=acc, organism=organism, method=name,
                             n_pred=len(res.islands), **m, **summarise_boundaries(be)))
            print(f"  {name:22s} P={m['precision']:6.2f} R={m['recall']:6.2f} "
                  f"F1={m['f1']:6.2f}", flush=True)
            be.insert(0, "method", name); be.insert(0, "accession", acc)
            brows.append(be)
            for s, e in res.islands:
                ivrows.append(dict(accession=acc, method=name, start=s + 1, end=e + 1))

    if not rows:
        print("nothing evaluated"); return
    pd.DataFrame(rows).to_csv(RESULTS / "case_studies.csv", index=False)
    pd.DataFrame(ivrows).to_csv(RESULTS / "case_study_intervals.csv", index=False)
    if brows:
        pd.concat(brows, ignore_index=True).to_csv(RESULTS / "case_study_boundaries.csv",
                                                   index=False)
    print("\nwrote results/case_studies.csv")


if __name__ == "__main__":
    main()
