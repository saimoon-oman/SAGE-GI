#!/usr/bin/env python
"""Are the predicted islands biologically meaningful, or just compositionally odd?

Accuracy against a curated island list says a method agrees with a previous method. It does
not say the regions matter. Genomic islands are interesting because they concentrate
virulence factors, antibiotic-resistance determinants and mobility machinery, so a useful
detector should enrich for those functions inside its calls.

For each method we take its predicted islands on *Salmonella* Typhi CT18, count annotated
genes whose product matches a curated keyword set, and compare the rate inside predictions
against the rate outside, with a Fisher exact test. Keywords are grouped so that mobility
genes -- which are the mechanistic signature of horizontal transfer -- can be reported
separately from virulence and resistance.

Annotation is the NCBI protein table for NC_003198.1, distributed with the GI-review
repository of Lu & Leong (github.com/icelu/GI_Prediction).

    python scripts/09_biological_relevance.py
"""
from __future__ import annotations
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from sagegi.paths import BENCH, RESULTS, WORK                       # noqa: E402

ACC = "NC_003198.1"
PTT = WORK / "raw" / "ct18" / "NC_003198.ptt"

#: Curated keyword sets. Deliberately conservative: these match product descriptions that are
#: unambiguous about function, rather than casting a wide net that would inflate every count.
CATEGORIES = {
    "mobility": [
        r"\btransposase\b", r"\bintegrase\b", r"\brecombinase\b", r"\btransposon\b",
        r"\binsertion sequence\b", r"\bIS\d{2,4}\b", r"\bphage\b", r"\bprophage\b",
        r"\bconjugal\b", r"\bconjugative\b", r"\bplasmid\b", r"\brelaxase\b",
        r"\bmobilization\b", r"\bexcisionase\b", r"\bresolvase\b",
    ],
    "virulence": [
        r"\bvirulence\b", r"\bpathogenicity\b", r"\binvasin\b", r"\bhemolysin\b",
        r"\bhaemolysin\b", r"\btoxin\b", r"\badhesin\b", r"\bfimbria", r"\bpilus\b",
        r"\bpili\b", r"\bsecretion system\b", r"\bflagell", r"\binvasion protein\b",
        r"\beffector\b", r"\bVi polysaccharide\b",
    ],
    "resistance": [
        r"\bresistance\b", r"\bbeta-lactamase\b", r"\bantibiotic\b", r"\befflux\b",
        r"\bmultidrug\b", r"\bchloramphenicol\b", r"\btetracycline\b", r"\bstreptomycin\b",
        r"\bsulfonamide\b", r"\baminoglycoside\b",
    ],
}
COMPILED = {k: re.compile("|".join(v), re.I) for k, v in CATEGORIES.items()}


def load_genes() -> pd.DataFrame:
    """Parse an NCBI .ptt protein table into start/end/product rows."""
    rows = []
    with open(PTT, encoding="utf-8", errors="replace") as fh:
        lines = fh.read().splitlines()
    header = next(i for i, l in enumerate(lines) if l.startswith("Location"))
    for line in lines[header + 1:]:
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        m = re.match(r"(\d+)\.\.(\d+)", parts[0])
        if not m:
            continue
        rows.append(dict(start=int(m.group(1)), end=int(m.group(2)),
                         gene=parts[4], product=parts[8]))
    df = pd.DataFrame(rows)
    for cat, rx in COMPILED.items():
        df[cat] = df["product"].fillna("").str.contains(rx)
    df["any_hgt"] = df[list(CATEGORIES)].any(axis=1)
    return df


def inside(genes: pd.DataFrame, intervals) -> np.ndarray:
    """Gene is 'inside' if its midpoint falls in a predicted island."""
    mid = ((genes.start + genes.end) / 2).values
    m = np.zeros(len(genes), dtype=bool)
    for s, e in intervals:
        m |= (mid >= s) & (mid <= e)
    return m


def enrichment(genes: pd.DataFrame, mask: np.ndarray, col: str):
    a = int((genes[col] & mask).sum())          # marker genes inside predictions
    b = int((~genes[col] & mask).sum())         # other genes inside
    c = int((genes[col] & ~mask).sum())         # marker genes outside
    d = int((~genes[col] & ~mask).sum())        # other genes outside
    odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
    rate_in = 100 * a / max(a + b, 1)
    rate_out = 100 * c / max(c + d, 1)
    return dict(n_in=a, n_out=c, rate_in=rate_in, rate_out=rate_out,
                odds_ratio=float(odds), p=float(p))


def main() -> None:
    if not PTT.exists():
        print(f"annotation not found: {PTT}\n"
              "fetch it from the GI-review repository first "
              "(GI_review/genome/NC_003198.ptt)")
        return
    genes = load_genes()
    print(f"{len(genes)} annotated genes on {ACC}")
    for cat in CATEGORIES:
        print(f"  {cat:12s} {int(genes[cat].sum()):4d} genes")
    print()

    iv = pd.read_csv(RESULTS / "case_study_intervals.csv")
    iv = iv[iv.accession == ACC]

    ref = pd.read_csv(BENCH / "case_study_islands.csv")
    ref = ref[ref.accession == ACC]
    sets = {"curated islands (reference)": list(zip(ref.start - 1, ref.end - 1))}
    for m in ["SSG-LUGIA-F", "Comp-only+center", "DNABERT-S-only", "SAGE-GI"]:
        g = iv[iv.method == m]
        if len(g):
            sets[m] = list(zip(g.start - 1, g.end - 1))

    rows = []
    for name, intervals in sets.items():
        mask = inside(genes, intervals)
        bp = sum(e - s + 1 for s, e in intervals)
        for cat in list(CATEGORIES) + ["any_hgt"]:
            e = enrichment(genes, mask, cat)
            rows.append(dict(method=name, category=cat, n_genes_called=int(mask.sum()),
                             pred_bp=bp, **e))
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "biological_relevance.csv", index=False)

    print(f"{'method':28s} {'category':11s} {'in%':>6s} {'out%':>6s} {'OR':>6s} {'p':>10s}")
    for r in out.itertuples():
        print(f"{r.method:28s} {r.category:11s} {r.rate_in:6.2f} {r.rate_out:6.2f} "
              f"{r.odds_ratio:6.2f} {r.p:10.2e}")
    print("\nwrote results/biological_relevance.csv")


if __name__ == "__main__":
    main()
