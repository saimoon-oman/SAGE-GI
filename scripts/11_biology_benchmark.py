#!/usr/bin/env python
"""Functional enrichment of predicted islands across the whole benchmark.

`scripts/09_biological_relevance.py` answers this for one chromosome in detail. This script
answers it for all 118, and adds the question that decides whether the method is worth using:

    the islands SAGE-GI finds that the baseline **misses** -- are they biologically coherent,
    or are they compositional noise?

If the extra calls are enriched for mobility and virulence machinery at rates comparable to
the calls both methods agree on, the additional recall is real signal. If they are not, the
extra recall is bought with junk and the F1 gain does not mean what it appears to.

Gene coordinates and product descriptions come from NCBI feature tables (`rettype=ft`),
downloaded once per accession and cached.

    python scripts/11_biology_benchmark.py --workers 3
"""
from __future__ import annotations
import argparse, json, os, re, sys, time, urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from _common import load_benchmark                                   # noqa: E402
from sagegi.paths import BENCH, RESULTS, WORK                        # noqa: E402

FT = WORK / "featuretables"
FT.mkdir(parents=True, exist_ok=True)
GLEN = json.loads((BENCH / "genome_lengths.json").read_text())
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Same conservative keyword sets as scripts/09; see that file for the rationale.
CATEGORIES = {
    "mobility": [r"\btransposase\b", r"\bintegrase\b", r"\brecombinase\b", r"\btransposon\b",
                 r"\binsertion sequence\b", r"\bIS\d{2,4}\b", r"\bphage\b", r"\bprophage\b",
                 r"\bconjugal\b", r"\bconjugative\b", r"\bplasmid\b", r"\brelaxase\b",
                 r"\bmobilization\b", r"\bexcisionase\b", r"\bresolvase\b"],
    "virulence": [r"\bvirulence\b", r"\bpathogenicity\b", r"\binvasin\b", r"\bhemolysin\b",
                  r"\bhaemolysin\b", r"\btoxin\b", r"\badhesin\b", r"\bfimbria", r"\bpilus\b",
                  r"\bpili\b", r"\bsecretion system\b", r"\bflagell", r"\binvasion protein\b",
                  r"\beffector\b"],
    "resistance": [r"\bresistance\b", r"\bbeta-lactamase\b", r"\bantibiotic\b", r"\befflux\b",
                   r"\bmultidrug\b", r"\bchloramphenicol\b", r"\btetracycline\b",
                   r"\bstreptomycin\b", r"\bsulfonamide\b", r"\baminoglycoside\b"],
}
COMPILED = {k: re.compile("|".join(v), re.I) for k, v in CATEGORIES.items()}
_COORD = re.compile(r"^[<>]?(\d+)\t[<>]?(\d+)\t(\w+)")

_gate, _last = __import__("threading").Semaphore(1), [0.0]


def fetch_ft(acc: str, retries: int = 4) -> Path:
    dest = FT / f"{acc}.ft"
    if dest.exists() and dest.stat().st_size > 500:
        return dest
    url = (f"{EUTILS}?db=nuccore&id={acc}&rettype=ft&retmode=text"
           f"&tool=SAGE-GI&email=1024052020@grad.cse.buet.ac.bd")
    for attempt in range(retries):
        try:
            with _gate:                                # respect the 3 req/s Entrez limit
                w = 0.4 - (time.time() - _last[0])
                if w > 0:
                    time.sleep(w)
                _last[0] = time.time()
            with urllib.request.urlopen(url, timeout=300) as r:
                txt = r.read().decode("utf-8", "replace")
            if len(txt) > 500:
                dest.write_text(txt, encoding="utf-8")
                return dest
            raise ValueError("short payload")
        except Exception:                              # noqa: BLE001
            time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"could not fetch feature table for {acc}")


def parse_ft(path: Path) -> pd.DataFrame:
    """Extract CDS intervals and their product description from an NCBI feature table."""
    rows, cur = [], None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        m = _COORD.match(line)
        if m:
            a, b, kind = int(m.group(1)), int(m.group(2)), m.group(3)
            if kind == "CDS":
                cur = dict(start=min(a, b), end=max(a, b), product="")
                rows.append(cur)
            else:
                cur = None
            continue
        if cur is not None and line.startswith("\t\t\t"):
            parts = line.strip().split("\t", 1)
            if len(parts) == 2 and parts[0] == "product":
                cur["product"] = parts[1]
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    for cat, rx in COMPILED.items():
        df[cat] = df["product"].fillna("").str.contains(rx)
    df["any_hgt"] = df[list(CATEGORIES)].any(axis=1)
    return df


def mask_for(genes: pd.DataFrame, intervals) -> np.ndarray:
    mid = ((genes.start + genes.end) / 2).values
    m = np.zeros(len(genes), dtype=bool)
    for s, e in intervals:
        m |= (mid >= s) & (mid <= e)
    return m


def subtract(a_iv, b_iv, genome_len):
    """Bases called by a but not by b, as intervals."""
    t = np.zeros(genome_len, dtype=bool)
    for s, e in a_iv:
        t[max(0, s):min(genome_len, e + 1)] = True
    for s, e in b_iv:
        t[max(0, s):min(genome_len, e + 1)] = False
    d = np.diff(np.concatenate(([0], t.view(np.uint8), [0])).astype(np.int8))
    return list(zip(np.flatnonzero(d == 1).tolist(), (np.flatnonzero(d == -1) - 1).tolist()))


MAIN_METHOD = "SAGE-GI"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=3)
    # Which configuration counts as "our method" is a question the evidence gets to settle,
    # so it is a flag rather than a hard-coded string.
    ap.add_argument("--main-method", default="SAGE-GI")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()
    global MAIN_METHOD
    MAIN_METHOD = args.main_method

    pos, _ = load_benchmark()
    sg = pd.read_csv(RESULTS / "sage_gi_intervals.csv")
    ssg = pd.read_csv(RESULTS / "ssg_lugia_intervals.csv")
    accs = sorted(set(sg.accession) & set(ssg.accession) & set(pos.accession))
    if args.limit:
        accs = accs[:args.limit]
    print(f"{len(accs)} genomes", flush=True)

    todo = [a for a in accs if not (FT / f"{a}.ft").exists()]
    if todo:
        print(f"downloading {len(todo)} feature tables...", flush=True)
        with ThreadPoolExecutor(args.workers) as ex:
            futs = {ex.submit(fetch_ft, a): a for a in todo}
            for i, f in enumerate(as_completed(futs), 1):
                try:
                    f.result()
                except Exception as exc:                 # noqa: BLE001
                    print(f"  {futs[f]}: {exc}")
                if i % 20 == 0 or i == len(todo):
                    print(f"  {i}/{len(todo)}", flush=True)

    # pooled contingency counts across genomes, per method and category
    counts: dict = {}

    def add(method, cat, a, b, c, d):
        k = (method, cat)
        p = counts.setdefault(k, [0, 0, 0, 0])
        p[0] += a; p[1] += b; p[2] += c; p[3] += d

    n_used = 0
    for acc in accs:
        f = FT / f"{acc}.ft"
        if not f.exists():
            continue
        genes = parse_ft(f)
        if genes.empty or len(genes) < 200:
            continue
        L = int(GLEN[acc])
        sets = {
            "curated islands": [(int(s) - 1, int(e) - 1) for s, e in
                                zip(pos[pos.accession == acc].start,
                                    pos[pos.accession == acc].end)],
        }
        for name, src in [("SSG-LUGIA-F", ssg), ("SSG-LUGIA-P", ssg),
                          (MAIN_METHOD, sg), ("SAGE-GI-P", sg)]:
            g = src[(src.accession == acc) & (src.method == name)]
            sets[name] = [(int(s) - 1, int(e) - 1) for s, e in zip(g.start, g.end)]
        if not sets[MAIN_METHOD] or not sets["SSG-LUGIA-F"]:
            continue
        sets["SAGE-GI only"] = subtract(sets[MAIN_METHOD], sets["SSG-LUGIA-F"], L)
        sets["both agree"] = [(s, e) for s, e in sets[MAIN_METHOD]
                              if any(min(e, be) >= max(s, bs)
                                     for bs, be in sets["SSG-LUGIA-F"])]
        n_used += 1
        for name, iv in sets.items():
            if not iv:
                continue
            m = mask_for(genes, iv)
            for cat in list(CATEGORIES) + ["any_hgt"]:
                col = genes[cat].values
                add(name, cat,
                    int((col & m).sum()), int((~col & m).sum()),
                    int((col & ~m).sum()), int((~col & ~m).sum()))

    rows = []
    for (method, cat), (a, b, c, d) in counts.items():
        odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append(dict(method=method, category=cat, n_marker_in=a, n_gene_in=a + b,
                         rate_in=100 * a / max(a + b, 1), rate_out=100 * c / max(c + d, 1),
                         odds_ratio=float(odds), p=float(p)))
    out = pd.DataFrame(rows)
    out.to_csv(RESULTS / "biology_benchmark.csv", index=False)

    print(f"\nPooled over {n_used} genomes with usable annotation\n")
    print(f"{'region set':18s} {'category':11s} {'in%':>6s} {'out%':>6s} {'OR':>6s} {'p':>10s}")
    order = ["curated islands", "SSG-LUGIA-P", "SSG-LUGIA-F", "SAGE-GI-P",
             MAIN_METHOD, "both agree", "SAGE-GI only"]
    for m in order:
        for cat in ["mobility", "virulence", "resistance", "any_hgt"]:
            r = out[(out.method == m) & (out.category == cat)]
            if len(r):
                r = r.iloc[0]
                print(f"{m:18s} {cat:11s} {r.rate_in:6.2f} {r.rate_out:6.2f} "
                      f"{r.odds_ratio:6.2f} {r.p:10.2e}")
    print("\nwrote results/biology_benchmark.csv")


if __name__ == "__main__":
    main()
