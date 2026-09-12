#!/usr/bin/env python
"""Download and parse the IslandPick benchmark into tidy CSV tables.

The IslandPick benchmark (Langille, Hsiao & Brinkman 2008, *BMC Bioinformatics* 9:329) is
distributed only as Excel additional files.  This script turns them into three tidy tables
that everything downstream consumes:

``data/benchmarks/islandpick_positive.csv``
    One row per curated genomic island: accession, start, end, length.
``data/benchmarks/islandpick_negative.csv``
    One row per curated non-island (backbone) region.
``data/benchmarks/islandpick_baselines.csv``
    Per-genome TP / FP / TN / FN for the six tools benchmarked in the original study
    (SIGI-HMM, Centroid, IslandPath-DIMOB, PAI-IDA, IslandPath-DINUC, AlienHunter).
    These are the numbers SSG-LUGIA used to build its Table 2, so reusing them keeps our
    comparison exactly commensurate with the published literature.
"""
from __future__ import annotations
import json, re, sys, urllib.request
from pathlib import Path

import pandas as pd
import xlrd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from sagegi.paths import BENCH, RAW

BASE = ("https://media.springernature.com/original/springer-static/esm/"
        "art%3A10.1186%2F1471-2105-9-329/MediaObjects/12859_2008_2314_MOESM{n}_ESM.xls")
FILES = {2: "positive", 4: "negative", 6: "accuracy"}
ACC_RE = re.compile(r"^(N[CZ]_\d+\.\d+|AC_\d+\.\d+)\s*,\s*(.*)$")


def ensure(n: int) -> Path:
    dest = RAW / f"langille_AF{n}.xls"
    if not dest.exists() or dest.stat().st_size < 10_000:
        req = urllib.request.Request(BASE.format(n=n), headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=300) as r:
            dest.write_bytes(r.read())
    return dest


def parse_regions(path: Path, kind: str) -> pd.DataFrame:
    """Rows are grouped under '<accession>, <description>' section headers."""
    sheet = xlrd.open_workbook(path).sheet_by_index(0)
    rows, acc, desc = [], None, None
    for r in range(sheet.nrows):
        v = sheet.row_values(r)
        first = str(v[0]).strip()
        if not first:
            continue
        m = ACC_RE.match(first)
        if m:
            acc, desc = m.group(1), m.group(2).strip()
            continue
        if first.lower().startswith("total") or first.lower().startswith("additional"):
            continue
        if first.lower().startswith("region_id"):
            continue
        if acc is None:
            continue
        try:
            start, end = int(float(v[1])), int(float(v[2]))
        except (ValueError, TypeError):
            continue
        if end < start:
            start, end = end, start
        rows.append(dict(accession=acc, organism=desc, region_id=first,
                         start=start, end=end, length=end - start + 1, kind=kind))
    return pd.DataFrame(rows)


def parse_accuracy(path: Path) -> pd.DataFrame:
    sheet = xlrd.open_workbook(path).sheet_by_index(0)
    hdr = None
    rows = []
    for r in range(sheet.nrows):
        v = sheet.row_values(r)
        if str(v[0]).strip() == "Tool":
            hdr = [str(x).strip() for x in v]
            continue
        if hdr is None or not str(v[0]).strip() or str(v[1]).strip() in ("", "Total"):
            continue
        rows.append(dict(zip(hdr, v)))
    df = pd.DataFrame(rows)
    num = ["True Positives", "False Positives", "Total Predicted", "True Negative",
           "False Negative", "Negative Dataset", "Positive Dataset"]
    for c in num:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.rename(columns={"Tool": "tool", "Accession": "accession", "Description": "organism",
                            "True Positives": "TP", "False Positives": "FP",
                            "True Negative": "TN", "False Negative": "FN"})
    return df[["tool", "accession", "organism", "TP", "FP", "TN", "FN"]]


def main() -> None:
    paths = {n: ensure(n) for n in FILES}
    pos = parse_regions(paths[2], "island")
    neg = parse_regions(paths[4], "backbone")
    acc = parse_accuracy(paths[6])
    acc = acc[acc.tool != "literature"].copy()

    pos.to_csv(BENCH / "islandpick_positive.csv", index=False)
    neg.to_csv(BENCH / "islandpick_negative.csv", index=False)
    acc.to_csv(BENCH / "islandpick_baselines.csv", index=False)

    accessions = sorted(set(pos.accession) | set(neg.accession) | set(acc.accession))
    (BENCH / "islandpick_accessions.json").write_text(json.dumps(accessions, indent=1))

    print(f"positive : {len(pos):5d} islands   over {pos.accession.nunique()} genomes"
          f"  ({pos.length.sum():,} bp)")
    print(f"negative : {len(neg):5d} regions   over {neg.accession.nunique()} genomes"
          f"  ({neg.length.sum():,} bp)")
    print(f"baselines: {len(acc):5d} rows      {sorted(acc.tool.unique())}")
    print(f"accessions written: {len(accessions)}")


if __name__ == "__main__":
    main()
