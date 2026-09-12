"""Shared helpers for the driver scripts."""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from sagegi.paths import BENCH, GENOMES        # noqa: E402
from sagegi.genome import encode, read_fasta   # noqa: E402


def load_benchmark():
    pos = pd.read_csv(BENCH / "islandpick_positive.csv")
    neg = pd.read_csv(BENCH / "islandpick_negative.csv")
    return pos, neg


def genome_codes(acc: str):
    seq = read_fasta(GENOMES / f"{acc}.fna")
    return encode(seq), len(seq)


def intervals_for(df: pd.DataFrame, acc: str):
    """0-based inclusive intervals for one accession (benchmark files are 1-based)."""
    sub = df[df.accession == acc]
    return [(int(s) - 1, int(e) - 1) for s, e in zip(sub.start, sub.end)]


def available_accessions():
    return sorted(p.stem for p in GENOMES.glob("*.fna") if p.stat().st_size > 1000)
