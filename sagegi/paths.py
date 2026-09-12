"""Central path configuration.

Heavy artefacts (genome FASTA, tile embeddings) live outside the repository so that the
repository stays small and cloud-sync friendly.  Override with the ``SAGEGI_WORK``
environment variable.
"""
from __future__ import annotations
import os
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
WORK = Path(os.environ.get("SAGEGI_WORK", r"C:/sagegi_work"))

DATA = REPO / "data"
BENCH = DATA / "benchmarks"
CONFIGS = REPO / "configs"
RESULTS = REPO / "results"
FIGURES = REPO / "figures"

GENOMES = WORK / "genomes"
EMB = WORK / "emb"
MODELS = WORK / "models"
RAW = WORK / "raw"
CACHE = WORK / "cache"

for _p in (DATA, BENCH, RESULTS, FIGURES, GENOMES, EMB, MODELS, RAW, CACHE):
    _p.mkdir(parents=True, exist_ok=True)
