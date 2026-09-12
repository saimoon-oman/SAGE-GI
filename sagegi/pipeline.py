"""End-to-end predictors.

``run_ssg_lugia``  – the compositional baseline, faithful to the published algorithm.
``run_sage_gi``    – the proposed pipeline (see ``sagegi.sage``).
"""
from __future__ import annotations

import numpy as np

from .anomaly import two_stage_envelope
from .features import extract_features
from .postprocess import (drop_short, median_smooth, project_to_nucleotides, runs)


def predict_from_scores(scores: np.ndarray, genome_len: int, params: dict) -> list[tuple[int, int]]:
    """Shared tail of every predictor: smooth → threshold → project → drop short islands."""
    smoothed = median_smooth(scores, params["median_filter_window_len"])
    alien = smoothed < 0
    track = project_to_nucleotides(alien, genome_len, params["w"], params["dw"],
                                   params.get("assign", "start"))
    return drop_short(runs(track), params["min_island_len"])


def run_ssg_lugia(code: np.ndarray, params: dict, features=None):
    """Predict genomic islands with the compositional SSG-LUGIA pipeline."""
    wf = features if features is not None else extract_features(code, params)
    _, scores = two_stage_envelope(wf.X, params)
    islands = predict_from_scores(scores, len(code), params)
    return islands, dict(scores=scores, features=wf)
