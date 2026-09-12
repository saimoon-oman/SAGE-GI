"""Pre-tuned parameter sets.

``SSG_LUGIA_*`` reproduce Table 1 of Ibtehaz et al. (2021) exactly.
``SAGE_GI_*`` are the settings tuned for this work (see ``scripts/40_tune.py``).
"""
from __future__ import annotations
from copy import deepcopy

_BASE = dict(w=10000, dw=100, karlin_mode="normalized", pca_dn=2, pca_amino_acid=2,
             pca_kmer4=2, entropy_features=True, support_fraction_model1=0.75,
             support_fraction_model2=0.9, median_filter_window_len=400,
             min_island_len=10000, assign="start")

SSG_LUGIA_F = {**_BASE, "contamination_model1": 0.15, "contamination_model2": 0.05}
SSG_LUGIA_R = {**_BASE, "contamination_model1": 0.20, "contamination_model2": 0.25}
SSG_LUGIA_P = {**_BASE, "contamination_model1": 0.075, "contamination_model2": 0.075}

SSG_LUGIA = {"SSG-LUGIA-F": SSG_LUGIA_F, "SSG-LUGIA-R": SSG_LUGIA_R, "SSG-LUGIA-P": SSG_LUGIA_P}


def get(name: str) -> dict:
    if name in SSG_LUGIA:
        return deepcopy(SSG_LUGIA[name])
    raise KeyError(name)
