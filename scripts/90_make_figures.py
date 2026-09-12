#!/usr/bin/env python
"""Build every figure in the paper from the tidy result tables.

Each figure is a function; missing inputs are skipped with a note rather than crashing, so
the script can be re-run at any point while the pipeline is still filling in results.

    python scripts/90_make_figures.py [--only fig3 fig4]
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from sagegi.paths import BENCH, EMB, FIGURES, RESULTS, WORK          # noqa: E402
from sagegi.viz import (BASELINE_COLOR, CAT4, GRID, INK, INK2, MUTED,  # noqa: E402
                        PRETTY, color_for, label_bars, marker_for, save,
                        use_paper_style)
from sagegi.evaluate import rates                                    # noqa: E402

MAIN4 = ["Comp-only", "DNABERT-2-only", "DNABERT-S-only", "SAGE-GI"]
PRETTY_MAIN = {"Comp-only": "Hand-crafted\ncomposition",
               "DNABERT-2-only": "DNABERT-2\nembedding",
               "DNABERT-S-only": "DNABERT-S\nembedding",
               "SAGE-GI": "SAGE-GI\n(fusion)"}


def read(name):
    p = RESULTS / name
    return pd.read_csv(p) if p.exists() else None


# ============================================================== Fig 1 - pipeline schematic
def fig1_schematic():
    fig, ax = plt.subplots(figsize=(7.4, 3.6))
    ax.set_xlim(0, 100); ax.set_ylim(-1, 50); ax.axis("off")

    def box(x, y, w, h, title, body, color):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.6,rounding_size=1.2",
                                    linewidth=1.2, edgecolor=color,
                                    facecolor=color + "14"))
        ax.text(x + w / 2, y + h - 2.6, title, ha="center", va="top", fontsize=7.6,
                weight="semibold", color=INK)
        ax.text(x + w / 2, y + h - 7.0, body, ha="center", va="top", fontsize=6.6,
                color=INK2, linespacing=1.45)

    def arrow(x0, y0, x1, y1):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=9,
                                     linewidth=1.0, color=MUTED, shrinkA=0, shrinkB=0))

    TOP, BOT, H = 26, 2, 19
    box(1, TOP, 22, H, "1 · Tile once",
        "genome split into\nnon-overlapping\n1 kb tiles; DNABERT-S\nembeds each tile", CAT4[2])
    box(26, TOP, 22, H, "2 · Host-reference",
        "L2-normalise, then PCA\nonto this chromosome's\nown principal subspace\n— the host defines its\nown frame of reference",
        CAT4[2])
    box(51, TOP, 22, H, "3 · Pool to windows",
        "overlap-weighted mean\nover the 10 kb / 100 bp\ngrid — cost stays linear\nin genome length", CAT4[0])
    box(76, TOP, 23, H, "4 · Composition\n(off by default)",
        "the 11 compositional\nfeatures can be\nconcatenated here \u2014 over\n118 genomes this earns\nnothing, so it is not used", MUTED)

    box(14, BOT, 23, H, "5 · Detect",
        "two-stage robust\nMinimum-Covariance-\nDeterminant anomaly\ndetection over all\nwindows", CAT4[1])
    box(41, BOT, 23, H, "6 · Localise",
        "median smoothing,\ncentred window →\nnucleotide projection,\nminimum island length", CAT4[1])
    box(68, BOT, 23, H, "7 · Sharpen edges",
        "CUSUM change-point\nsearch on the tile-\nresolution anomaly track\nfixes each boundary", CAT4[3])

    for x0, x1 in [(23, 26), (48, 51), (73, 76)]:
        arrow(x0, TOP + H / 2, x1, TOP + H / 2)
    for x0, x1 in [(37, 41), (64, 68)]:
        arrow(x0, BOT + H / 2, x1, BOT + H / 2)

    # wrap-around connector: down the right, back along the gutter, into box 5
    y_gut = (TOP + BOT + H) / 2
    ax.plot([87.5, 87.5], [TOP, y_gut], color=MUTED, linewidth=1.0, solid_capstyle="butt")
    ax.plot([87.5, 25.5], [y_gut, y_gut], color=MUTED, linewidth=1.0, solid_capstyle="butt")
    arrow(25.5, y_gut, 25.5, BOT + H)

    ax.text(1, 48.4, "frozen model, inference only — no fine-tuning, no annotation, "
                     "no reference genome", fontsize=6.8, color=MUTED, style="italic")
    save(fig, FIGURES / "fig1_pipeline")
    return "fig1_pipeline"


# ======================================================== Fig 2 - compute scaling argument
def fig2_compute():
    """Why tile-and-pool is what makes this feasible at all."""
    fig, ax = plt.subplots(figsize=(4.4, 3.2))
    lengths = np.array([1, 2, 4, 6, 8]) * 1e6
    w, dw, tile = 10000, 100, 1000
    naive = lengths / dw * w            # bases pushed through the model, per genome
    tiled = lengths                     # each base embedded exactly once
    ax.plot(lengths / 1e6, naive / 1e9, color=CAT4[1], marker=marker_for("Comp-only"),
            label="window-wise (naive)")
    ax.plot(lengths / 1e6, tiled / 1e9, color=CAT4[0], marker=marker_for("SAGE-GI"),
            label="tile-and-pool (SAGE-GI)")
    ax.set_yscale("log")
    ax.set_xlabel("chromosome length (Mbp)")
    ax.set_ylabel("bases pushed through the model (Gbp)")
    ax.set_title("Transformer work per genome", loc="left")
    ax.legend(loc="center right")
    ax.annotate(f"{w // dw}x less work\nat every genome size", xy=(6, tiled[3] / 1e9),
                xytext=(2.2, 0.25), fontsize=7.2, color=INK2,
                arrowprops=dict(arrowstyle="-|>", color=MUTED, linewidth=0.9))
    save(fig, FIGURES / "fig2_compute")
    return "fig2_compute"


# ================================================= Fig 3 - precision/recall on IslandPick
def _baseline_table():
    b = read("../data/benchmarks/islandpick_baselines.csv")
    if b is None:
        p = BENCH / "islandpick_baselines.csv"
        b = pd.read_csv(p) if p.exists() else None
    if b is None:
        return None
    r = pd.DataFrame([rates(dict(TP=t.TP, FP=t.FP, TN=t.TN, FN=t.FN)) for t in b.itertuples()])
    return pd.concat([b[["tool", "accession"]].reset_index(drop=True), r], axis=1)


def fig3_pr(restrict_to=None):
    ours = read("sage_gi_per_genome.csv")
    ssg = read("ssg_lugia_per_genome.csv")
    base = _baseline_table()
    if base is None or (ours is None and ssg is None):
        return None
    frames = []
    if ssg is not None:
        frames.append(ssg.assign(family="ssg")[["accession", "method", "precision", "recall", "f1"]])
    if ours is not None:
        frames.append(ours.assign(family="ours")[["accession", "method", "precision", "recall", "f1"]])
    mine = pd.concat(frames, ignore_index=True)
    accs = set(mine.accession)
    if restrict_to:
        accs &= set(restrict_to)
    base = base[base.accession.isin(accs)]
    mine = mine[mine.accession.isin(accs)]

    fig, ax = plt.subplots(figsize=(6.2, 4.4))

    # iso-F1 contours first, so everything else sits on top of them
    for f1 in (0.3, 0.4, 0.5, 0.6):
        r = np.linspace(max(f1 / (2 - f1) * 100 + 0.5, 15), 99, 300)
        p = f1 * r / (2 * r / 100 - f1)
        ok = (p > 0) & (p <= 100)
        ax.plot(r[ok], p[ok], color=GRID, linewidth=0.9, zorder=1)
        ax.annotate(f"F1={f1:.1f}", (r[ok][-1], p[ok][-1]), fontsize=6.4, color=MUTED,
                    textcoords="offset points", xytext=(-24, 4), zorder=1)

    # previously published tools: one neutral colour, labelled in place (they are well spread)
    OFF = {"SIGI-HMM": (7, -3), "IslandPath-DIMOB": (7, -3), "PAI-IDA": (7, -3),
           "Centroid": (7, -3), "IslandPath-DINUC": (7, -10), "AlienHunter": (-6, -12)}
    for tool, g in base.groupby("tool"):
        name = PRETTY.get(tool, tool)
        x, y = g.recall.mean(), g.precision.mean()
        ax.scatter(x, y, s=42, color=BASELINE_COLOR, marker="o", zorder=3, linewidths=0)
        ax.annotate(name, (x, y), textcoords="offset points",
                    xytext=OFF.get(name, (7, -3)), fontsize=7, color=INK2, zorder=3)

    # our methods: marker shape + colour + a legend, so labels never collide
    # Comp-only is our reproduction of SSG-LUGIA-F and plots on top of it; the fidelity
    # table makes that point properly, so it is omitted here.
    show = [m for m in ["SSG-LUGIA-P", "SSG-LUGIA-F", "SSG-LUGIA-R", "DNABERT-2-only",
                        "DNABERT-S-only", "SAGE-GI"] if m in set(mine.method)]
    for m in show:
        g = mine[mine.method == m]
        ax.scatter(g.recall.mean(), g.precision.mean(), s=95, color=color_for(m),
                   marker=marker_for(m), zorder=5, edgecolor="white", linewidths=1.2,
                   label=m)
    ax.set_xlabel("recall (%)"); ax.set_ylabel("precision (%)")
    ax.set_xlim(15, 100); ax.set_ylim(20, 100)
    ax.set_title(f"IslandPick benchmark ({len(accs)} genomes)", loc="left")
    ax.legend(loc="upper right", fontsize=7.2, handletextpad=0.3, borderpad=0.5,
              labelspacing=0.45, frameon=True, framealpha=0.94, edgecolor=GRID)
    save(fig, FIGURES / "fig3_precision_recall")
    return "fig3_precision_recall"


# ============================================================ Fig 4 - boundary resolution
def fig4_boundaries():
    frames = []
    for f in ("ssg_lugia_boundaries.csv", "sage_gi_boundaries.csv"):
        d = read(f)
        if d is not None:
            frames.append(d)
    if not frames:
        return None
    df = pd.concat(frames, ignore_index=True)
    df = df[df.matched]
    order = [m for m in ["SSG-LUGIA-F", "Comp-only", "DNABERT-2-only", "DNABERT-S-only",
                         "SAGE-GI/no-refine", "SAGE-GI"] if m in set(df.method)]
    if not order:
        return None

    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.4), gridspec_kw=dict(wspace=0.32))

    ax = axes[0]
    med = [np.median(np.concatenate([g.start_err.values, g.end_err.values])) / 1000
           for m in order for g in [df[df.method == m]]]
    bars = ax.bar(range(len(order)), med, width=0.62,
                  color=[color_for(m) for m in order], linewidth=0)
    label_bars(ax, bars, fmt="{:.1f}", dy=max(med) * 0.02)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([m.replace("SAGE-GI/", "SAGE-GI\n") for m in order],
                       rotation=30, ha="right", fontsize=7)
    ax.set_ylabel("median absolute boundary error (kb)")
    ax.set_title("How far off is each edge?", loc="left")

    ax = axes[1]
    for m in order:
        g = df[df.method == m]
        v = np.concatenate([g.signed_start_err.values, g.signed_end_err.values]) / 1000
        v = v[np.isfinite(v)]
        xs = np.linspace(-60, 60, 300)
        if len(v) > 5:
            from scipy.stats import gaussian_kde
            ax.plot(xs, gaussian_kde(np.clip(v, -60, 60))(xs), color=color_for(m),
                    linewidth=1.8, label=m)
    ax.axvline(0, color=MUTED, linewidth=0.9, linestyle="--")
    ax.set_xlabel("signed boundary error (kb)   predicted − true")
    ax.set_ylabel("density")
    ax.set_title("Is the error a bias or just noise?", loc="left")
    ax.legend(fontsize=6.6, loc="upper left")
    save(fig, FIGURES / "fig4_boundaries")
    return "fig4_boundaries"


# ==================================================== Fig 5 - synthetic insertion sweep
def fig5_synthetic():
    df = read("synthetic_insertions.csv")
    if df is None:
        return None

    #: Chimera embeddings are assembled from full-dimensional host and donor tiles in
    #: one coherent basis, so every channel is valid here; compositional features are
    #: additionally recomputed on the real chimeric sequence.
    VALID_ON_CHIMERA = {"Comp-only", "DNABERT-2-only", "DNABERT-S-only", "SAGE-GI"}
    df = df[df.method.isin(VALID_ON_CHIMERA)]
    if df.empty:
        return None
    methods = [m for m in MAIN4 if m in set(df.method)]
    df = df.copy()
    df["dist_bin"] = pd.qcut(df.tnf_euclidean, 2,
                             labels=["compositionally similar donor",
                                     "compositionally distant donor"])
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.3), sharey=True,
                             gridspec_kw=dict(wspace=0.12))
    for ax, (b, sub) in zip(axes, df.groupby("dist_bin", observed=True)):
        for m in methods:
            g = sub[sub.method == m].groupby("insert_len").detection_rate.mean()
            ax.plot(g.index / 1000, g.values, color=color_for(m), marker=marker_for(m),
                    label=PRETTY_MAIN.get(m, m).replace("\n", " "))
        ax.set_xscale("log")
        ax.set_xticks([5, 10, 20, 40, 80])
        ax.set_xticklabels(["5", "10", "20", "40", "80"])
        ax.set_xlabel("insert length (kb)")
        ax.set_title(str(b), loc="left", fontsize=8.5)
    axes[0].set_ylabel("inserts recovered (%)")
    axes[0].set_ylim(0, 102)
    axes[1].legend(fontsize=7, loc="lower right")
    save(fig, FIGURES / "fig5_synthetic")
    return "fig5_synthetic"


# ================================================================ Fig 6 - ablation study
def fig6_ablation():
    df = read("sage_gi_per_genome.csv")
    if df is None:
        return None
    order = [m for m in ["SAGE-GI", "SAGE-GI/no-refine", "SAGE-GI/no-detrend",
                         "SAGE-GI/start-assign", "DNABERT-S-only", "DNABERT-2-only",
                         "Comp-only"] if m in set(df.method)]
    if len(order) < 2:
        return None
    g = df[df.method.isin(order)].groupby("method")[["precision", "recall", "f1", "mabe"]].mean()
    g = g.loc[order]
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.4), gridspec_kw=dict(wspace=0.3))
    x = np.arange(len(order))
    ax = axes[0]
    for i, (metric, off) in enumerate([("precision", -0.26), ("recall", 0.0), ("f1", 0.26)]):
        bars = ax.bar(x + off, g[metric], width=0.24, linewidth=0,
                      color=[color_for(m) for m in order],
                      alpha=[1.0, 0.72, 0.48][i])
        if metric == "f1":
            label_bars(ax, bars, dy=1.0)
    ax.set_xticks(x); ax.set_xticklabels([m.replace("SAGE-GI/", "− ") for m in order],
                                         rotation=30, ha="right", fontsize=7)
    ax.set_ylabel("%")
    ax.set_title("Precision / recall / F1  (bar triplets, F1 labelled)", loc="left",
                 fontsize=8.2)
    ax = axes[1]
    bars = ax.bar(x, g["mabe"] / 1000, width=0.6, linewidth=0,
                  color=[color_for(m) for m in order])
    label_bars(ax, bars, fmt="{:.1f}", dy=(g['mabe'].max() / 1000) * 0.02)
    ax.set_xticks(x); ax.set_xticklabels([m.replace("SAGE-GI/", "− ") for m in order],
                                         rotation=30, ha="right", fontsize=7)
    ax.set_ylabel("median boundary error (kb)")
    ax.set_title("Boundary accuracy", loc="left", fontsize=8.2)
    save(fig, FIGURES / "fig6_ablation")
    return "fig6_ablation"


# ================================================= Fig 7 - CT18 genome map (case study)
def fig7_ct18_map(acc="NC_003198.1"):
    iv = read("case_study_intervals.csv")
    if iv is None:
        return None
    ref = pd.read_csv(BENCH / "case_study_islands.csv")
    ref = ref[ref.accession == acc]
    iv = iv[iv.accession == acc]
    if not len(ref) or not len(iv):
        return None
    glen = json.loads((BENCH / "genome_lengths.json").read_text())[acc]

    order = [m for m in ["SAGE-GI", "DNABERT-S-only", "SSG-LUGIA-F", "SSG-LUGIA-P",
                         "SSG-LUGIA-R"] if m in set(iv.method)]
    tracks = [("known islands", None)] + [(m, m) for m in order]

    fig, ax = plt.subplots(figsize=(7.4, 0.42 * len(tracks) + 1.25))
    h = 0.55
    for row, (label, method) in enumerate(tracks):
        y = len(tracks) - row - 1
        ax.add_patch(plt.Rectangle((0, y - h / 2), glen / 1e6, h, facecolor="#f2f1ee",
                                   edgecolor="none", zorder=1))
        if method is None:
            for r in ref.itertuples():
                ax.add_patch(plt.Rectangle((r.start / 1e6, y - h / 2),
                                           (r.end - r.start) / 1e6, h,
                                           facecolor=INK, edgecolor="none", zorder=2))
        else:
            g = iv[iv.method == method]
            for r in g.itertuples():
                ax.add_patch(plt.Rectangle((r.start / 1e6, y - h / 2),
                                           (r.end - r.start) / 1e6, h,
                                           facecolor=color_for(method), edgecolor="none",
                                           zorder=2))
        ax.text(-0.06 * glen / 1e6, y, label, ha="right", va="center", fontsize=7.6,
                color=INK if method is None else color_for(method),
                weight="semibold" if method in ("SAGE-GI", None) else "normal")

    ax.set_xlim(0, glen / 1e6)
    ax.set_ylim(-0.7, len(tracks) - 0.3)
    ax.set_yticks([])
    ax.set_xlabel("genome coordinate (Mbp)")
    ax.grid(axis="y", visible=False)
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)
    ax.set_title("Salmonella Typhi CT18: 19 curated islands and what each method calls",
                 loc="left", fontsize=9)
    save(fig, FIGURES / "fig7_ct18_map")
    return "fig7_ct18_map"



# --------------------------------------------------------------- Fig 8: representation
def fig8_representation():
    """Per-genome AUC for each representation, and the paired differences.

    Panel A shows that the separation is not carried by a few genomes; panel B shows the
    paired contrasts the text quotes, with the zero line making the win rate readable
    without consulting a table.
    """
    df = read("representation_auc.csv")
    if df is None:
        return None
    nod = df[df.detrend_bp == 0]
    if not len(nod):
        return None

    # Use each representation's own best dimensionality, exactly as Table 3 does, so the
    # figure and the table can never disagree.
    best_d = (nod.groupby(["representation", "n_components"]).auc.mean()
                 .reset_index().sort_values("auc", ascending=False)
                 .drop_duplicates("representation").set_index("representation")
                 .n_components.to_dict())
    sel = nod[[r.representation in best_d and r.n_components == best_d[r.representation]
               for r in nod.itertuples()]]
    piv = sel.pivot_table(index="accession", columns="representation", values="auc")

    NICE = {"DNABERT-S": "DNABERT-S\n(species-aware)",
            "compositional": "hand-crafted\ncomposition",
            "DNABERT-2": "DNABERT-2\n(general)",
            "NT-v2-50M": "Nucleotide Tr. v2\n(general)"}
    COL = {"DNABERT-S": CAT4[2], "compositional": CAT4[1],
           "DNABERT-2": CAT4[3], "NT-v2-50M": "#e87ba4"}
    order = [r for r in ["DNABERT-S", "compositional", "DNABERT-2", "NT-v2-50M"]
             if r in piv.columns]
    if len(order) < 2:
        return None

    rng = np.random.default_rng(0)
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.5),
                             gridspec_kw=dict(wspace=0.34, width_ratios=[1.15, 1.0]))

    # ---- panel A: per-genome AUC ------------------------------------------------------
    ax = axes[0]
    for i, rep in enumerate(order):
        v = piv[rep].dropna().values
        ax.scatter(v, np.full(len(v), i) + rng.uniform(-0.16, 0.16, len(v)),
                   s=7, color=COL[rep], alpha=0.35, linewidth=0, zorder=2)
        ax.scatter([v.mean()], [i], s=64, marker="|", color=INK, zorder=4, linewidth=1.6)
        ax.annotate(f"{v.mean():.3f}", (v.mean(), i - 0.34), ha="center", fontsize=7.2,
                    color=INK)
    ax.axvline(0.5, color=MUTED, lw=0.8, ls=(0, (3, 3)), zorder=1)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([NICE.get(r, r) for r in order], fontsize=7.4)
    ax.invert_yaxis()
    ax.set_xlabel("island-vs-backbone AUC (one point per chromosome)", fontsize=7.6)
    ax.set_title("A  Every representation, every genome", loc="left", fontsize=8.2)
    ax.set_ylim(len(order) - 0.45, -0.62)

    # ---- panel B: paired differences against DNABERT-S --------------------------------
    ax = axes[1]
    others = [r for r in order if r != "DNABERT-S"]
    for i, rep in enumerate(others):
        m = piv["DNABERT-S"].notna() & piv[rep].notna()
        d = (piv["DNABERT-S"][m] - piv[rep][m]).values
        ax.scatter(d, np.full(len(d), i) + rng.uniform(-0.16, 0.16, len(d)),
                   s=7, color=COL[rep], alpha=0.35, linewidth=0, zorder=2)
        ax.scatter([d.mean()], [i], s=64, marker="|", color=INK, zorder=4, linewidth=1.6)
        ax.annotate(f"{(d > 0).sum()}/{len(d)} genomes", (d.mean(), i - 0.34),
                    ha="center", fontsize=7.2, color=INK)
    ax.axvline(0.0, color=INK, lw=1.0, zorder=3)
    ax.set_yticks(range(len(others)))
    ax.set_yticklabels([NICE.get(r, r) for r in others], fontsize=7.4)
    ax.invert_yaxis()
    ax.set_xlabel("paired AUC advantage of DNABERT-S", fontsize=7.6)
    ax.set_title("B  What the species-aware stage buys", loc="left", fontsize=8.2)
    ax.set_ylim(len(others) - 0.45, -0.62)

    save(fig, FIGURES / "fig8_representation")
    return "fig8_representation"


FIGS = {"fig1": fig1_schematic, "fig2": fig2_compute, "fig3": fig3_pr,
        "fig4": fig4_boundaries, "fig5": fig5_synthetic, "fig6": fig6_ablation,
        "fig7": fig7_ct18_map, "fig8": fig8_representation}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*", default=None)
    args = ap.parse_args()
    use_paper_style()
    FIGURES.mkdir(parents=True, exist_ok=True)
    for key, fn in FIGS.items():
        if args.only and key not in args.only:
            continue
        try:
            out = fn()
        except Exception as exc:                              # noqa: BLE001
            print(f"  {key}: ERROR {type(exc).__name__}: {exc}")
            continue
        print(f"  {key}: {'wrote ' + out if out else 'skipped (inputs not ready)'}")


if __name__ == "__main__":
    main()
