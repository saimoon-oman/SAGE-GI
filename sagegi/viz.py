r"""Figure style and palette.

Colours come from a colour-vision-deficiency-validated categorical palette.  The eight-slot
order clears the adjacent-pair gates (worst CVD Delta-E 9.1 under protanopia, worst
normal-vision Delta-E 19.6, OKLab x100) and is used for bars and lines, where only
neighbouring series touch.  Scatter plots put every pair on screen at once, so they use the
four-slot subset blue/orange/aqua/violet, which clears the stricter all-pairs gates
(CVD 9.2, normal-vision 16.3) - and which happens to be exactly the four feature
representations compared in this study.

Three light-mode slots sit below 3:1 contrast against the page, so every chart that uses
them also carries direct value labels or a companion table, never colour alone.
"""
from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt

# eight-slot categorical order - adjacent-pair safe (bars, lines, stacks)
CAT8 = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
        "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
# four-slot subset - all-pairs safe (scatter, small multiples)
CAT4 = ["#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"]

INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8880"
GRID = "#e3e2dd"
SURFACE = "#ffffff"

#: Stable colour per method, so a figure that drops a series never repaints the others.
#: Colour encodes the *feature representation family*, never the rank of a bar:
#: blue = proposed fusion, orange = hand-crafted composition, aqua/violet = the two
#: foundation models on their own, neutral grey = previously published tools, which are
#: reference points rather than series of interest.  Within a family the variants are
#: separated by the categorical axis label (and by marker shape in scatter plots), not by
#: inventing tints - generated hues are what broke the colour-vision gates.
METHOD_COLOR = {
    "SAGE-GI":            CAT4[0],
    "SAGE-GI-P":          CAT4[0],
    "SAGE-GI-F":          CAT4[0],
    "SAGE-GI-R":          CAT4[0],
    "DNABERT-S-only":     CAT4[2],
    "DNABERT-2-only":     CAT4[3],
    "Comp-only":          CAT4[1],
    "SSG-LUGIA-F":        CAT4[1],
    "SSG-LUGIA-P":        CAT4[1],
    "SSG-LUGIA-R":        CAT4[1],
}
#: published tools are drawn in the muted neutral
BASELINE_COLOR = MUTED
#: secondary encoding for scatter plots, so identity never rests on colour alone
METHOD_MARKER = {
    "SAGE-GI": "o", "SAGE-GI-P": "o", "SAGE-GI-F": "P", "SAGE-GI-R": "X",
    "DNABERT-S-only": "s", "DNABERT-2-only": "^",
    "Comp-only": "D", "SSG-LUGIA-F": "D", "SSG-LUGIA-P": "v", "SSG-LUGIA-R": "<",
}

PRETTY = {
    "sigi_hmm": "SIGI-HMM", "centroid": "Centroid", "islandpath_dimob": "IslandPath-DIMOB",
    "pai_ida": "PAI-IDA", "islandpath_dinuc": "IslandPath-DINUC",
    "alien_hunter": "AlienHunter",
}


def use_paper_style(base: float = 8.5) -> None:
    """Journal-figure defaults: recessive axes, thin marks, no chartjunk."""
    mpl.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.facecolor": SURFACE,
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "font.family": "DejaVu Sans",
        "font.size": base,
        "axes.titlesize": base + 1.5,
        "axes.titleweight": "semibold",
        "axes.labelsize": base,
        "axes.labelcolor": INK2,
        "axes.edgecolor": GRID,
        "axes.linewidth": 0.8,
        "axes.grid": True,
        "axes.axisbelow": True,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "xtick.labelsize": base - 0.5,
        "ytick.labelsize": base - 0.5,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "legend.frameon": False,
        "legend.fontsize": base - 0.5,
        "lines.linewidth": 2.0,
        "lines.markersize": 5.5,
        "text.color": INK,
        "pdf.fonttype": 42,          # embed TrueType so the PDF is editable/searchable
        "ps.fonttype": 42,
    })


def color_for(method: str) -> str:
    """Colour by representation family; anything unrecognised is a published baseline."""
    return METHOD_COLOR.get(method, BASELINE_COLOR)


def marker_for(method: str) -> str:
    return METHOD_MARKER.get(method, "o")


def label_bars(ax, bars, fmt="{:.1f}", dy: float = 0.6, size: float = 7.0) -> None:
    """Direct value labels - required for the low-contrast palette slots."""
    for b in bars:
        h = b.get_height()
        ax.text(b.get_x() + b.get_width() / 2, h + dy, fmt.format(h),
                ha="center", va="bottom", fontsize=size, color=INK2)


def save(fig, path, also_png: bool = True) -> None:
    """Write the figure once, then mirror it where the paper and the website expect it.

    The manuscript includes figures from ``paper/figures/`` and the site serves them from
    ``docs/assets/``; copying here keeps all three in step without a manual step that is
    easy to forget after a re-run.
    """
    import shutil
    from pathlib import Path

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".pdf"))
    if also_png:
        fig.savefig(path.with_suffix(".png"))
    plt.close(fig)

    repo = path.parent.parent
    for dest, ext in ((repo / "paper" / "figures", ".pdf"), (repo / "docs" / "assets", ".png")):
        src = path.with_suffix(ext)
        if dest.is_dir() and src.exists():
            shutil.copy2(src, dest / src.name)
