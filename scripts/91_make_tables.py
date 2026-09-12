#!/usr/bin/env python
"""Generate every LaTeX table and every quoted number in the manuscript.

Nothing in the paper is typed by hand: `paper/numbers.tex` defines a macro for each figure
quoted in the prose, and the tables are written straight from `results/`. A value that has
not been computed yet is emitted as ``TBD``, so an unfinished pipeline shows up in the PDF
rather than silently going stale.

    python scripts/91_make_tables.py
"""
from __future__ import annotations
import json, sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, wilcoxon

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from sagegi.evaluate import rates                                  # noqa: E402
from sagegi.paths import BENCH, REPO, RESULTS                      # noqa: E402
from sagegi.viz import PRETTY                                      # noqa: E402

PAPER = REPO / "paper"
TABLES = PAPER / "tables"
TABLES.mkdir(parents=True, exist_ok=True)
MACROS: dict[str, str] = {}

#: published SSG-LUGIA numbers (Ibtehaz et al. 2021, Table 2) for the fidelity check
PUB_SSG = {"SSG-LUGIA-P": (60.280, 55.044, 52.376),
           "SSG-LUGIA-F": (55.105, 65.223, 55.079),
           "SSG-LUGIA-R": (40.011, 86.755, 49.062)}

#: SAGE-GI is the embedding-only pipeline. The fused variant is identical except for the
#: compositional block, which the ablation shows earns nothing (see DISPLAY below).
MAIN = "DNABERT-S-only"
FUSED = "SAGE-GI"
BASE = "SSG-LUGIA-F"

#: Configuration names are internal; these are what the paper calls them.
DISPLAY = {"DNABERT-S-only": "SAGE-GI",
           "SAGE-GI": "SAGE-GI + composition",
           "SAGE-GI-P": "SAGE-GI + composition (precision)",
           "SAGE-GI-R": "SAGE-GI + composition (recall)"}


def disp(name: str) -> str:
    return DISPLAY.get(name, name)


def _check(name):
    r"""LaTeX macro names may contain letters only -- a digit silently turns e.g.
    \d2fused into the accent command \d followed by the text '2fused', which leaks into
    the preamble and derails the whole compile."""
    if not name.isalpha():
        raise ValueError(f"invalid LaTeX macro name {name!r}: letters only")
    return name


def mac(name, value, fmt="{:.2f}", suffix=""):
    _check(name)
    bad = value is None or (isinstance(value, float) and not np.isfinite(value))
    MACROS[name] = "TBD" if bad else fmt.format(value) + suffix


def pval(name, p):
    """Format a p-value the way a journal expects."""
    _check(name)
    if p is None or not np.isfinite(p):
        MACROS[name] = "TBD"
    elif p < 1e-15:
        MACROS[name] = r"<10^{-15}"
    else:
        e = int(np.floor(np.log10(p)))
        MACROS[name] = rf"{p/10**e:.1f}\times 10^{{{e}}}"


def read(name):
    p = RESULTS / name
    return pd.read_csv(p) if p.exists() else None


def baseline_rates():
    p = BENCH / "islandpick_baselines.csv"
    if not p.exists():
        return None
    b = pd.read_csv(p)
    r = pd.DataFrame([rates(dict(TP=t.TP, FP=t.FP, TN=t.TN, FN=t.FN)) for t in b.itertuples()])
    return pd.concat([b[["tool", "accession"]].reset_index(drop=True), r], axis=1)


def paired(piv, a, b, col=None):
    """Mean difference and Wilcoxon p for two methods over the genomes both cover."""
    if a not in piv or b not in piv:
        return np.nan, np.nan, 0, 0
    x, y = piv[a], piv[b]
    m = x.notna() & y.notna()
    if m.sum() < 5:
        return np.nan, np.nan, 0, 0
    d = (x[m] - y[m])
    try:
        p = wilcoxon(x[m], y[m], zero_method="zsplit").pvalue
    except ValueError:
        p = np.nan
    return float(d.mean()), float(p), int((d > 0).sum()), int(m.sum())


# ------------------------------------------------------------------------ benchmark facts
def benchmark_facts():
    pos = pd.read_csv(BENCH / "islandpick_positive.csv")
    neg = pd.read_csv(BENCH / "islandpick_negative.csv")
    mac("nIslands", len(pos), "{:,d}")
    mac("nBackbone", len(neg), "{:,d}")
    mac("nIslandBp", pos.length.sum() / 1e6, "{:.1f}", r"\,Mbp")
    mac("nBackboneBp", neg.length.sum() / 1e6, "{:.1f}", r"\,Mbp")
    mac("nGenomes", pos.accession.nunique(), "{:d}")


# ------------------------------------------------------------- Table 1: main benchmark
ROW_PUB = ["sigi_hmm", "centroid", "pai_ida", "islandpath_dimob", "islandpath_dinuc",
           "alien_hunter"]
ROW_OURS = ["SSG-LUGIA-P", "SSG-LUGIA-F", "SSG-LUGIA-R", "DNABERT-2-only",
            "NT-v2-only", "SAGE-GI", "SAGE-GI-P", "SAGE-GI-R", "DNABERT-S-only"]


def table1():
    base = baseline_rates()
    ssg = read("ssg_lugia_per_genome.csv")
    ours = read("sage_gi_per_genome.csv")
    if base is None:
        return None
    have = set()
    for d in (ssg, ours):
        if d is not None:
            have |= set(d.accession)
    base = base[base.accession.isin(have or set(base.accession))]

    rows = []
    for t in ROW_PUB:
        g = base[base.tool == t]
        if len(g):
            rows.append((PRETTY.get(t, t), g.precision.mean(), g.recall.mean(),
                         g.f1.mean(), np.nan, "published"))
    ti = read("treasureisland_per_genome.csv")
    if ti is not None and len(ti):
        rows.append(("TreasureIsland", ti.precision.mean(), ti.recall.mean(),
                     ti.f1.mean(), ti.mabe.mean() / 1000 if "mabe" in ti else np.nan,
                     "published"))
    for src in (ssg, ours):
        if src is None:
            continue
        for m in ROW_OURS:
            g = src[src.method == m]
            if len(g):
                rows.append((m, g.precision.mean(), g.recall.mean(), g.f1.mean(),
                             g.mabe.mean() / 1000 if "mabe" in g else np.nan, "this work"))
    df = pd.DataFrame(rows, columns=["Method", "P", "R", "F1", "MABE", "src"])
    df = df.drop_duplicates("Method", keep="last")
    df.to_csv(RESULTS / "table1_main.csv", index=False)

    best = {c: df[c].max() for c in ("P", "R", "F1")}
    best["MABE"] = df.MABE.min()

    def cell(v, col):
        if not np.isfinite(v):
            return "--"
        s = f"{v:.2f}"
        return rf"\textbf{{{s}}}" if np.isclose(v, best[col]) else s

    L = [r"\begin{table}[t]",
         r"\caption{Genomic island prediction on the IslandPick benchmark",
         r"(\nGenomes{} chromosomes). Metrics are averaged over genomes, not over pooled",
         r"nucleotides. MABE is the median absolute boundary error in kb, which the",
         r"published studies did not report. Best value per column in bold.}\label{tab:main}",
         r"\begin{tabular*}{\columnwidth}{@{\extracolsep{\fill}}lrrrr@{}}", r"\toprule",
         r"Method & P (\%) & R (\%) & F1 (\%) & MABE (kb) \\", r"\midrule"]
    prev = None
    for r_ in df.itertuples():
        if prev == "published" and r_.src == "this work":
            L.append(r"\midrule")
        prev = r_.src
        L.append(f"{disp(r_.Method)} & {cell(r_.P,'P')} & {cell(r_.R,'R')} & "
                 f"{cell(r_.F1,'F1')} & {cell(r_.MABE,'MABE')} \\\\")
    L += [r"\botrule", r"\end{tabular*}",
          r"\begin{tablenotes}\item The first block is quoted from Additional File 6 of",
          r"Langille \emph{et al.}~\cite{langille2008islandpick}; the second was produced by",
          r"this study on the same chromosomes under the same protocol.\end{tablenotes}",
          r"\end{table}"]
    (TABLES / "table1_main.tex").write_text("\n".join(L), encoding="utf-8")

    for key, m in [("ssg", BASE), ("sg", MAIN), ("embS", "DNABERT-S-only"),
                   ("embD", "DNABERT-2-only"), ("ntOnly", "NT-v2-only"),
                   ("sgNT", "SAGE-GI(NT)")]:
        g = df[df.Method == m]
        if not len(g) and ours is not None:
            # not every configuration is shown in Table 1, but its numbers are still quoted
            h = ours[ours.method == m]
            if len(h):
                g = pd.DataFrame([{"P": h.precision.mean(), "R": h.recall.mean(),
                                   "F1": h.f1.mean(),
                                   "MABE": h.mabe.mean() / 1000 if "mabe" in h else np.nan}])
        for suf, col in [("P", "P"), ("R", "R"), ("F", "F1")]:
            mac(key + suf, g[col].iloc[0] if len(g) else None)
        mac(key + "MABE", g.MABE.iloc[0] if len(g) else None, "{:.1f}", r"\,kb")
    return df


# ------------------------------------------------- fidelity of the SSG-LUGIA re-implementation
def fidelity():
    ssg = read("ssg_lugia_per_genome.csv")
    if ssg is None:
        return
    worst, L = 0.0, []
    for m, (p, r, f) in PUB_SSG.items():
        g = ssg[ssg.method == m]
        if not len(g):
            continue
        d = g.f1.mean() - f
        worst = max(worst, abs(d))
        L.append(f"{m} & {g.precision.mean():.2f} & {g.recall.mean():.2f} & "
                 f"{g.f1.mean():.2f} & {p:.2f} & {r:.2f} & {f:.2f} & {d:+.2f} \\\\")
    mac("reproF", worst, "{:.2f}")
    head = [r"\begin{table}[t]",
            r"\caption{Fidelity of our SSG-LUGIA re-implementation. Left: this study, over",
            r"all \nGenomes{} chromosomes. Right: the values published in Table 2",
            r"of~\cite{ibtehaz2021ssglugia}.}\label{tab:fidelity}",
            r"\begin{tabular*}{\columnwidth}{@{\extracolsep{\fill}}lrrrrrrr@{}}", r"\toprule",
            r"& \multicolumn{3}{c}{this work} & \multicolumn{3}{c}{published} & \\",
            r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}",
            r"Variant & P & R & F1 & P & R & F1 & $\Delta$F1 \\", r"\midrule"]
    (TABLES / "table_fidelity.tex").write_text(
        "\n".join(head + L + [r"\botrule", r"\end{tabular*}", r"\end{table}"]),
        encoding="utf-8")


# ------------------------------------------------------------------ Table 2: ablations
ABL = [(MAIN, "full model (fusion, no detrending)"),
       ("SAGE-GI/no-refine", "without CUSUM boundary refinement"),
       ("SAGE-GI/start-assign", "leading-edge window projection"),
       ("SAGE-GI/detrend=250kb", "with positional detrending"),
       ("SAGE-GI(D2)", "DNABERT-2 substituted for DNABERT-S"),
       ("DNABERT-S-only", "embedding channel only"),
       ("DNABERT-2-only", "general-purpose embedding only"),
       ("SAGE-GI", "with the compositional channel added"),
       ("SAGE-GI(NT)", "Nucleotide Transformer v2 substituted for DNABERT-S"),
       ("NT-v2-only", "Nucleotide Transformer v2 embedding only"),
       ("Comp-only+center", "compositional channel only")]


def table2():
    df = read("sage_gi_per_genome.csv")
    if df is None:
        return
    piv = df.pivot_table(index="accession", columns="method", values="f1")
    rows = []
    for m, desc in ABL:
        g = df[df.method == m]
        if not len(g):
            continue
        d, p, wins, n = paired(piv, MAIN, m)
        rows.append((m, desc, g.precision.mean(), g.recall.mean(), g.f1.mean(),
                     g.mabe.mean() / 1000, d, p, wins, n))
    if not rows:
        return
    out = pd.DataFrame(rows, columns=["Method", "Description", "P", "R", "F1", "MABE",
                                      "dF1", "p", "wins", "n"])
    out.to_csv(RESULTS / "table2_ablation.csv", index=False)

    def fp(p):
        if not np.isfinite(p):
            return "--"
        if p < 1e-15:
            return r"$<10^{-15}$"
        e = int(np.floor(np.log10(p)))
        return rf"${p/10**e:.0f}\times 10^{{{e}}}$"

    L = [r"\begin{table*}[t]",
         r"\caption{Ablation of SAGE-GI on the IslandPick benchmark. $\Delta$F1 and its",
         r"$p$-value are from a paired Wilcoxon signed-rank test against the full model over",
         r"the \nGenomes{} chromosomes; \emph{wins} counts genomes on which the full model",
         r"scores higher.}\label{tab:ablation}",
         r"\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}lrrrrrrr@{}}", r"\toprule",
         r"Configuration & P (\%) & R (\%) & F1 (\%) & MABE (kb) & $\Delta$F1 & $p$ & wins \\",
         r"\midrule"]
    for r_ in out.itertuples():
        dd = "--" if r_.Method == MAIN else f"{r_.dF1:+.2f}"
        pp = "--" if r_.Method == MAIN else fp(r_.p)
        ww = "--" if r_.Method == MAIN else f"{r_.wins}/{r_.n}"
        L.append(f"{r_.Description} & {r_.P:.2f} & {r_.R:.2f} & {r_.F1:.2f} & "
                 f"{r_.MABE:.1f} & {dd} & {pp} & {ww} \\\\")
    L += [r"\botrule", r"\end{tabular*}", r"\end{table*}"]
    (TABLES / "table2_ablation.tex").write_text("\n".join(L), encoding="utf-8")

    for nm, m in [("noRefine", "SAGE-GI/no-refine"), ("detrend", "SAGE-GI/detrend=250kb"),
                  ("startAssign", "SAGE-GI/start-assign"), ("dtwoFused", "SAGE-GI(D2)")]:
        g = out[out.Method == m]
        mac(nm + "F", g.F1.iloc[0] if len(g) else None)
        mac(nm + "MABE", g.MABE.iloc[0] if len(g) else None, "{:.1f}", r"\,kb")
        if len(g):
            mac(nm + "Delta", -g.dF1.iloc[0])
            pval("p" + nm, g.p.iloc[0])
    g = out[out.Method == "SAGE-GI/detrend=250kb"]
    mac("detrendCost", g.dF1.iloc[0] if len(g) else None, "{:.1f}")

    # headline gains against the published baseline and the compositional channel
    ssg = read("ssg_lugia_per_genome.csv")
    if ssg is not None:
        both = piv.join(ssg.pivot_table(index="accession", columns="method",
                                        values="f1")[[BASE]])
        d, p, wins, n = paired(both, MAIN, BASE)
        mac("gainF", d, "{:.1f}")
        pval("pGain", p)
        mac("gainWins", wins, "{:d}")
        mac("gainN", n, "{:d}")
    d, p, _, _ = paired(piv, MAIN, "SAGE-GI(D2)")
    mac("dFsd", d, "+{:.1f}")
    pval("pFsd", p)
    d, p, _, _ = paired(piv, "DNABERT-S-only", "DNABERT-2-only")
    mac("dFsdSolo", d, "+{:.1f}")
    # same contrast against the third model, which is what generalises the claim beyond
    # the DNABERT pair
    d, p, _, _ = paired(piv, MAIN, "SAGE-GI(NT)")
    mac("dFsnt", d, "+{:.1f}")
    pval("pFsnt", p)
    # What the compositional channel is worth on top of the embedding. The two configurations
    # differ in exactly one setting, so this is the cleanest available test of whether the
    # hand-crafted block earns its place at all. Orientation is FUSED - MAIN, matching the
    # bootstrap interval in confidence_intervals(), so a negative value means the extra block
    # costs F1 and compWins counts the genomes on which the fusion is ahead.
    d, p, wins, n = paired(piv, FUSED, MAIN)
    mac("dFcomp", d, "{:+.2f}")
    pval("pFcomp", p)
    mac("compWins", wins, "{:d}")
    mac("compN", n, "{:d}")


# ------------------------------------------------------- Table 3: representation quality
def table3_representation():
    df = read("representation_auc.csv")
    if df is None:
        return
    nod = df[df.detrend_bp == 0]
    best = (nod.groupby(["representation", "n_components"])[["auc", "contrast"]]
              .mean().reset_index())
    pick = {}
    for rep in best.representation.unique():
        g = best[best.representation == rep].sort_values("auc", ascending=False)
        pick[rep] = (int(g.n_components.iloc[0]), g.auc.iloc[0], g.contrast.iloc[0])

    piv = nod.pivot_table(index="accession", columns=["representation", "n_components"],
                          values="auc")

    def col(rep):
        return piv[(rep, pick[rep][0])]

    L = [r"\begin{table}[t]",
         r"\caption{How well each representation separates curated islands from curated",
         r"backbone, before any detector is involved: Mann--Whitney AUC of a fixed robust",
         r"anomaly statistic, averaged over the \nGenomes{} chromosomes. \emph{Contrast} is",
         r"the ratio of the mean statistic inside islands to outside.}\label{tab:repr}",
         r"\begin{tabular*}{\columnwidth}{@{\extracolsep{\fill}}lrrr@{}}", r"\toprule",
         r"Representation & dims & AUC & contrast \\", r"\midrule"]
    order = ["DNABERT-S", "compositional", "DNABERT-2", "NT-v2-50M"]
    nice = {"DNABERT-S": "DNABERT-S embedding (species-aware)",
            "compositional": "hand-crafted composition",
            "DNABERT-2": "DNABERT-2 embedding (general-purpose)",
            "NT-v2-50M": "Nucleotide Transformer v2 (general-purpose)"}
    for rep in order:
        if rep not in pick:
            continue
        d, a, c = pick[rep]
        bold = r"\textbf{%s}" if rep == "DNABERT-S" else "%s"
        L.append(f"{nice[rep]} & {d} & {bold % f'{a:.3f}'} & {c:.2f} \\\\")
    L += [r"\botrule", r"\end{tabular*}", r"\end{table}"]
    (TABLES / "table3_representation.tex").write_text("\n".join(L), encoding="utf-8")

    if "DNABERT-S" in pick:
        mac("aucS", pick["DNABERT-S"][1], "{:.3f}")
        mac("contrastS", pick["DNABERT-S"][2], "{:.2f}")
        mac("dimS", pick["DNABERT-S"][0], "{:d}")
    if "compositional" in pick:
        mac("aucC", pick["compositional"][1], "{:.3f}")
        mac("contrastC", pick["compositional"][2], "{:.2f}")
    if "DNABERT-2" in pick:
        mac("aucD", pick["DNABERT-2"][1], "{:.3f}")
        mac("contrastD", pick["DNABERT-2"][2], "{:.2f}")
    if "NT-v2-50M" in pick:
        mac("aucNT", pick["NT-v2-50M"][1], "{:.3f}")
        mac("contrastNT", pick["NT-v2-50M"][2], "{:.2f}")

    def cmp(a, b):
        x, y = col(a), col(b)
        m = x.notna() & y.notna()
        return float((x[m] - y[m]).mean()), wilcoxon(x[m], y[m]).pvalue, \
            int((x[m] > y[m]).sum()), int(m.sum())

    if {"DNABERT-S", "DNABERT-2"} <= set(pick):
        d, p, w, n = cmp("DNABERT-S", "DNABERT-2")
        mac("dAUCsd", d, "+{:.3f}"); pval("pAUCsd", p)
        mac("winsSD", w, "{:d}"); mac("nRepr", n, "{:d}")
    if {"DNABERT-S", "compositional"} <= set(pick):
        d, p, w, n = cmp("DNABERT-S", "compositional")
        mac("dAUCsc", d, "+{:.3f}"); pval("pAUCsc", p)
    if {"compositional", "DNABERT-2"} <= set(pick):
        d, p, _, _ = cmp("compositional", "DNABERT-2")
        mac("dAUCcd", d, "+{:.3f}"); pval("pAUCcd", p)
    # The third model is what turns "DNABERT-S beats DNABERT-2" into a claim about
    # species-aware pre-training in general, so its comparison needs the same treatment.
    if {"DNABERT-S", "NT-v2-50M"} <= set(pick):
        d, p, w, n = cmp("DNABERT-S", "NT-v2-50M")
        mac("dAUCsnt", d, "+{:.3f}"); pval("pAUCsnt", p)
        mac("winsSNT", w, "{:d}")
    if {"compositional", "NT-v2-50M"} <= set(pick):
        d, p, _, _ = cmp("compositional", "NT-v2-50M")
        mac("dAUCcnt", d, "+{:.3f}"); pval("pAUCcnt", p)


# ---------------------------------------------------------------- Table 4: case studies
def table4_case_studies():
    df = read("case_studies.csv")
    if df is None:
        return
    from importlib import import_module
    pub = import_module("07_case_studies").CT18_PUBLISHED if False else {
        "GI-SVM": (89.5, 44.6, 59.6), "EGID": (77.9, 53.5, 63.4),
        "SIGI-HMM": (24.1, 55.6, 33.7), "IslandViewer": (65.4, 67.0, 66.2),
        "GIHunter": (82.7, 67.6, 74.4), "IslandPath-DIMOB": (55.3, 78.8, 65.0),
        "tRNAcc": (28.6, 99.3, 44.4), "IslandPick": (6.0, 100.0, 11.4)}
    ct = df[df.accession == "NC_003198.1"]
    L = [r"\begin{table}[t]",
         r"\caption{\emph{Salmonella} Typhi CT18, evaluated genome-wide against the 19",
         r"curated islands. The upper block is quoted from Table 3",
         r"of~\cite{ibtehaz2021ssglugia}, which collected it",
         r"from~\cite{lu2016computational}.}\label{tab:ct18}",
         r"\begin{tabular*}{\columnwidth}{@{\extracolsep{\fill}}lrrr@{}}", r"\toprule",
         r"Method & P (\%) & R (\%) & F1 (\%) \\", r"\midrule"]
    for k, (r_, p_, f_) in sorted(pub.items(), key=lambda kv: -kv[1][2]):
        L.append(f"{k} & {p_:.1f} & {r_:.1f} & {f_:.1f} \\\\")
    L.append(r"\midrule")
    for m in ["SSG-LUGIA-P", "SSG-LUGIA-F", "DNABERT-S-only", "SAGE-GI"]:
        g = ct[ct.method == m]
        if len(g):
            L.append(f"{disp(m)} & {g.precision.iloc[0]:.1f} & {g.recall.iloc[0]:.1f} & "
                     f"{g.f1.iloc[0]:.1f} \\\\")
    L += [r"\botrule", r"\end{tabular*}", r"\end{table}"]
    (TABLES / "table4_ct18.tex").write_text("\n".join(L), encoding="utf-8")

    g = ct[ct.method == MAIN]
    mac("ctEighteenF", g.f1.iloc[0] if len(g) else None, "{:.1f}")
    mac("ctEighteenP", g.precision.iloc[0] if len(g) else None, "{:.1f}")
    mac("ctEighteenR", g.recall.iloc[0] if len(g) else None, "{:.1f}")
    for acc, key in [("NC_002935.2", "diph"), ("NC_011770.1", "pae")]:
        s = df[(df.accession == acc) & (df.method == MAIN)]
        b = df[(df.accession == acc) & (df.method == "SSG-LUGIA-F")]
        mac(key + "F", s.f1.iloc[0] if len(s) else None, "{:.1f}")
        mac(key + "Base", b.f1.iloc[0] if len(b) else None, "{:.1f}")


# ------------------------------------------------------------ Table 5: cross-validation
def table5_cv():
    p = RESULTS / "tuning_sage_DNABERT-S_t1000.json"
    if not p.exists():
        return
    rep = json.loads(p.read_text())
    L = [r"\begin{table}[t]",
         r"\caption{In-sample versus grouped 5-fold cross-validated performance. The",
         r"in-sample column follows the tuning protocol of~\cite{ibtehaz2021ssglugia}; the",
         r"cross-validated column selects the configuration on four fifths of the",
         r"chromosomes and scores it on the held-out fifth.}\label{tab:cv}",
         r"\begin{tabular*}{\columnwidth}{@{\extracolsep{\fill}}lrrrrrr@{}}", r"\toprule",
         r"& \multicolumn{3}{c}{in-sample} & \multicolumn{3}{c}{5-fold held out} \\",
         r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}",
         r"Variant & P & R & F1 & P & R & F1 \\", r"\midrule"]
    name = {"precision": "SAGE-GI-P", "f1": "SAGE-GI-F", "recall": "SAGE-GI-R"}
    for r_ in rep:
        a, b = r_["in_sample"], r_["cv_heldout"]
        L.append(f"{name.get(r_['objective'], r_['objective'])} & {a['precision']:.2f} & "
                 f"{a['recall']:.2f} & {a['f1']:.2f} & {b['precision']:.2f} & "
                 f"{b['recall']:.2f} & {b['f1']:.2f} \\\\")
        if r_["objective"] == "f1":
            mac("cvF", b["f1"]); mac("cvP", b["precision"]); mac("cvR", b["recall"])
            mac("insF", a["f1"])
    L += [r"\botrule", r"\end{tabular*}", r"\end{table}"]
    (TABLES / "table5_cv.tex").write_text("\n".join(L), encoding="utf-8")


# ------------------------------------------------------------ Table 6: biological relevance
def table6_biology():
    df = read("biological_relevance.csv")
    if df is None:
        return
    order = ["curated islands (reference)", "SSG-LUGIA-F", "Comp-only+center", MAIN]
    cats = [("mobility", "mobility genes"), ("virulence", "virulence factors"),
            ("any_hgt", "any HGT marker")]
    L = [r"\begin{table}[t]",
         r"\caption{Functional enrichment inside predicted islands on \emph{Salmonella} Typhi",
         r"CT18. Odds ratio and one-sided Fisher exact $p$ for annotated genes whose product",
         r"matches each keyword set, inside predictions versus outside. Resistance genes are",
         r"omitted: in CT18 they sit on the IncHI1 plasmid pHCM1 rather than the chromosome, so",
         r"neither the curated islands nor any method enriches for them.}\label{tab:biology}",
         r"\begin{tabular*}{\columnwidth}{@{\extracolsep{\fill}}llrrr@{}}",
         r"\toprule",
         r"Method & Category & in (\%) & out (\%) & OR \\",
         r"\midrule"]
    for m in order:
        g = df[df.method == m]
        if not len(g):
            continue
        for i, (key, label) in enumerate(cats):
            r_ = g[g.category == key]
            if not len(r_):
                continue
            r_ = r_.iloc[0]
            L.append(f"{m if i == 0 else ''} & {label} & {r_.rate_in:.2f} & "
                     f"{r_.rate_out:.2f} & {r_.odds_ratio:.2f} " + r"\\")
        L.append(r"\addlinespace")
    L += [r"\botrule", r"\end{tabular*}", r"\end{table}"]
    (TABLES / "table6_biology.tex").write_text("\n".join(L), encoding="utf-8")

    for key, cat in [("bioMob", "mobility"), ("bioVir", "virulence"), ("bioAny", "any_hgt")]:
        for who, m in [("SG", MAIN), ("SSG", "SSG-LUGIA-F"),
                       ("Ref", "curated islands (reference)")]:
            g = df[(df.method == m) & (df.category == cat)]
            mac(key + who, g.odds_ratio.iloc[0] if len(g) else None, "{:.1f}")
        g = df[(df.method == MAIN) & (df.category == cat)]
        if len(g):
            pval("p" + key, g.p.iloc[0])


# --------------------------------------------- Table 7: benchmark-wide functional enrichment
def paired_ci(piv, a, b, n_boot=10000, seed=0):
    """Bootstrap 95% CI for the mean paired difference over genomes."""
    if a not in piv or b not in piv:
        return None, None, None
    x, y = piv[a], piv[b]
    m = x.notna() & y.notna()
    d = (x[m] - y[m]).to_numpy()
    if len(d) < 5:
        return None, None, None
    rng = np.random.default_rng(seed)
    boot = rng.choice(d, size=(n_boot, len(d)), replace=True).mean(axis=1)
    return float(d.mean()), float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))


def start_assign_effect():
    """What crediting a window to its centre rather than its leading bases actually buys.

    It removes a systematic leftward displacement of start coordinates. It does *not*
    measurably change F1 or median absolute boundary error, and the paper should not imply
    otherwise.
    """
    b = read("sage_gi_boundaries.csv")
    if b is None or "signed_start_err" not in b:
        return
    b = b[b.matched.astype(bool)]
    piv = b.pivot_table(index=["accession", "true_start"], columns="method",
                        values="signed_start_err")
    if MAIN not in piv or "SAGE-GI/start-assign" not in piv:
        return
    a, c = piv[MAIN], piv["SAGE-GI/start-assign"]
    m = a.notna() & c.notna()
    if m.sum() < 20:
        return
    shift = (a[m] - c[m]).mean()
    try:
        p = wilcoxon(a[m], c[m]).pvalue
    except ValueError:
        p = np.nan
    mac("startShift", shift / 1000, "{:.1f}", r"\,kb")
    pval("pStartShift", p)
    mac("startShiftN", int(m.sum()), "{:,d}")
    mac("signedStartMain", a[m].mean() / 1000, "{:.1f}", r"\,kb")


def refine_effect():
    r"""What CUSUM boundary refinement actually buys, measured edge by edge.

    It does *not* lower the median absolute boundary error --- it raises it slightly, and the
    manuscript said the opposite until we measured it. What it does is place far more edges
    very close to the truth: a change-point search commits, which pays off where the local
    signal is clear and costs where it is not. The honest summary is therefore the fraction of
    edges within a tight tolerance, tested edge-wise with an exact McNemar test on the paired
    binary outcomes (equivalently a two-sided sign test on the discordant pairs).
    """
    b = read("sage_gi_boundaries.csv")
    if b is None or "start_err" not in b:
        return
    b = b[b.matched.astype(bool)]
    ps = b.pivot_table(index=["accession", "true_start"], columns="method", values="start_err")
    pe = b.pivot_table(index=["accession", "true_start"], columns="method", values="end_err")
    NR = "SAGE-GI/no-refine"
    if MAIN not in ps or NR not in ps:
        return
    m = ps[MAIN].notna() & ps[NR].notna() & pe[MAIN].notna() & pe[NR].notna()
    if m.sum() < 20:
        return
    # both edges of every matched island, refined against unrefined
    fine = np.concatenate([ps.loc[m, MAIN].to_numpy(), pe.loc[m, MAIN].to_numpy()])
    coarse = np.concatenate([ps.loc[m, NR].to_numpy(), pe.loc[m, NR].to_numpy()])
    mac("refineEdgeN", len(fine), "{:,d}")

    for tol, key in [(1000, "One"), (2000, "Two")]:
        x, y = fine <= tol, coarse <= tol
        gained, lost = int((x & ~y).sum()), int((~x & y).sum())
        # exact McNemar == two-sided binomial test on the discordant pairs
        p = binomtest(min(gained, lost), gained + lost, 0.5).pvalue if gained + lost else np.nan
        mac("refine" + key + "Kb", 100 * x.mean(), "{:.1f}")
        mac("noRefine" + key + "Kb", 100 * y.mean(), "{:.1f}")
        pval("pRefine" + key + "Kb", p)
        mac("refine" + key + "Gained", gained, "{:d}")
        mac("refine" + key + "Lost", lost, "{:d}")

    # ...and the metric it does not improve, reported so the trade-off is explicit
    try:
        pval("pRefineMedian", wilcoxon(fine, coarse).pvalue)
    except ValueError:
        pass
    mac("refineMedian", np.median(fine) / 1000, "{:.1f}", r"\,kb")
    mac("noRefineMedian", np.median(coarse) / 1000, "{:.1f}", r"\,kb")


def boundary_vs_baseline():
    """Is SAGE-GI's boundary error actually lower than the baseline's? Paired, over the
    genomes both methods cover.

    Table 1 shows 18.0\\,kb against 18.8\\,kb, but those are unpaired means over slightly
    different genome sets. Paired on the chromosomes both cover the gap is much smaller and
    not significant, and the abstract must not claim it.
    """
    ours, ssg = read("sage_gi_per_genome.csv"), read("ssg_lugia_per_genome.csv")
    if ours is None or ssg is None or "mabe" not in ours:
        return
    a = ours[ours.method == MAIN].set_index("accession").mabe
    b = ssg[ssg.method == BASE].set_index("accession").mabe
    j = a.index.intersection(b.index)
    a, b = a[j], b[j]
    m = a.notna() & b.notna()
    if m.sum() < 5:
        return
    mac("mabeGain", (b[m] - a[m]).mean() / 1000, "{:.2f}", r"\,kb")
    pval("pMabeGain", wilcoxon(a[m], b[m]).pvalue)
    mac("mabeGainWins", int((a[m] < b[m]).sum()), "{:d}")
    mac("mabeGainN", int(m.sum()), "{:d}")


def edges_vs_baseline():
    r"""Edge placement against the published baseline, on the islands both methods match.

    This is where the "boundary resolution" claim in the title has to be earned or dropped,
    and the answer is genuinely two-sided: SAGE-GI places roughly three times as many edges
    within a kilobase of the truth, yet its *median* edge error is higher. A change-point
    search commits --- it lands almost exactly where the local signal is clear and further
    away where it is not, which sharpens the head of the error distribution and thickens its
    tail. Both halves are reported.
    """
    sg, ssg = read("sage_gi_boundaries.csv"), read("ssg_lugia_boundaries.csv")
    if sg is None or ssg is None or "start_err" not in sg:
        return
    key = ["accession", "true_start"]
    a = sg[sg.matched.astype(bool) & (sg.method == MAIN)].set_index(key)[["start_err", "end_err"]]
    b = ssg[ssg.matched.astype(bool) & (ssg.method == BASE)].set_index(key)[["start_err", "end_err"]]
    j = a.index.intersection(b.index)
    if len(j) < 20:
        return
    a, b = a.loc[j], b.loc[j]
    fa = np.concatenate([a.start_err.to_numpy(), a.end_err.to_numpy()])
    fb = np.concatenate([b.start_err.to_numpy(), b.end_err.to_numpy()])
    ok = ~(np.isnan(fa) | np.isnan(fb))
    fa, fb = fa[ok], fb[ok]
    mac("edgeN", len(fa), "{:,d}")
    mac("edgeIslandN", len(j), "{:,d}")
    for tol, k in [(1000, "One"), (2000, "Two")]:
        x, y = fa <= tol, fb <= tol
        g, l = int((x & ~y).sum()), int((~x & y).sum())
        mac("edge" + k + "Sg", 100 * x.mean(), "{:.1f}")
        mac("edge" + k + "Ssg", 100 * y.mean(), "{:.1f}")
        pval("pEdge" + k, binomtest(min(g, l), g + l, 0.5).pvalue if g + l else np.nan)
    # the half of the picture that does not favour us
    mac("edgeMedSg", np.median(fa) / 1000, "{:.1f}", r"\,kb")
    mac("edgeMedSsg", np.median(fb) / 1000, "{:.1f}", r"\,kb")
    pval("pEdgeMed", wilcoxon(fa, fb).pvalue)


def audit_self_comparisons():
    """Guard against the defect this script actually shipped.

    ``dFcomp`` was computed as ``paired(piv, MAIN, "DNABERT-S-only")`` while ``MAIN`` *was*
    ``"DNABERT-S-only"``, so the manuscript reported the compositional ablation as exactly
    +0.00 with p = 1.0 on 0 of 118 chromosomes --- a comparison of a configuration with
    itself, printed in the abstract. A difference of exactly zero next to a p-value of exactly
    one is the signature; fail loudly rather than emit it.
    """
    bad = []
    for name, val in MACROS.items():
        if val not in ("+0.00", "0.00", "-0.00"):
            continue
        stem = name[1:] if name.startswith("d") else name
        pk = next((k for k in ("p" + stem, "p" + stem[0].upper() + stem[1:]) if k in MACROS), None)
        if pk and MACROS[pk] == r"1.0\times 10^{0}":
            bad.append(f"{name}={val} with {pk}=1.0")
    if bad:
        raise SystemExit("self-comparison detected (a method compared with itself): "
                         + "; ".join(bad))


def confidence_intervals():
    """Intervals for the comparisons the abstract and contributions rest on."""
    ours = read("sage_gi_per_genome.csv")
    ssg = read("ssg_lugia_per_genome.csv")
    rep = read("representation_auc.csv")
    if ours is None:
        return
    piv = ours.pivot_table(index="accession", columns="method", values="f1")
    if ssg is not None:
        piv = piv.join(ssg.pivot_table(index="accession", columns="method", values="f1")[[BASE]])

    for key, a, b in [("gain", MAIN, BASE), ("comp", FUSED, MAIN),
                      ("fsd", MAIN, "DNABERT-2-only"), ("fsnt", MAIN, "NT-v2-only")]:
        mean, lo, hi = paired_ci(piv, a, b)
        if mean is not None:
            mac("ci" + key.capitalize() + "Lo", lo, "{:.1f}")
            mac("ci" + key.capitalize() + "Hi", hi, "{:.1f}")

    if rep is not None:
        nod = rep[rep.detrend_bp == 0]
        best = (nod.groupby(["representation", "n_components"]).auc.mean().reset_index()
                   .sort_values("auc", ascending=False).drop_duplicates("representation")
                   .set_index("representation").n_components.to_dict())
        sel = nod[[r.representation in best and r.n_components == best[r.representation]
                   for r in nod.itertuples()]]
        rp = sel.pivot_table(index="accession", columns="representation", values="auc")
        for key, a, b in [("Aucsd", "DNABERT-S", "DNABERT-2"),
                          ("Aucsnt", "DNABERT-S", "NT-v2-50M"),
                          ("Aucsc", "DNABERT-S", "compositional")]:
            mean, lo, hi = paired_ci(rp, a, b)
            if mean is not None:
                mac("ci" + key + "Lo", lo, "{:.3f}")
                mac("ci" + key + "Hi", hi, "{:.3f}")


def treasureisland_split():
    """Decompose TreasureIsland's score by whether it recognises the genome.

    Its classifier is a supervised SVM trained on the Benbow compilation, whose main
    component is the revised version of this benchmark, so a large part of our test set is
    related to its training set. The package reports a per-genome similarity check; splitting
    on it measures the effect instead of arguing about it.
    """
    ti = read("treasureisland_per_genome.csv")
    ours = read("sage_gi_per_genome.csv")
    if ti is None or "out_of_distribution" not in ti:
        return
    ood = ti.out_of_distribution.astype(bool)
    mac("tiFamF", ti[~ood].f1.mean() if (~ood).any() else None)
    mac("tiOodF", ti[ood].f1.mean() if ood.any() else None)
    mac("tiFamN", int((~ood).sum()), "{:d}")
    mac("tiOodN", int(ood.sum()), "{:d}")
    if ours is not None:
        g = ours[ours.method == MAIN].set_index("accession").f1
        fam = ti.loc[~ood, "accession"]
        oodacc = ti.loc[ood, "accession"]
        mac("sgFamF", g.reindex(fam).mean())
        mac("sgOodF", g.reindex(oodacc).mean())


def table8_protocol_sensitivity():
    """How much of each method's reported precision rests on unscored sequence.

    The official protocol grades ~13% of a chromosome. Recomputing precision with every base
    outside a curated island counted as a negative gives a lower bound; the *gap* between the
    two is what the benchmark cannot see. It is a property of the benchmark, not of any
    method, but methods differ enormously in how exposed they are to it.
    """
    df = read("protocol_sensitivity.csv")
    if df is None:
        return
    agg = (df.groupby("method")[["precision_official", "precision_strict", "called_frac"]]
             .mean())
    agg["drop"] = agg.precision_official - agg.precision_strict
    order = [m for m in ["TreasureIsland", MAIN, FUSED, "SAGE-GI-P", BASE, "SSG-LUGIA-P",
                         "Comp-only+center"] if m in agg.index]
    L = [r"\begin{table}[t]",
         r"\caption{How far each method's reported precision depends on sequence the protocol",
         r"never scores. \emph{Official} is the published protocol (curated regions only);",
         r"\emph{strict} counts every base outside a curated island as a negative, a lower",
         r"bound, since the unscored majority certainly contains uncurated islands. The",
         r"\emph{gap} is the informative quantity.}\label{tab:protocol}",
         r"\begin{tabular*}{\columnwidth}{@{\extracolsep{\fill}}lrrrr@{}}", r"\toprule",
         r"Method & called (\%) & P official & P strict & gap \\", r"\midrule"]
    for m in order:
        r_ = agg.loc[m]
        L.append(f"{disp(m)} & {r_.called_frac:.1f} & {r_.precision_official:.1f} & "
                 f"{r_.precision_strict:.1f} & {r_['drop']:.1f} \\\\")
    L += [r"\botrule", r"\end{tabular*}", r"\end{table}"]
    (TABLES / "table8_protocol.tex").write_text("\n".join(L), encoding="utf-8")

    for key, m in [("sens" + "Ti", "TreasureIsland"), ("sens" + "Sg", MAIN),
                   ("sens" + "Ssg", BASE)]:
        if m in agg.index:
            mac(key + "Drop", agg.loc[m, "drop"], "{:.0f}")
            mac(key + "Strict", agg.loc[m, "precision_strict"], "{:.1f}")


def table7_biology_benchmark():
    df = read("biology_benchmark.csv")
    if df is None:
        return
    order = [("curated islands", "curated islands (reference)"),
             ("SSG-LUGIA-P", "SSG-LUGIA-P"),
             ("SSG-LUGIA-F", "SSG-LUGIA-F"),
             ("SAGE-GI-P", "SAGE-GI (precision-oriented)"),
             (MAIN, "SAGE-GI"),
             ("both agree", "regions both methods call"),
             ("SAGE-GI only", "regions only SAGE-GI calls")]
    cats = [("mobility", "mobility"), ("virulence", "virulence"), ("any_hgt", "any marker")]

    L = [r"\begin{table}[t]",
         r"\caption{Functional enrichment pooled over the benchmark. Odds ratio for annotated",
         r"genes matching each keyword set, inside a region set versus outside it. The last two",
         r"rows separate the regions the two methods agree on from the regions only SAGE-GI",
         r"calls -- the ones responsible for its higher recall.}\label{tab:biologybench}",
         r"\begin{tabular*}{\columnwidth}{@{\extracolsep{\fill}}lrrr@{}}",
         r"\toprule",
         r"Region set & mobility & virulence & any marker \\",
         r"\midrule"]
    for key, label in order:
        g = df[df.method == key]
        if not len(g):
            continue
        vals = []
        for cat, _ in cats:
            r = g[g.category == cat]
            vals.append(f"{r.odds_ratio.iloc[0]:.2f}" if len(r) else "--")
        if key == "both agree":
            L.append(r"\midrule")
        L.append(f"{label} & {vals[0]} & {vals[1]} & {vals[2]} " + r"\\")
    L += [r"\botrule", r"\end{tabular*}", r"\end{table}"]
    (TABLES / "table7_biology_benchmark.tex").write_text("\n".join(L), encoding="utf-8")

    def orat(method, cat):
        r = df[(df.method == method) & (df.category == cat)]
        return r.odds_ratio.iloc[0] if len(r) else None

    mac("bioBenchRefMob", orat("curated islands", "mobility"))
    mac("bioBenchSSGMob", orat("SSG-LUGIA-F", "mobility"))
    mac("bioBenchSGMob", orat(MAIN, "mobility"))
    mac("bioBenchSGPMob", orat("SAGE-GI-P", "mobility"))
    mac("bioBenchBothMob", orat("both agree", "mobility"))
    mac("bioBenchOnlyMob", orat("SAGE-GI only", "mobility"))
    r = df[(df.method == "SAGE-GI only") & (df.category == "mobility")]
    if len(r):
        pval("pbioBenchOnly", r.p.iloc[0])
        mac("bioBenchOnlyN", r.n_marker_in.iloc[0], "{:,.0f}")


# ------------------------------------ how much of each chromosome does a method call?
def called_fractions():
    """The protocol scores only the curated regions; everything else is unpenalised."""
    pos = pd.read_csv(BENCH / "islandpick_positive.csv")
    neg = pd.read_csv(BENCH / "islandpick_negative.csv")
    glen = json.loads((BENCH / "genome_lengths.json").read_text())
    accs = sorted(set(pos.accession))
    pf = np.mean([pos[pos.accession == a].length.sum() / glen[a] for a in accs]) * 100
    nf = np.mean([neg[neg.accession == a].length.sum() / glen[a] for a in accs]) * 100
    mac("posFrac", pf, "{:.1f}")
    mac("negFrac", nf, "{:.1f}")
    mac("scoredFrac", pf + nf, "{:.1f}")

    def frac(df, method):
        g = df[df.method == method]
        if not len(g) or "pred_bp" not in g:
            return None
        return 100 * (g.pred_bp / g.genome_len).mean()

    sg = read("sage_gi_per_genome.csv")
    ssg = read("ssg_lugia_per_genome.csv")
    ti = read("treasureisland_per_genome.csv")
    if sg is not None:
        mac("sgFrac", frac(sg, MAIN), "{:.1f}")
    if ssg is not None:
        mac("ssgFrac", frac(ssg, BASE), "{:.1f}")
    if ti is not None:
        mac("tiF", ti.f1.mean())
        mac("tiP", ti.precision.mean())
        mac("tiR", ti.recall.mean())
        mac("tiFrac", 100 * (ti.pred_bp / ti.genome_len).mean(), "{:.1f}")
        mac("tiN", len(ti), "{:d}")
        mac("tiOOD", int(ti.out_of_distribution.sum()) if "out_of_distribution" in ti else None,
            "{:d}")


# ------------------------------------------------------------------- synthetic insertions
def synthetic_numbers():
    df = read("synthetic_insertions.csv")
    if df is None:
        return
    for key, m in [("synSG", MAIN), ("synComp", "Comp-only"),
                   ("synS", "DNABERT-S-only"), ("synDtwo", "DNABERT-2-only")]:
        g = df[df.method == m]
        mac(key, g.detection_rate.mean() if len(g) else None, "{:.1f}", r"\%")
        s = g[g.insert_len <= 10000]
        mac(key + "Short", s.detection_rate.mean() if len(s) else None, "{:.1f}", r"\%")
        mac(key + "MABE", g.mabe.mean() / 1000 if len(g) else None, "{:.1f}", r"\,kb")


#: every macro the manuscript may quote. Pre-registering them as TBD means a missing result
#: shows up in the PDF as "TBD" instead of breaking the compile with an undefined control
#: sequence, and makes it obvious which pipeline stage has not been run.
EXPECTED = [
    "nIslands", "nBackbone", "nIslandBp", "nBackboneBp", "nGenomes", "reproF",
    "aucS", "aucC", "aucD", "contrastS", "contrastC", "contrastD", "dimS",
    "dAUCsd", "pAUCsd", "dAUCsc", "pAUCsc", "dAUCcd", "pAUCcd", "winsSD", "nRepr",
    "gainF", "pGain", "gainWins", "gainN", "dFsd", "pFsd", "dFsdSolo",
    "detrendCost", "pdetrend", "detrendF", "detrendMABE",
    "noRefineF", "noRefineMABE", "noRefineDelta", "pnoRefine",
    "startAssignF", "startAssignMABE", "startAssignDelta", "pstartAssign",
    "dtwoFusedF", "dtwoFusedMABE", "dtwoFusedDelta", "pdtwoFused",
    "cvF", "cvP", "cvR", "insF",
    "ctEighteenF", "ctEighteenP", "ctEighteenR", "diphF", "diphBase", "paeF", "paeBase",
    "synSG", "synComp", "synS", "synDtwo",
    "synSGShort", "synCompShort", "synSShort", "synDtwoShort",
    "synSGMABE", "synCompMABE", "synSMABE", "synDtwoMABE",
    "bioMobSG", "bioMobSSG", "bioMobRef", "pbioMob",
    "bioVirSG", "bioVirSSG", "bioVirRef", "pbioVir",
    "bioAnySG", "bioAnySSG", "bioAnyRef", "pbioAny",
    "aucNT", "contrastNT", "dAUCsnt", "pAUCsnt", "ntOnlyF", "ntOnlyP", "ntOnlyR",
    "ntOnlyMABE", "sgNTF", "sgNTP", "sgNTR", "sgNTMABE", "winsSNT", "dAUCcnt", "pAUCcnt",
    "dFsnt", "pFsnt", "dFcomp", "pFcomp", "compWins", "compN",
    "tiF", "tiP", "tiR", "tiFrac", "sgFrac", "ssgFrac", "scoredFrac",
    "bioBenchRefMob", "bioBenchSSGMob", "bioBenchSGMob", "bioBenchSGPMob",
    "bioBenchBothMob", "bioBenchOnlyMob", "pbioBenchOnly", "bioBenchOnlyN",
    "posFrac", "negFrac", "tiN", "tiOOD",
    "startShift", "pStartShift", "startShiftN", "signedStartMain",
    "refineEdgeN", "refineOneKb", "noRefineOneKb", "pRefineOneKb",
    "refineOneGained", "refineOneLost",
    "refineTwoKb", "noRefineTwoKb", "pRefineTwoKb", "refineTwoGained", "refineTwoLost",
    "refineMedian", "noRefineMedian", "pRefineMedian",
    "mabeGain", "pMabeGain", "mabeGainWins", "mabeGainN",
    "edgeN", "edgeIslandN", "edgeOneSg", "edgeOneSsg", "pEdgeOne",
    "edgeTwoSg", "edgeTwoSsg", "pEdgeTwo",
    "edgeMedSg", "edgeMedSsg", "pEdgeMed",
    "ciGainLo", "ciGainHi", "ciCompLo", "ciCompHi",
    "ciFsdLo", "ciFsdHi", "ciFsntLo", "ciFsntHi",
    "ciAucsdLo", "ciAucsdHi", "ciAucsntLo", "ciAucsntHi",
    "ciAucscLo", "ciAucscHi",
    "sensTiDrop", "sensTiStrict", "sensSgDrop", "sensSgStrict",
    "sensSsgDrop", "sensSsgStrict",
    "tiFamF", "tiOodF", "tiFamN", "tiOodN", "sgFamF", "sgOodF",
]


def main():
    for _k in EXPECTED:
        MACROS[_k] = "TBD"
    benchmark_facts()
    table1()
    fidelity()
    table2()
    table3_representation()
    table4_case_studies()
    table5_cv()
    table6_biology()
    start_assign_effect()
    refine_effect()
    boundary_vs_baseline()
    edges_vs_baseline()
    confidence_intervals()
    treasureisland_split()
    table7_biology_benchmark()
    table8_protocol_sensitivity()
    called_fractions()
    synthetic_numbers()
    audit_self_comparisons()

    # Always emit every table file, so an unfinished stage shows a visible placeholder
    # instead of breaking the compile with "File not found".
    ALL_TABLES = {
        "table1_main.tex": "IslandPick benchmark results",
        "table_fidelity.tex": "SSG-LUGIA re-implementation fidelity",
        "table2_ablation.tex": "component ablation",
        "table3_representation.tex": "representation quality",
        "table4_ct18.tex": "case study: Salmonella Typhi CT18",
        "table5_cv.tex": "cross-validated performance",
        "table6_biology.tex": "functional enrichment of predicted islands",
        "table7_biology_benchmark.tex": "benchmark-wide functional enrichment",
    }
    for fname, what in ALL_TABLES.items():
        f = TABLES / fname
        if not f.exists():
            f.write_text("\n".join([
                r"\begin{table}[t]",
                r"\caption{" + what + r" --- pipeline stage not yet run.}",
                r"\centering\emph{TBD}",
                r"\end{table}", ""]), encoding="utf-8")
            print(f"  placeholder: {fname}")

    body = ["% Auto-generated by scripts/91_make_tables.py --- do not edit by hand.",
            r"\makeatletter"]
    for k, v in sorted(MACROS.items()):
        body.append(rf"\providecommand{{\{k}}}{{}}\renewcommand{{\{k}}}{{{v}}}")
    body.append(r"\makeatother")
    (PAPER / "numbers.tex").write_text("\n".join(body) + "\n", encoding="utf-8")

    ready = sum(1 for v in MACROS.values() if v != "TBD")
    print(f"numbers.tex: {ready}/{len(MACROS)} macros resolved")
    missing = sorted(k for k, v in MACROS.items() if v == "TBD")
    if missing:
        print("  still TBD:", ", ".join(missing))
    for f in sorted(TABLES.glob("*.tex")):
        print("  table:", f.name)


if __name__ == "__main__":
    main()
