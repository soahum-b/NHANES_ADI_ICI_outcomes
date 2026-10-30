#!/usr/bin/env python3
"""
fig_bolte3_panels.py — Bolte-parameterisation version of the hPDI/uPDI x SES panel
figure: healthful (top row) and unhealthful (bottom row) plant-based indices against
PIR, education and family income, with the design-based significance of each gradient.

Scores are hPDI_b / uPDI_b from 23_pdi_bolte.R (Bolte et al. 2023 grouping: nuts+legumes
collapsed, meat+poultry collapsed, 16 of their 17 groups constructible in NHANES).
Reads results/bolte/{pd_main,trend_tests}.csv  (29_bolte_panels.R).
Output: figures/bolte/bolte3_pdi_vs_ses.{png,pdf,svg}
"""
import pandas as pd, matplotlib.pyplot as plt, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figsave import save_fig

# set by the DHQ wrappers; defaults keep the FPED-component figures unchanged
RESULTS_DIR = 'results/bolte'
SCORES = ['hPDI_b', 'uPDI_b']   # column names in pd_main; wrappers override
OUT_PNG = 'figures/bolte/bolte3_pdi_vs_ses.png'
GROUPING = 'Bolte et al. 2023 food grouping'


COL = {"hPDI_b": "#7B4FA8", "uPDI_b": "#C0603A"}
TITLE = {"hPDI_b": "Healthful Plant-Based (hPDI)", "uPDI_b": "Unhealthful Plant-Based (uPDI)"}
ORDER_PIR = ["0-1","1-2","2-3","3-4","4-5"]
ORDER_EDU = ["<9th","9-11th","HS/GED","SomeCollege","CollegeGrad"]
LAB_EDU   = ["<9th","9-11th","HS\nGED","Some\ncollege","College\ngrad"]
ORDER_INC = ["<5k","5-10k","10-15k","15-20k","20-25k","25-35k",
             "35-45k","45-55k","55-65k","65-75k","75-100k","100k+"]


def _panel(ax, pd_main, ann, score, spec, ylim, letter=None,
           show_xlab=True, show_ylab=True, ylab=None):
    """Draw ONE score x SES panel. Factored out of main() so the composite figure and
    the standalone per-panel figures are drawn by identical code and cannot drift."""
    by, order, xlab, xticklab, rot = spec
    c = COL[score]
    d = pd_main[(pd_main.score == score) & (pd_main["by"] == by)].copy()
    d["group"] = pd.Categorical(d["group"], categories=order, ordered=True)
    d = d.sort_values("group")
    x = range(len(d))
    ax.errorbar(x, d["mean"], yerr=[d["mean"]-d.ci_l, d.ci_u-d["mean"]],
                fmt="o-", color=c, capsize=3, lw=1.6, ms=5)
    ax.set_xticks(list(x))
    ax.set_xticklabels(xticklab if xticklab else order,
                       rotation=rot, ha="right" if rot else "center", fontsize=8)
    if show_xlab:
        ax.set_xlabel(xlab, fontsize=9.5, fontweight="bold")
    if show_ylab:
        ax.set_ylabel(ylab, fontsize=10, fontweight="bold")
    # x-only margins: a bare margins() would re-expand y and undo the shared ylim
    ax.margins(x=0.08)
    ax.set_ylim(*ylim)
    # Put the annotation where the line ISN'T. A RISING series starts low on the left,
    # vacating the TOP-left; a FALLING series starts high, vacating the BOTTOM-left.
    rising = d["mean"].iloc[-1] >= d["mean"].iloc[0]
    ax.text(0.03, 0.95 if rising else 0.05, ann(score, by),
            transform=ax.transAxes, va="top" if rising else "bottom",
            ha="left", fontsize=8.5, fontweight="bold")
    if letter:
        ax.text(-0.02, 1.06, letter, transform=ax.transAxes,
                fontsize=13, fontweight="bold", va="bottom", ha="right")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    return d


def main():
    pd_main = pd.read_csv(f"{RESULTS_DIR}/pd_main.csv")
    tt = pd.read_csv(f"{RESULTS_DIR}/trend_tests.csv")

    def ann(score, by):
        r = tt[(tt.score == score) & (tt["by"] == by)]
        if r.empty:
            return ""
        r = r.iloc[0]
        kind = "Trend" if r.test == "linear trend" else "Assoc."
        return f"{kind} FDR p {r.p_fdr_disp}" if str(r.p_fdr_disp).startswith("<") \
               else f"{kind} FDR p = {r.p_fdr_disp}"

    specs = [("pir_bin", ORDER_PIR, "PIR Ratio", None, 0),
             ("educ_f", ORDER_EDU, "Education", LAB_EDU, 0),
             ("income_lab", ORDER_INC, "Annual Family Income, USD", None, 55)]

    fig, axes = plt.subplots(2, 3, figsize=(13.2, 7.0))

    # ONE shared y-range across ALL SIX panels, not per panel and not per row.
    # Without this matplotlib autoscales each panel independently, so six different
    # ranges make the gradients look steeper or flatter than each other purely from
    # axis scaling. Legitimate here because hPDI and uPDI always share ONE scale --
    # 16-80 for the 16-group default, 17-85 when the 17-group wrappers override
    # SCORES -- and their CI extremes sit in one narrow band either way (47.9-53.6
    # for the 16-group indices, 52.6-62.1 for the 17-group), so a common axis
    # costs almost nothing in resolution and makes the two indices directly
    # comparable. Derived from the CI extremes over every panel, padded for the
    # significance annotation which sits top-left or bottom-left.
    plotted = pd_main[pd_main.score.isin(SCORES)]
    ylo, yhi = plotted.ci_l.min(), plotted.ci_u.max()
    pad = 0.10 * (yhi - ylo)
    YLIM = (ylo - pad, yhi + pad)

    letters = iter("ABCDEF")
    for row, score in enumerate(SCORES):
        for col, spec in enumerate(specs):
            _panel(axes[row, col], pd_main, ann, score, spec, YLIM, next(letters),
                   show_xlab=(row == 1), show_ylab=(col == 0),
                   ylab=f"{'hPDI' if row == 0 else 'uPDI'} Score")
        axes[row, 1].set_title(TITLE[score], fontsize=10.5, fontweight="bold",
                               color=COL[score], pad=12)

    fig.suptitle("Plant-based diet indices and socioeconomic position \u2014 " + GROUPING,
                 y=1.00, fontsize=11.5, fontweight="bold")
    fig.tight_layout()
    os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
    save_fig(fig, OUT_PNG)
    print("saved", OUT_PNG)


if __name__ == "__main__":
    main()
