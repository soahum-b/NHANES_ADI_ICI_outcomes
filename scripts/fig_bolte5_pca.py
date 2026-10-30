#!/usr/bin/env python3
"""
fig_bolte5_pca.py — PCA biplots of diet composition in the poorest vs richest PIR
quartile, on the BOLTE grouping (16 groups).

Reuses _biplot from fig8_pca_biplot.py so the Bolte panels and the primary ones are
drawn by the same code and cannot drift apart; only the input tables differ
(results/bolte/pca_*.csv from 29_bolte_panels.R).

Q4 is the income top-code group (PIR >= 5.0): NHANES censors PIR at 5, and the
survey-weighted 75th percentile IS exactly 5.0, so the richest quartile is defined by
the ceiling rather than by a quantile cut.

Output: figures/bolte/bolte5_pca.{png,pdf,svg}
"""
import pandas as pd, matplotlib.pyplot as plt, os, sys
from matplotlib.lines import Line2D
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figsave import save_fig

# set by the DHQ wrappers; defaults keep the FPED-component figures unchanged
RESULTS_DIR = 'results/bolte'
OUT_PNG = 'figures/bolte/bolte5_pca.png'
GROUPING = 'Bolte et al. 2023 grouping'

from fig8_pca_biplot import _biplot, BLUE, GREY, RED, POOR, RICH, POOR_BG, RICH_BG


def main():
    load = pd.read_csv(f"{RESULTS_DIR}/pca_loadings.csv")
    var = pd.read_csv(f"{RESULTS_DIR}/pca_variance.csv").set_index("group")
    sco = pd.read_csv(f"{RESULTS_DIR}/pca_scores.csv")

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.2))
    panels = [(axes[0], "Q1", "PIR Q1 (poorest)", POOR, POOR_BG),
              (axes[1], "Q4", "PIR Q4 (richest, top-coded)", RICH, RICH_BG)]
    for ax, g, title, fc, bg in panels:
        L = load[load.group == g].sort_values("num")
        n_g = int(var.loc[g, "n"]) if "n" in var.columns else None
        ttl = (f"PCA of diet composition \u2014 {title}"
               + (f", n = {n_g:,}" if n_g is not None else ""))
        _biplot(ax, L, var.loc[g], sco[sco.group == g], ttl, frame_color=fc, frame_bg=bg)

    LA = load[load.group == "Q1"].sort_values("num")
    fig.text(0.995, 0.46, "\n".join(f"{int(r.num)}. {r.label}" for _, r in LA.iterrows()),
             fontsize=8, va="center", ha="right", family="DejaVu Sans")
    handles = [Line2D([0], [0], color=BLUE, lw=2, marker=">", label="PC1 Spearman |\u03c1| > 0.5"),
               Line2D([0], [0], color=GREY, lw=2, marker=">", label="PC2 Spearman |\u03c1| > 0.5"),
               Line2D([0], [0], color=RED,  lw=2, marker=">", label="Weak correlation")]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.995, 0.99),
               frameon=False, fontsize=8, title="Spearman correlation", title_fontsize=8.5)
    fig.suptitle("Diet-component structure within income groups \u2014 " + GROUPING, y=1.02, fontsize=11, fontweight="bold")
    # The background cloud is a subsample (file size); every loading, variance figure and
    # arrow uses the full group. Without this a reader sees "n = 6,244" above a few hundred
    # points and reads the density as real.
    if "n_plotted" in var.columns:
        npl = int(var["n_plotted"].max())
        fig.text(0.01, 0.015,
                 f"Background points are a random subsample of {npl:,} participants per panel, "
                 f"shown for legibility; component loadings, variance explained and arrow "
                 f"directions are computed on the full group (n in panel titles).",
                 fontsize=7.6, style="italic", color="#555", ha="left", va="bottom")
    fig.tight_layout(rect=[0, 0, 0.82, 1])
    os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
    save_fig(fig, OUT_PNG)
    print("saved", OUT_PNG)


if __name__ == "__main__":
    main()
