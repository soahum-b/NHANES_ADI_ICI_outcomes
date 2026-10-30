#!/usr/bin/env python3
"""
fig8_pca_biplot.py / (shared builder) — PCA biplots of the 18 plant-based diet
components within income extreme groups, in the style of the example: numbered
arrows coloured by Spearman correlation with the PCs (blue |rho(PC1)|>0.5, grey
|rho(PC2)|>0.5, red weak), a numbered food legend, and a light individual-score
scatter behind. Reads results/pca_{loadings,variance,scores}.csv (from 12_pca_pdi.R).

Running this file builds BOTH deliverables:
  figures/fig8_pca_quartiles.png  (Q1 vs Q4)
  figures/fig9_pca_tertiles.png   (T1 vs T3)
Run from project root.
"""
import pandas as pd, numpy as np, matplotlib.pyplot as plt, os
from matplotlib.lines import Line2D
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figsave import save_fig

BLUE, GREY, RED = "#1f6fd0", "#8a8a8a", "#d62728"
CMAP = {"PC1": BLUE, "PC2": GREY, "weak": RED}
# income-group frame colours (match the project palette: poorest orange, richest blue)
POOR, RICH = "#c0603a", "#1f5fa8"
POOR_BG, RICH_BG = "#fdfbfa", "#fafcfe"   # near-white group tints (see note above)
FRAME_LW = 1.4

def _biplot(ax, load, var, scores, title, frame_color=None, frame_bg=None):
    # coloured background tint + boundary frame to distinguish the income groups
    if frame_bg is not None:
        ax.set_facecolor(frame_bg)
    if frame_color is not None:
        for sp in ax.spines.values():
            sp.set_visible(True); sp.set_edgecolor(frame_color); sp.set_linewidth(FRAME_LW)
    # background individual scores
    ax.scatter(scores["pc1"], scores["pc2"], s=9, color="#9a9a9a", alpha=0.5, zorder=1, edgecolor="none")
    ax.axhline(0, color="#d9d9d9", lw=0.8, zorder=0); ax.axvline(0, color="#d9d9d9", lw=0.8, zorder=0)
    # arrows: loading direction = (rho_pc1, rho_pc2), scaled to fill the score cloud
    R = 2.6
    for _, r in load.iterrows():
        x, y = r["rho_pc1"] * R, r["rho_pc2"] * R
        ax.annotate("", xy=(x, y), xytext=(0, 0),
                    arrowprops=dict(arrowstyle="-|>", color=CMAP[r["color"]], lw=1.6, alpha=0.9), zorder=3)
        # number label just past the arrow tip
        lx, ly = x * 1.12, y * 1.12
        ax.text(lx, ly, str(int(r["num"])), fontsize=8, fontweight="bold", color="#222",
                ha="center", va="center", zorder=4)
    lim = R * 1.28
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
    ax.set_xlabel(f"PC1 ({var['pc1_var']:.1f}%)", fontsize=9)
    ax.set_ylabel(f"PC2 ({var['pc2_var']:.1f}%)", fontsize=9)
    ax.set_title(title, fontsize=10, fontweight="bold",
                 color=(frame_color if frame_color is not None else "#000"))

def _panel_pair(gA, gB, titleA, titleB, outfile, suptitle):
    load = pd.read_csv("results/pca_loadings.csv")
    var = pd.read_csv("results/pca_variance.csv").set_index("group")
    sco = pd.read_csv("results/pca_scores.csv")
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.2))
    # left panel = poorest (orange frame), right = richest (blue frame)
    frames = [(POOR, POOR_BG), (RICH, RICH_BG)]
    for (ax, g, title), (fc, bg) in zip([(axes[0], gA, titleA), (axes[1], gB, titleB)], frames):
        L = load[load.group == g].sort_values("num")
        _biplot(ax, L, var.loc[g], sco[sco.group == g], title, frame_color=fc, frame_bg=bg)
    # shared numbered legend (food names) down the right, from group A's ordering
    LA = load[load.group == gA].sort_values("num")
    legend_txt = "\n".join(f"{int(r.num)}. {r.label}" for _, r in LA.iterrows())
    fig.text(0.995, 0.5, legend_txt, fontsize=8, va="center", ha="right", family="DejaVu Sans")
    # colour key for the arrows
    handles = [Line2D([0],[0], color=BLUE, lw=2, marker=">", label="PC1 Spearman |ρ| > 0.5"),
               Line2D([0],[0], color=GREY, lw=2, marker=">", label="PC2 Spearman |ρ| > 0.5"),
               Line2D([0],[0], color=RED,  lw=2, marker=">", label="Weak correlation")]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.995, 0.99),
               frameon=False, fontsize=8, title="Spearman correlation", title_fontsize=8.5)
    fig.suptitle(suptitle, y=1.02, fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 0.82, 1])
    os.makedirs("figures", exist_ok=True)
    save_fig(fig, outfile)
    print(f"saved {outfile}")

def main():
    _panel_pair("Q1", "Q4",
                "PCA biplot of plant-based diet composition — PIR Q1 (poorest)",
                "PCA biplot of plant-based diet composition — PIR Q4 (richest, top-coded)",
                "figures/fig8_pca_quartiles.png",
                "PDI components in Q1 vs Q4: PCA of within-group plant-based diet composition")
    _panel_pair("T1", "T3",
                "PCA biplot of plant-based diet composition — PIR T1 (poorest third)",
                "PCA biplot of plant-based diet composition — PIR T3 (richest third)",
                "figures/fig9_pca_tertiles.png",
                "PDI components in T1 vs T3: PCA of within-group plant-based diet composition")

if __name__ == "__main__":
    main()
