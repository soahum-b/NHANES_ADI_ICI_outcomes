#!/usr/bin/env python3
"""
fig_bolte4_food_odds.py — which foods carry the SES gradient, on the BOLTE grouping.

Odds of being a high consumer (above the survey-weighted 67th percentile of that food)
per +1 unit of income-to-poverty ratio, adjusted for age, sex and race/ethnicity.
PIR enters as a CONTINUOUS predictor -- this is a per-unit slope, not a group contrast.

16 Bolte groups. ALCOHOL IS ABSENT because it is not one of Bolte's categories, so the
BH-FDR family here is 16 tests rather than the 19 of our own food-odds figure.
Reads results/bolte/food_odds.csv (29_bolte_panels.R).
Output: figures/bolte/bolte4_food_odds.{png,pdf,svg}
"""
import pandas as pd, numpy as np, matplotlib.pyplot as plt, os, sys
import matplotlib.ticker as ticker
from matplotlib.transforms import blended_transform_factory
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figsave import save_fig

# set by the DHQ wrappers; defaults keep the FPED-component figures unchanged
RESULTS_DIR = 'results/bolte'
OUT_PNG = 'figures/bolte/bolte4_food_odds.png'
GROUPING = 'Bolte et al. 2023 grouping (16 groups, no alcohol)'


BLUE, ORANGE, GREY = "#1f5fa8", "#c0603a", "#9aa0a6"


def main():
    d = pd.read_csv(f"{RESULTS_DIR}/food_odds.csv").sort_values("OR").reset_index(drop=True)
    y = np.arange(len(d))
    cols = [GREY if (lo < 1 < hi) else (BLUE if o > 1 else ORANGE)
            for o, lo, hi in zip(d.OR, d.OR_lo, d.OR_hi)]

    fig, ax = plt.subplots(figsize=(9.4, 6.6))
    ax.axvline(1, color="#444", lw=1, ls="--", zorder=1)
    for i, r in d.iterrows():
        ax.plot([r.OR_lo, r.OR_hi], [i, i], color=cols[i], lw=2.0, zorder=2)
        ax.scatter(r.OR, i, color=cols[i], s=34, zorder=3)
    ax.set_yticks(y); ax.set_yticklabels(d.label, fontsize=9)
    ax.set_xscale("log")
    ticks = [0.6, 0.8, 1.0, 1.2, 1.5]
    ax.set_xticks(ticks); ax.set_xticklabels([str(t) for t in ticks])
    # a log axis also draws MINOR ticks with their own labels (7x10^-1, 9x10^-1),
    # which collide with the explicit major labels above. Silence them.
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_xlabel("Odds ratio per +1 PIR unit  (high consumer vs not)",
                  fontsize=10, fontweight="bold")
    ax.set_ylim(-0.8, len(d) - 0.2)

    trans = blended_transform_factory(ax.transAxes, ax.transData)
    ax.text(1.04, len(d) - 0.4, "OR (95% CI)", transform=trans, ha="left", va="bottom",
            fontsize=8.6, fontweight="bold", color="#333")
    for i, r in d.iterrows():
        ax.text(1.04, i, f"{r.OR:.2f} ({r.OR_lo:.2f}\u2013{r.OR_hi:.2f})", transform=trans,
                ha="left", va="center", fontsize=8.2, color=cols[i])

    ax.text(0.0, -0.105, "\u2190 Less likely at higher SES", transform=ax.transAxes,
            ha="left", fontsize=9, color=ORANGE)
    ax.text(0.62, -0.105, "More likely at higher SES \u2192", transform=ax.transAxes,
            ha="left", fontsize=9, color=BLUE)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)

    fig.suptitle("Which foods carry the SES gradient? \u2014 " + GROUPING,
                 y=0.99, fontsize=10.5, fontweight="bold")
    fig.tight_layout()
    fig.subplots_adjust(right=0.70)
    os.makedirs(os.path.dirname(OUT_PNG), exist_ok=True)
    save_fig(fig, OUT_PNG)
    print("saved", OUT_PNG)


if __name__ == "__main__":
    main()
