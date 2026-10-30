#!/usr/bin/env python3
"""
fig_dhq_split_panels.py — panels A-F of the 17-group (DHQ whole-item) SES figure as SIX
SEPARATE standalone figures, one per score x SES measure.

  bolte17_a_hpdi_pir.png        A   hPDI vs income-to-poverty ratio
  bolte17_b_hpdi_education.png  B   hPDI vs education
  bolte17_c_hpdi_income.png     C   hPDI vs annual family income
  bolte17_d_updi_pir.png        D   uPDI vs income-to-poverty ratio
  bolte17_e_updi_education.png  E   uPDI vs education
  bolte17_f_updi_income.png     F   uPDI vs annual family income

Reuses _panel() from fig_bolte3_panels.py, so these and the combined figure are drawn by
identical code and cannot drift apart; only the layout differs.

THE Y-RANGE IS SHARED PER INDEX: the three hPDI panels share one range and the three
uPDI panels another, rather than autoscaling individually. Standalone figures are the case where per-panel
autoscaling misleads most: a reader comparing two separate images has no shared axis to
anchor on, so a gradient can look twice as steep purely from scaling.
Each panel carries its own axis labels and a descriptive title, plus the composite's
letter so a figure can be traced back to the panel it came from.

Reads results/dhq/{pd_main,trend_tests}.csv (36_dhq_scores_panels.R).
Run from the project root.
"""
import os, sys
import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig_bolte3_panels as B
from figsave import save_fig

B.RESULTS_DIR = "results/dhq"
B.SCORES = ["hPDI", "uPDI"]
B.COL = {"hPDI": "#7B4FA8", "uPDI": "#C0603A"}
B.TITLE = {"hPDI": "Healthful Plant-Based (hPDI)", "uPDI": "Unhealthful Plant-Based (uPDI)"}

OUTDIR = "figures/dhq/panels"
NAMES = {("hPDI", "pir_bin"):    ("A", "bolte17_a_hpdi_pir"),
         ("hPDI", "educ_f"):     ("B", "bolte17_b_hpdi_education"),
         ("hPDI", "income_lab"): ("C", "bolte17_c_hpdi_income"),
         ("uPDI", "pir_bin"):    ("D", "bolte17_d_updi_pir"),
         ("uPDI", "educ_f"):     ("E", "bolte17_e_updi_education"),
         ("uPDI", "income_lab"): ("F", "bolte17_f_updi_income")}


def main():
    pd_main = pd.read_csv(f"{B.RESULTS_DIR}/pd_main.csv")
    tt = pd.read_csv(f"{B.RESULTS_DIR}/trend_tests.csv")

    def ann(score, by):
        r = tt[(tt.score == score) & (tt["by"] == by)]
        if r.empty:
            return ""
        r = r.iloc[0]
        kind = "Trend" if r.test == "linear trend" else "Assoc."
        return (f"{kind} FDR p {r.p_fdr_disp}" if str(r.p_fdr_disp).startswith("<")
                else f"{kind} FDR p = {r.p_fdr_disp}")

    specs = {"pir_bin":    ("pir_bin", B.ORDER_PIR, "PIR Ratio", None, 0),
             "educ_f":     ("educ_f", B.ORDER_EDU, "Education", B.LAB_EDU, 0),
             "income_lab": ("income_lab", B.ORDER_INC, "Annual Family Income, USD", None, 55)}

    # ONE range PER INDEX, not one across all six. This follows the standing instruction
    # for these panels ("a constant start and stop for the hPDI y axes, and similarly for
    # uPDI"): the three SES panels of an index share an axis so their gradients are
    # directly comparable, while hPDI and uPDI keep their own. Forcing a single range over
    # both indices here would be worse than on the combined figure, because hPDI sits
    # ~53-58 and uPDI ~57-62: a common 52-63 axis leaves each standalone panel using under
    # half its height, and standalone figures have no neighbouring panel to justify it.
    YLIM = {}
    for sc in B.SCORES:
        p_ = pd_main[pd_main.score == sc]
        pad = 0.10 * (p_.ci_u.max() - p_.ci_l.min())
        YLIM[sc] = (p_.ci_l.min() - pad, p_.ci_u.max() + pad)

    os.makedirs(OUTDIR, exist_ok=True)
    for (score, by), (letter, stem) in NAMES.items():
        wide = by == "income_lab"          # 12 rotated tick labels need the extra width
        fig, ax = plt.subplots(figsize=(5.6 if wide else 4.6, 4.1))
        B._panel(ax, pd_main, ann, score, specs[by], YLIM[score], letter,
                 show_xlab=True, show_ylab=True, ylab=f"{score} Score")
        ax.set_title(B.TITLE[score], fontsize=10, fontweight="bold",
                     color=B.COL[score], pad=10)
        fig.tight_layout()
        save_fig(fig, f"{OUTDIR}/{stem}.png")
        plt.close(fig)
        print(f"saved {OUTDIR}/{stem}.png  ({letter}: {score} vs {specs[by][2]})")


if __name__ == "__main__":
    main()
