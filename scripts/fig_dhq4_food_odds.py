#!/usr/bin/env python3
"""Same figure as fig_bolte4_food_odds.py, built on the DHQ whole-item classification
(36_dhq_scores_panels.R -> results/dhq/) instead of FPED components. Imports that
script's panel code so the two versions cannot drift apart; only the inputs differ."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fig_bolte4_food_odds as M
M.RESULTS_DIR = "results/dhq"
M.OUT_PNG = "figures/dhq/bolte17_food_odds.png"
M.GROUPING = "Bolte et al. 2023 grouping, all 17 groups (whole-item, DHQ dictionary)"
if __name__ == "__main__":
    os.makedirs("figures/dhq", exist_ok=True)
    M.main()
