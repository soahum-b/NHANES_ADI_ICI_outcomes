#!/usr/bin/env python3
"""
37_dhq_pca_loadings_table.py — PC1/PC2 loadings table for the DHQ whole-item PCA,
counterpart to the primary Table 5 and the Bolte table from 33_.

build_pca_loadings_table() below is INLINED VERBATIM from 20_refresh_whitepaper_tables.py
(the project's shared table builder) so this bundle stands alone. It is copied rather than
reimplemented specifically so the SIGN ALIGNMENT matches the primary and 16-group tables:
a component's sign is arbitrary, so raw per-group loadings cannot be compared column by
column until each group's axes are anchored. Any change here must be mirrored there.

Emits both the |rho| >= 0.30 display cut and the complete 17-row table, because the cut
is a presentation convention and every group is in the PCA regardless.
Reads results/dhq/pca_loadings.csv. Run from the project root.
"""
import pandas as pd, numpy as np


def build_pca_loadings_table(src="results/pca_loadings.csv",
                             out="results/pca_loadings_table.csv",
                             groups=("Q1", "Q4"),
                             anchors=("animal_fat", "fruits")):
    """PC1/PC2 loadings per income group, SIGN-ALIGNED so groups are comparable.

    A principal component is sign-indeterminate: the whole axis can flip with no
    change in meaning, so raw per-group loadings cannot be compared column-by-column
    (Q4's PC2 came out sign-flipped relative to Q1's). We therefore orient every
    group's PC1 so animal fat is positive and PC2 so fruit is positive, then only
    genuine differences in magnitude remain. Keeps foods reaching |rho| >= 0.30.
    """
    import numpy as np
    L = pd.read_csv(src)
    al = []
    for g, s in L.groupby("group"):
        s = s.set_index("food").copy()
        s["rho_pc1"] *= np.sign(s.loc[anchors[0], "rho_pc1"])
        s["rho_pc2"] *= np.sign(s.loc[anchors[1], "rho_pc2"])
        al.append(s.assign(group=g).reset_index())
    A = pd.concat(al)
    strong = A.groupby("food")[["rho_pc1", "rho_pc2"]].apply(lambda d: d.abs().max().max()) >= 0.30
    foods = strong[strong].index.tolist()
    piv = A[A.food.isin(foods)].pivot_table(index=["food", "label"], columns="group",
                                            values=["rho_pc1", "rho_pc2"])
    piv = piv.loc[:, [(pc, g) for pc in ("rho_pc1", "rho_pc2") for g in groups]]
    piv.columns = [f"\u03c1({pc}) {g}" for pc in ("PC1", "PC2") for g in groups]
    piv = (piv.reset_index().drop(columns="food").rename(columns={"label": "Food group"})
              .sort_values(f"\u03c1(PC1) {groups[0]}", key=lambda c: c.abs(), ascending=False).round(2))
    piv.to_csv(out, index=False)
    return piv

SRC = "results/dhq/pca_loadings.csv"
OUT = "results/dhq/bolte17_pca_loadings_table.csv"
FULL = OUT.replace(".csv", "_full.csv")   # derived, so it cannot diverge from OUT

ANCHORS = ("Animal_Fat", "Fruits")
tbl = build_pca_loadings_table(src=SRC, out=OUT, groups=("Q1", "Q4"), anchors=ANCHORS)

raw = pd.read_csv(SRC)
al = []
for g, s in raw.groupby("group"):
    s = s.copy()
    si = s.set_index("food")
    s["rho_pc1"] *= np.sign(si.loc[ANCHORS[0], "rho_pc1"])
    s["rho_pc2"] *= np.sign(si.loc[ANCHORS[1], "rho_pc2"])
    al.append(s)
A = pd.concat(al)
piv = (A.pivot_table(index="label", columns="group", values=["rho_pc1", "rho_pc2"])
         .round(2))
piv.columns = [f"\u03c1({'PC1' if a=='rho_pc1' else 'PC2'}) {b}" for a, b in piv.columns]
piv["Max |\u03c1|"] = piv.abs().max(axis=1).round(2)
piv["In main table"] = np.where(piv["Max |\u03c1|"] >= 0.30, "yes", "no")
piv.sort_values("Max |\u03c1|", ascending=False).to_csv(FULL)

print(f"wrote {OUT} — {len(tbl)} of {raw.food.nunique()} groups reach |rho| >= 0.30")
print(f"wrote {FULL} — all {len(piv)} groups")
print(piv.sort_values("Max |\u03c1|", ascending=False).to_string())
