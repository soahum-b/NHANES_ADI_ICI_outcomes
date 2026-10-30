#!/usr/bin/env python3
"""
35_dhq_bolte_reclassify.py — reclassify NHANES food codes into the Bolte plant-based
index categories, using DHQ_Dictionary_V2.csv's Plant_Based_Index_Bolte column as the
labelled reference set.

WHY THIS IS A MATCH AND NOT A JOIN: the dictionary keys on DHQ questionnaire item
numbers (FoodID 1.0-1176.5); NHANES foods are 8-digit FNDDS codes (11000000-95342000).
No key links them, so each NHANES food code is matched to its nearest DHQ item by
description and inherits that item's category. Matching is IDF-weighted cosine over
description tokens, so a shared rare word ("tangerine") counts far more than a shared
common one ("fresh").

TEXT ALONE IS NOT TRUSTED. The first digit of an FNDDS code is a major food group --
verified empirically against the data, not assumed:
    1 dairy · 2 meat/poultry/fish · 3 eggs · 4 legumes/nuts/seeds
    5 grains · 6 fruit · 7 vegetables · 8 fats/oils · 9 beverages/sugars
Each Bolte category's PLAUSIBLE leading digits are derived from where its own
high-confidence matches land, and any assignment landing outside its category's
plausible set is flagged for review rather than silently accepted. Review tiers:
    accept  high cosine AND plausible FNDDS group
    review  everything else, ranked worst-first for adjudication

Outputs:
  results/bolte/dhq_food_classification.csv   every NHANES code, category, tier, evidence
  results/bolte/dhq_review_queue.csv          only the rows needing adjudication
  results/bolte/dhq_category_summary.csv      per-category counts and diagnostics
Run from the project root, after 34b_nhanes_food_codes.R.

RESTRICTED INPUTS: the DHQ dictionary and the manual overrides live in data_restricted/,
which is never committed (see README). Both are REQUIRED: the 351 overrides are part of
the method, so their absence is an error, not a fallback -- an earlier version skipped
them silently when the file was missing and regenerated a different classification.
"""
import pandas as pd, numpy as np, re, sys, os
from collections import Counter

DICT = sys.argv[1] if len(sys.argv) > 1 else "data_restricted/raw/dhq/DHQ_Dictionary_V2.csv"
NH   = "data/interim/nhanes_codes_all.csv"
INTAKE = "data/interim/nhanes_code_intake.csv"   # weighted grams per food code
OVERRIDES = "data_restricted/interim/dhq_manual_overrides.csv"   # hand adjudications; win over the matcher
for _p, _what in ((DICT, "DHQ dictionary"), (OVERRIDES, "manual overrides")):
    if not os.path.exists(_p):
        sys.exit(f"{_what} not found at {_p}. It is a restricted input (not in the public repo); "
                 f"obtain it out of band or use the shipped data/interim/dhq_food_classification_public.csv "
                 f"with 36_dhq_scores_panels.R instead of re-running this script.")
HI_COS = 0.45          # head cosine above which the IDENTITY match is strong
FULL_FLOOR = 0.30      # full-description cosine the QUALIFIERS must also clear

STOP = {"and","the","with","from","not","nfs","reg","part","sandwiches","added","other",
        "type","made","prepared","eaten","only","all","any","etc","ns","as","to","or","in",
        "on","of","for","no","cooked","raw","fresh"}
def _stem(w):
    """Crude singularisation. Without it 'Egg omelet' never matches the dictionary item
    'Eggs', which collapsed the Eggs category to a single food code on the first pass."""
    for suf, rep in (("ies", "ie"), ("oes", "o"), ("ses", "s"), ("xes", "x"),
                     ("ches", "ch"), ("shes", "sh"), ("s", "")):
        if w.endswith(suf) and len(w) - len(suf) + len(rep) >= 3:
            return w[: len(w) - len(suf)] + rep
    return w

# DHQ names food GROUPS where FNDDS names specific foods, so the two never share a token
# for whole families: 'Watermelon, raw' vs 'Other melons', 'Ginger ale' vs 'Soda or pop'.
# Those were left UNMATCHED and read as dictionary gaps until the auditor checked -- the
# dictionary covers them, the matcher could not reach them. Map specific -> generic.
HYPERNYM = {
    "watermelon": "melon", "cantaloupe": "melon", "honeydew": "melon",
    "gingerale": "soda", "cola": "soda", "pop": "soda", "softdrink": "soda",
    "seltzer": "water", "sparkling": "water",
    "clementine": "orange", "tangerine": "orange", "mandarin": "orange", "tangelo": "orange",
    "scallion": "onion", "shallot": "onion", "zucchini": "squash", "courgette": "squash",
    "garbanzo": "bean", "chickpea": "bean", "lentil": "bean", "edamame": "soybean",
    "arugula": "lettuce", "romaine": "lettuce", "endive": "lettuce", "escarole": "lettuce",
    "pollock": "fish", "tilapia": "fish", "halibut": "fish", "mahi": "fish",
    "flounder": "fish", "snapper": "fish", "haddock": "fish", "swordfish": "fish",
}

QUALIFIER = {
    "unsweet": "unsweetened", "unsweeten": "unsweetened", "unsweetened": "unsweetened",
    "sugarfree": "artsweetened", "sugarless": "artsweetened", "diet": "artsweetened",
    "sweet": "sweetened", "sweeten": "sweetened", "sweetened": "sweetened",
    "presweetened": "sweetened", "sweetener": "sweetened",
    "decaf": "decaffeinated", "decaffeinated": "decaffeinated",
}
def _qual(w):
    return QUALIFIER.get(w, HYPERNYM.get(w, w))

def toks(s):
    s = re.sub(r"sugar[- ]free", "sugarfree", str(s).lower())
    s = re.sub(r"no cr\s*/?\s*sug", "unsweetened", s)      # DHQ shorthand for unsweetened
    s = re.sub(r"ginger ale", "gingerale", s)
    s = re.sub(r"soft drink", "softdrink", s)
    s = re.sub(r"pre[- ]sweetened", "presweetened", s)
    return [_qual(_stem(w)) for w in re.sub(r"[^a-z ]", " ", s).split()
            if len(w) > 2 and w not in STOP]

def head(s):
    """The food's IDENTITY, not how it was prepared.

    FNDDS descriptions are 'identity, qualifier, qualifier': 'Fennel bulb, cooked, made
    with oil'. Matching the whole string lets a preparation clause decide the category --
    measured on the first pass, 87% of codes assigned to a fat category got there because
    oil/margarine/butter appeared in a 'made with ...' clause, not because the food is a
    fat. That reproduces the very defect this reclassification is meant to fix, so only
    the leading clause is matched. DHQ items are split on '/' for the same reason
    ('Butter/ reg/ on bread,' -> 'Butter').
    """
    s = str(s).split(",")[0].split("/")[0]
    return re.sub(r"\b(cooked|raw|fresh|from|made|with|NS|NFS)\b.*", "", s, flags=re.I)

dq = pd.read_csv(DICT).dropna(subset=["FoodDescription", "Plant_Based_Index_Bolte"])
nh = pd.read_csv(NH, dtype={"code": str}).dropna(subset=["DESCRIPTION"])
nh = nh.drop_duplicates("code").reset_index(drop=True)

# ---- IDF over the NHANES corpus (the side we are classifying) ----------------
docs_nh = [toks(head(x)) for x in nh.DESCRIPTION]
df_counts = Counter(t for d in docs_nh for t in set(d))
N = len(docs_nh)
idf = {t: np.log(1 + N / (1 + c)) for t, c in df_counts.items()}
def vec(ts):
    v = {}
    for t in ts:
        v[t] = v.get(t, 0.0) + idf.get(t, np.log(1 + N))
    n = np.sqrt(sum(x * x for x in v.values())) or 1.0
    return {t: x / n for t, x in v.items()}

vq_head = [vec(toks(head(x))) for x in dq.FoodDescription]
vq_full = [vec(toks(x)) for x in dq.FoodDescription]
cats = dq.Plant_Based_Index_Bolte.to_numpy()
dsc  = dq.FoodDescription.to_numpy()
docs_full = [toks(x) for x in nh.DESCRIPTION]

best_i, best_s, best_f = np.zeros(N, int), np.zeros(N), np.zeros(N)
for i, (dh, df) in enumerate(zip(docs_nh, docs_full)):
    vh, vf = vec(dh), vec(df)
    sh = np.array([sum(vh.get(t, 0.0) * w for t, w in q.items()) for q in vq_head])
    top = sh.max()
    cand = np.flatnonzero(sh >= top - 1e-9)                 # all items tying on identity
    sf = np.array([sum(vf.get(t, 0.0) * w for t, w in vq_full[j].items()) for j in cand])
    j = int(cand[int(np.argmax(sf))])                       # tie-break on the qualifiers
    best_i[i], best_s[i], best_f[i] = j, top, sf.max()

unmatched = best_s <= 0.0
out = nh.assign(category=np.where(unmatched, "UNMATCHED", cats[best_i]),
                matched_dhq_item=np.where(unmatched, "", dsc[best_i]),
                cosine=np.round(best_s, 3),
                cosine_full=np.round(best_f, 3),
                fndds_group=nh.code.str[0])

# ---- plausible FNDDS groups per category, LEARNED from strong matches --------
plaus = {
    "Dairy": set("1"), "Meat": set("2"), "Fish_or_Seafood": set("2"), "Eggs": set("3"),
    "Nuts_and_Legumes": set("4"), "Whole_Grains": set("5"), "Refined_Grains": set("5"),
    "Sweets_and_desserts": set("59"),            # baked goods sit in 5, candy/sugars in 9
    "Fruits": set("6"), "Fruit_Juices": set("69"),
    "Vegetables": set("7"), "Potatoes": set("7"),
    "Vegetable_Oils": set("8"), "MARGARINE": set("8"), "Animal_Fat": set("189"),
    "Tea_and_Coffee": set("9"), "Sugar_Sweetened_Beverages": set("9"),
    "Alcoholic_Beverages": set("9"), "WATER": set("9"), "ART_SWEET": set("9"),
    "PROTEIN": set("9"),
    # mixed dishes: FPED cannot build this group at all, and in FNDDS its members are
    # spread across dairy, meat, grain-based dishes and dressings
    "Miscellaneous_animal_based_foods": set("1258"),
    "UNMATCHED": set(),
}
out["fndds_plausible"] = [g in plaus.get(c, set()) for c, g in zip(out.category, out.fndds_group)]
out["tier"] = np.where((out.cosine >= HI_COS) & (out.cosine_full >= FULL_FLOOR)
                       & out.fndds_plausible, "accept", "review")
out["review_reason"] = np.where(
    out.tier == "accept", "",
    np.where(out.cosine < HI_COS,
             np.where(out.fndds_plausible, "weak text match", "weak match + implausible food group"),
             np.where(~out.fndds_plausible, "strong text match but implausible FNDDS group",
                      "identity matches but qualifiers disagree")))

# ---- whole vs refined is decided by the QUALIFIER, not the identity ---------
# Head-anchoring (needed so "cooked with oil" cannot make a vegetable a fat) strips exactly
# the word that separates Whole_Grains from Refined_Grains -- and that is the HEALTHFUL vs
# LESS-HEALTHFUL plant boundary, the most load-bearing distinction in the index. Measured
# before this rule: 157 codes said "whole grain"/"whole wheat"/"brown rice" in their
# description and NONE were in Whole_Grains. The grain block therefore reads the full
# description, overriding the text match.
WHOLE_RX = re.compile(r"whole.{0,3}(grain|wheat)|brown rice|wild rice|oatmeal|\boat\b|"
                      r"\bbulgur\b|\bquinoa\b|\bbarley\b|\bbuckwheat\b|\bmillet\b|"
                      r"shredded wheat|bran\b", re.I)
GRAIN_CATS = {"Whole_Grains", "Refined_Grains"}
is_grainish = out.fndds_group.eq("5") | out.category.isin(GRAIN_CATS)
says_whole = out.DESCRIPTION.str.contains(WHOLE_RX, na=False)
# only move codes the matcher already read as grain, or FNDDS grain codes it mis-binned
promote = is_grainish & says_whole & ~out.category.isin(["UNMATCHED"]) & out.category.ne("Whole_Grains")
demote = out.category.eq("Whole_Grains") & ~says_whole
out.loc[promote, "category"] = "Whole_Grains"
out.loc[demote, "category"] = "Refined_Grains"
print(f"grain rule: {int(promote.sum())} code(s) -> Whole_Grains, "
      f"{int(demote.sum())} -> Refined_Grains (qualifier overrides the text match)")

# ---- whole fruit vs juice: again decided by the qualifier -------------------
JUICE_RX = re.compile(r"\bjuice\b|\bnectar\b|\blemonade\b|\blimeade\b|\bcider\b", re.I)
fruitish = out.category.isin(["Fruits", "Fruit_Juices"])
is_juice = out.DESCRIPTION.str.contains(JUICE_RX, na=False)
to_juice = fruitish & is_juice & out.category.eq("Fruits")
to_fruit = out.category.eq("Fruit_Juices") & ~is_juice
out.loc[to_juice, "category"] = "Fruit_Juices"
out.loc[to_fruit, "category"] = "Fruits"
print(f"juice rule: {int(to_juice.sum())} code(s) -> Fruit_Juices, {int(to_fruit.sum())} -> Fruits")

# ---- manual adjudications override the matcher entirely --------------------
# Existence was checked at the top; this block is unconditional on purpose.
if True:
    ov = pd.read_csv(OVERRIDES, dtype={"code": str})
    unknown = set(ov.code) - set(out.code)
    if unknown:
        sys.exit(f"override codes absent from the NHANES code list: {sorted(unknown)}")
    bad = set(ov.category) - set(plaus) - {"UNMATCHED"}
    if bad:
        sys.exit(f"override categories not in the dictionary scheme: {sorted(bad)}")
    m = out.code.isin(set(ov.code))
    out.loc[m, "category"] = out.loc[m, "code"].map(dict(zip(ov.code, ov.category)))
    out.loc[m, "tier"] = "manual"
    out.loc[m, "review_reason"] = "manually adjudicated"
    out["adjudication_note"] = out.code.map(
        dict(zip(ov.code, ov.get("note", pd.Series(index=ov.index, dtype=str)).fillna("")))).fillna("")
    print(f"applied {m.sum()} manual adjudication(s) from {OVERRIDES}")

os.makedirs("results/bolte", exist_ok=True)
intake = pd.read_csv(INTAKE, dtype={"code": str})
out = out.merge(intake, on="code", how="left")
out["grams_wt"] = out.grams_wt.fillna(0.0)
out["pct_of_total_intake"] = (100 * out.grams_wt / out.grams_wt.sum()).round(4)

cols = ["code","DESCRIPTION","category","tier","cosine","cosine_full","fndds_group",
        "matched_dhq_item","review_reason","pct_of_total_intake"]
if "adjudication_note" in out.columns:
    cols = cols + ["adjudication_note"]
(out[cols].sort_values(["tier","pct_of_total_intake"], ascending=[True, False])
          .to_csv("results/bolte/dhq_food_classification.csv", index=False))
# Committed, DHQ-free cut: only what 36_ needs (code, DESCRIPTION, category, tier). The
# matcher diagnostics, including matched_dhq_item which quotes DHQ item text verbatim,
# are deliberately absent. Written here so the public file is a script output, not a
# hand-made artefact that could drift from the full classification.
PUBLIC = "data/interim/dhq_food_classification_public.csv"
(out.sort_values(["tier","pct_of_total_intake"], ascending=[True, False])
    [["code","DESCRIPTION","category","tier"]]
    .to_csv(PUBLIC, index=False))
print(f"wrote {PUBLIC} (public four-column cut, no DHQ text)")
rq = (out[out.tier == "review"].sort_values("grams_wt", ascending=False)
         .assign(intake_rank=lambda d: range(1, len(d) + 1)))
rq[["intake_rank"] + cols].to_csv("results/bolte/dhq_review_queue.csv", index=False)

cum = rq.grams_wt.cumsum() / rq.grams_wt.sum()
cover = {p: int((cum < p).sum()) + 1 for p in (0.5, 0.8)}

summ = (out.groupby("category")
          .agg(n_codes=("code","size"),
               n_accept=("tier", lambda s: (s=="accept").sum()),
               median_cosine=("cosine","median"))
          .assign(pct_accept=lambda d: (100*d.n_accept/d.n_codes).round(0))
          .sort_values("n_codes", ascending=False).reset_index())
summ.to_csv("results/bolte/dhq_category_summary.csv", index=False)

top = out.sort_values("pct_of_total_intake", ascending=False).head(50)
top[cols].to_csv("results/bolte/dhq_top_intake_audit.csv", index=False)
print("TOP 12 FOOD CODES BY INTAKE -- always inspect these, whatever tier they land in:")
for _, r in top.head(12).iterrows():
    print(f"  {r.pct_of_total_intake:6.2f}%  {str(r.DESCRIPTION)[:44]:44s} -> "
          f"{str(r.category)[:26]:26s} [{r.tier}]")
print(f"  (top 50 = {top.pct_of_total_intake.sum():.1f}% of all grams; "
      f"full list in results/bolte/dhq_top_intake_audit.csv)\n")
print(f"NHANES food codes classified: {len(out):,}")
# ---- reader-facing export: every NHANES food code and its category --------------
# dhq_food_classification.csv carries the matcher's working columns (cosine, fndds_group,
# matched item, review reason). This is the same rows stripped to what a reader of the
# classification needs, plus an explicit statement of whether each category feeds the
# plant-based index: 5 of the 22 categories do NOT (WATER, PROTEIN, ART_SWEET, MARGARINE,
# Alcoholic_Beverages), and UNMATCHED means no category was assigned at all.
HEALTHY = ["Whole_Grains","Fruits","Vegetables","Nuts_and_Legumes","Vegetable_Oils","Tea_and_Coffee"]
LESSHEALTHY = ["Fruit_Juices","Refined_Grains","Potatoes","Sugar_Sweetened_Beverages","Sweets_and_desserts"]
ANIMAL = ["Animal_Fat","Dairy","Eggs","Fish_or_Seafood","Meat","Miscellaneous_animal_based_foods"]
BLOCK = {**{c: "healthful plant" for c in HEALTHY},
         **{c: "less-healthful plant" for c in LESSHEALTHY},
         **{c: "animal" for c in ANIMAL}}

_exp = out.copy()
_exp["index_block"] = _exp["category"].map(BLOCK)
_exp.loc[_exp["index_block"].isna() & _exp["category"].eq("UNMATCHED"), "index_block"] = "unclassified"
_exp["index_block"] = _exp["index_block"].fillna("not in index")
_exp["settled"] = _exp["tier"].isin(["accept", "manual"]).map({True: "yes", False: "needs review"})
_exp = (_exp[["code","DESCRIPTION","category","index_block","settled","pct_of_total_intake"]]
        .rename(columns={"code":"NHANES food code","DESCRIPTION":"Food description",
                         "category":"Assigned category","index_block":"Index block",
                         "settled":"Assignment settled","pct_of_total_intake":"% of total grams consumed"})
        .sort_values(["Index block","Assigned category","% of total grams consumed"],
                     ascending=[True, True, False]))
assert len(_exp) == out.code.nunique(), "export lost rows"
_exp.to_csv("results/bolte/dhq_food_categories_for_analyses.csv", index=False)
print(f"wrote results/bolte/dhq_food_categories_for_analyses.csv - all {len(_exp):,} NHANES food codes, "
      f"{_exp['Assigned category'].nunique()} categories, "
      f"{_exp['Index block'].isin(['healthful plant','less-healthful plant','animal']).sum():,} "
      f"in one of the three index blocks")

print(f"  intake covered by accept tier: {out.loc[out.tier=='accept','pct_of_total_intake'].sum():.1f}%"
      f" | review tier: {out.loc[out.tier=='review','pct_of_total_intake'].sum():.1f}%")
print(f"  review queue ranked by consumption: top {cover[0.5]} codes = 50% of review-tier intake,"
      f" top {cover[0.8]} = 80%")
print(f"  accept : {(out.tier=='accept').sum():,}  ({100*(out.tier=='accept').mean():.0f}%)")
print(f"  review : {(out.tier=='review').sum():,}  ({100*(out.tier=='review').mean():.0f}%)")
print("\nreview reasons:")
print(out[out.tier=="review"].review_reason.value_counts().to_string())
print("\nper-category:")
print(summ.to_string(index=False))
