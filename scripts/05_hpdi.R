#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 05_hpdi.R  —  Overall (oPDI/PDI), Healthful (hPDI) and Unhealthful (uPDI)
#               Plant-Based Diet Indices
# Satija et al. 2017 (J Am Coll Cardiol 70:411), NHANES 2011-2018 day-1 recalls.
#
# Method: 18 food groups collapsed into 3 categories —
#   HEALTHY PLANT   : whole grains, fruits, vegetables (non-starchy), nuts,
#                     legumes, vegetable oils, tea/coffee
#   LESS-HEALTHY PLANT: fruit juices, refined grains, potatoes/starchy veg,
#                     sugar-sweetened beverages, sweets/desserts
#   ANIMAL          : animal fat, dairy, eggs, fish/seafood, meat (red+processed),
#                     POULTRY  <-- see note
# NOTE ON THE 18th ANIMAL GROUP: this header previously read "misc animal" while the
#   code below has always used poultry (DR1I_PF_POULT). The CODE is what ran; this
#   comment was stale. Satija's published animal category ends in a "miscellaneous
#   animal foods" group; we use poultry as the 18th group instead, which narrows it.
#   The JACC 2017 full text is paywalled from this environment, so the original's exact
#   group table could not be re-verified here -- treat poultry-vs-misc-animal as a
#   KNOWN DEVIATION of our implementation, not as a match to the published index.
# NO ENERGY ADJUSTMENT: quintiles are taken on raw FPED servings/day, and the SES models
#   adjust for age, sex and race/ethnicity only -- not total energy. So a person who eats
#   more of everything ranks higher in every group. Documented, not corrected.
# Each group -> population quintiles (1..5) within the analytic sample, weighted.
#   oPDI: ALL PLANT groups POSITIVE (healthy + less-healthy), animal REVERSE.
#         Satija's "overall" PDI (often written PDI) — makes no quality
#         distinction among plant foods; measures plant-vs-animal only.
#   hPDI: healthy-plant scored POSITIVE (q->score), less-healthy-plant & animal REVERSE.
#   uPDI: less-healthy-plant POSITIVE, healthy-plant & animal REVERSE.
# Per-person index = sum over 18 groups. Higher hPDI = more healthful plant-based.
# The three indices share identical group definitions and quintile cut-points;
# they differ ONLY in which groups are reverse-scored, so all are on the same
# 18-90 scale and are directly comparable.
#
# Inputs: food-level FPED (fped_dr1iff_*.sas7bdat) for plant/animal group servings,
#   the SAME FPED per-food file for added-sugar drinks and food-code level detail.
#   NOTE: this script does NOT read CDC's DR1IFF_*.xpt. Everything food-level comes from
#   the USDA FPED fped_dr1iff_*.sas7bdat, which already carries the food code, the gram
#   weight and the pattern components. The shared DR1IFF name stem is coincidence.
# Output: data/interim/hpdi_scores.csv (SEQN, cycle, oPDI, hPDI, uPDI + 18 group servings)
# ------------------------------------------------------------------------------
suppressMessages({library(haven); library(dplyr); library(tidyr); library(survey)})

cyc <- tibble::tribble(
  ~suf, ~yy,   ~cycle,
  "G",  "1112","2011-2012",
  "H",  "1314","2013-2014",
  "I",  "1516","2015-2016",
  "J",  "1718","2017-2018")

# ---- Build per-food 18-group servings from food-level FPED ------------------
# FPED components are in cup-eq (fruit, veg, dairy), oz-eq (grains, protein),
# tsp-eq (added sugars), grams (oils, solid fats). We keep each group in its
# native FPED unit; quantiles are within-group so units need only be consistent.
read_food <- function(suf, yy, cycle) {
  # Food-level FPED is self-contained: food code (DR1IFDCD), grams (DR1IGRMS),
  # and all food-pattern components per food record. No DR1IFF join needed.
  read_sas(sprintf("data/raw/fped/fped_dr1iff_%s.sas7bdat", yy)) %>% mutate(cycle=cycle)
}
foods <- purrr::pmap_dfr(cyc, read_food)
cat("Food records loaded:", nrow(foods), " people:", length(unique(foods$SEQN)), "\n")

# Sugar-sweetened beverages: FNDDS food codes beginning 924 (soft drinks),
# 926 (fruit drinks/-ades), plus sweetened coffees/teas — flagged via food code.
# Tea/coffee (healthy plant group): codes 92100000-92309000 unsweetened.
fc <- foods$DR1IFDCD
foods <- foods %>% mutate(
  fdc4 = as.integer(substr(sprintf("%.0f", DR1IFDCD),1,4)),  # %.0f avoids sci-notation (e.g. 1.1e+07) -> NA
  is_ssb   = fdc4 %in% c(9240,9241,9242,9250,9260,9261,9262,9264,9266) ,
  is_coffeetea = fdc4 %in% c(9210,9220,9230))

# ---- Collapse to 18 groups (per food record) --------------------------------
g <- foods %>% transmute(
  SEQN, cycle,
  # healthy plant
  whole_grains = DR1I_G_WHOLE,
  fruits       = DR1I_F_TOTAL - DR1I_F_JUICE,
  veg_healthy  = DR1I_V_DRKGR + DR1I_V_REDOR_TOTAL + DR1I_V_OTHER,   # non-starchy veg
  nuts         = DR1I_PF_NUTSDS,
  legumes      = DR1I_V_LEGUMES + DR1I_PF_LEGUMES + DR1I_PF_SOY,
  veg_oils     = DR1I_OILS,
  tea_coffee   = ifelse(is_coffeetea, DR1IGRMS/240, 0),             # cup-eq approx
  # less-healthy plant
  fruit_juice  = DR1I_F_JUICE,
  refined_grains = DR1I_G_REFINED,
  potatoes     = DR1I_V_STARCHY_TOTAL,
  ssb          = ifelse(is_ssb, DR1IGRMS/240, 0),                   # cup-eq approx
  sweets       = DR1I_ADD_SUGARS,                                   # tsp-eq added sugar as sweets proxy
  # animal
  animal_fat   = DR1I_SOLID_FATS,
  dairy        = DR1I_D_TOTAL,
  eggs         = DR1I_PF_EGGS,
  fish         = DR1I_PF_SEAFD_HI + DR1I_PF_SEAFD_LOW,
  # Satija separates red/processed meat from poultry; keep them as distinct animal groups
  meat         = DR1I_PF_MEAT + DR1I_PF_CUREDMEAT + DR1I_PF_ORGAN,   # red + processed + organ
  poultry      = DR1I_PF_POULT)                                       # 18th group: poultry (was placeholder)

groups <- c("whole_grains","fruits","veg_healthy","nuts","legumes","veg_oils","tea_coffee",
            "fruit_juice","refined_grains","potatoes","ssb","sweets",
            "animal_fat","dairy","eggs","fish","meat","poultry")
healthy_plant <- c("whole_grains","fruits","veg_healthy","nuts","legumes","veg_oils","tea_coffee")
less_healthy  <- c("fruit_juice","refined_grains","potatoes","ssb","sweets")
animal_grp    <- c("animal_fat","dairy","eggs","fish","meat","poultry")

# ---- Aggregate to per-person daily servings ---------------------------------
per <- g %>% group_by(SEQN, cycle) %>% summarise(across(all_of(groups), ~sum(.x, na.rm=TRUE)), .groups="drop")
cat("Per-person rows:", nrow(per), "\n")

# ---- Restrict to analytic sample (align with main dataset) ------------------
ana <- readRDS("data/interim/nhanes_analytic.rds")   # adults 20+, reliable recall
per <- per %>% semi_join(ana, by="SEQN")
per <- per %>% left_join(ana %>% select(SEQN, SDMVSTRA, SDMVPSU, WTDRD1_POOL), by="SEQN")
cat("Per-person rows in analytic sample:", nrow(per), "\n")

# ---- Weighted quintile scoring per group ------------------------------------
# Satija: zero-consumers get the lowest score. Non-zero split into quintiles by
# weighted cutpoints; each group contributes 1..5. Reverse for the negatively
# scored categories. Using survey-weighted quantiles for representativeness.
des <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTDRD1_POOL, nest=TRUE, data=per)
options(survey.lonely.psu="adjust")

# inline weighted quantile (type-7-ish on the weighted ECDF)
wtd_quantile <- function(x, w, probs) {
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w)/sum(w)
  sapply(probs, function(p) x[which(cw >= p)[1]])
}
score_group <- function(x, w, reverse=FALSE) {
  # weighted quintile breaks on positive consumers; zeros -> quintile 1
  pos <- x>0
  s <- rep(1L, length(x))
  if (sum(pos)>2) {
    qs <- unique(wtd_quantile(x[pos], w[pos], probs=seq(.2,.8,.2)))
    br <- unique(c(-Inf, qs, Inf))
    s[pos] <- as.integer(cut(x[pos], breaks=br, labels=FALSE, include.lowest=TRUE))
    # map so lowest positive bin starts above the zero-consumers where possible
    s[pos] <- s[pos] + (5L - (length(br)-1L))
  }
  s <- pmin(5L, pmax(1L, s))
  if (reverse) s <- 6L - s
  s
}
w <- per$WTDRD1_POOL
score_mat <- sapply(groups, function(gn) {
  rev <- !(gn %in% healthy_plant)   # for hPDI: healthy-plant positive, others reverse
  score_group(per[[gn]], w, reverse=rev)
})
# hPDI = sum with healthy-plant positive, less-healthy & animal reversed (done above)
per$hPDI <- rowSums(score_mat)

# uPDI: less-healthy-plant positive; healthy-plant & animal reversed
score_mat_u <- sapply(groups, function(gn) {
  rev <- !(gn %in% less_healthy)
  score_group(per[[gn]], w, reverse=rev)
})
per$uPDI <- rowSums(score_mat_u)

# oPDI (overall PDI): ALL plant groups positive (healthy AND less-healthy),
# animal groups reversed. No quality distinction among plant foods.
plant_all <- c(healthy_plant, less_healthy)
score_mat_o <- sapply(groups, function(gn) {
  rev <- !(gn %in% plant_all)          # only animal groups reversed
  score_group(per[[gn]], w, reverse=rev)
})
per$oPDI <- rowSums(score_mat_o)

# --- structural check: the 3 indices must share cut-points, so for any person
# the animal-group contribution is identical in oPDI and hPDI, and each index
# must lie on the same 18-90 scale.
stopifnot(all(per$oPDI >= 18 & per$oPDI <= 90),
          all(per$hPDI >= 18 & per$hPDI <= 90),
          all(per$uPDI >= 18 & per$uPDI <= 90))
# oPDI and hPDI differ only in the sign of the 5 less-healthy-plant groups, so
# their difference must equal that block's score flip: hPDI - oPDI = sum(6-2s) over less_healthy
lh_scores <- sapply(less_healthy, function(gn) score_group(per[[gn]], w, reverse=FALSE))
stopifnot(all(abs((per$hPDI - per$oPDI) - rowSums(6L - 2L*lh_scores)) < 1e-8))

cat("\n=== oPDI summary ===\n"); print(summary(per$oPDI))
cat("=== hPDI summary ===\n"); print(summary(per$hPDI))
cat("=== uPDI summary ===\n"); print(summary(per$uPDI))
cat("\ncorrelations:  oPDI-hPDI", round(cor(per$oPDI, per$hPDI),3),
    " oPDI-uPDI", round(cor(per$oPDI, per$uPDI),3),
    " hPDI-uPDI", round(cor(per$hPDI, per$uPDI),3), "\n")

out <- per %>% select(SEQN, cycle, oPDI, hPDI, uPDI, all_of(groups))
write.csv(out, "data/interim/hpdi_scores.csv", row.names=FALSE)
saveRDS(out, "data/interim/hpdi_scores.rds")
cat("Saved data/interim/hpdi_scores.{csv,rds}\n")
