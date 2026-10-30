#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 36_dhq_scores_panels.R — plant-based indices and every figure input, built on the
#   DHQ-dictionary WHOLE-ITEM classification (35_dhq_bolte_reclassify.py) instead of
#   FPED's decomposed components.
#
# WHY THIS DIFFERS FROM 29_bolte_panels.R: FPED splits each food into pattern components,
# so one food contributes fractionally to several groups -- which is why french fries
# counted toward the HEALTHFUL "vegetable oils" group. Here every food code belongs to
# exactly ONE category, so a food counts once, where it belongs.
#
# This also constructs Bolte's 17th group, Miscellaneous_animal_based_foods, which FPED
# cannot build at all (it disaggregates mixed dishes rather than pooling them).
#
# Group servings are grams/day per category. Units need not be comparable ACROSS groups
# because each group is rank-scored 1-5 separately before being summed.
#
# EXCLUDED from the index, consistent with the published grouping: WATER, ART_SWEET,
# PROTEIN, MARGARINE and Alcoholic_Beverages are dictionary categories that are not PDI
# components. Codes still in the review tier are kept but flagged in the coverage report.
#
# Output: data/interim/dhq_scores.rds and results/dhq/{pd_main,trend_tests,food_odds,
#         pca_loadings,pca_variance,pca_scores,coverage}.csv
# ------------------------------------------------------------------------------
suppressMessages({library(haven); library(dplyr); library(tidyr); library(purrr)
                  library(tibble); library(survey); library(readr)})
options(survey.lonely.psu = "adjust")
set.seed(1)
dir.create("results/dhq", showWarnings = FALSE, recursive = TRUE)

HEALTHY <- c("Whole_Grains","Fruits","Vegetables","Nuts_and_Legumes","Vegetable_Oils","Tea_and_Coffee")
LESSHEALTHY <- c("Fruit_Juices","Refined_Grains","Potatoes","Sugar_Sweetened_Beverages","Sweets_and_desserts")
ANIMAL <- c("Animal_Fat","Dairy","Eggs","Fish_or_Seafood","Meat","Miscellaneous_animal_based_foods")
GROUPS <- c(HEALTHY, LESSHEALTHY, ANIMAL)
stopifnot(length(GROUPS) == 17)

# Classification input. The full file (results/bolte/) exists only after 35_ has been run,
# which needs the restricted DHQ inputs. The committed four-column cut in data/interim/ is
# the same rows with the matcher diagnostics stripped; this script uses only code,
# category and tier, so the two are interchangeable here. Prefer the full file when
# present (it is the fresher of the two after a re-run of 35_), else fall back.
CLS_FULL   <- "results/bolte/dhq_food_classification.csv"
CLS_PUBLIC <- "data/interim/dhq_food_classification_public.csv"
# $DHQ_CLASSIFICATION wins over both, for a classification held outside the tree. If it
# is set it must exist: falling through to another file on a typo would silently score
# the indices on a different classification than the one asked for.
CLS_ENV  <- Sys.getenv("DHQ_CLASSIFICATION", unset = "")
if (nzchar(CLS_ENV) && !file.exists(CLS_ENV))
  stop("DHQ_CLASSIFICATION is set but the file does not exist: ", CLS_ENV)
cls_path <- if (nzchar(CLS_ENV)) CLS_ENV else
            if (file.exists(CLS_FULL)) CLS_FULL else CLS_PUBLIC
if (!file.exists(cls_path))
  stop("no food classification found: run scripts/35_dhq_bolte_reclassify.py (needs data_restricted/) ",
       "or restore the committed ", CLS_PUBLIC)
cat("classification source:", cls_path, "\n")
cls <- read_csv(cls_path, col_types = cols(code = col_character(), .default = col_guess()))
stopifnot(all(c("code", "category", "tier") %in% names(cls)))
ana <- readRDS("data/interim/nhanes_analytic.rds")

foods <- pmap_dfr(tibble(yy = c("1112","1314","1516","1718")), function(yy)
  read_sas(sprintf("data/raw/fped/fped_dr1iff_%s.sas7bdat", yy)) %>%
    select(SEQN, DR1IFDCD, DR1IGRMS)) %>%
  mutate(code = sprintf("%.0f", DR1IFDCD)) %>%
  semi_join(ana %>% select(SEQN), by = "SEQN") %>%
  left_join(cls %>% select(code, category, tier), by = "code")

# coverage: what share of grams is classified into an index group at all
cov <- foods %>% left_join(ana %>% select(SEQN, WTDRD1_POOL), by = "SEQN") %>%
  mutate(g = DR1IGRMS * WTDRD1_POOL,
         status = case_when(category %in% GROUPS & tier %in% c("accept","manual") ~ "in index (settled)",
                            category %in% GROUPS ~ "in index (review tier)",
                            is.na(category) | category == "UNMATCHED" ~ "unclassified",
                            TRUE ~ "outside index (water/alcohol/etc)")) %>%
  group_by(status) %>% summarise(pct_of_grams = round(100*sum(g, na.rm=TRUE)/sum(.$g, na.rm=TRUE), 1),
                                 .groups = "drop")
write_csv(cov, "results/dhq/coverage.csv")
cat("=== gram coverage of the index ===\n"); print(as.data.frame(cov), row.names = FALSE)

serv <- foods %>% filter(category %in% GROUPS) %>%
  group_by(SEQN, category) %>% summarise(g = sum(DR1IGRMS, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = category, values_from = g, values_fill = 0)
missing_grp <- setdiff(GROUPS, names(serv))
for (m in missing_grp) serv[[m]] <- 0
stopifnot(length(setdiff(GROUPS, names(serv))) == 0)

dat <- ana %>% left_join(serv, by = "SEQN") %>%
  mutate(across(all_of(GROUPS), ~ coalesce(.x, 0)))

# ---- survey-weighted rank scoring, 1-5 per group -----------------------------
# NOT the same rule as 05_hpdi.R, despite both producing a 1-5 score per group:
#   here   non-consumers -> 1; consumers -> weighted QUARTILES c(.25,.5,.75) scored
#          1-4 then +1, i.e. 2-5. Non-consumers are strictly separated from consumers.
#   05_    non-consumers -> 1; consumers -> weighted QUINTILES seq(.2,.8,.2), with the
#          bins shifted up by (5 - nbins), so the lowest consumer bin CAN land on 1 and
#          share it with non-consumers.
# Both give 17-85 / 18-90 totals, but the within-group distributions differ. Describe
# each index by its own script -- do not call this one "quintile-scored".
wq <- function(x, w, p) { o <- order(x); x <- x[o]; w <- w[o]
                          cw <- cumsum(w)/sum(w); sapply(p, function(q) x[which(cw >= q)[1]]) }
sg <- function(x, w, reverse = FALSE) {
  pos <- x > 0; s <- rep(1L, length(x))
  if (sum(pos) > 0) { cuts <- wq(x[pos], w[pos], c(.25,.5,.75,1))
                      s[pos] <- as.integer(cut(x[pos], c(-Inf, cuts[1:3], Inf), labels = FALSE)) + 1L }
  if (reverse) s <- 6L - s; s }
w <- dat$WTDRD1_POOL
build <- function(pos, neg) {
  P <- sapply(pos, function(v) sg(dat[[v]], w, FALSE))
  N <- sapply(neg, function(v) sg(dat[[v]], w, TRUE))
  rowSums(P) + rowSums(N) }
dat$hPDI_d <- build(HEALTHY, c(LESSHEALTHY, ANIMAL))
dat$uPDI_d <- build(LESSHEALTHY, c(HEALTHY, ANIMAL))
dat$oPDI_d <- build(c(HEALTHY, LESSHEALTHY), ANIMAL)
stopifnot(all(dat$hPDI_d >= 17, dat$hPDI_d <= 85))
saveRDS(dat, "data/interim/dhq_scores.rds")

des <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTDRD1_POOL, nest = TRUE, data = dat)
SC <- c(hPDI_d = "hPDI", uPDI_d = "uPDI")

# ---- pd_main: weighted means by each SES measure ------------------------------
bins <- list(pir_bin = "pir_bin", educ_f = "educ_f", income_lab = "income_lab")
dat2 <- dat %>% mutate(
  pir_bin = cut(PIR, breaks = c(0,1,2,3,4,5), include.lowest = TRUE,
                labels = c("0-1","1-2","2-3","3-4","4-5")),
  income_lab = factor(income_cat, levels = 1:12,
    labels = c("<5k","5-10k","10-15k","15-20k","20-25k","25-35k",
               "35-45k","45-55k","55-65k","65-75k","75-100k","100k+")))
des2 <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTDRD1_POOL, nest=TRUE, data=dat2)
pd <- map_dfr(names(SC), function(v) map_dfr(names(bins), function(b) {
  d <- subset(des2, !is.na(get(b)))
  m <- svyby(as.formula(paste0("~", v)), as.formula(paste0("~", b)), d, svymean, na.rm = TRUE)
  tibble(score = SC[[v]], by = b, group = as.character(m[[1]]),
         mean = m[[2]], se = m$se, ci_l = m[[2]] - 1.96*m$se, ci_u = m[[2]] + 1.96*m$se) }))
write_csv(pd, "results/dhq/pd_main.csv")

# ---- trend / association tests, mirroring 14_trend_tests.R ---------------------
tt <- map_dfr(names(SC), function(v) map_dfr(c("PIR","income_cat","educ_f"), function(x) {
  d <- subset(des2, !is.na(get(x)))
  f <- svyglm(as.formula(paste(v, "~", x)), design = d)
  if (x == "educ_f") {
    rt <- survey::regTermTest(f, as.formula(paste0("~", x)))
    tibble(score = SC[[v]], ses = x, test = "joint Wald", effect = NA_real_, p = as.numeric(rt$p))
  } else {
    cf <- summary(f)$coefficients
    tibble(score = SC[[v]], ses = x, test = "linear trend",
           effect = cf[x, "Estimate"], p = cf[x, "Pr(>|t|)"]) } }))
tt$p_fdr_BH <- p.adjust(tt$p, "BH")
fmtp <- function(p) ifelse(p < 1e-4, "<0.0001", sprintf("%.4f", p))
tt <- tt %>%
  mutate(by = recode(ses, PIR = "pir_bin", income_cat = "income_lab", educ_f = "educ_f"),
         ses = recode(ses, income_cat = "income", educ_f = "education"),
         p_disp = fmtp(p), p_fdr_disp = fmtp(p_fdr_BH)) %>%
  select(score, ses, by, test, effect, p, p_fdr_BH, p_disp, p_fdr_disp)
write_csv(tt, "results/dhq/trend_tests.csv")

# ---- food odds: high consumer of each group vs PIR -----------------------------
fo <- map_dfr(GROUPS, function(g) {
  d0 <- dat %>% mutate(.v = .data[[g]])
  des0 <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTDRD1_POOL, nest=TRUE, data=d0)
  thr <- tryCatch(as.numeric(coef(svyquantile(~.v, subset(des0, .v > 0), 2/3, na.rm=TRUE))),
                  error = function(e) NA_real_)
  if (is.na(thr)) return(NULL)
  d0$.hi <- as.integer(d0$.v >= thr)
  des0 <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTDRD1_POOL, nest=TRUE, data=d0)
  f <- tryCatch(svyglm(.hi ~ PIR + age + sex + raceeth, design = subset(des0, !is.na(PIR)),
                       family = quasibinomial()), error = function(e) NULL)
  if (is.null(f)) return(NULL)
  cf <- summary(f)$coefficients
  tibble(food = g, label = gsub("_", " ", g), OR = exp(cf["PIR","Estimate"]),
         OR_lo = exp(cf["PIR","Estimate"] - 1.96*cf["PIR","Std. Error"]),
         OR_hi = exp(cf["PIR","Estimate"] + 1.96*cf["PIR","Std. Error"]),
         p = cf["PIR","Pr(>|t|)"]) })
fo$p_fdr_BH <- p.adjust(fo$p, "BH")
fo <- fo %>% arrange(desc(OR)) %>% mutate(num = row_number())
write_csv(fo, "results/dhq/food_odds.csv")

# ---- PCA within the extreme PIR quartiles, mirroring 12_pca_pdi.R --------------
pdat <- dat %>% filter(!is.na(PIR))
desp <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTDRD1_POOL, nest=TRUE, data=pdat)
q1c <- as.numeric(coef(svyquantile(~PIR, desp, quantiles = .25, na.rm = TRUE)))
grp <- list(Q1 = pdat$PIR < q1c, Q4 = pdat$PIR >= 5.0)
stopifnot(all(vapply(grp, sum, integer(1)) > 0))
pca_one <- function(mask, gname) {
  g <- pdat[mask, ]; X <- log1p(as.matrix(g[, GROUPS])); wt <- g$WTDRD1_POOL
  cw <- cov.wt(X, wt = wt, cor = TRUE); ev <- eigen(cw$cor)
  ve <- ev$values/sum(ev$values)*100; load <- ev$vectors[, 1:2]
  Xs <- scale(X, center = cw$center, scale = sqrt(diag(cw$cov))); sc <- Xs %*% load
  r1 <- apply(g[, GROUPS], 2, function(v) cor(v, sc[,1], method="spearman"))
  r2 <- apply(g[, GROUPS], 2, function(v) cor(v, sc[,2], method="spearman"))
  if (r1[which.max(abs(r1))] < 0) { r1 <- -r1; sc[,1] <- -sc[,1] }
  if (r2[which.max(abs(r2))] < 0) { r2 <- -r2; sc[,2] <- -sc[,2] }
  list(loadings = tibble(group = gname, num = seq_along(GROUPS), food = GROUPS,
                         label = gsub("_", " ", GROUPS), rho_pc1 = r1, rho_pc2 = r2,
                         color = ifelse(abs(r1) > .5, "PC1", ifelse(abs(r2) > .5, "PC2", "weak"))),
       variance = tibble(group = gname, n = nrow(g), pc1_var = ve[1], pc2_var = ve[2]),
       scores = tibble(group = gname, pc1 = sc[,1], pc2 = sc[,2])) }
res <- imap(grp, ~pca_one(.x, .y))
write_csv(bind_rows(map(res, "loadings")), "results/dhq/pca_loadings.csv")
write_csv(bind_rows(map(res, "variance")), "results/dhq/pca_variance.csv")
write_csv(bind_rows(map(res, "scores")),   "results/dhq/pca_scores.csv")

cat("\n=== score distributions (17 groups, range 17-85) ===\n")
print(as.data.frame(dat %>% summarise(across(c(hPDI_d, uPDI_d, oPDI_d),
      list(mean = ~round(mean(.x),1), sd = ~round(sd(.x),1), min = min, max = max)))), row.names = FALSE)
cat("\n=== PCA variance explained ===\n")
print(as.data.frame(bind_rows(map(res, "variance"))), row.names = FALSE)
cat("\n=== trend tests ===\n")
print(as.data.frame(tt), row.names = FALSE)
cat("\nWrote data/interim/dhq_scores.rds and results/dhq/*.csv\n")
