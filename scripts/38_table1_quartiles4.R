#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 38_table1_quartiles4.R — survey-weighted "Table 1" across ALL FOUR PIR quartiles.
#
# 13_table1.R reports only the EXTREME groups (Q1 vs Q4, T1 vs T3) with a two-group
# survey t-test. This script keeps the same weighting, the same variables and the
# same formatting but shows Q1/Q2/Q3/Q4 side by side, which needs different tests:
#
#   continuous   omnibus  Wald test on the 3 df of pir_q in svyglm(y ~ pir_q)
#                         (a survey t-test is defined for two groups only)
#                P-trend  linear term of the ordinal quartile rank 1..4, as in
#                         16_pir_quartile_trend.R
#   categorical  omnibus  Rao-Scott chi-square over the 4 x levels table
#                P-trend  not reported: these are nominal, so a linear term over
#                         the levels is not meaningful
#
# Why BOTH p columns: the omnibus p answers only "do the four groups differ at
# all" and is satisfied by any one group sitting apart; P-trend answers "is the
# pattern graded across increasing income", which is the question these quartiles
# are actually asked to settle. They can disagree in both directions.
#
# QUARTILE CUTS are the project's established, top-code-aware definition (as in
# 16_pir_quartile_trend.R): PIR is censored at 5.0 and its weighted 75th
# percentile IS 5.0, so the Q3/Q4 boundary is placed AT the ceiling and Q4 is the
# tied-at-5.0 group. Do not replace this with a naive quantile cut -- that yields
# an empty top group.
#
# SMD is the Q4-vs-Q1 standardised difference, carried over because at this sample
# size every p is <0.0001 and only the SMD distinguishes a large gap from a trivial
# one. The two helpers are copied verbatim from 13_table1.R; their output was checked
# against that script's Q4-vs-Q1 column (11/11 cells identical) when this was written.
#
# Set BOLTE17=1 to report the plant-based indices on the Bolte et al. 2023
# 17-group parameterisation instead of the primary 18-group FPED one.
# CATEGORICAL CELLS are "n (weighted %)": the count is UNWEIGHTED (people observed)
# and the percentage is SURVEY-WEIGHTED (population share). Counts therefore do not
# reconcile against the weighted percentages, and for variables with item-missing
# values they sum to less than the group n in the column header.
# Writes results/[bolte17_][median_]table1_4quartiles.csv
# ------------------------------------------------------------------------------
suppressMessages({library(dplyr); library(survey)})
options(survey.lonely.psu = "adjust")

BOLTE17 <- nzchar(Sys.getenv("BOLTE17"))
# MEDIAN=1 reports the survey-weighted MEDIAN (IQR) instead of mean (SD) for the
# continuous rows, and pairs it with a test of the matching quantity: a DESIGN-BASED
# Kruskal-Wallis (svyranktest) rather than the mean-based Wald test. This fills the one
# gap in the project's Table 1 family -- 17_table1_demographics_quartile.R already gives
# UNWEIGHTED median (range) with an unweighted Kruskal-Wallis, and the default here gives
# WEIGHTED mean (SD); neither gives a weighted median.
MEDIAN  <- nzchar(Sys.getenv("MEDIAN"))
SUF     <- paste0(if (BOLTE17) "bolte17_" else "", if (MEDIAN) "median_" else "")

a <- readRDS("data/interim/nhanes_analytic.rds")
if (BOLTE17) {
  d <- readRDS("data/interim/dhq_scores.rds")
  stopifnot(all(c("hPDI_d","uPDI_d","oPDI_d") %in% names(d)))
  # drop the primary columns BEFORE joining: joining while present suffixes them
  # away and the table silently reports the wrong parameterisation.
  a <- a %>% select(-any_of(c("oPDI","hPDI","uPDI"))) %>%
    left_join(d %>% select(SEQN, oPDI = oPDI_d, hPDI = hPDI_d, uPDI = uPDI_d), by = "SEQN")
  stopifnot(sum(!is.na(a$hPDI)) > 0)
}

a <- a %>% mutate(
  bmi_cat = cut(bmi, c(-Inf, 25, 30, Inf), right = FALSE,
                labels = c("< 25 (normal/under)", "25-< 30 (overweight)", ">= 30 (obese)")),
  smoke_f = factor(ifelse(is.na(eversmoke), NA, ifelse(eversmoke == 1, "Yes", "No")),
                   levels = c("No", "Yes")),
  sex_f   = factor(sex, levels = c("Male", "Female")),
  race_f  = factor(raceeth, levels = c("NHWhite","NHBlack","MexAmerican","OtherHispanic","OtherRace"),
                   labels = c("Non-Hispanic White","Non-Hispanic Black","Mexican American",
                              "Other Hispanic","Other race/multi")))

des0 <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTDRD1_POOL, data=a, nest=TRUE)
qc <- as.numeric(coef(svyquantile(~PIR, subset(des0, !is.na(PIR)),
                                  quantiles=c(.25,.5,.75), na.rm=TRUE)))
cat(sprintf("Weighted PIR quartile cuts: %.3f / %.3f / %.3f (ceiling 5.0)\n", qc[1], qc[2], qc[3]))
a$pir_q <- cut(a$PIR, breaks=c(-Inf, qc[1], qc[2], 5.0, Inf), right=FALSE,
               labels=c("Q1","Q2","Q3","Q4"))
a$qrank <- as.integer(a$pir_q)

des <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTDRD1_POOL, data=a, nest=TRUE)
dq  <- subset(des, !is.na(pir_q))          # subpopulation, not row-dropping
GRPS <- c("Q1","Q2","Q3","Q4")
ns   <- sapply(GRPS, function(g) sum(a$pir_q == g, na.rm = TRUE))

cont_vars <- list(age="Age, years", bmi="BMI, kg/m2", PIR="PIR (grouping variable)",
                  HEI2020_TOTAL="HEI-2020", NRF93="NRF9.3", oPDI="oPDI", hPDI="hPDI", uPDI="uPDI")
cat_vars  <- list(sex_f="Sex, n (weighted %)", race_f="Race / ethnicity, n (weighted %)", educ_f="Education, n (weighted %)",
                  bmi_cat="BMI category, n (weighted %)", smoke_f="Ever smoker, n (weighted %)")

# --- survey-weighted SMD helpers, identical to 13_table1.R ---------------------
smd_cont <- function(dg, v, grpvar, g1, g2) {
  m  <- svyby(as.formula(paste0("~",v)), as.formula(paste0("~",grpvar)), dg, svymean, na.rm=TRUE)
  vv <- svyby(as.formula(paste0("~",v)), as.formula(paste0("~",grpvar)), dg, svyvar,  na.rm=TRUE)
  rownames(m) <- m[[1]]; rownames(vv) <- vv[[1]]
  abs((m[g1,2]-m[g2,2]) / sqrt((vv[g1,2]+vv[g2,2])/2))
}
smd_cat <- function(des, grpvar, v, g1, g2) {
  p1 <- as.numeric(svymean(as.formula(paste0("~",v)), subset(des, base::get(grpvar)==g1 & !is.na(base::get(v))), na.rm=TRUE))
  p2 <- as.numeric(svymean(as.formula(paste0("~",v)), subset(des, base::get(grpvar)==g2 & !is.na(base::get(v))), na.rm=TRUE))
  k <- length(p1)
  if (k == 2) return(abs((p1[2]-p2[2]) / sqrt((p1[2]*(1-p1[2]) + p2[2]*(1-p2[2]))/2)))
  p1 <- p1[-k]; p2 <- p2[-k]
  S <- ((diag(p1)-p1%*%t(p1)) + (diag(p2)-p2%*%t(p2)))/2
  d <- p1 - p2
  as.numeric(sqrt(t(d) %*% solve(S) %*% d))
}

fmt_cont <- function(dsub, v) {
  if (MEDIAN) {
    # survey-weighted quartiles of the variable within this group -> median (IQR)
    q <- coef(svyquantile(as.formula(paste0("~",v)), dsub, quantiles=c(.25,.5,.75), na.rm=TRUE))
    return(sprintf("%.1f (%.1f-%.1f)", q[2], q[1], q[3]))
  }
  m <- svymean(as.formula(paste0("~",v)), dsub, na.rm=TRUE)
  s <- sqrt(coef(svyvar(as.formula(paste0("~",v)), dsub, na.rm=TRUE))[1])
  sprintf("%.1f (%.1f)", coef(m)[1], s)
}
STAT_LAB <- if (MEDIAN) "median (IQR)" else "mean (SD)"
pfmt <- function(p) if (is.na(p)) "" else if (p < 1e-4) "<0.0001" else sprintf("%.4f", p)

rows <- list()
add <- function(ch, lv, ov, g, p, ptr, smd, test="")
  rows[[length(rows)+1]] <<- data.frame(Characteristic=ch, Level=lv, Overall=ov,
    Q1=g[1], Q2=g[2], Q3=g[3], Q4=g[4], p=p, p_trend=ptr, smd=smd, test=test,
    stringsAsFactors=FALSE)

for (v in names(cont_vars)) {
  ov <- fmt_cont(des, v)
  g  <- sapply(GRPS, function(gg) fmt_cont(subset(dq, pir_q == gg), v))
  if (v == "PIR") { add(cont_vars[[v]], STAT_LAB, ov, g, NA, NA, NA, "not tested"); next }
  if (MEDIAN) {
    # rank-based omnibus, so the test matches the displayed statistic
    pom <- svyranktest(as.formula(paste(v, "~ pir_q")), dq, test="KruskalWallis")$p.value
  } else {
    f   <- svyglm(as.formula(paste(v, "~ pir_q")), design = dq)
    pom <- regTermTest(f, "pir_q")$p[1]
  }
  ptr <- summary(svyglm(as.formula(paste(v, "~ qrank")), design = dq))$coefficients["qrank","Pr(>|t|)"]
  sd  <- tryCatch(smd_cont(dq, v, "pir_q", "Q1", "Q4"), error=function(e) NA)
  # P-trend stays the linear term of the ordinal rank, which is a MEAN-based test even
  # when medians are displayed. Labelled so, rather than left to look rank-based.
  add(cont_vars[[v]], STAT_LAB, ov, g, pom, ptr, sd,
      if (MEDIAN) "Design KruskalWallis / P-trend (mean-based)" else "Wald (3 df) / P-trend")
}

for (v in names(cat_vars)) {
  ov <- svymean(as.formula(paste0("~",v)), subset(des, !is.na(base::get(v))), na.rm=TRUE)
  levs <- levels(a[[v]])
  pm <- sapply(GRPS, function(gg)
    as.numeric(svymean(as.formula(paste0("~",v)),
                       subset(dq, pir_q == gg & !is.na(base::get(v))), na.rm=TRUE)))
  pom <- tryCatch(svychisq(as.formula(paste0("~", v, "+pir_q")), dq, statistic="F")$p.value,
                  error=function(e) NA)
  sd  <- tryCatch(smd_cat(dq, "pir_q", v, "Q1", "Q4"), error=function(e) NA)
  n_ov  <- table(a[[v]][!is.na(a[[v]])])
  n_grp <- sapply(GRPS, function(gg) table(a[[v]][a$pir_q == gg & !is.na(a[[v]])]))
  cell  <- function(n, pct) sprintf("%s (%.1f)", format(n, big.mark=",", trim=TRUE), 100*pct)
  for (i in seq_along(levs))
    add(if (i == 1) cat_vars[[v]] else "", levs[i],
        cell(n_ov[i], as.numeric(ov)[i]),
        cell(n_grp[i, ], pm[i, ]),
        if (i == 1) pom else NA, NA, if (i == 1) sd else NA,
        if (i == 1) "Rao-Scott chi-sq" else "")
}

tab <- bind_rows(rows)
isp <- !is.na(tab$p)
tab$fdr <- NA_real_; tab$fdr[isp] <- p.adjust(tab$p[isp], method="BH")
tab$p       <- sapply(tab$p, pfmt)
tab$p_trend <- sapply(tab$p_trend, pfmt)
tab$fdr     <- sapply(tab$fdr, pfmt)
tab$smd     <- ifelse(is.na(tab$smd), "", sprintf("%.2f", tab$smd))
tab <- tab[, c("Characteristic","Level","Overall","Q1","Q2","Q3","Q4",
               "p","p_trend","fdr","smd","test")]
names(tab) <- c("Characteristic","Level",sprintf("Overall (N=%s)", format(nrow(a), big.mark=",")),
                sprintf("Q1 (n=%s)", format(ns[1], big.mark=",")),
                sprintf("Q2 (n=%s)", format(ns[2], big.mark=",")),
                sprintf("Q3 (n=%s)", format(ns[3], big.mark=",")),
                sprintf("Q4 (n=%s)", format(ns[4], big.mark=",")),
                "p (omnibus, 4 groups)","P for trend","FDR p","SMD (Q4 vs Q1)","Test")
out <- paste0("results/", SUF, "table1_4quartiles.csv")
write.csv(tab, out, row.names=FALSE)
cat(sprintf("\nWrote %s  (%d rows; groups %s)\n", out, nrow(tab),
            paste(sprintf("%s=%s", GRPS, format(ns, big.mark=",")), collapse=" ")))
