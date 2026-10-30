#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 13_table1.R — survey-weighted "Table 1" of baseline cohort characteristics for
# the analytic sample (n = 19,605 adults >=20), overall and compared across income
# extreme groups: quartiles Q1 (poorest) vs Q4 (top-coded ceiling) and tertiles
# T1 vs T3. Continuous vars: survey-weighted mean (SD). Categorical: survey-
# weighted %. Between-group p-values are DESIGN-BASED: svyttest for continuous,
# svychisq (Rao-Scott F) for categorical, comparing the two extreme groups only.
# PIR is the grouping variable, shown for transparency (differences are by
# construction). Writes results/table1_quartile.csv and results/table1_tertile.csv.
# ------------------------------------------------------------------------------
suppressMessages({library(haven); library(dplyr); library(survey)})
options(survey.lonely.psu = "adjust")

# BOLTE17 switch. Set BOLTE17=1 to report the plant-based indices computed on the
# Bolte et al. 2023 17-group parameterisation (hPDI_d/uPDI_d/oPDI_d from
# 36_dhq_scores_panels.R) INSTEAD of the primary 18-group FPED ones. Everything else --
# design, groups, tests, formatting -- is identical, so the two tables are directly
# comparable and the primary outputs are untouched. HEI-2020 and NRF9.3 are not
# plant-index scores and carry over unchanged.
BOLTE17 <- nzchar(Sys.getenv("BOLTE17"))
SUF     <- if (BOLTE17) "bolte17_" else ""

a <- readRDS("data/interim/nhanes_analytic.rds")
if (BOLTE17) {
  d <- readRDS("data/interim/dhq_scores.rds")
  stopifnot(all(c("hPDI_d","uPDI_d","oPDI_d") %in% names(d)))
  # drop the primary PDI columns first: joining while they are present would suffix
  # them away and leave the table silently reporting the wrong parameterisation.
  a <- a %>% select(-any_of(c("oPDI","hPDI","uPDI"))) %>%
    left_join(d %>% select(SEQN, oPDI = oPDI_d, hPDI = hPDI_d, uPDI = uPDI_d), by = "SEQN")
  stopifnot(sum(!is.na(a$hPDI)) > 0)
  cat("BOLTE17: reporting the 17-group plant-based indices (n non-missing hPDI:",
      sum(!is.na(a$hPDI)), ")\n")
}

# --- derived categorical: BMI category, ever-smoker as Yes/No ------------------
a <- a %>% mutate(
  bmi_cat = cut(bmi, c(-Inf, 25, 30, Inf), right = FALSE,
                labels = c("< 25 (normal/under)", "25–< 30 (overweight)", ">= 30 (obese)")),
  smoke_f = factor(ifelse(is.na(eversmoke), NA, ifelse(eversmoke == 1, "Yes", "No")),
                   levels = c("No", "Yes")),
  sex_f   = factor(sex, levels = c("Male", "Female")),
  race_f  = factor(raceeth, levels = c("NHWhite","NHBlack","MexAmerican","OtherHispanic","OtherRace"),
                    labels = c("Non-Hispanic White","Non-Hispanic Black","Mexican American",
                               "Other Hispanic","Other race/multi")))

# --- income groups (same weighted cuts as 11_pir_groups.R) ---------------------
des0 <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTDRD1_POOL, data=a, nest=TRUE)
tc  <- as.numeric(coef(svyquantile(~PIR, subset(des0,!is.na(PIR)), quantiles=c(1/3,2/3), na.rm=TRUE)))
q1c <- as.numeric(coef(svyquantile(~PIR, subset(des0,!is.na(PIR)), quantiles=.25, na.rm=TRUE)))
a$grp_q <- factor(ifelse(is.na(a$PIR), NA, ifelse(a$PIR <  q1c, "Q1",
                   ifelse(a$PIR >= 5.0, "Q4", NA))), levels = c("Q1","Q4"))
a$grp_t <- factor(ifelse(is.na(a$PIR), NA, ifelse(a$PIR <  tc[1], "T1",
                   ifelse(a$PIR >= tc[2], "T3", NA))), levels = c("T1","T3"))

des <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTDRD1_POOL, data=a, nest=TRUE)

# --- variable plan -------------------------------------------------------------
cont_vars <- list(age="Age, years", bmi="BMI, kg/m²", PIR="PIR (grouping variable)",
                  HEI2020_TOTAL="HEI-2020", NRF93="NRF9.3", oPDI="oPDI", hPDI="hPDI", uPDI="uPDI")
cat_vars  <- list(sex_f="Sex", race_f="Race / ethnicity", educ_f="Education",
                  bmi_cat="BMI category", smoke_f="Ever smoker")

fmt_cont <- function(d, v) {
  f <- as.formula(paste0("~", v))
  m <- as.numeric(svymean(f, d, na.rm=TRUE))
  s <- sqrt(as.numeric(svyvar(f, d, na.rm=TRUE)))
  sprintf("%.1f (%.1f)", m, s)
}
pct_levels <- function(d, v) {              # weighted % per factor level
  f <- as.formula(paste0("~", v))
  m <- svymean(f, subset(d, !is.na(base::get(v))), na.rm=TRUE)
  setNames(as.numeric(m)*100, sub(v, "", names(m)))
}
# Standardized mean difference (sample-size-independent effect size; |SMD|>0.1 = meaningful
# imbalance by convention). Continuous: (mean1-mean2)/sqrt(mean of the two weighted variances).
# Categorical: multinomial SMD (Yang & Dalton 2012) via the pooled covariance of the level
# proportions (last level dropped for full rank; binary reduces to the usual single-proportion SMD).
smd_cont <- function(d12g, v, grpvar, g1, g2) {
  m  <- svyby(as.formula(paste0("~",v)), as.formula(paste0("~",grpvar)), d12g, svymean, na.rm=TRUE)
  vv <- svyby(as.formula(paste0("~",v)), as.formula(paste0("~",grpvar)), d12g, svyvar,  na.rm=TRUE)
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

build_table <- function(grpvar, g1, g2, col_overall="Overall") {
  dO <- des; d1 <- subset(des, base::get(grpvar)==g1); d2 <- subset(des, base::get(grpvar)==g2)
  nO <- nrow(a); n1 <- sum(a[[grpvar]]==g1, na.rm=TRUE); n2 <- sum(a[[grpvar]]==g2, na.rm=TRUE)
  d12 <- subset(des, base::get(grpvar) %in% c(g1,g2))
  rows <- list()
  # `test` names the statistical test behind the p-value on that row, so the
  # printed table is self-documenting (reviewers should not have to read the
  # methods to know which test produced a p).
  TEST_CONT <- "Survey t-test"
  TEST_CAT  <- "Rao-Scott chi-sq"
  add <- function(variable, level, ov, c1, c2, p, smd, test="") rows[[length(rows)+1]] <<-
    data.frame(Characteristic=variable, Level=level, Overall=ov, G1=c1, G2=c2,
               p=p, smd=smd, test=test, stringsAsFactors=FALSE)
  # continuous
  for (v in names(cont_vars)) {
    if (v == "PIR") {   # grouping variable — differences are by construction, don't balance-test
      add(cont_vars[[v]], "mean (SD)", fmt_cont(dO,v), fmt_cont(d1,v), fmt_cont(d2,v), NA, NA,
          "not tested")
      next
    }
    tt <- tryCatch(svyttest(as.formula(paste0(v," ~ ",grpvar)), d12), error=function(e) NULL)
    p  <- if (is.null(tt)) NA else tt$p.value
    sd <- tryCatch(smd_cont(d12, v, grpvar, g1, g2), error=function(e) NA)
    add(cont_vars[[v]], "mean (SD)", fmt_cont(dO,v), fmt_cont(d1,v), fmt_cont(d2,v), p, sd, TEST_CONT)
  }
  # categorical
  for (v in names(cat_vars)) {
    ch <- tryCatch(svychisq(as.formula(paste0("~",v," + ",grpvar)), d12, statistic="F"),
                   error=function(e) NULL)
    p  <- if (is.null(ch)) NA else as.numeric(ch$p.value)
    sd <- tryCatch(smd_cat(des, grpvar, v, g1, g2), error=function(e) NA)
    oO <- pct_levels(dO,v); o1 <- pct_levels(d1,v); o2 <- pct_levels(d2,v)
    levs <- names(oO); first <- TRUE
    for (lv in levs) {
      add(if (first) paste0(cat_vars[[v]], ", %") else "", lv,
          sprintf("%.1f", oO[lv]), sprintf("%.1f", o1[lv]), sprintf("%.1f", o2[lv]),
          if (first) p else NA, if (first) sd else NA, if (first) TEST_CAT else "")
      first <- FALSE
    }
  }
  tab <- bind_rows(rows)
  # BH-FDR across the tested variables (one test per variable = rows with non-NA p)
  praw <- tab$p
  tab$fdr <- NA_real_
  tab$fdr[!is.na(praw)] <- p.adjust(praw[!is.na(praw)], "BH")
  # pretty formatting
  fp <- function(x) ifelse(is.na(x), "", ifelse(x < 1e-4, "<0.0001", sprintf("%.4f", x)))
  tab$p   <- fp(tab$p)
  tab$fdr <- fp(tab$fdr)
  tab$smd <- ifelse(is.na(tab$smd), "", sprintf("%.2f", tab$smd))
  # reorder to Characteristic, Level, Overall, G1, G2, p, FDR, SMD, Test
  tab <- tab[, c("Characteristic","Level","Overall","G1","G2","p","fdr","smd","test")]
  names(tab) <- c("Characteristic","Level",
                  sprintf("%s (N=%s)", col_overall, format(nO, big.mark=",")),
                  sprintf("%s (n=%s)", g1, format(n1, big.mark=",")),
                  sprintf("%s (n=%s)", g2, format(n2, big.mark=",")),
                  "p (G1 vs G2)", "FDR p", "SMD", "Test")
  tab
}

t1q <- build_table("grp_q", "Q1", "Q4")
t1t <- build_table("grp_t", "T1", "T3")
dir.create("results", showWarnings = FALSE)
write.csv(t1q, paste0("results/", SUF, "table1_quartile.csv"), row.names=FALSE)
write.csv(t1t, paste0("results/", SUF, "table1_tertile.csv"),  row.names=FALSE)

cat(sprintf("Income groups: Q1 PIR<%.2f | Q4 PIR>=5.0 | T1 PIR<%.2f | T3 PIR>=%.2f\n\n",
            q1c, tc[1], tc[2]))
cat("=== TABLE 1 — cohort characteristics, Overall vs income QUARTILES (Q1 vs Q4) ===\n")
print(as.data.frame(t1q), row.names=FALSE)
cat("\n=== TABLE 1 — cohort characteristics, Overall vs income TERTILES (T1 vs T3) ===\n")
print(as.data.frame(t1t), row.names=FALSE)
cat(sprintf("\nWrote results/%stable1_quartile.csv, results/%stable1_tertile.csv\n", SUF, SUF))
