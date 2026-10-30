#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 39_table1_missing.R — supplement to the four-quartile Table 1: how much data is
# MISSING for each row of that table, overall and within each PIR quartile.
#
# Table 1 reports statistics over non-missing values only. For continuous rows the
# median/mean is computed on whoever had a measurement; for categorical rows the
# percentages are shares of the non-missing within the group. Neither states how
# many were left out. This table does.
#
# THE TEST IS THE POINT. A count of missing values is only reassuring if it is
# similar across the groups being compared. Missingness that CONCENTRATES in one
# income group biases the comparison in that row, so each variable gets a
# design-based Rao-Scott chi-square of {missing, observed} x {Q1..Q4}. A small p
# means the row's missingness is differential, not that there is a lot of it.
#
# ALSO COUNTED, and larger than any item-level gap: participants with missing PIR,
# who cannot be placed in a quartile at all and so are absent from every quartile
# column of Table 1. They are retained in the analytic sample and appear in the
# Overall column, which is why the quartile ns do not sum to the overall N.
#
# Set BOLTE17=1 to assess the 17-group plant-based indices instead of the primary
# 18-group ones. Writes results/[bolte17_]table1_missing_supplement.csv
# ------------------------------------------------------------------------------
suppressMessages({library(dplyr); library(survey)})
options(survey.lonely.psu = "adjust")

BOLTE17 <- nzchar(Sys.getenv("BOLTE17"))
SUF     <- if (BOLTE17) "bolte17_" else ""

a <- readRDS("data/interim/nhanes_analytic.rds")
if (BOLTE17) {
  d <- readRDS("data/interim/dhq_scores.rds")
  stopifnot(all(c("hPDI_d","uPDI_d","oPDI_d") %in% names(d)))
  a <- a %>% select(-any_of(c("oPDI","hPDI","uPDI"))) %>%
    left_join(d %>% select(SEQN, oPDI = oPDI_d, hPDI = hPDI_d, uPDI = uPDI_d), by = "SEQN")
}

a <- a %>% mutate(
  bmi_cat = cut(bmi, c(-Inf, 25, 30, Inf), right = FALSE,
                labels = c("< 25", "25-< 30", ">= 30")),
  smoke_f = factor(ifelse(is.na(eversmoke), NA, ifelse(eversmoke == 1, "Yes", "No")),
                   levels = c("No", "Yes")),
  sex_f   = factor(sex), race_f = factor(raceeth))

des0 <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTDRD1_POOL, data=a, nest=TRUE)
qc <- as.numeric(coef(svyquantile(~PIR, subset(des0, !is.na(PIR)),
                                  quantiles=c(.25,.5,.75), na.rm=TRUE)))
a$pir_q <- cut(a$PIR, breaks=c(-Inf, qc[1], qc[2], 5.0, Inf), right=FALSE,
               labels=c("Q1","Q2","Q3","Q4"))
des <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTDRD1_POOL, data=a, nest=TRUE)
dq  <- subset(des, !is.na(pir_q))
GRPS <- c("Q1","Q2","Q3","Q4")
ns   <- sapply(GRPS, function(g) sum(a$pir_q == g, na.rm=TRUE))
N    <- nrow(a)

# rows in the same order as Table 1
vars <- c("age"="Age, years", "bmi"="BMI, kg/m2", "PIR"="PIR (grouping variable)",
          "HEI2020_TOTAL"="HEI-2020", "NRF93"="NRF9.3", "oPDI"="oPDI", "hPDI"="hPDI",
          "uPDI"="uPDI", "sex_f"="Sex", "race_f"="Race / ethnicity",
          "educ_f"="Education", "bmi_cat"="BMI category", "smoke_f"="Ever smoker")

cellfmt <- function(n, d) sprintf("%s (%.1f)", format(n, big.mark=",", trim=TRUE), 100*n/d)
rows <- list()
for (v in names(vars)) {
  miss  <- is.na(a[[v]])
  n_ov  <- sum(miss)
  n_grp <- sapply(GRPS, function(g) sum(miss & a$pir_q == g, na.rm=TRUE))
  # design-based test of differential missingness; undefined when nothing is missing
  p <- if (n_ov == 0 || all(n_grp == 0)) NA_real_ else {
    a$.miss <- factor(miss, levels=c(FALSE, TRUE), labels=c("observed","missing"))
    dqm <- subset(svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTDRD1_POOL,
                            data=a, nest=TRUE), !is.na(pir_q))
    tryCatch(svychisq(~.miss + pir_q, dqm, statistic="F")$p.value, error=function(e) NA_real_)
  }
  rows[[length(rows)+1]] <- data.frame(
    Characteristic = vars[[v]],
    Overall = cellfmt(n_ov, N),
    Q1 = cellfmt(n_grp[1], ns[1]), Q2 = cellfmt(n_grp[2], ns[2]),
    Q3 = cellfmt(n_grp[3], ns[3]), Q4 = cellfmt(n_grp[4], ns[4]),
    p = p, stringsAsFactors = FALSE)
}
tab <- bind_rows(rows)
tab$test <- ifelse(is.na(tab$p), "no missing values", "Rao-Scott chi-sq")
tab$test[tab$Characteristic == "PIR (grouping variable)"] <-
  "0 within quartiles by construction; see final row"
isp <- !is.na(tab$p)
tab$fdr <- NA_real_; tab$fdr[isp] <- p.adjust(tab$p[isp], method="BH")
pf <- function(p) ifelse(is.na(p), "", ifelse(p < 1e-4, "<0.0001", sprintf("%.4f", p)))
tab$p <- pf(tab$p); tab$fdr <- pf(tab$fdr)

# structural exclusion: missing PIR cannot be assigned a quartile at all
npir <- sum(is.na(a$PIR))
tab <- bind_rows(tab, data.frame(
  Characteristic = "PIR missing (no quartile assignable)",
  Overall = cellfmt(npir, N), Q1="-", Q2="-", Q3="-", Q4="-",
  p = "", test = "structural - excluded from all quartile columns", fdr = ""))

tab <- tab[, c("Characteristic","Overall","Q1","Q2","Q3","Q4","p","fdr","test")]
names(tab) <- c("Characteristic", sprintf("Overall (N=%s)", format(N, big.mark=",")),
                sprintf("Q1 (n=%s)", format(ns[1], big.mark=",")),
                sprintf("Q2 (n=%s)", format(ns[2], big.mark=",")),
                sprintf("Q3 (n=%s)", format(ns[3], big.mark=",")),
                sprintf("Q4 (n=%s)", format(ns[4], big.mark=",")),
                "p (differential across quartiles)","FDR p","Test")
out <- paste0("results/", SUF, "table1_missing_supplement.csv")
write.csv(tab, out, row.names=FALSE)
cat(sprintf("\nWrote %s  (cells are n missing (%% of that column's n))\n", out))
print(as.data.frame(tab), row.names=FALSE)
