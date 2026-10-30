#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 34b_nhanes_food_codes.R — the two inputs 35_dhq_bolte_reclassify.py needs.
#
# WHY THIS SCRIPT EXISTS: both files were originally produced by ad-hoc one-liners
# typed at the shell, so `35_` read inputs that NO script in the repository
# produced. Anyone cloning this project could not reproduce them. Numbered 34b so
# it sorts between 34_ and 35_, i.e. immediately before its consumer.
#
#   data/interim/nhanes_codes_all.csv   one row per distinct food code x description
#                                       across the four cycles. The code is written
#                                       with sprintf("%.0f") because read_sas returns
#                                       it as a double and as.character() renders
#                                       round values in scientific notation.
#   data/interim/nhanes_code_intake.csv survey-weighted grams consumed per food code,
#                                       restricted to the analytic sample, used to
#                                       rank the adjudication queue by how much of the
#                                       diet each code actually accounts for.
#
# Reads data/raw/fped/fped_dr1iff_{1112,1314,1516,1718}.sas7bdat (from 01_download.sh)
# and data/interim/nhanes_analytic.rds (from 02_harmonize_score.R).
# ------------------------------------------------------------------------------
suppressMessages({library(haven); library(dplyr); library(purrr); library(readr)})

YY <- c("1112", "1314", "1516", "1718")
f <- map_dfr(YY, function(yy)
  read_sas(sprintf("data/raw/fped/fped_dr1iff_%s.sas7bdat", yy)) %>%
    select(SEQN, DR1IFDCD, DR1IGRMS, DESCRIPTION)) %>%
  mutate(code = sprintf("%.0f", DR1IFDCD))

dir.create("data/interim", showWarnings = FALSE, recursive = TRUE)

codes <- f %>% distinct(code, DESCRIPTION) %>% arrange(code)
write_csv(codes, "data/interim/nhanes_codes_all.csv")

ana <- readRDS("data/interim/nhanes_analytic.rds") %>% select(SEQN, WTDRD1_POOL)
intake <- f %>%
  semi_join(ana, by = "SEQN") %>%
  left_join(ana, by = "SEQN") %>%
  group_by(code) %>%
  summarise(grams_wt = sum(DR1IGRMS * WTDRD1_POOL, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(grams_wt))
write_csv(intake, "data/interim/nhanes_code_intake.csv")

cat(sprintf("Wrote data/interim/nhanes_codes_all.csv   (%s code x description rows, %s distinct codes)\n",
            format(nrow(codes), big.mark = ","), format(length(unique(codes$code)), big.mark = ",")))
cat(sprintf("Wrote data/interim/nhanes_code_intake.csv (%s codes with weighted intake)\n",
            format(nrow(intake), big.mark = ",")))
