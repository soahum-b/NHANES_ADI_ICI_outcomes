#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 02_harmonize_score.R
# Harmonize NHANES 2011-2018 (cycles G,H,I,J), compute HEI-2020 (primary) and the
# published NRF9.3 Nutrient-Rich Foods index (sensitivity), build SES measures +
# covariates, and write the analytic dataset with full survey-design columns.
#
# HEI-2020 note: for ages 2+, HEI-2020 components and scoring standards are identical
# to HEI-2015 (USDA/NCI 2023 review; the rename reflects alignment with the 2020-2025
# Dietary Guidelines). The scoring code below is therefore unchanged; only the label
# is updated. See docs/METHODS.md.
#
# NRF9.3 (Fulgoni, Keast & Drewnowski, J Nutr 2009;139:1549): 9 nutrients to encourage
# (protein, fiber, vit A, C, E, calcium, iron, potassium, magnesium; each as %DV per
# 100 kcal, capped at 100%) minus 3 to limit (saturated fat, added sugars, sodium; as
# %MRV per 100 kcal). Replaces the earlier bespoke nutrient-density index.
#
# Diet score inputs:
#   - USDA FPED per-person day-1 (.sas7bdat): 37 food-pattern components + design vars
#   - NHANES DR1TOT (.xpt): energy, saturated fat, sodium, MUFA, PUFA, fiber
#   - NHANES DEMO (.xpt): SES + demographics ; BMX: BMI ; SMQ: smoking
#
# Pooled weights: 4 two-year cycles -> day-1 dietary weight = WTDRD1 / 4  (NHANES rule)
# Analytic restriction: adults 20+, reliable day-1 recall (DR1DRSTZ==1), WTDRD1>0
# ------------------------------------------------------------------------------
suppressMessages({library(haven); library(dplyr); library(tidyr)})

cyc <- tibble::tribble(
  ~suf, ~yy,   ~cycle,
  "G",  "1112","2011-2012",
  "H",  "1314","2013-2014",
  "I",  "1516","2015-2016",
  "J",  "1718","2017-2018")

read_cycle <- function(suf, yy, cycle) {
  demo <- read_xpt(sprintf("data/raw/cdc/DEMO_%s.xpt", suf))
  dr1  <- read_xpt(sprintf("data/raw/cdc/DR1TOT_%s.xpt", suf))
  bmx  <- read_xpt(sprintf("data/raw/cdc/BMX_%s.xpt", suf))
  smq  <- read_xpt(sprintf("data/raw/cdc/SMQ_%s.xpt", suf))
  fped <- read_sas(sprintf("data/raw/fped/fped_dr1tot_%s.sas7bdat", yy))

  demo2 <- demo %>% transmute(
    SEQN, SDMVSTRA, SDMVPSU, WTINT2YR, WTMEC2YR,
    RIDAGEYR, RIAGENDR, RIDRETH1,
    DMDEDUC2,                    # adult education (20+)
    INDFMPIR,                    # family poverty-income ratio
    INDFMIN2)                    # family income category (harmonized across G-J)
  dr1n <- dr1 %>% transmute(
    SEQN, WTDRD1, DR1DRSTZ, DR1TKCAL,
    DR1TSFAT, DR1TSODI, DR1TMFAT, DR1TPFAT, DR1TFIBE,
    # NRF9.3 'encourage' nutrients (protein/fiber already/here + 7 micros)
    DR1TPROT, DR1TVARA, DR1TVC, DR1TATOC, DR1TCALC, DR1TIRON, DR1TPOTA, DR1TMAGN)
  bmxn <- bmx %>% transmute(SEQN, BMXBMI)
  smqn <- smq %>% transmute(SEQN, SMQ020)
  fpedn <- fped %>% select(SEQN, starts_with("DR1T_"))

  demo2 %>%
    inner_join(dr1n, by="SEQN") %>%
    left_join(fpedn, by="SEQN") %>%
    left_join(bmxn,  by="SEQN") %>%
    left_join(smqn,  by="SEQN") %>%
    mutate(cycle=cycle, suffix=suf)
}

raw <- purrr::pmap_dfr(cyc, read_cycle)
cat("Stacked rows (all ages, day-1 recall present):", nrow(raw), "\n")

# ---- HEI-2020 (density method, per-person day-1; identical scoring to HEI-2015) ----
# Legume dual-count (standard simplified per-person method): legumes added to
# vegetables (cup eq) AND protein (oz eq).
kcal_min <- 1  # avoid divide-by-zero; real filter is DR1TKCAL>0 below
hei <- raw %>% mutate(
  KCAL = DR1TKCAL,
  # derived component intakes
  FWHOLEFRT   = DR1T_F_TOTAL - DR1T_F_JUICE,
  VTOTAL_LEG  = DR1T_V_TOTAL + DR1T_V_LEGUMES,
  GREEN_BEAN  = DR1T_V_DRKGR + DR1T_V_LEGUMES,
  PROT_TOTAL  = DR1T_PF_MPS_TOTAL + DR1T_PF_EGGS + DR1T_PF_NUTSDS + DR1T_PF_SOY + DR1T_PF_LEGUMES,
  SEAPLANT    = DR1T_PF_SEAFD_HI + DR1T_PF_SEAFD_LOW + DR1T_PF_SOY + DR1T_PF_NUTSDS + DR1T_PF_LEGUMES,
  MONOPOLY    = DR1TMFAT + DR1TPFAT,
  # densities per 1000 kcal
  d_totfruit = DR1T_F_TOTAL / (KCAL/1000),
  d_whfruit  = FWHOLEFRT    / (KCAL/1000),
  d_totveg   = VTOTAL_LEG   / (KCAL/1000),
  d_grebean  = GREEN_BEAN   / (KCAL/1000),
  d_whgrain  = DR1T_G_WHOLE / (KCAL/1000),
  d_dairy    = DR1T_D_TOTAL / (KCAL/1000),
  d_totprot  = PROT_TOTAL   / (KCAL/1000),
  d_seaplant = SEAPLANT     / (KCAL/1000),
  fa_ratio   = ifelse(DR1TSFAT>0, MONOPOLY/DR1TSFAT, NA_real_),
  d_refgrain = DR1T_G_REFINED / (KCAL/1000),
  d_sodium   = (DR1TSODI/1000) / (KCAL/1000),                 # g per 1000 kcal
  pct_addsug = 100 * (DR1T_ADD_SUGARS*16) / KCAL,             # % energy (16 kcal/tsp eq)
  pct_satfat = 100 * (DR1TSFAT*9) / KCAL)                     # % energy

# scoring helpers: adequacy (higher better), moderation (lower better)
sc_adeq <- function(x, maxpt, std_max) pmin(maxpt, pmax(0, maxpt * x / std_max))
sc_mod  <- function(x, maxpt, std_min, std_max)  # full pts at <=std_min, 0 at >=std_max
  pmin(maxpt, pmax(0, maxpt * (std_max - x)/(std_max - std_min)))
sc_fa   <- function(r, lo=1.2, hi=2.5) ifelse(is.na(r), 0, pmin(10, pmax(0, 10*(r-lo)/(hi-lo))))

hei <- hei %>% mutate(
  HEI_totfruit = sc_adeq(d_totfruit, 5, 0.8),
  HEI_whfruit  = sc_adeq(d_whfruit,  5, 0.4),
  HEI_totveg   = sc_adeq(d_totveg,   5, 1.1),
  HEI_grebean  = sc_adeq(d_grebean,  5, 0.2),
  HEI_whgrain  = sc_adeq(d_whgrain, 10, 1.5),
  HEI_dairy    = sc_adeq(d_dairy,   10, 1.3),
  HEI_totprot  = sc_adeq(d_totprot,  5, 2.5),
  HEI_seaplant = sc_adeq(d_seaplant, 5, 0.8),
  HEI_fatacid  = sc_fa(fa_ratio),
  HEI_refgrain = sc_mod(d_refgrain, 10, 1.8, 4.3),
  HEI_sodium   = sc_mod(d_sodium,   10, 1.1, 2.0),
  HEI_addsug   = sc_mod(pct_addsug, 10, 6.5, 26),
  HEI_satfat   = sc_mod(pct_satfat, 10, 8,   16))

hei_comp <- c("HEI_totfruit","HEI_whfruit","HEI_totveg","HEI_grebean","HEI_whgrain",
              "HEI_dairy","HEI_totprot","HEI_seaplant","HEI_fatacid",
              "HEI_refgrain","HEI_sodium","HEI_addsug","HEI_satfat")
hei <- hei %>% mutate(HEI2020_TOTAL = rowSums(across(all_of(hei_comp))))

# ---- SES measures (reported separately) ------------------------------------
# Education (DMDEDUC2): 1<9th,2 9-11,3 HS/GED,4 some college/AA,5 college grad+
# Recode 7(refused)/9(don't know) to NA. Ordered factor.
# Income (INDFMIN2): harmonized family income category G-J (1-15 usable, 77/99 NA,
#   12='20k+',13='<20k' are non-ordinal -> NA for ordinal use). PIR is the continuous SES.
hei <- hei %>% mutate(
  educ = case_when(DMDEDUC2 %in% 1:5 ~ DMDEDUC2, TRUE ~ NA_real_),
  educ_f = factor(educ, levels=1:5,
                  labels=c("<9th","9-11th","HS/GED","SomeCollege","CollegeGrad")),
  # INDFMIN2 ordinal income: codes 1-10 = <5k..65-75k; 14 = 75-100k; 15 = 100k+.
  # Codes 12 ('20k+') & 13 ('<20k') are non-ordinal overlaps -> NA; 77/99 -> NA.
  income_cat = case_when(
    INDFMIN2 %in% 1:10 ~ INDFMIN2,
    INDFMIN2 == 14 ~ 11,   # $75,000-$99,999
    INDFMIN2 == 15 ~ 12,   # $100,000+
    TRUE ~ NA_real_),
  PIR = INDFMPIR,
  age = RIDAGEYR,
  sex = factor(RIAGENDR, levels=c(1,2), labels=c("Male","Female")),
  raceeth = factor(RIDRETH1, levels=1:5,
             labels=c("MexAmerican","OtherHispanic","NHWhite","NHBlack","OtherRace")),
  bmi = BMXBMI,
  eversmoke = case_when(SMQ020==1 ~ 1, SMQ020==2 ~ 0, TRUE ~ NA_real_))

# ---- NRF9.3 Nutrient-Rich Foods index (published sensitivity index) ---------
# Fulgoni, Keast & Drewnowski, J Nutr 2009;139:1549. 9 nutrients to encourage as %DV per
# 100 kcal (each capped at 100%) minus 3 to limit as %MRV per 100 kcal.
# Reference amounts (standard NRF9.3 usage):
DV  <- c(prot=50, fiber=28, vitA=900, vitC=90, vitE=15,   # protein g, fiber g, A ugRAE, C mg, E mg
         calc=1300, iron=18, pota=4700, magn=420)          # Ca mg, Fe mg, K mg, Mg mg
MRV <- c(satfat=20, addsug=50, sodium=2300)                # g, g, mg
pct_cap <- function(amt, kcal, dv) pmin(100, 100 * (amt*(100/kcal)) / dv)  # %DV per 100 kcal, capped
pct_lim <- function(amt, kcal, mrv) 100 * (amt*(100/kcal)) / mrv           # %MRV per 100 kcal
nrf <- hei %>% mutate(
  addsug_g = DR1T_ADD_SUGARS * 4.2,                        # tsp-eq -> grams (~4.2 g/tsp)
  nrf_encourage =
    pct_cap(DR1TPROT, KCAL, DV["prot"])  + pct_cap(DR1TFIBE, KCAL, DV["fiber"]) +
    pct_cap(DR1TVARA, KCAL, DV["vitA"])  + pct_cap(DR1TVC,   KCAL, DV["vitC"])  +
    pct_cap(DR1TATOC, KCAL, DV["vitE"])  + pct_cap(DR1TCALC, KCAL, DV["calc"])  +
    pct_cap(DR1TIRON, KCAL, DV["iron"])  + pct_cap(DR1TPOTA, KCAL, DV["pota"])  +
    pct_cap(DR1TMAGN, KCAL, DV["magn"]),
  nrf_limit =
    pct_lim(DR1TSFAT,  KCAL, MRV["satfat"]) +
    pct_lim(addsug_g,  KCAL, MRV["addsug"]) +
    pct_lim(DR1TSODI,  KCAL, MRV["sodium"]),
  NRF93 = nrf_encourage - nrf_limit)

# ---- Analytic restriction --------------------------------------------------
ana <- nrf %>% filter(RIDAGEYR>=20, DR1DRSTZ==1, WTDRD1>0, DR1TKCAL>0)
cat("Analytic sample (adults 20+, reliable recall, wt>0):", nrow(ana), "\n")

# pooled 4-cycle day-1 weight
ana <- ana %>% mutate(WTDRD1_POOL = WTDRD1/4)

# ---- Save analytic dataset --------------------------------------------------
keep <- ana %>% transmute(
  SEQN, cycle, suffix,
  SDMVSTRA, SDMVPSU, WTDRD1, WTDRD1_POOL, WTINT2YR, WTMEC2YR,
  age, sex, raceeth, bmi, eversmoke,
  PIR, educ, educ_f, income_cat,
  HEI2020_TOTAL, !!!rlang::syms(hei_comp),
  NRF93, nrf_encourage, nrf_limit,
  KCAL=DR1TKCAL)

dir.create("data/interim", showWarnings=FALSE, recursive=TRUE)

# ---- Merge oPDI/hPDI/uPDI if already computed (scripts/05_hpdi.R) ------------
# The PDI family needs the individual-food files; computed separately and merged
# here so the analytic dataset carries all diet scores. Run 05_hpdi.R first;
# if its output is absent this merge is skipped (HEI-2020 + NRF9.3 still written).
# oPDI is tolerated as absent for backward compatibility with pre-oPDI outputs.
if (file.exists("data/interim/hpdi_scores.rds")) {
  hp <- readRDS("data/interim/hpdi_scores.rds")
  pdi_cols <- intersect(c("oPDI","hPDI","uPDI"), names(hp))
  hp <- hp %>% select(SEQN, all_of(pdi_cols))
  keep <- keep %>% left_join(hp, by="SEQN")
  cat("Merged PDI (", paste(pdi_cols, collapse=", "), "): non-missing",
      sum(!is.na(keep$hPDI)), "of", nrow(keep), "\n")
  if (!"oPDI" %in% pdi_cols)
    cat("  NOTE: oPDI absent from hpdi_scores.rds — re-run scripts/05_hpdi.R.\n")
}

write.csv(keep, "data/interim/nhanes_analytic.csv", row.names=FALSE)
saveRDS(keep, "data/interim/nhanes_analytic.rds")

cat("\n=== HEI-2020 total (unweighted summary) ===\n"); print(summary(keep$HEI2020_TOTAL))
cat("\n=== NRF9.3 (unweighted summary) ===\n"); print(summary(keep$NRF93))
cat("\n=== per-cycle analytic n ===\n"); print(table(keep$cycle))
cat("\n=== PIR non-missing:", sum(!is.na(keep$PIR)),
    " educ non-missing:", sum(!is.na(keep$educ)),
    " income non-missing:", sum(!is.na(keep$income_cat)), "\n")
cat("Saved data/interim/nhanes_analytic.{csv,rds}\n")
