# Bolte 17-group plant-based dietary index in NHANES 2011–2018

Replication code for the healthful, unhealthful and overall plant-based dietary indices
(hPDI, uPDI, oPDI) built on the **17 food groups of Bolte et al. 2023**
(*JAMA Oncol* 9(5):705–709, doi:10.1001/jamaoncol.2022.7753), which adapt the 18 groups of
Satija et al. 2017 by collapsing nuts and legumes into one group. Indices are computed on
NHANES 2011–2018 day-1 dietary recalls for adults aged 20 and over and related to
socioeconomic position (income-to-poverty ratio, education, family income).

All estimation is survey-weighted and design-based. R produces every statistic; Python
only draws figures and renders tables. Neither side calls the other, so each can be
checked independently.

## Contents

1. [Requirements](#requirements)
2. [Quick start](#quick-start)
3. [Full pipeline](#full-pipeline)
4. [Repository layout and data provenance](#repository-layout-and-data-provenance)
5. [Methods notes](#methods-notes)
6. [Citation](#citation)
7. [Revision notes](#revision-notes)

## Requirements

- **R 4.5.3** with `survey` 4.5, `emmeans` 2.0.4, `haven`, `dplyr`, `tidyr`, `tibble`,
  `purrr`, `readr`, `broom` (versions in `docs/environment_r.txt`).
- **Python 3.11** with `matplotlib`, `pandas`, `numpy` (versions in
  `docs/environment_python.txt`; the other packages listed there belong to the parent
  project and are not imported by any script in this repository).
- `curl` and `bsdtar` (`libarchive-tools` on Debian/Ubuntu; included in macOS) for
  `01_download.sh`, which unpacks the USDA self-extracting `.exe` archives.
- Several hundred MB of disk for the raw downloads.

Both environment files pin direct dependencies only. There is no lock file, so transitive
versions are not fixed.

## Quick start

Everything is run **from the repository root**, never from `scripts/`.

The repository ships the finished food-code classification
(`data/interim/dhq_food_classification_public.csv`), so the full set of 17-group results
can be reproduced without the restricted inputs described below:

```bash
bash    scripts/01_download.sh                     # CDC XPT files + USDA FPED (network)
Rscript scripts/02_harmonize_score.R               # analytic sample, HEI-2020, NRF9.3
Rscript scripts/05_hpdi.R                          # 18-group indices (Satija)
Rscript scripts/02_harmonize_score.R               # second pass: attaches the 18-group indices
Rscript scripts/36_dhq_scores_panels.R             # 17-group indices + all figure inputs
python  scripts/37_dhq_pca_loadings_table.py       # PCA loadings tables
python  scripts/fig_dhq3_pdi_vs_ses.py             # figures (see below for the rest)
```

`36_` prints `classification source: data/interim/dhq_food_classification_public.csv`
when it is using the shipped file.

## Full pipeline

### Scripts, inputs and outputs

| Script | Reads | Writes |
|---|---|---|
| `01_download.sh` | network | `data/raw/cdc/*.xpt`, `data/raw/fped/*.sas7bdat` |
| `02_harmonize_score.R` | raw CDC + FPED per-person files; `data/interim/hpdi_scores.rds` if present | `data/interim/nhanes_analytic.{csv,rds}` |
| `05_hpdi.R` | FPED per-food files; `nhanes_analytic.rds` | `data/interim/hpdi_scores.{csv,rds}` (oPDI, hPDI, uPDI, 18 groups) |
| `34b_nhanes_food_codes.R` | FPED per-food files; `nhanes_analytic.rds` | `data/interim/nhanes_codes_all.csv`, `nhanes_code_intake.csv` |
| `35_dhq_bolte_reclassify.py` | the two files above; **`data_restricted/`** | `results/bolte/dhq_*.csv` (full classification, review queue, summaries); `data/interim/dhq_food_classification_public.csv` |
| `36_dhq_scores_panels.R` | FPED per-food files; `nhanes_analytic.rds`; a classification file | `data/interim/dhq_scores.rds`; `results/dhq/{coverage,pd_main,trend_tests,food_odds,pca_loadings,pca_variance,pca_scores}.csv` |
| `37_dhq_pca_loadings_table.py` | `results/dhq/pca_loadings.csv` | `results/dhq/bolte17_pca_loadings_table{,_full}.csv` |
| `13_table1.R` | `nhanes_analytic.rds` (+ `dhq_scores.rds` with `BOLTE17=1`) | `results/[bolte17_]table1_{quartile,tertile}.csv` |
| `38_table1_quartiles4.R` | same | `results/[bolte17_][median_]table1_4quartiles.csv` |
| `39_table1_missing.R` | same | `results/[bolte17_]table1_missing_supplement.csv` |
| `fig_dhq*.py` | `results/dhq/*.csv` | `figures/dhq/*.png` with vector twins in `figures/dhq/{pdf,svg}/` |

### Order, including the classification step

```bash
bash    scripts/01_download.sh
Rscript scripts/02_harmonize_score.R
Rscript scripts/05_hpdi.R
Rscript scripts/02_harmonize_score.R               # second pass, see below
Rscript scripts/34b_nhanes_food_codes.R            # only needed before 35_
python  scripts/35_dhq_bolte_reclassify.py         # needs data_restricted/
Rscript scripts/36_dhq_scores_panels.R
python  scripts/37_dhq_pca_loadings_table.py
```

**`02_` runs twice by design.** `05_hpdi.R` reads the analytic dataset that `02_` writes,
and `02_` merges the 18-group index columns into that dataset only if `05_` has already
produced them. The first pass builds the dataset without the indices; the second attaches
them. Skipping the second pass leaves `13_`, `38_` and `39_` without `oPDI`/`hPDI`/`uPDI`
unless `BOLTE17=1` is set.

**Which classification `36_` uses.** In order of precedence:

1. `$DHQ_CLASSIFICATION`, if set. The file must exist and carry `code`, `category` and
   `tier` columns; the script stops rather than falling back to another file.
2. `results/bolte/dhq_food_classification.csv`, the full output of `35_`, if present.
3. `data/interim/dhq_food_classification_public.csv`, the committed four-column cut.

The script prints which one it used.

### Tables

The table scripts report the **18-group** indices by default. Set `BOLTE17=1` for the
17-group indices from `36_`; everything else (design, groups, tests, formatting) is
identical, so the two variants are directly comparable.

```bash
BOLTE17=1           Rscript scripts/13_table1.R             # cohort table, Q1 vs Q4 and T1 vs T3
BOLTE17=1           Rscript scripts/38_table1_quartiles4.R  # all four PIR quartiles, mean (SD)
BOLTE17=1 MEDIAN=1  Rscript scripts/38_table1_quartiles4.R  # same, median (IQR)
BOLTE17=1           Rscript scripts/39_table1_missing.R     # missingness supplement
```

`MEDIAN=1` switches the continuous rows of `38_` to survey-weighted median (IQR) and
replaces the omnibus Wald test with a design-based Kruskal–Wallis test (`svyranktest`),
so the test matches the statistic displayed. P-for-trend remains the linear term of the
quartile rank, which is mean-based, and is labelled as such in the output.

### Figures

```bash
python scripts/fig_dhq3_pdi_vs_ses.py     # hPDI/uPDI vs PIR, education, income (6 panels)
python scripts/fig_dhq4_food_odds.py      # odds of high consumption per +1 PIR, by food group
python scripts/fig_dhq5_pca.py            # PCA biplots, poorest vs richest PIR quartile
python scripts/fig_dhq_split_panels.py    # the 6 panels as standalone figures
```

The `fig_dhq*` scripts are thin wrappers that set the input directory, output path and
grouping label on the `fig_bolte*` modules, so the 16- and 17-group versions of each
figure are drawn by one code path. `figsave.py` writes PNG (300 dpi), PDF and SVG together,
with TrueType font embedding. `fig8_pca_biplot.py` is present only because
`fig_bolte5_pca.py` imports its `_biplot` helper and colour constants; its own `main()`
renders 18-group figures from files this repository does not produce. Do not run it
directly.

## Repository layout and data provenance

| Path | Committed | Contents |
|---|---|---|
| `scripts/` | yes | all code |
| `docs/` | yes | environment listings |
| `data/interim/dhq_food_classification_public.csv` | yes | NHANES food code → Bolte category, four columns (see below) |
| `data_restricted/raw/dhq/DHQ_Dictionary_V2.csv` | no | annotated DHQ item list |
| `data_restricted/interim/dhq_manual_overrides.csv` | no | 351 hand adjudications with notes |
| `data/raw/cdc/`, `data/raw/fped/` | no | public downloads, re-fetched by `01_download.sh` |
| `data/interim/` (other files), `results/`, `figures/` | no | generated |

`.gitignore` enforces this, including excluding `data_restricted/` as a whole.

### The restricted inputs

`35_dhq_bolte_reclassify.py` is the only script that reads either restricted file, and it
refuses to run without both. The 351 manual overrides are part of the method; a run that
silently skipped them would produce a different classification from the committed one.
Collaborators who need to re-run `35_` obtain `data_restricted/` out of band.

**`DHQ_Dictionary_V2.csv`** is 453 items from NCI's Diet History Questionnaire carrying
five hand-built classification columns: `NOVA`, `FoodGroup_Bolte`,
`Plant_Based_Index_Bolte`, `aMED_Groups` and `Plant_Based_Index`. This pipeline uses
`Plant_Based_Index_Bolte` only. The DHQ itself is public, and NCI distributes the DHQ III
nutrient and food-group database (<https://epi.grants.cancer.gov/dhq3/>), but NCI's food
groups follow the USDA pyramid-servings scheme; the five classification columns are
research work product and are not that database. The category definitions they encode are
published separately (Bolte et al. 2023, eTable 2; Satija et al. 2017), so an equivalent
mapping can be rebuilt from those sources.

**`dhq_manual_overrides.csv`** holds code, category, adjudicator and a free-text note per
override. Several notes quote DHQ item descriptions verbatim, which is why the file is
restricted rather than committed.

### The committed classification

`data/interim/dhq_food_classification_public.csv` is written by `35_` and carries four
columns for each of the 7,418 NHANES food codes:

| Column | Source |
|---|---|
| `code` | 8-digit FNDDS food code (CDC/USDA, public) |
| `DESCRIPTION` | NHANES per-food description from the FPED file (public) |
| `category` | one of the 22 dictionary category labels, or `UNMATCHED` |
| `tier` | `accept`, `manual` or `review`; this project's own label |

The full classification in `results/bolte/` additionally carries the matcher diagnostics,
including `matched_dhq_item`, which quotes the DHQ item each code was matched to on 5,708
rows. Those columns are stripped from the committed file, so no DHQ row is reproduced in
it. Plain food names shared by both lists ("Applesauce", "Cream cheese") do occur in the
NHANES descriptions; that is shared vocabulary, not copied content, since NCI builds the
DHQ database from NHANES recall data.

The code-to-category mapping is nonetheless *derived from* the DHQ annotations and remains
a derivative work. Stripping columns addresses verbatim redistribution, not attribution;
cite the annotations' authors regardless.

## Methods notes

Points where a Methods paragraph written from memory, or from the 18-group scripts, would
be wrong for the 17-group index.

- **Group membership comes from the DHQ dictionary classification, not from FPED's
  food-pattern components.** The per-food FPED file supplies the food code and gram
  weight; every food code belongs to exactly one category, so a food counts once. The
  FPED-component build is a different, 16-group index that cannot construct
  *miscellaneous animal-based foods*, because FPED disaggregates mixed dishes.
- **Food codes are assigned by automated matching, not by hand.** IDF-weighted cosine
  similarity on the food's leading description clause, with the FNDDS code hierarchy as a
  plausibility check, plus rules for the whole/refined grain and whole fruit/juice
  boundaries, which identity-only matching would otherwise discard. Of 7,418 codes, 2,617
  are auto-accepted, 351 manually adjudicated and 4,450 remain in a review tier, covering
  73.0%, 1.8% and 25.2% of consumed grams; a further 4.6% of grams is unclassified. Report
  this coverage; do not describe the assignment as manually verified.
- **Five dictionary categories are excluded** from the index: water, alcoholic beverages,
  margarine, artificial sweeteners and protein supplements.
- **Per-group scoring differs between `36_` and `05_`.** Both give a 1–5 score per group.
  In `36_` non-consumers score 1 and consumers are cut at the weighted 25th, 50th and 75th
  percentiles among consumers into scores 2–5. In `05_` non-consumers score 1 and
  consumers are cut into weighted quintiles, so the lowest consumer quintile can also
  score 1. The 17-group totals run 17–85 and the 18-group totals 18–90.
- **"High consumer" also differs.** In `36_` it is the weighted 2/3 quantile among
  consumers of the group, coded `>=`. The 16- and 18-group scripts use the weighted 2/3
  quantile over all adults with a `>0` fallback, because the among-consumers quantile is
  unstable for rare foods.
- **No energy adjustment** is applied to the indices or the models. Someone who eats more
  of everything ranks higher in every group.
- **Education is an unordered five-level factor** with indicator contrasts against
  `<9th grade`, deliberately not an ordered or linear term, because the association is
  non-monotonic. Its p-value is an omnibus joint Wald test on 4 df: it supports "the levels
  differ" and nothing about any specific pair. Per-level contrasts are not computed.
- **PIR is top-coded at 5.0.** 17.6% of sampled adults (about a quarter by survey weight)
  sit at the ceiling, so the weighted 75th percentile is 5.0 and the top quartile is
  defined by the ceiling, not by a quantile cut. A naive `quantile()` split yields an empty
  top group; the scripts use a top-code-aware cut.
- **Item-missing PIR (about 9%) is handled as a subpopulation.** The survey design is built
  once on the full analytic sample and then `subset()`; rows are never dropped before
  `svydesign`, so the variance estimator keeps the full stratum/PSU structure.
- `options(survey.lonely.psu = "adjust")` is set throughout, because repeated subsetting
  can leave a stratum with a single PSU.

## Citation

- Bolte LA, et al. *JAMA Oncol.* 2023;9(5):705–709. doi:10.1001/jamaoncol.2022.7753
- Satija A, et al. *J Am Coll Cardiol.* 2017;70(4):411–422. doi:10.1016/j.jacc.2017.05.047
- Lumley T. *J Stat Softw.* 2004;9(8):1–19. doi:10.18637/jss.v009.i08
- R Core Team (2026). *R: A Language and Environment for Statistical Computing.*

## Revision notes

### 2026-10-30

The restricted-tier split (DHQ dictionary and overrides moved to `data_restricted/`, a
four-column public classification committed) was introduced in the previous revision
without updating the scripts. Two defects, both confirmed by running the code:

- `35_` looked for the overrides at their old path, and because the lookup was guarded by
  `os.path.exists()` it ran without them and reported nothing: zero codes in the `manual`
  tier instead of 351.
- `36_` read only `results/bolte/dhq_food_classification.csv`, so the committed public
  file was never used and the public tier did not run.

Fixed in this revision:

- `35_` reads both restricted inputs from `data_restricted/` and exits if either is
  missing; it also writes the public four-column file, so that file is a script output.
- `36_` resolves its classification input in the order described above and stops if
  `$DHQ_CLASSIFICATION` is set to a missing path.
- `37_` carries an inlined copy of `build_pca_loadings_table()`, previously imported from
  a script outside this repository.
- `.gitignore` no longer whitelists the overrides under `data/interim/`.
- Header comments in `05_` and `36_` corrected: `05_` reads the USDA FPED per-food file,
  not CDC's `DR1IFF_*.xpt`, and `36_` is not "identical in form" to `05_` (see Methods
  notes).

Known and left as is: `01_download.sh` fetches two CDC files no script here reads,
`DR2TOT_*.xpt` (kept for the usual-intake sensitivity analysis in the parent project) and
`DR1IFF_*.xpt`. The food-level inputs the scripts use are the USDA `fped_dr1iff_*.sas7bdat`
files from the download script's third block; the shared name stem is coincidental.
