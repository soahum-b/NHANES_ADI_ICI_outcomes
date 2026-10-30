# Bolte 17-group plant-based dietary index in NHANES 2011–2018

Replication code for the healthful / unhealthful / overall plant-based dietary indices
(hPDI, uPDI, oPDI) built on the **17 food groups of Bolte et al. 2023**
(*JAMA Oncol* 9(5):705–709, doi:10.1001/jamaoncol.2022.7753), which adapt the Satija
et al. 2017 categories by collapsing nuts and legumes into one group.

All estimation is survey-weighted and design-based. **R produces every statistic;
Python only draws figures and renders tables** — neither crosses over, so the two sides
can be checked independently.

---

## Repository layout — what is committed and what is not

The bundle is split so that a public push carries no third-party item text.

| Path | Committed? | Why |
|---|---|---|
| `scripts/` | yes | all code |
| `data/interim/dhq_food_classification_public.csv` | **yes** | derived NHANES-code→category mapping, four columns, no DHQ item text (checked — see below) |
| `data_restricted/raw/dhq/DHQ_Dictionary_V2.csv` | **no** | the annotated DHQ item list |
| `data_restricted/interim/dhq_manual_overrides.csv` | **no** | 351 hand adjudications; 5 of the notes quote DHQ item text verbatim |
| `data/raw/cdc/`, `data/raw/fped/` | no | several hundred MB, re-fetchable with `01_download.sh` |
| `results/`, `figures/` | no | generated |

`.gitignore` enforces all of this, including excluding `data_restricted/` wholesale.

**The split is exact, not approximate:** `35_dhq_bolte_reclassify.py` is the only script
that reads either restricted file, and it is also the only script that needs them. It
refuses to run without them (the 351 overrides are part of the method, so a silent skip
would regenerate a different classification). `36_dhq_scores_panels.R` reads the full
`results/bolte/dhq_food_classification.csv` when `35_` has been run and otherwise falls
back to the committed `data/interim/dhq_food_classification_public.csv`, which `35_`
itself writes. So the committed tier runs end to end from `36_` onward, skipping `34b_`
and `35_`.

### Is there a public source for the DHQ dictionary?

**Public.** NCI's Diet History Questionnaire is a freely available FFQ for adults, and NCI
distributes the DHQ III nutrient and food-group database as a downloadable CSV built from
NHANES 24-hour recall data (<https://epi.grants.cancer.gov/dhq3/>).

**Not public.** `DHQ_Dictionary_V2.csv` is not that database. It is 453 DHQ items carrying
five hand-built classification columns — `NOVA`, `FoodGroup_Bolte`,
`Plant_Based_Index_Bolte`, `aMED_Groups` and `Plant_Based_Index`. NCI publishes food groups
on the USDA pyramid-servings scheme, not these. The mappings are research work product;
this pipeline uses one of the five columns.

**The category definitions are separately published:** Bolte et al.'s eTable 2 lists the 17
groups and their constituent items, and Satija et al. 2017 gives the parent 18-group
scheme. An equivalent mapping can be rebuilt from those.

### What was checked, and one correction

An earlier version of this README claimed the derived classification "contains no DHQ
content". **That was wrong** — the full `dhq_food_classification.csv` carries verbatim DHQ
item descriptions in `matched_dhq_item` on 5,708 of 7,418 rows. The committed file is a
four-column cut (`code`, `DESCRIPTION`, `category`, `tier`) with all seven
matcher-diagnostic columns removed, including both DHQ-bearing ones.

Provenance of what remains: `code` is an FNDDS 8-digit food code (public, CDC/USDA);
`DESCRIPTION` is NHANES per-food FPED text (public); `category` is a Bolte group label from
their published eTable 2; `tier` is this project's own matcher label.

Eleven distinctive DHQ strings do appear in the committed file, and all of them were
checked against the raw NHANES food descriptions: each occurs **inside a genuine NHANES
FPED description** (e.g. "Macaroni and cheese", "Cream cheese,"). The overlap is expected
rather than suspicious — NCI builds the DHQ database from NHANES recall data, so the two
food lists share vocabulary and formatting. By contrast, "Avocado/ guacamole," and
"Catsup or ketchup" appear nowhere in NHANES text, and those are in the **restricted**
override notes, not the committed file.

**One caveat the column surgery does not remove:** the code-to-category mapping is
*derived from* the DHQ annotations, so it remains a derivative work. Stripping columns
addresses verbatim redistribution, not attribution — cite whoever built the annotations
regardless.

## Run order

Run every script **from the project root**, not from `scripts/`.

```bash
bash   scripts/01_download.sh              # CDC XPT + USDA FPED (needs network + bsdtar/libarchive)
Rscript scripts/02_harmonize_score.R       # analytic sample; HEI-2020 + NRF9.3
Rscript scripts/05_hpdi.R                  # 18-group plant indices
Rscript scripts/02_harmonize_score.R       # RE-RUN: see the note below
Rscript scripts/34b_nhanes_food_codes.R    # food codes + weighted intake per code
python  scripts/35_dhq_bolte_reclassify.py # assign codes to the 17 Bolte groups (needs data_restricted/)
Rscript scripts/36_dhq_scores_panels.R     # 17-group scores + all figure inputs
python  scripts/37_dhq_pca_loadings_table.py
```

### `02_` is run twice, and that is not a mistake

`02_harmonize_score.R` merges the plant-index columns **if** `05_hpdi.R` has already
written them, while `05_hpdi.R` reads the analytic dataset that `02_` produces. The order
is therefore **02 → 05 → 02**. The first pass builds the dataset without the indices; the
second attaches them.

### Running the committed tier only (no DHQ dictionary)

`35_` writes the committed cut, so on a clone that has it you skip both `34b_` and `35_`:

```bash
bash    scripts/01_download.sh
Rscript scripts/02_harmonize_score.R && Rscript scripts/05_hpdi.R && Rscript scripts/02_harmonize_score.R
Rscript scripts/36_dhq_scores_panels.R      # prints "classification source: data/interim/..."
python  scripts/37_dhq_pca_loadings_table.py
```

### Tables — the 17-group variants are environment switches

```bash
BOLTE17=1              Rscript scripts/13_table1.R            # cohort table, extreme groups
BOLTE17=1              Rscript scripts/38_table1_quartiles4.R # all four PIR quartiles, mean (SD)
BOLTE17=1 MEDIAN=1     Rscript scripts/38_table1_quartiles4.R # same, median (IQR)
BOLTE17=1              Rscript scripts/39_table1_missing.R    # missingness supplement
```

Without `BOLTE17=1` these scripts report the **18-group** indices instead.

`36_dhq_scores_panels.R` takes one more switch: `DHQ_CLASSIFICATION=/path/to/file.csv`
points it at a classification held outside the tree (it needs the `code`, `category` and
`tier` columns). When set, the file must exist; the script stops rather than falling back. `MEDIAN=1`
switches the continuous rows to survey-weighted median (IQR) and swaps the omnibus test
from a Wald test to a **design-based Kruskal–Wallis** (`svyranktest`), so the test matches
the statistic being displayed.

### Figures

```bash
python scripts/fig_dhq3_pdi_vs_ses.py    # hPDI/uPDI vs PIR, education, income (6 panels)
python scripts/fig_dhq4_food_odds.py     # food-group odds per +1 PIR
python scripts/fig_dhq5_pca.py           # PCA biplots, poorest vs richest quartile
python scripts/fig_dhq_split_panels.py   # the 6 panels as standalone figures
```

The `fig_dhq*` scripts are thin wrappers that parameterise `fig_bolte*` by input
directory, output path and grouping label, so the 16- and 17-group versions are drawn by
the same code and cannot drift apart. `figsave.py` writes PNG, SVG and PDF together.

`fig8_pca_biplot.py` is included **only as the shared biplot helper** that `fig_dhq5_pca.py` imports (`_biplot` and the colour constants). Importing it needs nothing on disk — verified — but its own `main()` renders the *18-group* figures and reads `results/pca_*.csv`, which this bundle does not produce. Do not run it directly; it is a library here.

---

## Methods notes that are easy to get wrong

- **The 17 groups come from the DHQ dictionary classification, not from FPED's food-pattern
  components.** The per-food FPED file supplies only the food code and the gram weight.
  The FPED-component build is a *different*, 16-group index, which cannot construct
  *miscellaneous animal-based foods* because FPED disaggregates mixed dishes.
- **Food codes are assigned by automated matching**, not by hand: cosine similarity on food
  descriptions, with the FNDDS code hierarchy used as a validity check, plus rules for the
  whole/refined-grain and whole-fruit/juice boundaries (those qualifiers decide the
  healthful/less-healthful split and identity-only matching discards them). Of 7,418 codes,
  2,617 are auto-accepted, **351 manually adjudicated**, and 4,450 sit in a review tier —
  covering 73.0%, 1.8% and 25.2% of consumed grams respectively. A further **4.6% of grams
  is unclassified**. Report this coverage; do not describe the assignment as manually
  verified.
- **Five dictionary categories are excluded** from the index — water, alcoholic beverages,
  margarine, artificial sweetener and protein supplements — not two.
- **No energy adjustment** is applied, to either the index or the models. Someone who eats
  more of everything ranks higher in every group.
- **Education is an UNORDERED five-level factor** with indicator contrasts against
  `<9th grade`, deliberately not an ordered or linear term, because the association is
  non-monotonic. Its p-value is an **omnibus joint Wald test on 4 df** — it licenses "the
  levels differ" but no claim about any specific pair. Per-level education contrasts are
  **not** computed on this scheme.
- **PIR is top-coded at 5.0.** 17.6% of sampled adults (about a quarter by survey weight)
  sit at the ceiling, so the weighted 75th percentile *is* 5.0 and the top quartile is
  defined by the ceiling rather than by a quantile cut. Scripts use a top-code-aware cut;
  a naive `quantile()` split yields an empty top group.
- **"High consumer" is defined differently in the two schemes.** Here (`36_`) it is the
  weighted 2/3 quantile **among consumers** of that group, coded `>=`. The 16-/18-group
  script uses the weighted 2/3 quantile **over all adults** with a `>0` fallback, because
  the among-consumers quantile is unstable for rare foods. A Methods sentence written for
  one is wrong for the other.
- Item-missing PIR (~9%) is handled by **subpopulation (domain) analysis** — the design is
  built once on the full sample and then `subset()`, never row-filtered before
  `svydesign`, so the variance estimator keeps the full stratum/PSU structure.
- `survey.lonely.psu = "adjust"` is set, because repeated subsetting can leave a stratum
  with one PSU.

## Environment

See `docs/environment_r.txt` and `docs/environment_python.txt`. R 4.5.3 with `survey` 4.5
and `emmeans` 2.0.4 does the estimation; CPython 3.11.15 with `matplotlib`, `pandas` and
`numpy` draws the figures. These pin **direct** dependencies only — there is no lock file,
so transitive versions are not fixed.

## Citation

- Bolte LA, et al. *JAMA Oncol.* 2023;9(5):705–709. doi:10.1001/jamaoncol.2022.7753
- Satija A, et al. *J Am Coll Cardiol.* 2017;70(4):411–422. doi:10.1016/j.jacc.2017.05.047
- Lumley T. *J Stat Softw.* 2004;9(8):1–19. doi:10.18637/jss.v009.i08
- R Core Team (2026). *R: A Language and Environment for Statistical Computing.*

## Revision notes

### 2026-10-30: restricted-tier split made consistent with the code

The previous revision moved the DHQ dictionary and the manual overrides into
`data_restricted/` and shipped a four-column public classification, but the scripts
were not updated to match. Two consequences, both verified by running the scripts:

- `35_dhq_bolte_reclassify.py` still looked for the overrides at their old path and,
  because that lookup was guarded by `os.path.exists()`, **ran without them and said
  nothing**: zero codes landed in the `manual` tier instead of 351. Anyone re-running
  `35_` would have regenerated a classification that differs from the committed one.
- `36_dhq_scores_panels.R` only read `results/bolte/dhq_food_classification.csv`, which
  `35_` produces, so the committed public file was never read by any script and the
  claim that the public tier "runs from `36_` onward" did not hold.

What changed:

- `35_` reads both restricted inputs from `data_restricted/` and exits with an error if
  either is missing. The overrides are part of the method, not an optional refinement.
- `35_` now writes `data/interim/dhq_food_classification_public.csv` itself, so the
  committed file is a script output and cannot drift from the full classification.
- `36_` uses the full classification when present and otherwise falls back to the public
  file, printing which it used.
- `37_dhq_pca_loadings_table.py` carries an inlined copy of `build_pca_loadings_table()`
  (previously imported from a script outside this bundle).
- `.gitignore` no longer whitelists the overrides under `data/interim/`.

Not changed, and worth knowing: the per-group rank scoring in `36_` is not the same rule
as `05_hpdi.R`, despite both producing a 1-5 score per group. `05_` cuts consumers into weighted quintiles
(the lowest consumer bin can share score 1 with non-consumers); `36_` cuts consumers
into weighted quartiles scored 2 to 5, with non-consumers at 1. Both are on a 1 to 5
scale per group but they are different rules, so describe each index by its own script.
`01_download.sh` also fetches two CDC files that no script in this bundle reads: `DR2TOT`
(it supports the usual-intake sensitivity analysis, which is outside this replication set)
and the CDC `DR1IFF_*.xpt` individual-foods file. The food-level inputs that `05_`, `34b_`
and `36_` actually read are the USDA FPED files `fped_dr1iff_*.sas7bdat`, fetched in the
third block of the download script; the shared `DR1IFF` name stem is a coincidence of
naming, not a dependency.
