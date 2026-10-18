# mortselect

[![Forecast explorer](https://img.shields.io/badge/Explore-mortality_forecasts-247466?style=flat-square)](https://hg-kaibai.github.io/mortselect/)
[![Release downloads](https://img.shields.io/github/downloads/HG-KaiBai/mortselect/total?style=flat-square&label=release%20downloads)](https://github.com/HG-KaiBai/mortselect/releases)
[![Page views](https://hits.dwyl.com/HG-KaiBai/mortselect.svg?style=flat-square)](https://hits.dwyl.com/HG-KaiBai/mortselect)
[![GitHub stars](https://img.shields.io/github/stars/HG-KaiBai/mortselect?style=flat-square)](https://github.com/HG-KaiBai/mortselect/stargazers)
[![Last commit](https://img.shields.io/github/last-commit/HG-KaiBai/mortselect?style=flat-square)](https://github.com/HG-KaiBai/mortselect/commits/master)
[![License: GPL v3+](https://img.shields.io/badge/license-GPL_%E2%89%A5_3-blue?style=flat-square)](https://www.gnu.org/licenses/gpl-3.0.html)

## Interactive Forecast Explorer

**[Open the Mortality Forecast Explorer →](https://hg-kaibai.github.io/mortselect/)**

Explore the saved paper results by age, population, or year. The explorer has
interactive charts, method-specific filters, a log scale, shareable views,
and CSV/SVG downloads. It runs entirely in the browser on GitHub Pages, so it
does not require an active R server or a running local computer.

The supplied results cover **2011–2025**. CBD forecasts cover ages **55–90**;
the other methods cover **0–100**. These are saved point forecasts, not observed
mortality rates or newly fitted projections. All displayed values are labelled
as **mortality rate × 100**; CSV downloads also include the unscaled rates.

To run the redesigned native Shiny version from a clone of this repository:

```r
install.packages(c("shiny", "ggplot2", "jsonlite"))
shiny::runApp("shiny")
```

See [the website documentation](docs/README.md) for data export and local preview
instructions. The release-download badge counts GitHub release assets, not
R package installs. The page-view badge counts requests to the README badge,
not unique visitors or GitHub's private traffic analytics.

## About the Package

`mortselect` is an R package for mortality forecasting with data-driven
information borrowing. It implements workflows based on three published
methodological papers:

- deletion-substitution-addition (DSA) population selection;
- age-specific band selection (ASB) for single- and multi-population
  mortality forecasting;
- bivariate model-based ensemble (BMBE) forecasts with ACF- and CBD-type base
  learners.

The package works with mortality datasets represented as lists of
`demography::demogdata` objects. The public functions expose age ranges, year
windows, population exclusions, and forecast horizons, so the workflows can be
used beyond the original Human Mortality Database study periods.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("HG-KaiBai/mortselect")
```

Alternatively, with `pak`:

```r
install.packages("pak")
pak::pak("HG-KaiBai/mortselect")
```

On Windows, source-package installation may require Rtools if a dependency
needs compilation.

## Main Workflows

The exported forecasting functions are:

- `forecast_dsa()` fits the DSA population-selection workflow.
- `forecast_asb_single()` fits single-population age-specific band forecasts.
- `forecast_asb_multi()` fits multi-population age-specific band forecasts.
- `forecast_bmbe_acf()` fits a BMBE workflow using ACF-type base learners.
- `forecast_bmbe_cbd()` fits a BMBE workflow using CBD-type base learners.

Supporting data helpers include:

- `prepare_mortselect_data()`
- `validate_demogdata_list()`
- `validate_prepared_mortselect_data()`
- `create_population_age_dictionary()`
- `bind_training_rates()`

## Data Preparation

```r
library(mortselect)

prepared <- prepare_mortselect_data(
  dataset,
  ages = 0:100,
  series = c("female", "male", "total"),
  learning_years = 1980:2004,
  validation_years = 2005:2012,
  test_years = 2013:2020,
  exclude_labels = c("Example population to omit")
)

prepared_dsa <- prepare_mortselect_data(
  dataset,
  ages = 0:100,
  learning_years = 1980:2004,
  validation_years = 2005:2012,
  test_years = 2013:2020,
  exclude_labels = c("Example population to omit"),
  relabel_sex = TRUE
)

prepared_cbd <- prepare_mortselect_data(
  dataset,
  ages = 55:90,
  series = "female",
  learning_years = 1980:2004,
  validation_years = 2005:2012,
  test_years = 2013:2020,
  combine.upper = TRUE
)
```

`prepare_mortselect_data()` validates the requested mortality series, filters
populations with sufficient year coverage, applies optional exclusions, repairs
zero or missing mortality rates, and returns `learning`, `validation`, `test`,
`training`, and `full` datasets.

If a population has an age row with no positive finite mortality rates available
for interpolation, that population is skipped consistently across all prepared
splits and reported in `prepared$metadata$skipped`.

The forecasting functions expect this prepared object and check its format
before model fitting. They do not call `prepare_mortselect_data()` internally.

## Forecasting Examples

```r
dsa_forecast <- forecast_dsa(
  prepared_dsa,
  COG = 10,
  N = 3,
  LOSS = "SSE",
  t = 10
)

single_asb <- forecast_asb_single(
  prepared,
  Series = "female",
  t = 10
)

multi_asb <- forecast_asb_multi(
  prepared,
  Series = "male",
  t = 10
)

acf_ensemble <- forecast_bmbe_acf(
  prepared,
  Series = "female",
  t = 10
)

cbd_ensemble <- forecast_bmbe_cbd(
  prepared_cbd,
  Series = "female",
  t = 10
)
```

Each workflow returns its forecast object. To also write an `.RData` result,
pass `save_file = "result.RData"`.

## Parallel Execution

All forecasting workflows can run independent target populations or bivariate
model pairs in parallel. Use `parallel = TRUE` and optionally set `workers`:

```r
multi_asb <- forecast_asb_multi(
  prepared,
  Series = "female",
  t = 10,
  parallel = TRUE,
  workers = 4
)

cbd_ensemble <- forecast_bmbe_cbd(
  prepared_cbd,
  Series = "female",
  t = 10,
  parallel = TRUE,
  workers = 4
)
```

When `workers` is not supplied, the package uses one fewer than the number of
physical CPU cores, capped at the number of independent tasks.

## Default Study Windows

The defaults reproduce the empirical setup in the source scripts and papers:

- learning years: `1970:1994`
- validation years: `1995:2002`
- testing years: `2003:2010`
- final training years: learning plus validation
- ASB/DSA age range: `0:100`
- CBD age range: `55:90`

For a new mortality dataset, set the year windows explicitly and prefer
`exclude_labels` over numeric exclusions.

## References

Please cite the relevant methodological paper when using a workflow from this
package.

1. Diao, L., Meng, Y., & Weng, C. (2021). A DSA algorithm for mortality
   forecasting. *North American Actuarial Journal, 25*(3), 438-458.
   <https://doi.org/10.1080/10920277.2020.1806884>

2. Meng, Y., Diao, L., & Weng, C. (2026). Mortality prediction via
   age-specific band selection. *Scandinavian Actuarial Journal, 2026*(4),
   443-467. <https://doi.org/10.1080/03461238.2025.2537925>

3. Diao, L., Meng, Y., Weng, C., & Wirjanto, T. (2023). Enhancing mortality
   forecasting through bivariate model-based ensemble. *North American
   Actuarial Journal, 27*(4), 751-770.
   <https://doi.org/10.1080/10920277.2023.2167832>

Method mapping:

- `forecast_dsa()` follows the DSA algorithm paper.
- `forecast_asb_single()` and `forecast_asb_multi()` follow the age-specific
  band selection paper.
- `forecast_bmbe_acf()` and `forecast_bmbe_cbd()` follow the BMBE paper.
