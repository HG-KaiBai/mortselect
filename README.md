# mortselect

`mortselect` is an R package for mortality forecasting with data-driven
information borrowing. It implements workflows based on three published
methodological papers:

- deletion-substitution-addition (DSA) population selection;
- age-specific band selection (ASS/ASB) for single- and multi-population
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
- `forecast_ass_single()` fits single-population age-specific band forecasts.
- `forecast_ass_multi()` fits multi-population age-specific band forecasts.
- `forecast_bmbe_acf()` fits a BMBE workflow using ACF-type base learners.
- `forecast_bmbe_cbd()` fits a BMBE workflow using CBD-type base learners.

Supporting data helpers include:

- `prepare_mortselect_data()`
- `validate_demogdata_list()`
- `create_population_age_dictionary()`
- `bind_training_rates()`

## Data Preparation

```r
library(mortselect)

prepared <- prepare_mortselect_data(
  dataset,
  ages = 0:100,
  learning_years = 1980:2004,
  validation_years = 2005:2012,
  test_years = 2013:2020,
  exclude_labels = c("Example population to omit")
)
```

`prepare_mortselect_data()` validates the requested mortality series, filters
populations with sufficient year coverage, applies optional exclusions, repairs
zero mortality rates, and returns `learning`, `validation`, `test`, `training`,
and `full` datasets.

## Forecasting Examples

```r
dsa_forecast <- forecast_dsa(
  dataset,
  COG = 10,
  N = 3,
  LOSS = "SSE",
  t = 10,
  learning_years = 1980:2004,
  validation_years = 2005:2012,
  test_years = 2013:2020
)

single_ass <- forecast_ass_single(
  dataset,
  Series = "female",
  t = 10,
  ages = 0:100,
  learning_years = 1980:2004,
  validation_years = 2005:2012,
  test_years = 2013:2020
)

multi_ass <- forecast_ass_multi(
  dataset,
  Series = "male",
  t = 10,
  ages = 0:100,
  learning_years = 1980:2004,
  validation_years = 2005:2012,
  test_years = 2013:2020
)

acf_ensemble <- forecast_bmbe_acf(
  dataset,
  Series = "female",
  t = 10
)

cbd_ensemble <- forecast_bmbe_cbd(
  dataset,
  Series = "female",
  t = 10,
  ages = 55:90
)
```

Each workflow returns its forecast object. To also write an `.RData` result,
pass `save_file = "result.RData"`.

## Default Study Windows

The defaults reproduce the empirical setup in the source scripts and papers:

- learning years: `1970:1994`
- validation years: `1995:2002`
- testing years: `2003:2010`
- final training years: learning plus validation
- ASS/DSA age range: `0:100`
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
- `forecast_ass_single()` and `forecast_ass_multi()` follow the age-specific
  band selection paper.
- `forecast_bmbe_acf()` and `forecast_bmbe_cbd()` follow the BMBE paper.


