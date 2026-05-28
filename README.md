# mortselect

`mortselect` packages mortality-forecasting workflows from the papers in
`papers/`:

- DSA population selection for augmented common factor forecasts.
- Age-specific band selection for single- and multi-population forecasts.
- Bivariate model-based ensembles using ACF- and CBD-type base learners.

The code expects mortality data as a list of `demography::demogdata` objects,
such as objects read from the Human Mortality Database with the `demography`
package. The public functions expose age and year windows as arguments, so the
methods are not tied to the original 1970-2010 study period.

## Data preparation

```r
library(mortselect)

prepared <- prepare_mortselect_data(
  dataset,
  ages = 0:100,
  learning_years = 1980:2004,
  validation_years = 2005:2012,
  test_years = 2013:2020,
  exclude_labels = c("Example region to omit")
)
```

`prepare_mortselect_data()` validates that each population has the requested
series, filters populations with sufficient year coverage, optionally removes
named populations, repairs zero mortality rates, and returns prepared
`learning`, `validation`, `test`, `training`, and `full` datasets.

## Forecasting examples

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

single_asb <- forecast_ass_single(
  dataset,
  Series = "female",
  t = 10,
  ages = 0:100,
  learning_years = 1980:2004,
  validation_years = 2005:2012,
  test_years = 2013:2020
)

multi_asb <- forecast_ass_multi(
  dataset,
  Series = "male",
  t = 10,
  ages = 0:100,
  learning_years = 1980:2004,
  validation_years = 2005:2012,
  test_years = 2013:2020
)

acf_ensemble <- forecast_bmbe_acf(dataset, Series = "female", t = 10)
cbd_ensemble <- forecast_bmbe_cbd(dataset, Series = "female", t = 10, ages = 55:90)
```

Each top-level workflow returns its forecast object. To reproduce the older
script behavior of writing an `.RData` file, pass `save_file = "result.RData"`.

## Notes

The default year windows reproduce the empirical setup used in the source
scripts and papers:

- learning: 1970-1994
- validation: 1995-2002
- testing: 2003-2010
- final training: learning plus validation

For a different mortality dataset, supply the windows explicitly and use
`exclude_labels` rather than numeric exclusions whenever possible.
