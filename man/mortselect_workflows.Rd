\name{forecast_dsa}
\alias{forecast_dsa}
\alias{forecast_bmbe_acf}
\alias{forecast_bmbe_cbd}
\alias{forecast_ass_single}
\alias{forecast_ass_multi}
\title{Mortality forecasting workflows}
\usage{
forecast_dsa(Dataset_raw, COG = 30, N = 5, LOSS = "SSE", t = 15,
  ages = 0:100, learning_years = 1970:1994,
  validation_years = 1995:2002, test_years = 2003:2010,
  training_years = NULL, fit_years = NULL, exclude_labels = NULL,
  exclude_indices = NULL, save_file = NULL)

forecast_bmbe_acf(Dataset_raw, Series = "female", t = 15, ages = 0:100,
  learning_years = 1970:1994, validation_years = 1995:2002,
  test_years = 2003:2010, training_years = NULL, fit_years = NULL,
  exclude_labels = NULL, exclude_indices = NULL, save_file = NULL)

forecast_bmbe_cbd(Dataset_raw, Series = "female", t = 15, ages = 55:90,
  learning_years = 1970:1994, validation_years = 1995:2002,
  test_years = 2003:2010, training_years = NULL, fit_years = NULL,
  exclude_labels = NULL, exclude_indices = NULL, save_file = NULL)

forecast_ass_single(Dataset_raw, Series = "female", t = 15, ages = 0:100,
  learning_years = 1970:1994, validation_years = 1995:2002,
  test_years = 2003:2010, training_years = NULL, fit_years = NULL,
  exclude_labels = NULL, exclude_indices = NULL, smoothing_window = 5,
  save_file = NULL)

forecast_ass_multi(Dataset_raw, Series = "female", t = 15, ages = 0:100,
  learning_years = 1970:1994, validation_years = 1995:2002,
  test_years = 2003:2010, training_years = NULL, fit_years = NULL,
  exclude_labels = NULL, exclude_indices = NULL, smoothing_window = 9,
  save_file = NULL)
}
\arguments{
\item{Dataset_raw}{A list of \code{demography::demogdata} mortality objects.}
\item{COG}{Maximum selected group size for DSA.}
\item{N}{Maximum number of country-specific components for DSA.}
\item{LOSS}{Selection loss: \code{"SSE"}, \code{"MAE"}, \code{"MSPE"}, or \code{"MAPE"}.}
\item{t}{Number of years to forecast.}
\item{ages}{Ages to retain.}
\item{learning_years}{Candidate-learning years.}
\item{validation_years}{Validation years.}
\item{test_years}{Held-out testing years.}
\item{training_years}{Final training years. Defaults to learning plus validation years.}
\item{fit_years}{Full data years retained before splitting.}
\item{exclude_labels}{Optional population labels to remove.}
\item{exclude_indices}{Optional numeric population positions to remove after filtering.}
\item{save_file}{Optional \code{.RData} path. If \code{NULL}, no file is written.}
\item{Series}{Mortality series to forecast.}
\item{smoothing_window}{Window size for smoothing selected age-band bandwidths.}
}
\value{
A named list of forecast objects. Passing \code{save_file} also writes the
returned object to disk.
}
\description{
Top-level workflows for the DSA algorithm, bivariate ACF and CBD ensembles,
and single- or multi-population age-specific band selection.
}
\details{
The default year windows reproduce the empirical setup in the source scripts.
For other mortality datasets, set \code{learning_years},
\code{validation_years}, \code{test_years}, and optionally
\code{training_years} and \code{fit_years}.
}
