\name{prepare_mortselect_data}
\alias{prepare_mortselect_data}
\alias{validate_demogdata_list}
\alias{create_population_age_dictionary}
\alias{bind_training_rates}
\title{Prepare mortality data for mortselect workflows}
\usage{
prepare_mortselect_data(dataset, ages = 0:100,
  learning_years = 1970:1994, validation_years = 1995:2002,
  test_years = 2003:2010, training_years = NULL, fit_years = NULL,
  series = c("female", "male", "total"), exclude_labels = NULL,
  exclude_indices = NULL, fix_zero = TRUE, relabel_sex = FALSE,
  combine.upper = FALSE)

validate_demogdata_list(dataset, series = c("female", "male", "total"))

create_population_age_dictionary(dataset, ages = dataset[[1]]$age)

bind_training_rates(dataset, series = "female", ages = dataset[[1]]$age)
}
\arguments{
\item{dataset}{A list of \code{demography::demogdata} mortality objects.}
\item{ages}{Integer vector of ages to retain.}
\item{learning_years}{Years used to construct the candidate learning pool.}
\item{validation_years}{Years used for validation or model selection.}
\item{test_years}{Years used as the held-out testing period.}
\item{training_years}{Years used to fit final models. Defaults to learning plus validation years.}
\item{fit_years}{Years kept in the full fitting data. Defaults to training plus testing years.}
\item{series}{Mortality series to validate or bind, usually \code{"female"}, \code{"male"}, or \code{"total"}.}
\item{exclude_labels}{Optional population labels to remove.}
\item{exclude_indices}{Optional numeric positions to remove after year filtering.}
\item{fix_zero}{If \code{TRUE}, replace zero rates with adjacent-year averages.}
\item{relabel_sex}{If \code{TRUE}, split each population into female and male records.}
\item{combine.upper}{Passed to \code{demography::extract.ages()}.}
}
\value{
\code{prepare_mortselect_data()} returns a list with \code{full},
\code{learning}, \code{validation}, \code{test}, \code{training}, and
\code{metadata}. The other helpers return validation status, a dictionary data
frame, or a population-age training-rate matrix.
}
\description{
Validates, filters, splits, labels, and optionally relabels mortality data for
the forecasting workflows.
}
