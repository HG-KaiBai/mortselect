\name{DSA_V5}
\alias{DSA_V5}
\title{Deletion-substitution-addition population selection}
\usage{
DSA_V5(target, Pool, COG, N, LossFunction)
}
\arguments{
\item{target}{Target \code{demogdata} population.}
\item{Pool}{Named candidate-pool list of \code{demogdata} populations.}
\item{COG}{Maximum group size to consider.}
\item{N}{Maximum number of country-specific components.}
\item{LossFunction}{Selection loss: \code{"SSE"}, \code{"MAE"},
\code{"MSPE"}, or \code{"MAPE"}.}
}
\value{
A list containing the target, tuning settings, candidate names, selected groups,
and loss matrix.
}
\description{
Runs the DSA algorithm used by \code{\link{forecast_dsa}} to identify auxiliary
populations for augmented common factor mortality forecasting.
}
