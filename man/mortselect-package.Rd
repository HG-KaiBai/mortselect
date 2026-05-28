\docType{package}
\name{mortselect-package}
\alias{mortselect}
\alias{mortselect-package}
\title{Mortality forecasting with selected information borrowing}
\description{
Research implementations for deletion-substitution-addition population
selection, age-specific band selection, and bivariate model-based ensemble
mortality forecasting.
}
\details{
Input data are expected to be lists of \code{demography::demogdata} mortality
objects. Use \code{\link{prepare_mortselect_data}} to validate and split a
mortality dataset before fitting models.
}
\references{
Diao, L., Meng, Y., and Weng, C. (2020). A DSA Algorithm for Mortality
Forecasting. North American Actuarial Journal.

Meng, Y., Diao, L., and Weng, C. Mortality Prediction via Age-Specific Band
Selection.

Diao, L., Meng, Y., Weng, C., and Wirjanto, T. (2023). Enhancing Mortality
Forecasting through Bivariate Model-Based Ensemble. North American Actuarial
Journal.
}
