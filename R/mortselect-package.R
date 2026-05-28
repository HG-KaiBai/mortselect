#' mortselect: mortality forecasting with selected information borrowing
#'
#' The package collects research implementations for three related mortality
#' forecasting ideas:
#'
#' * deletion-substitution-addition (DSA) selection for choosing auxiliary
#'   populations before fitting an augmented common factor model;
#' * age-specific band selection (single- and multi-population versions) for
#'   borrowing information from nearby ages or similar population-age series;
#' * bivariate model-based ensembles (BMBE) with augmented common factor and
#'   CBD-type base learners.
#'
#' Input data are expected to be `demography::demogdata` mortality objects or
#' lists of those objects. Use [prepare_mortselect_data()] to adapt an arbitrary
#' mortality dataset to the age and year windows needed by the modeling
#' workflows.
#'
#' @references
#' Diao, L., Meng, Y., and Weng, C. (2020). "A DSA Algorithm for Mortality
#' Forecasting." North American Actuarial Journal.
#'
#' Meng, Y., Diao, L., and Weng, C. "Mortality Prediction via Age-Specific Band
#' Selection."
#'
#' Diao, L., Meng, Y., Weng, C., and Wirjanto, T. (2023). "Enhancing Mortality
#' Forecasting through Bivariate Model-Based Ensemble." North American
#' Actuarial Journal.
#'
#' @keywords internal
"_PACKAGE"
