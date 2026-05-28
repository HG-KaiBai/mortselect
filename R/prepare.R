#' Validate a list of mortality data objects
#'
#' Checks that an input object is a list of `demogdata`-like mortality
#' records. Each record must contain `year`, `age`, `rate`, `pop`, and
#' `label` fields, and the requested series must be available in both
#' `rate` and `pop`.
#'
#' @param dataset A list of `demography::demogdata` mortality objects.
#' @param series Character vector of mortality series to check. Common
#'   values are `"female"`, `"male"`, and `"total"`.
#' @return Invisibly returns `TRUE` if validation succeeds.
#' @export
validate_demogdata_list <- function(dataset, series = c("female", "male", "total")) {
  if (!is.list(dataset)) {
    stop("`dataset` must be a list of demogdata-like objects.", call. = FALSE)
  }
  if (length(dataset) == 0) {
    stop("`dataset` must contain at least one population.", call. = FALSE)
  }

  required_fields <- c("year", "age", "rate", "pop", "label")
  for (i in seq_along(dataset)) {
    obj <- dataset[[i]]
    missing_fields <- setdiff(required_fields, names(obj))
    if (length(missing_fields) > 0) {
      stop(
        "Population at position ", i, " is missing fields: ",
        paste(missing_fields, collapse = ", "),
        call. = FALSE
      )
    }

    missing_rates <- setdiff(series, names(obj$rate))
    missing_pop <- setdiff(series, names(obj$pop))
    if (length(missing_rates) > 0 || length(missing_pop) > 0) {
      stop(
        "Population `", obj$label, "` does not contain the requested series ",
        "in both `rate` and `pop`.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

#' Prepare mortality data for mortselect workflows
#'
#' Splits a list of `demography::demogdata` objects into the learning,
#' validation, testing, and training windows used by the DSA, BMBE, ACF/CBD,
#' and age-specific band-selection workflows. The defaults reproduce the
#' empirical windows used by the accompanying papers, while all age and year
#' ranges can be supplied for another mortality dataset.
#'
#' @param dataset A list of `demography::demogdata` mortality objects.
#' @param ages Integer vector of ages to retain.
#' @param learning_years Years used to construct the candidate learning pool.
#' @param validation_years Years used for validation or model selection.
#' @param test_years Years used as the held-out testing period.
#' @param training_years Years used to fit final models. Defaults to
#'   `c(learning_years, validation_years)`.
#' @param fit_years Years kept in the full fitting data. Defaults to all
#'   supplied learning, validation, and testing years.
#' @param series Character vector of series to validate.
#' @param exclude_labels Optional population labels to remove.
#' @param exclude_indices Optional numeric positions to remove after
#'   filtering for required years. This is provided only for reproducing
#'   legacy studies that excluded populations by position.
#' @param fix_zero Logical; if `TRUE`, replace zero mortality rates by the
#'   average of adjacent years through [fix_mortality()].
#' @param relabel_sex Logical; if `TRUE`, split each population into female
#'   and male `demogdata` records named `<label>_female` and `<label>_male`.
#' @param combine.upper Logical passed to `demography::extract.ages()`.
#' @return A list with `full`, `learning`, `validation`, `test`, `training`,
#'   and `metadata` elements.
#' @export
prepare_mortselect_data <- function(dataset,
                                    ages = 0:100,
                                    learning_years = 1970:1994,
                                    validation_years = 1995:2002,
                                    test_years = 2003:2010,
                                    training_years = NULL,
                                    fit_years = NULL,
                                    series = c("female", "male", "total"),
                                    exclude_labels = NULL,
                                    exclude_indices = NULL,
                                    fix_zero = TRUE,
                                    relabel_sex = FALSE,
                                    combine.upper = FALSE) {
  training_years <- training_years %||% sort(unique(c(learning_years, validation_years)))
  fit_years <- fit_years %||% sort(unique(c(training_years, test_years)))

  dataset <- dataset[!vapply(dataset, is.null, logical(1))]
  validate_demogdata_list(dataset, series = series)

  required_years <- sort(unique(c(fit_years, learning_years, validation_years, test_years, training_years)))
  has_required_years <- vapply(
    dataset,
    function(obj) min(required_years) >= min(obj$year) && max(required_years) <= max(obj$year),
    logical(1)
  )
  dataset <- dataset[has_required_years]

  if (length(dataset) == 0) {
    stop("No populations contain the requested year range.", call. = FALSE)
  }

  if (!is.null(exclude_labels)) {
    labels <- vapply(dataset, function(obj) obj$label, character(1))
    dataset <- dataset[!labels %in% exclude_labels]
  }
  if (!is.null(exclude_indices)) {
    keep <- setdiff(seq_along(dataset), exclude_indices)
    dataset <- dataset[keep]
  }
  if (length(dataset) == 0) {
    stop("No populations remain after exclusions.", call. = FALSE)
  }

  subset_ages <- function(obj) {
    demography::extract.ages(obj, ages = ages, combine.upper = combine.upper)
  }
  subset_years <- function(obj, years) {
    demography::extract.years(obj, years = years)
  }
  label_data <- function(x) {
    names(x) <- vapply(x, GetLabel, character(1))
    x
  }
  finalize <- function(x) {
    x <- label_data(x)
    if (fix_zero) {
      x <- fix_mortality(x)
    }
    if (relabel_sex) {
      x <- relabel(x)
    }
    x
  }

  aged <- lapply(dataset, subset_ages)
  full <- finalize(lapply(aged, subset_years, years = fit_years))
  learning <- finalize(lapply(aged, subset_years, years = learning_years))
  validation <- finalize(lapply(aged, subset_years, years = validation_years))
  test <- finalize(lapply(aged, subset_years, years = test_years))
  training <- finalize(lapply(aged, subset_years, years = training_years))

  list(
    full = full,
    learning = learning,
    validation = validation,
    test = test,
    training = training,
    metadata = list(
      ages = ages,
      learning_years = learning_years,
      validation_years = validation_years,
      test_years = test_years,
      training_years = training_years,
      fit_years = fit_years,
      series = series,
      excluded_labels = exclude_labels,
      excluded_indices = exclude_indices,
      relabel_sex = relabel_sex
    )
  )
}

#' Build a population-age dictionary
#'
#' @param dataset A named list of prepared `demogdata` objects.
#' @param ages Ages to include.
#' @return A data frame with `population`, `age`, and `target` columns.
#' @export
create_population_age_dictionary <- function(dataset, ages = dataset[[1]]$age) {
  if (length(dataset) == 0) {
    stop("`dataset` must contain at least one population.", call. = FALSE)
  }
  out <- do.call(
    rbind,
    lapply(names(dataset), function(population) {
      data.frame(
        population = population,
        age = ages,
        target = paste(population, ages, sep = "."),
        stringsAsFactors = FALSE
      )
    })
  )
  row.names(out) <- NULL
  out
}

#' Bind training rates into a population-age matrix
#'
#' @param dataset A named list of prepared `demogdata` objects.
#' @param series Mortality series to bind.
#' @param ages Ages to include.
#' @return A matrix with rows named `<population>.<age>` and columns
#'   corresponding to the training years.
#' @export
bind_training_rates <- function(dataset, series = "female", ages = dataset[[1]]$age) {
  if (length(dataset) == 0) {
    stop("`dataset` must contain at least one population.", call. = FALSE)
  }

  pieces <- lapply(names(dataset), function(population) {
    rates <- dataset[[population]]$rate[[series]]
    rates <- rates[as.character(ages), , drop = FALSE]
    rownames(rates) <- paste(population, ages, sep = ".")
    rates
  })
  do.call(rbind, pieces)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
