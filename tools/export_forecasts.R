# Run from the repository: Rscript tools/export_forecasts.R ../Shiny docs/data
# Only forecast matrices and their dimension labels are exported.
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1L) args[1] else "../Shiny"
output <- if (length(args) >= 2L) args[2] else "docs/data"
stopifnot(requireNamespace("jsonlite", quietly = TRUE))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
specs <- list(
  DSA = c("2025_DSA_Prediction.RData", "DSA", "DSA", "Deletion-substitution-addition", "Population selection", "2021"),
  SingleASS = c("2025_Single_ASS_Prediction.RData", "SingleASS", "Single-population ASB", "Age-specific band selection", "Single-population selection", "2026"),
  MultiASS = c("2025_Multi_ASS_Prediction.RData", "MultiASS", "Multi-population ASB", "Age-specific band selection", "Multi-population selection", "2026"),
  ACF = c("2025_ACF_Prediction.RData", "ACF", "ACF ensemble", "Bivariate model-based ensemble", "ACF-type base learners", "2023"),
  CBD = c("2025_CBD_55-90_Prediction.RData", "CBD", "CBD ensemble", "Bivariate model-based ensemble", "CBD-type base learners", "2023")
)
manifest <- list(schema_version = 1L, rate_unit = "mortality rate", methods = list())
all_populations <- character()
for (id in names(specs)) {
  spec <- specs[[id]]
  env <- new.env(parent = emptyenv())
  load(file.path(input, spec[1]), envir = env)
  entries <- env[[spec[2]]]
  first <- entries[[1]]$forecast
  ages <- as.integer(rownames(first))
  years <- as.integer(colnames(first))
  populations <- sort(unique(vapply(entries, `[[`, "", "target")))
  series <- lapply(entries, function(entry) {
    m <- entry$forecast
    stopifnot(identical(dimnames(m), dimnames(first)), all(is.finite(m)), all(m >= 0))
    list(population = entry$target, sex = entry$Series, values = unname(m))
  })
  payload <- list(method = id, ages = ages, years = years, series = unname(series))
  file <- paste0(id, ".json")
  jsonlite::write_json(payload, file.path(output, file), auto_unbox = TRUE, digits = NA, matrix = "rowmajor")
  # Round-trip every number, age and year before publishing.
  restored <- jsonlite::read_json(file.path(output, file), simplifyVector = TRUE)
  stopifnot(identical(restored$ages, ages), identical(restored$years, years))
  for (i in seq_along(entries)) {
    stopifnot(isTRUE(all.equal(unname(entries[[i]]$forecast), restored$series$values[[i]], tolerance = 1e-14)))
  }
  manifest$methods[[length(manifest$methods) + 1L]] <- list(
    id = id, label = spec[3], description = spec[4], family = spec[5], paper_year = spec[6],
    file = file, source_file = spec[1], source_md5 = unname(tools::md5sum(file.path(input, spec[1]))),
    ages = ages, years = years, populations = populations,
    sexes = sort(unique(vapply(entries, `[[`, "", "Series"))), series_count = length(series)
  )
  all_populations <- union(all_populations, populations)
  cat(id, ":", length(series), "series;", length(ages), "ages;", length(years), "years verified\n")
}
manifest$population_count <- length(all_populations)
jsonlite::write_json(manifest, file.path(output, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
