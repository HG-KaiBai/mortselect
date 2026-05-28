.mort_lapply <- function(X, FUN, ..., parallel = FALSE, workers = NULL) {
  if (!isTRUE(parallel) || length(X) <= 1) {
    return(lapply(X, FUN, ...))
  }

  workers <- workers %||% max(1, parallel::detectCores(logical = FALSE) - 1)
  workers <- min(as.integer(workers), length(X))
  if (is.na(workers) || workers < 1) {
    workers <- 1
  }
  if (workers == 1) {
    return(lapply(X, FUN, ...))
  }

  cl <- parallel::makeCluster(workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(
    cl,
    suppressPackageStartupMessages(requireNamespace("mortselect", quietly = TRUE))
  )
  parallel::parLapply(cl, X, FUN, ...)
}
