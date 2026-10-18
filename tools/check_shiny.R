# Run from the repository root: Rscript tools/check_shiny.R
setwd("shiny")
env <- new.env()
sys.source("app.R", envir=env)
shiny::testServer(env$server, {
  for (method in unname(env$methods)) {
    d <- env$forecasts[[method]]
    population <- "Canada"
    row <- which(d$series$population==population & d$series$sex=="female")
    ages_available <- d$ages
    # The age-zero case specifically detects numeric-row off-by-one errors.
    for (age in c(min(ages_available),60,max(ages_available))) {
      session$setInputs(method=method,mode="age",populations=population,sex="female",age=age,log=FALSE)
      actual <- rows()
      stopifnot(nrow(actual)==15, all(actual$age==age),
        isTRUE(all.equal(actual$mortality_rate, unname(d$series$values[[row]][match(age,d$ages),]))))
      ggplot2::ggplot_build(chart())
    }
    chosen <- c(min(ages_available),60,max(ages_available))
    session$setInputs(mode="population",ages=as.character(chosen),log=TRUE)
    stopifnot(nrow(rows())==45, all(sort(unique(rows()$age))==sort(chosen)))
    ggplot2::ggplot_build(chart())
    session$setInputs(mode="year",year=2025,sex=c("female","male"),log=FALSE)
    stopifnot(nrow(rows())==length(ages_available)*2, all(rows()$year==2025))
    ggplot2::ggplot_build(chart())
    cat(method, ": all views, age boundaries, exact rates and chart generation passed\n")
  }
  session$setInputs(populations=character())
  err <- tryCatch({rows();NULL}, error=identity)
  stopifnot(inherits(err,"shiny.silent.error"))
  cat("Empty selection validation passed\n")
})
