#' Step1: extract ages to 0:100
#'
#' @param Obj a data set
#'
#' @return a data set with 0 < age < 100
myfun1=function(Obj)
{
  return(extract.ages(Obj,ages=c(0:100)))
}

#' Step1: extract ages to 55:90
#'
#' @param Obj a data set
#'
#' @return a data set with 55 < age < 90
myfun1_CBD=function(Obj)
{
  return(extract.ages(Obj,ages=c(55:90)))
}


#' Step1: extract ages to 0:100
#'
#' @param Obj a data set
#'
#' @return a data set with 0 < age < 100
myfun1_ASS=function(Obj)
{
  return(extract.ages(Obj,ages=c(0:100),combine.upper = FALSE))
}


#' Title get label (name of the countries)
#'
#' @param Obj a data set
#'
#' @return the label from data set
GetLabel<-function(Obj)
{
  return(Obj$label)
}


#' Step2: extract time period to 1970:2010, TrainingSet:1970-2002 TestSet:2003-2010
#'
#' @param Obj the data set
#'
#' @return a data set contains year 1970 to 2010
myfun2=function(Obj)
{
  return(Obj$year[1]<=1970 && Obj$year[length(Obj$year)]>=2010)
}


#' Step2: extract time period to 1970:2010, TrainingSet:1970-2002 TestSet:2003-2010
#'
#' @param Obj the data set
#'
#' @return a data set contains year 1970 to 2010
myfun2_5=function(Obj)
{
  return(extract.years(Obj,years=c(1970:2010)))
}


#' Learning:1970-1994
#'
#' @param Obj a data set
#'
#' @return data set from year 1970 to 1994
myfun3=function(Obj)
{
  return(extract.years(Obj,years=c(1970:1994)))
}

#' Validating:1995-2002
#'
#' @param Obj a data set
#'
#' @return data set from year 1995 to 2002
myfun4=function(Obj)
{
  return(extract.years(Obj,years=c(1995:2002)))
}

#' testing:2003-2010
#'
#' @param Obj a data set
#'
#' @return data set from year 2003 to 2010
myfun5=function(Obj)
{
  return(extract.years(Obj,years=c(2003:2010)))
}

#' training=learning+Validating:1970-2002
#'
#' @param Obj a data set
#'
#' @return data set from year 1970 to 2002
myfun6=function(Obj)
{
  return(extract.years(Obj,years=c(1970:2002)))
}

#' clean for learning and training set #(countries)=30
#'
#' @param Obj a data set
#'
#' @return a cleaned set
clean_v1=function(Obj)
{
  temp = Obj[-c(12,32,34)]  #Delete subregional data
  temp = temp[-c(9,13,19,33)] #Delete regional data contains 0 death rate
  temp = temp[-12] #remove Ireland
  names(temp)=sapply(temp,GetLabel)
  return(temp)
}

#' clean for validating and testing set #(countries)=35
#'
#' @param Obj a data set
#'
#' @return a cleaned set
clean_v2=function(Obj)
{
  temp = Obj[-c(12,32,34)]  #Delete subregional data
  names(temp)=sapply(temp,GetLabel)
  return(temp)
}


#' Interpolate invalid mortality rates
#'
#' @param Obj the data set
#'
#' @return Mortality data with zero, missing, and infinite rates repaired.
#'   Populations with unrecoverable mortality rows are skipped and recorded in
#'   the `"skipped"` attribute.
fix_mortality = function(Obj, warn = TRUE){
  repaired <- list()
  skipped <- data.frame(
    population = character(),
    series = character(),
    age = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
  for(i in 1: (length(Obj)))
  {
    message("Current country is ",names(Obj)[i])
    skip_population <- FALSE
    skip_record <- NULL
    for (gender in intersect(c("female", "male", "total"), names(Obj[[i]]$rate)))
    {
      MortalityMatrix=Obj[[i]]$rate[[gender]]
      repaired_matrix <- tryCatch(
        repair_mortality_matrix(MortalityMatrix, names(Obj)[i], gender),
        error = function(e) {
          skip_population <<- TRUE
          skip_record <<- data.frame(
            population = names(Obj)[i],
            series = gender,
            age = {
              age <- attr(e, "age", exact = TRUE)
              if (is.null(age)) NA_character_ else age
            },
            reason = conditionMessage(e),
            stringsAsFactors = FALSE
          )
          NULL
        }
      )
      if (skip_population) {
        break
      }
      MortalityMatrix=repaired_matrix
      Obj[[i]]$rate[[gender]]=MortalityMatrix
    }
    if (skip_population) {
      skipped <- rbind(skipped, skip_record)
      next
    }
    repaired[[length(repaired) + 1]] <- Obj[[i]]
    names(repaired)[length(repaired)] <- names(Obj)[i]
  }
  if (nrow(skipped) > 0 && warn) {
    warning(
      "Skipped ", nrow(skipped), " population(s) with unrecoverable mortality data: ",
      paste(paste0(skipped$population, " / ", skipped$series, " / age ", skipped$age), collapse = "; "),
      call. = FALSE
    )
  }
  attr(repaired, "skipped") <- skipped
  return(repaired)
}

repair_mortality_matrix <- function(MortalityMatrix, population, series) {
  invalid <- !is.finite(MortalityMatrix) | MortalityMatrix <= 0
  if (!any(invalid, na.rm = TRUE)) {
    return(MortalityMatrix)
  }

  for (rowIndex in seq_len(nrow(MortalityMatrix))) {
    row <- as.numeric(MortalityMatrix[rowIndex, ])
    invalid_row <- !is.finite(row) | row <= 0
    if (!any(invalid_row, na.rm = TRUE)) {
      next
    }

    valid_idx <- which(!invalid_row)
    if (length(valid_idx) == 0) {
      age_label <- rownames(MortalityMatrix)[rowIndex]
      if (is.null(age_label) || is.na(age_label)) {
        age_label <- rowIndex
      }
      err <- simpleError(paste0(
        "Population `", population, "`, series `", series,
        "`, age `", age_label, "` has no positive finite mortality rates ",
        "available for interpolation."
      ))
      attr(err, "age") <- as.character(age_label)
      stop(err)
    }

    if (length(valid_idx) == 1) {
      row[invalid_row] <- row[valid_idx]
    } else {
      row[invalid_row] <- stats::approx(
        x = valid_idx,
        y = row[valid_idx],
        xout = which(invalid_row),
        rule = 2
      )$y
    }
    MortalityMatrix[rowIndex, ] <- row
  }

  MortalityMatrix
}

#' Relabel by separate the Male and Female
#'
#' @param Obj the data set
#'
#' @return Separate Gender Data set
#' @import demography
relabel=function(Obj){
  res = list()
  s = 1
  for(i in 1:length(Obj))
  {
    for(gender in c("female", "male")){
      # female
      label=paste(Obj[[i]]$label,"_",gender,sep="")
      Ages=Obj[[i]]$age
      Year=Obj[[i]]$year
      rate=Obj[[i]]$rate[gender]
      pop=Obj[[i]]$pop[gender]
      Type="mortality"
      lambda=0
      obj <- list(year = Year, age = Ages, rate = rate,
                  pop = pop, type = Type, label = label,
                  lambda = lambda)
      temp=structure(obj, class = "demogdata")
      res[[s]]=temp
      s=s+1
    }
  }
  names(res)=sapply(res,GetLabel)
  return(res)
}
