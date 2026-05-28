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


#' Interpolate 0-mortality rate
#'
#' @param Obj the data set
#'
#' @return 0-mortality rate fixed
fix_mortality = function(Obj){
  for(i in 1: (length(Obj)))
  {
    message("Current country is ",names(Obj)[i])
    for (gender in c("female", "male", "total"))
    {
      MortalityMatrix=Obj[[i]]$rate[[gender]]
      if(sum(MortalityMatrix==0)>0)
      {
        Index=which(MortalityMatrix==0,arr.ind=TRUE)
        for(count in 1:nrow(Index))
        {
          rowIndex=Index[count,1]
          colIndex=Index[count,2]
          if(colIndex==1)
          {
            ApproxValue=MortalityMatrix[rowIndex,colIndex+1]
          }
          else if (colIndex==ncol(MortalityMatrix))
          {
            ApproxValue=MortalityMatrix[rowIndex,colIndex-1]
          }
          else
          {
            ApproxValue=mean(c(MortalityMatrix[rowIndex,colIndex-1],MortalityMatrix[rowIndex,colIndex+1]))
          }
          MortalityMatrix[rowIndex,colIndex]=ApproxValue
        }
      }
      Obj[[i]]$rate[[gender]]=MortalityMatrix
    }
  }
  return(Obj)
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
