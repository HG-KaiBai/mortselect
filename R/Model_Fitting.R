#' Grouping demographic data from different countries with same time period
#'
#' @param obj1 one demogdata
#' @param obj2 one demogdata
#' @param startyear year strart to merge
#' @param endyear year end to merge
#'
#' @return The merged demogdata
#' @import demography
#' @export
#'
Groupdemogdata_V2<-function(obj1,obj2,startyear,endyear)
  ## V2: female and male are equally treated as new populations
{
  GroupMember<-c(obj1$label,obj2$label) #Used to record members in the group
  Ages=(0:obj1$age[length(obj1$age)])
  Year=c(startyear:endyear)
  obj1_r=extract.years(obj1,Year)
  obj2_r=extract.years(obj2,Year)
  Pop=Map("+",obj1_r$pop,obj2_r$pop)
  #return(Pop)
  weight1 = mapply(FUN = `/`, obj1_r$pop, Pop, SIMPLIFY = FALSE)
  weight2 = mapply(FUN = `/`, obj2_r$pop, Pop, SIMPLIFY = FALSE)
  Data=Map("+",mapply("*", obj1_r$rate, weight1, SIMPLIFY = FALSE),mapply("*",obj2_r$rate, weight2, SIMPLIFY = FALSE))
  names(Data)="Combined"
  Type="mortality"
  Label=GroupMember
  lambda=0
  obj <- list(year = Year, age = Ages, rate = Data,
              pop = Pop, type = Type, label = Label,
              lambda = lambda)
  return(structure(obj, class = "demogdata"))
}


#' ## Function of prefiting demographic data and calculate in-sample SSE for the target country
#' Fit the augmented common factor models,Obj is demogdata after grouping
#' ACF model_with more country specific components (N, number of components, is a tuning parameter)
#'
#' @param Obj one demogdata
#' @param Target target country
#' @param N number of country-specific components
#' @param Series "female"/"male"/"total"
#' @param Training The Training set
#'
#' @return The error
#' @export
#'
prefit=function(Obj,Target,N,Series, Training)
{
  ## common factor
  Bx_i=list()
  Kt_i=list()
  for(i in 1:length(Obj$label))
  {
    demogdata=Training[[which(names(Training)==Obj$label[i])]]
    Bx_i[[i]]=lca(demogdata,series=Series, max.age=100,interpolate=TRUE)$bx
    Kt_i[[i]]=lca(demogdata,series=Series, max.age=100,interpolate=TRUE)$kt
  }
  Bx=rowMeans(sapply(Bx_i,unlist))
  Kt=ts(rowMeans(sapply(Kt_i,unlist)),start=Target$year[1])

  ## mean
  ax=lca(Target,series=Series,max.age=100,interpolate=TRUE)$ax
  clogratesfit <- outer(Bx, Kt)
  logratesfit<- sweep(clogratesfit,1,ax,"+") # ax + B_x * Kt
  log_mx=log(Target$rate[[which(names(Target$rate)==Series)]])
  Residual_c=log_mx-logratesfit
  ##N, number of country-specific components
  svd.mx=svd(Residual_c)
  for (i in 1:N)
  {
    bx=svd.mx$u[,i]/sum(svd.mx$u[,i])
    kt=ts(svd.mx$d[i]*svd.mx$v[,i]*sum(svd.mx$u[,i]),start=start(Kt),deltat=1)
    logratesfit=logratesfit+outer(bx,kt) # a_x + B_x*K_t + sum(bx*kt)
  }
  ## Fitted logrates
  LogFit=logratesfit # log fitted value
  return(list(LogFit, log_mx))
}


#' Similar as the prefit, only treet female and male as populations
#'
#' @param Obj one demogdata
#' @param Target target function
#' @param N number of country-specific components
#' @param Training The Training set
#'
#' @return The error
#'
#' @import forecast
#' @export
#'
prefit_v2=function(Obj,Target,N, Training, h = length(Target$year))
  ## Function of prefiting demographic data and calculate in-sample SSE for the target country
  ## Fit the augmented common factor models,Obj is demogdata after grouping
  ## ACF model_with more country specific components (N, number of components, is a tuning parameter)
  ## V2: treat female and male as populations
{
  ## common factor
  Bx_i=list()
  Kt_i=list()
  for(i in 1:length(Obj$label))
  { #-8.7279142 -9.0970117 -9.3611633 -8.2742470 -8.8182983 -8.2586825  -9.2305431  -8.5122056 -8.8810366 -8.8250780
    demogdata=Training[[which(names(Training)==Obj$label[i])]]
    LogRates_i=log(demogdata$rate[[1]])
    ax_i=apply(LogRates_i, 1, mean, na.rm = TRUE)
    clogrates_i <- sweep(LogRates_i, 1, ax_i)
    svd.mx_i <- svd(clogrates_i)
    sumu_i <- sum(svd.mx_i$u[, 1])
    Bx_i[[i]]=svd.mx_i$u[, 1]/sumu_i
    Kt_i[[i]]=svd.mx_i$d[1] * svd.mx_i$v[, 1] * sumu_i
  }
  Bx=rowMeans(sapply(Bx_i,unlist))
  Kt=ts(rowMeans(sapply(Kt_i,unlist)),start=Target$year[1])
  Kt_forecast = forecast(Kt, h = h)[[2]] ############## h = 8 is from 2002 - 1994
  ## mean
  logrates=log(Training[[which(names(Training)==Target$label)]]$rate[[1]]) ######log(Target$rate[[1]])
  ax=apply(logrates, 1, mean, na.rm = TRUE)

  clogratesfit <- outer(Bx, Kt) ###################

  logratesfit<- sweep(clogratesfit,1,ax,"+")
  forecast_fit = sweep(outer(Bx, Kt_forecast),1,ax,"+") # The forecast fitted
  log_mx=log(Training[[which(names(Training)==Target$label)]]$rate[[1]]) #log(Target$rate[[1]])    ##############
  Residual_c=log_mx-logratesfit
  ##N, number of country-specific components
  svd.mx=svd(Residual_c)
  for (i in 1:N)
  {
    bx=svd.mx$u[,i]/sum(svd.mx$u[,i])
    kt=ts(svd.mx$d[i]*svd.mx$v[,i]*sum(svd.mx$u[,i]),start=start(Kt),deltat=1)

    kt_forecast = forecast(kt, h = h)[[2]] #################
    forecast_fit=forecast_fit+outer(bx,kt_forecast) ################
    logratesfit=logratesfit+outer(bx,kt)
  }
  ## Fitted logrates
  LogFit=logratesfit
  ForecastFit = forecast_fit
  return(list(ForecastFit, log(Target$rate[[1]])))
}


loss=function(Fit_true)
{
  ## Residual
  Residual=Fit_true[[2]]-Fit_true[[1]]
  ## in sample SSE
  SSE=sum(Residual^2)
  ## in sample MAE
  MAE=mean(abs(Residual))
  ## in sample MSPE
  MSPE=mean((Residual/Fit_true[[2]])^2)
  ## in sample MAPE
  MAPE=mean(abs(Residual/Fit_true[[2]]))
  Out=data.frame(SSE=SSE,MAE=MAE,MSPE=MSPE,MAPE=MAPE)
  return(Out)
}





