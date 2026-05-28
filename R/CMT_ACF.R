#' Grouping demographic data from different countries with same time period, obj1 and obj2 are two demogdata
#' We should extract demogdata w.r.t time by startyear and endyear
#' delta is the shift in time from obj2 to obj1
#'
#' @param obj1 pop1
#' @param obj2 pop2
#' @param startyear start
#' @param endyear end
#' @param delta age bad
#' @import demography
#'
#' @return demogdata
#'
GroupdemogdataWithDrift<-function(obj1,obj2,startyear,endyear,delta){
  GroupMember<-c(obj1$label,obj2$label) #Used to record members in the group
  Ages=(0:obj1$age[length(obj1$age)])
  Year=c(startyear:endyear)
  obj1_r=extract.years(obj1,Year)
  obj2_r=extract.years(obj2,Year)

  obj2_shift_rate=list()
  obj2_shift_pop=list()
  obj1_shift_rate=list()
  obj1_shift_pop=list()
  if(delta>0){
    YearShifted=c((startyear-delta):endyear)
  }else if (delta<0){
    YearShifted=c(startyear:(endyear-delta))
  }else if (delta==0){
    YearShifted=Year
  }
  for(i in 1:3){
    if(delta==0){
      # Obj2
      obj2_shift_rate[[i]]=obj2_r$rate[[i]]
      obj2_shift_pop[[i]]=obj2_r$pop[[i]]
      # Obj1
      obj1_shift_rate[[i]]=obj1_r$rate[[i]]
      obj1_shift_pop[[i]]=obj1_r$pop[[i]]
    } else {
      # Obj2
      temp_rate=matrix(0,ncol=abs(delta),nrow=nrow(obj2_r$rate[[i]]))
      temp_pop=matrix(0,ncol=abs(delta),nrow=nrow(obj2_r$rate[[i]]))
      obj2_shift_rate[[i]]=cbind(obj2_r$rate[[i]],temp_rate)
      obj2_shift_pop[[i]]=cbind(obj2_r$pop[[i]],temp_pop)
      # Obj1
      temp_rate=matrix(0,ncol=abs(delta),nrow=nrow(obj1_r$rate[[i]]))
      temp_pop=matrix(0,ncol=abs(delta),nrow=nrow(obj1_r$rate[[i]]))
      obj1_shift_rate[[i]]=cbind(temp_rate,obj1_r$rate[[i]])
      obj1_shift_pop[[i]]=cbind(temp_pop,obj1_r$pop[[i]])
    }
    colnames(obj2_shift_rate[[i]])=YearShifted
    colnames(obj2_shift_pop[[i]])=YearShifted
    colnames(obj1_shift_rate[[i]])=YearShifted
    colnames(obj1_shift_pop[[i]])=YearShifted
  }
  names(obj2_shift_rate)= c("female","male","total")
  names(obj2_shift_pop)= c("female","male","total")
  names(obj1_shift_rate)= c("female","male","total")
  names(obj1_shift_pop)= c("female","male","total")

  obj2_shift <- list(year = YearShifted, age = Ages, rate = obj2_shift_rate,
                     pop = obj2_shift_pop, type = "mortality", label = paste(obj2$label,"_shift"),
                     lambda = 0)
  obj2_shift=structure(obj2_shift, class = "demogdata")

  obj1_shift <- list(year = YearShifted, age = Ages, rate = obj1_shift_rate,
                     pop = obj1_shift_pop, type = "mortality", label = paste(obj1$label,"_shift"),
                     lambda = 0)
  obj1_shift=structure(obj1_shift, class = "demogdata")

  Pop=Map("+",obj1_shift$pop,obj2_shift$pop)
  #return(Pop)
  weight1 = mapply(FUN = `/`, obj1_shift$pop, Pop, SIMPLIFY = FALSE)
  weight2 = mapply(FUN = `/`, obj2_shift$pop, Pop, SIMPLIFY = FALSE)
  Data=Map("+",mapply("*", obj1_shift$rate, weight1, SIMPLIFY = FALSE),mapply("*",obj2_shift$rate, weight2, SIMPLIFY = FALSE))
  Type="mortality"
  Label=GroupMember
  lambda=0
  obj <- list(year = YearShifted, age = Ages, rate = Data,
              pop = Pop, type = Type, label = Label,
              lambda = lambda)
  return(structure(obj, class = "demogdata"))
}


#' Model ACF on Candidate set
#'
#' @param obj1 pop1
#' @param obj2 pop2
#' @param delta age band
#' @param Series gender
#' @param Candidate candidate set
#' @import demography
#' @import forecast
#'
#' @return Fit ACF on candidate set
#'
ModelingACF=function(obj1,obj2,delta,Series, Candidate){
  target_label1=obj1$label
  data.real1=Candidate[[which(names(Candidate)==target_label1)]]
  year.real=data.real1$year

  target_label2=obj2$label
  data.real2=Candidate[[which(names(Candidate)==target_label2)]]

  ## CommonFactor
  trainingdata=GroupdemogdataWithDrift(obj1,obj2,obj1$year[1],obj1$year[length(obj1$year)],delta)
  Bx=lca(trainingdata,series=Series,max.age=100,interpolate=TRUE,adjust = "none")$bx
  Kt=lca(trainingdata,series=Series,max.age=100,interpolate=TRUE,adjust = "none")$kt
  time=c(start(Kt)[1]:end(Kt)[1])

  ## obj1
  ax1=rowSums(log(obj1$rate[[which(names(obj1$rate)==Series)]]))/length(trainingdata)
  #ax1=lca(obj1,series=Series,max.age=100,interpolate=TRUE,adjust = "none")$ax
  clogratesfit <- outer(Bx, Kt)
  logratesfit1<- sweep(clogratesfit,1,ax1,"+")
  logratesfit_used1=logratesfit1[,c(which(time==obj1$year[1]):which(time==obj1$year[length(obj1$year)]))]

  MortalityMatrix1=data.real1$rate[[which(names(obj1$rate)==Series)]]
  log_mx1=log(MortalityMatrix1)
  Residual_c1=log_mx1-logratesfit_used1
  SSE_common1=sum(Residual_c1^2)
  ## country-specific factor
  svd.mx1=svd(Residual_c1)
  bx1=svd.mx1$u[,1]/sum(svd.mx1$u[,1])
  kt1=ts(svd.mx1$d[1]*svd.mx1$v[,1]*sum(svd.mx1$u[,1]),start=start(Kt),deltat=1)
  logrates_fcast1=logratesfit_used1+outer(bx1,kt1)

  SSE_total1=sum((log_mx1-logrates_fcast1)^2)

  ## obj2
  ax2=rowSums(log(obj2$rate[[which(names(obj2$rate)==Series)]]))/length(trainingdata)
  #ax2=lca(obj2,series=Series,max.age=100,interpolate=TRUE,adjust = "none")$ax
  clogratesfit <- outer(Bx, Kt)
  logratesfit2<- sweep(clogratesfit,1,ax2,"+")
  logratesfit_used2=logratesfit2[,c(which(time==obj2$year[1]):which(time==obj2$year[length(obj2$year)]))]

  MortalityMatrix2=data.real2$rate[[which(names(obj2$rate)==Series)]]
  log_mx2=log(MortalityMatrix2)
  Residual_c2=log_mx2-logratesfit_used2
  SSE_common2=sum(Residual_c2^2)
  ## country-specific factor
  svd.mx2=svd(Residual_c2)
  bx2=svd.mx2$u[,1]/sum(svd.mx2$u[,1])
  kt2=ts(svd.mx2$d[1]*svd.mx2$v[,1]*sum(svd.mx2$u[,1]),start=start(Kt),deltat=1)
  logrates_fcast2=logratesfit_used2+outer(bx2,kt2)

  SSE_total2=sum((log_mx2-logrates_fcast2)^2)

  Out=list(obj1=target_label1,obj2=target_label2,Series=Series,delta=delta,
           SSE_common1=SSE_common1,SSE_total1=SSE_total1,
           SSE_common2=SSE_common2,SSE_total2=SSE_total2)
  return(Out)
}



#' Model ACF on Validate set
#'
#' @param obj1 pop1
#' @param obj2 pop2
#' @param delta age band
#' @param Series gender
#' @param Dataset_Validate Validate set
#' @import demography
#' @import forecast
#'
#' @return Fit ACF on Validate set
#'
ValidateACF=function(obj1,obj2,delta,Series, Dataset_Validate, h = NULL){
  target_label1=obj1$label
  data.real1=Dataset_Validate[[which(names(Dataset_Validate)==target_label1)]]
  year.real1=data.real1$year

  target_label2=obj2$label
  data.real2=Dataset_Validate[[which(names(Dataset_Validate)==target_label2)]]
  h <- h %||% ncol(data.real1$rate[[Series]])

  ## Forecasting common factor
  trainingdata=GroupdemogdataWithDrift(obj1,obj2,obj1$year[1],obj1$year[length(obj1$year)],delta)
  Bx=lca(trainingdata,series=Series,max.age=100,interpolate=TRUE,adjust = "none")$bx
  Kt=lca(trainingdata,series=Series,max.age=100,interpolate=TRUE,adjust = "none")$kt
  time=c(start(Kt)[1]:end(Kt)[1])
  if(delta > 0 & delta<h){
    temp1=auto.arima(Kt) %>% forecast(h=h)
    temp2=auto.arima(Kt) %>% forecast(h=(h-delta))
    ComKT1_f=temp1$mean
    ComKT2_f=c(Kt[(length(Kt)-delta+1):length(Kt)],temp2$mean)
  }else if (delta<0 & delta>-h){
    temp2=auto.arima(Kt) %>% forecast(h=h)
    temp1=auto.arima(Kt) %>% forecast(h=(h-abs(delta)))
    ComKT2_f=temp2$mean
    ComKT1_f=c(Kt[(length(Kt)-abs(delta)+1):length(Kt)],temp1$mean)
  }else if(delta==0){
    temp1=auto.arima(Kt) %>% forecast(h=h)
    temp2=auto.arima(Kt) %>% forecast(h=h)
    ComKT1_f=temp1$mean
    ComKT2_f=temp2$mean
  }else if(delta>=h){
    temp1=auto.arima(Kt) %>% forecast(h=h)
    ComKT1_f=temp1$mean
    ComKT2_f=Kt[(length(Kt)-h+1):length(Kt)]
  }else if (delta<=-h){
    temp2=auto.arima(Kt) %>% forecast(h=h)
    ComKT2_f=temp2$mean
    ComKT1_f=Kt[(length(Kt)-h+1):length(Kt)]
  }

  ## obj1
  ax1=rowSums(log(obj1$rate[[which(names(obj1$rate)==Series)]]))/length(trainingdata)
  #ax1=lca(obj1,series=Series,max.age=100,interpolate=TRUE,adjust = "none")$ax
  clogratesfit <- outer(Bx, Kt)
  logratesfit1<- sweep(clogratesfit,1,ax1,"+")
  logratesfit1_used=logratesfit1[,c(which(time==obj1$year[1]):which(time==obj1$year[length(obj1$year)]))]

  MortalityMatrix1=obj1$rate[[which(names(obj1$rate)==Series)]]
  MortalityMatrix.validate1=data.real1$rate[[which(names(data.real1$rate)==Series)]]
  log_mx1=log(MortalityMatrix1)
  Residual_c1=log_mx1-logratesfit1_used
  logrates_fcast1=ax1+outer(Bx,ComKT1_f)
  ## country-specific factor
  svd.mx1=svd(Residual_c1)
  bx1=svd.mx1$u[,1]/sum(svd.mx1$u[,1])
  kt1=ts(svd.mx1$d[1]*svd.mx1$v[,1]*sum(svd.mx1$u[,1]),start=start(Kt),deltat=1)
  kt1_f=auto.arima(kt1) %>% forecast(h=h)
  logrates_fcast1=logrates_fcast1+outer(bx1,kt1_f$mean)

  SSE_V1=sum((log(MortalityMatrix.validate1)-logrates_fcast1)^2)

  ## obj2
  ax2=rowSums(log(obj2$rate[[which(names(obj2$rate)==Series)]]))/length(trainingdata)
  #ax2=lca(obj2,series=Series,max.age=100,interpolate=TRUE,adjust = "none")$ax
  clogratesfit <- outer(Bx, Kt)
  logratesfit2<- sweep(clogratesfit,1,ax2,"+")
  year1=obj2$year[1]-delta
  year2=obj2$year[length(obj2$year)]-delta
  logratesfit2_used=logratesfit2[,c(which(time==year1):which(time==year2))]

  MortalityMatrix2=obj2$rate[[which(names(obj2$rate)==Series)]]
  MortalityMatrix.validate2=data.real2$rate[[which(names(data.real2$rate)==Series)]]
  log_mx2=log(MortalityMatrix2)
  Residual_c2=log_mx2-logratesfit2_used
  logrates_fcast2=ax2+outer(Bx,ComKT2_f)
  ## country-specific factor
  svd.mx2=svd(Residual_c2)
  bx2=svd.mx2$u[,1]/sum(svd.mx2$u[,1])
  kt2=ts(svd.mx2$d[1]*svd.mx2$v[,1]*sum(svd.mx2$u[,1]),start=start(Kt),deltat=1)
  kt2_f=auto.arima(kt2) %>% forecast(h=h)
  logrates_fcast2=logrates_fcast2+outer(bx2,kt2_f$mean)

  SSE_V2=sum((log(MortalityMatrix.validate2)-logrates_fcast2)^2)



  Out=list(obj1=target_label1,obj2=target_label2,Series=Series,delta=delta,
           SSE_V1=SSE_V1,SSE_V2=SSE_V2)
  return(Out)
}


#' BM ACF on both candidate and validate
#'
#' @param target target
#' @param reference ref
#' @param Series gender
#' @param Candidate candidate set
#' @param Dataset_Validate validate set
#'
#' @return the best delta
#'
JointModelACF=function(target,reference,Series, Candidate, Dataset_Validate, h = NULL){
  target_label=target$label
  reference_label=reference$label
  Out=data.frame()
  for(delta in -10:10)
  {
    M=ModelingACF(target,reference,delta,Series, Candidate)
    V=ValidateACF(target,reference,delta,Series, Dataset_Validate, h = h)
    temp=data.frame(target=target_label,reference=reference_label,delta=delta,
                    SSE_common1=M$SSE_common1,SSE_total1=M$SSE_total1,
                    SSE_common2=M$SSE_common2,SSE_total2=M$SSE_total2,
                    SSE_V1=V$SSE_V1,SSE_V2=V$SSE_V2)
    Out=rbind(Out,temp)
  }
  return(Out)
}



#' precdiction on test set
#'
#' @param Obj target
#' @param Series gender
#' @param Dataset_training traing set
#' @param Dataset_Test test set
#' @param h period to forecast
#'
#' @import demography
#' @import forecast
#' @return the prediction
Predict=function(Obj,Series, Dataset_training, Dataset_Test, h= 8){
  target=as.character(Obj$target)
  reference=as.character(Obj$reference)
  delta=Obj$delta
  Series=Series

  #training+validate set
  data.target.trainvalidate=Dataset_training[[which(names(Dataset_training)==target)]]
  data.ref.trainvalidate=Dataset_training[[which(names(Dataset_training)==reference)]]
  trainingdata=data.target.trainvalidate


  #test set
  data.test=Dataset_Test[[which(names(Dataset_Test)==target)]]
  year.real = c((trainingdata$year[length(trainingdata$year)]+1):
                  (trainingdata$year[length(trainingdata$year)]+h))
  #year.real=data.test$year
  # merge data with delta
  GroupedData=GroupdemogdataWithDrift(data.target.trainvalidate,data.ref.trainvalidate,
                                      data.target.trainvalidate$year[1],data.target.trainvalidate$year[length(data.target.trainvalidate$year)],
                                      delta)
  ## Forecasting with ACF model
  Bx=lca(GroupedData,series=Series,max.age=100,interpolate=TRUE,adjust = "none")$bx
  Kt=lca(GroupedData,series=Series,max.age=100,interpolate=TRUE,adjust = "none")$kt
  time=c(start(Kt)[1]:end(Kt)[1])
  ## Common Kt forcast
  if(delta >= 0){
    temp=auto.arima(Kt) %>% forecast(h=h)
    ComKT_f=temp$mean
  }else{ ####################
    temp1=Kt[c(which(time==year.real[1]):length(Kt))]
    if(abs(delta)< h){
      temp2=auto.arima(Kt) %>% forecast(h=(h-abs(delta)))
      ComKT_f=c(temp1,temp2$mean)
    }else{
      ComKT_f=temp1[c((length(temp1)-(h-1)):length(temp1))]
    }
  }
  ## mean
  ax=rowSums(log(trainingdata$rate[[which(names(trainingdata$rate)==Series)]]))/length(GroupedData)
  #ax=lca(trainingdata,series=Series,max.age=100,interpolate=TRUE,adjust="none")$ax
  clogratesfit <- outer(Bx, Kt)
  logratesfit<- sweep(clogratesfit,1,ax,"+")
  logratesfit_used=logratesfit[,c(which(time==trainingdata$year[1]):which(time==trainingdata$year[length(trainingdata$year)]))]

  MortalityMatrix=data.target.trainvalidate$rate[[which(names(data.target.trainvalidate$rate)==Series)]]
  log_mx=log(MortalityMatrix)
  Residual_c=log_mx-logratesfit_used
  logrates_fcast=ax+outer(Bx,ComKT_f)
  ## country-specific factor
  svd.mx=svd(Residual_c)
  bx=svd.mx$u[,1]/sum(svd.mx$u[,1])
  kt=ts(svd.mx$d[1]*svd.mx$v[,1]*sum(svd.mx$u[,1]),start=start(Kt),deltat=1)
  kt_f=auto.arima(kt) %>% forecast(h=h)
  logrates_fcast=logrates_fcast+outer(bx,kt_f$mean)
  colnames(logrates_fcast) = year.real
  #SSE=sum((log(data.test$rate[[which(names(data.test$rate)==Series)]])-logrates_fcast)^2)
  Out=list(target=target,reference=reference,delta=delta,Series=Series,lograte_f=logrates_fcast,Kt=Kt,kt=kt)
  #Out=list(target=target,reference=reference,delta=delta,Series=Series,lograte_f=logrates_fcast,SSE=SSE,Kt=Kt,kt=kt)
  return(Out)
}



#' Validation with simple avg
#'
#' @param Obj target
#' @param Series gender
#' @param Dataset_Test test set
#' @param Result_DeltaChoice_SumSSE_V best delta
#' @param LogRatePrediction prediction
#'
#' @return the sse
#'
PerformanceMeasure_ACF_SimAvg=function(Obj,Series, Dataset_Test,Result_DeltaChoice_SumSSE_V,LogRatePrediction){
  target=Obj$label
  data.test=Dataset_Test[[which(names(Dataset_Test)==target)]]

  Index=which(Result_DeltaChoice_SumSSE_V$target==target)
  logrates_fcast=0
  for(i in 1:length(Index)){
    logrates_fcast=logrates_fcast+LogRatePrediction[[Index[i]]]$lograte_f
  }
  logrates_fcast=logrates_fcast/length(Index)
  pred = exp(logrates_fcast)
  SSE=sum((log(data.test$rate[[which(names(data.test$rate)==Series)]])-logrates_fcast)^2)
  Out=list(target=target,Series=Series, forecast=pred ,SSE=SSE)
  return(Out)
}



#' Validate_Forecast
#'
#' @param target country
#' @param GroupedData groupdata
#' @param Series female/male
#' @param Dataset_Validate valid set
#'
#' @return Validate_Forecast
Validate_Forecast=function(target,GroupedData,Series,Dataset_Validate, h = NULL)
{
  Series=Series
  target_label=target$label
  data.real=Dataset_Validate[[which(names(Dataset_Validate)==target_label)]]
  year.real=data.real$year
  h <- h %||% ncol(data.real$rate[[Series]])
  trainingdata=GroupedData
  year.used=GroupedData$year

  ## Forecasting with ACF model
  Bx=lca(trainingdata,series=Series,max.age=100,interpolate=TRUE,adjust = "none")$bx
  Kt=lca(trainingdata,series=Series,max.age=100,interpolate=TRUE,adjust = "none")$kt
  time=c(start(Kt)[1]:end(Kt)[1])
  if(year.used[length(year.used)]<=(year.real[1]-1))
  {
    temp=auto.arima(Kt) %>% forecast(h=h)
    ComKT_f=temp$mean
  }
  else
  {
    temp1=Kt[c(which(time==year.real[1]):length(Kt))]
    if(length(temp1)<h)
    {
      temp2=auto.arima(Kt) %>% forecast(h=(h-length(temp1)))
      ComKT_f=c(temp1,temp2$mean)
    }
    else
    {
      ComKT_f=temp1[c((length(temp1)-h+1):length(temp1))]
    }
  }

  ## mean
  ax=rowSums(log(target$rate[[which(names(target$rate)==Series)]]))/length(GroupedData)
  #ax=lca(target,series=Series,max.age=100,interpolate=TRUE,adjust = "none")$ax
  clogratesfit <- outer(Bx, Kt)
  logratesfit<- sweep(clogratesfit,1,ax,"+")
  logratesfit_used=logratesfit[,c(which(time==target$year[1]):which(time==target$year[length(target$year)]))]

  MortalityMatrix=target$rate[[which(names(target$rate)==Series)]]
  MortalityMatrix.validate=data.real$rate[[which(names(data.real$rate)==Series)]]
  log_mx=log(MortalityMatrix)
  Residual_c=log_mx-logratesfit_used
  logrates_fcast=ax+outer(Bx,ComKT_f)
  ## country-specific factor
  svd.mx=svd(Residual_c)
  bx=svd.mx$u[,1]/sum(svd.mx$u[,1])
  kt=ts(svd.mx$d[1]*svd.mx$v[,1]*sum(svd.mx$u[,1]),start=start(Kt),deltat=1)
  kt_f=auto.arima(kt) %>% forecast(h=h)
  logrates_fcast=logrates_fcast+outer(bx,kt_f$mean)

  SSE=sum((log(MortalityMatrix.validate)-logrates_fcast)^2)
  Out=list(target=target_label,Series=Series,SSE=SSE,logrates_fcast=logrates_fcast)
  return(Out)
}


#' RankAndAvg
#'
#' @param Obj Obj
#' @param Series female / male
#' @param Result_DeltaChoice_SumSSE_V best delta
#' @param Dataset_Validate valid set
#' @param Validate_Forecast_Result forecast
#'
#' @return forecast result
RankAndAvg=function(Obj,Series, Result_DeltaChoice_SumSSE_V, Dataset_Validate,
                    Validate_Forecast_Result)
{
  Series=Series
  target_label=Obj$label
  temp=Result_DeltaChoice_SumSSE_V[which(Result_DeltaChoice_SumSSE_V$target==target_label),]
  tempSorted=temp[order(temp$SSE),]

  data.real=Dataset_Validate[[which(names(Dataset_Validate)==target_label)]]
  Result=data.frame()
  max_groups <- max(1, nrow(tempSorted))
  for(Groupsize in seq_len(max_groups))
  {
    CandidatePool=as.character(tempSorted$reference[c(1:Groupsize)])
    Index=which(Result_DeltaChoice_SumSSE_V$target==target_label & Result_DeltaChoice_SumSSE_V$reference %in% CandidatePool)

    logrates_fcast=0
    for(i in 1:length(Index))
    {
      logrates_fcast=logrates_fcast+Validate_Forecast_Result[[Index[i]]]$logrates_fcast
    }
    logrates_fcast=logrates_fcast/length(Index)

    SSE=sum((log(data.real$rate[[which(names(data.real$rate)==Series)]])-logrates_fcast)^2)
    temp=data.frame(target=target_label,Groupsize=Groupsize,SSE=SSE)
    Result=rbind(Result,temp)
  }
  GroupsizeChoice=which(Result$SSE==min(Result$SSE))
  ReferencePool=as.character(tempSorted$reference[c(1:GroupsizeChoice)])
  Out=list(target=target_label,GroupsizeChoice=GroupsizeChoice,ReferencePool=ReferencePool,Result=Result)
  return(Out)
}


#' PerformanceMeasure_RankAvg
#'
#' @param Obj Obj
#' @param Series female / male
#' @param Dataset_Test test set
#' @param PredictionRule_RankAvg rank avg
#' @param Result_DeltaChoice_SumSSE_V best delta
#' @param LogRatePrediction log pred
#'
#' @return PerformanceMeasure_RankAvg
PerformanceMeasure_RankAvg=function(Obj,Series,Dataset_Test,PredictionRule_RankAvg,
                                    Result_DeltaChoice_SumSSE_V,LogRatePrediction){
  Series=Series
  target=Obj$label
  data.test=Dataset_Test[[which(names(Dataset_Test)==target)]]

  ReferencePool=PredictionRule_RankAvg[[which(names(PredictionRule_RankAvg)==target)]]$ReferencePool
  Index=which(Result_DeltaChoice_SumSSE_V$target==target & Result_DeltaChoice_SumSSE_V$reference %in% ReferencePool)

  logrates_fcast=0
  for(i in 1:length(Index))
  {
    logrates_fcast=logrates_fcast+LogRatePrediction[[Index[i]]]$lograte_f
  }
  logrates_fcast=logrates_fcast/length(Index)
  pred = exp(logrates_fcast)
  #SSE=sum((log(data.test$rate[[which(names(data.test$rate)==Series)]])-logrates_fcast)^2)
  Out=list(target=target,Series=Series, forecast = pred)
  #Out=list(target=target,Series=Series,SSE=SSE, forecast = pred)
  return(Out)
}





















