#' compute binomial likelihood
#'
#' @param obs observations
#' @param fit fitted
#' @param exposure exposure
#' @param weight the w
#'
#'
computeLogLikBinomial=function(obs, fit, exposure, weight){
  ind <- (weight > 0)
  res <- array(NA, dim = length(weight))
  res[ind] <- weight[ind] * (exposure[ind] * (obs[ind] * log(fit[ind]) +
                                                (1 - obs[ind]) * log(1 - fit[ind])) +
                                                lchoose(round(exposure[ind]),
                                                round(obs[ind] * exposure[ind])))
  sum(res, na.rm = TRUE)
}



#' forecast multivariate rwd
#'
#' @param object
#' @param h
#' @param level
#' @param fan
#'
#'
forecast.mrwd=function(object, h = 10, level = c(80, 95), fan = FALSE, ...){
  x <- object$x
  nn <- 1:h
  nYear <- ncol(x)
  N <- nrow(x)
  yearsFor <- (as.numeric(colnames(x)[nYear]) + 1):(as.numeric(colnames(x)[nYear]) +
                                                      h)
  mean <- x[, nYear] + t(array(nn, c(h, N))) * array(object$drift,
                                                     c(N, h))
  rownames(mean) <- rownames(x)
  colnames(mean) <- yearsFor
  if (fan)
    level <- seq(51, 99, by = 3)
  else {
    if (min(level) > 0 & max(level) < 1)
      level <- 100 * level
    else if (min(level) < 0 | max(level) > 99.99)
      stop("Confidence limit out of range")
  }
  nn <- 1:h
  se <- sqrt(t(array(nn, c(h, N))) * array(diag(object$sigma),
                                           c(N, h)))
  nconf <- length(level)
  z <- qnorm(0.5 + level/200)
  lower <- upper <- array(NA, c(N, h, nconf), dimnames = list(rownames(x),
                                                              yearsFor, paste(level, "%", sep = "")))
  for (i in 1:nconf) {
    lower[, , i] <- mean - z[i] * se
    upper[, , i] <- mean + z[i] * se
  }
  structure(list(model = object, level = level, mean = mean,
                 lower = lower, upper = upper), class = "mrwdForecast")
}


#' Benchmark models, Fit single-population CBD1 model
#'
#' @param Dataset_training
#' @param Series
#' @param Dataset_Test
#'
#' @import StMoMo
#' @import forecast
#' @return
#'
single_CBD1 = function(Dataset_training, Series, Dataset_Test){
  CBD1model=cbd()
  CBD1result=data.frame()
  CBD1Error=list()
  for(index in 1:length(Dataset_training))
  {
    fitData=StMoMoData(Dataset_training[[index]], series = Series, type = "central")
    target.label=names(Dataset_training)[index]
    data.real.central=StMoMoData(Dataset_Test[[which(names(Dataset_Test)==target.label)]],series=Series,type="central")
    data.real.initial=central2initial(data.real.central)
    CBD1fit <- fit(CBD1model, data = central2initial(fitData))
    Forecast=forecast(CBD1fit, h = 8, level = c(55, 90))
    qxt=Forecast$rates
    logit_qxt=log(qxt/(1-qxt))
    qxt.real=data.real.initial$Dxt/data.real.initial$Ext
    logit_qxt.real=log(qxt.real/(1-qxt.real))
    #performance meaasure
    SSE=sum((logit_qxt.real-logit_qxt)^2)
    temp=data.frame(target=target.label,SSE=SSE)
    CBD1result=rbind(CBD1result,temp)
    Residual=logit_qxt.real-logit_qxt
    Error_CBD1=(logit_qxt.real-logit_qxt)^2
    CBD1Error[[index]]=list(target=target.label,Residual=Residual,SquaredError=Error_CBD1)
    names(CBD1Error)[index]=target.label
  }
  write.csv(CBD1result,file="CBD1result_logitqxt_female.csv")
  save(CBD1Error,file="CBD1Error_female.RData")
}


#' Help with the BMBE
#'
#' @param fitData
#' @param delta
#' @param index
#' @import reshape2
#'
#' @return
model_pop = function(fitData, delta, index = 1){
  ## Population####
  ## Death number and Exposure
  Dxt=fitData$Dxt
  Ext=fitData$Ext
  ModelingDataD <- (reshape2::melt(Dxt, value.name = "D",
                                    varnames = c("x", "t")))
  ModelingDataE <- (reshape2::melt(Ext, value.name = "E",
                                    varnames = c("x", "t")))
  ModelingData <- merge(ModelingDataD, ModelingDataE)
  ModelingData$j= index
  ## Cohort
  ModelingData <- transform(ModelingData, c = t - x)
  ## Age, year and optimal matrix
  ages=fitData$ages
  years=fitData$years
  nAges <- length(ages)
  nYears <- length(years)
  oxt <- matrix(0, nrow = nAges, ncol = nYears)
  rownames(oxt) <- ages
  colnames(oxt) <- years
  ModelingData$o <- 0
  ## Weight
  wxt <- matrix(1, nrow = nAges, ncol = nYears)
  rownames(wxt) <- ages
  colnames(wxt) <- years
  fitDataW <- (reshape2::melt(wxt, value.name = "w",
                              varnames = c("x", "t")))
  ModelingData <- merge(ModelingData, fitDataW)
  ## CBD B_i term
  ModelingData$B2 <- ModelingData$x-mean(ages)
  ## Year_shfted term
  if(index == 1){
    ModelingData$t_shifted=ModelingData$t
  } else{
    ModelingData$t_shifted=ModelingData$t-delta
  }
  return(ModelingData)
}




#' FitBivariateCBD1withDelta_modeling
#'
#' @param index1
#' @param index2
#' @param Series
#' @param delta
#' @param Candidate
#'
#' @import reshape2
#' @import gnm
#'
FitBivariateCBD1withDelta_modeling=function(index1,index2,Series,delta, Candidate){
  fitData1=StMoMoData(Candidate[[index1]], series = Series, type = "initial")
  fitData2=StMoMoData(Candidate[[index2]], series = Series, type = "initial")
  target.label1=names(Candidate)[index1]
  target.label2=names(Candidate)[index2]
  ages=fitData2$ages
  years=fitData2$years
  ModelingData1 = model_pop(fitData1, delta, 1)
  ModelingData2 = model_pop(fitData2, delta, 2)
  ## Aggregate Data ####
  ModelingData=rbind(ModelingData1,ModelingData2)
  ## Model fitting ####
  Formula="D/E ~ -1 + offset(o) + factor(t_shifted) + B2:factor(t):factor(j)"
  fittingModel <- gnm(formula = as.formula(Formula),
                      data = ModelingData, family = binomial(link = "logit"),
                      weights = ModelingData$E * ModelingData$w,
                      verbose = verbose)
  ## loglikelihood and BIC ####
  qhatxt=fittingModel$fitted.values
  loglik <- computeLogLikBinomial(obs = ModelingData$D/ModelingData$E, fit = qhatxt,
                                  exposure = ModelingData$E, weight = ModelingData$w)
  npar=length(fittingModel$coefficients)
  nobs=length(fittingModel$fitted.values)
  BIC=0.5*npar*log(nobs)-loglik
  years.shifted=sort(unique(ModelingData$t_shifted))
  ## Extract series ####
  Kt=fittingModel$coefficients[1:length(years.shifted)];names(Kt)=years.shifted
  kt_1=fittingModel$coefficients[(length(years.shifted)+1):(length(years.shifted)+length(years))];names(kt_1)=years
  kt_2=fittingModel$coefficients[(length(years.shifted)+1+length(years)):(length(years.shifted)+2*length(years))];names(kt_2)=years
  Out=list(label1=target.label1,label2=target.label2,ages=ages,years=years,Deltat=delta,years.shifted=years.shifted,
           loglik=loglik,BIC=BIC,
           Kt=Kt,kt_1=kt_1,kt_2=kt_2)
  return(Out)
}

#' Extrapolate those bivariate-population CBD-1 with Delta model
#'
#' @param Obj kt forecast
#' @param h period to forecast
#'
#' @import StMoMo
#'
ExtrapolateBivarateCBD1withDelta=function(Obj, h = 8){
  label1=Obj$label1;label2=Obj$label2;ages=Obj$ages;years=Obj$years
  years.shifted=Obj$years.shifted
  delta=Obj$Deltat
  Kt=Obj$Kt;kt_1=Obj$kt_1;kt_2=Obj$kt_2
  B2=ages-mean(ages)
  # kt_1 and kt_2
  MRWD=mrwd(rbind(kt_1,kt_2))
  F=forecast.mrwd(MRWD,h=h,level = c(55, 90), fan = FALSE)
  kt_1f=F$mean[1,]
  kt_2f=F$mean[2,]
  # Kt
  Kt_f=rwf(Kt,h=h,drift = TRUE)$mean
  names(Kt_f)=seq_len(h)+years.shifted[length(years.shifted)]
  Kt.value=c(Kt,Kt_f)
  # depending on the value of delta, Kt would use different value
  if(delta>=0) ###############################
  {
    index2=c((length(Kt.value)-delta-(h-1)):(length(Kt.value)-delta))
    Kt1.value=Kt_f
    Kt2.value=Kt.value[index2]
  }else  ####################
  {
    index1=c((length(Kt.value)-abs(delta)-(h-1)):(length(Kt.value)-abs(delta)))
    Kt1.value=Kt.value[index1]
    Kt2.value=Kt_f
  }
  temp1=sweep(outer(B2,kt_1f),2,Kt1.value,"+")
  qxt1=exp(temp1)/(exp(temp1)+1)
  rownames(qxt1)=ages
  mxt1=-log(1-qxt1)
  temp2=sweep(outer(B2,kt_2f),2,Kt2.value,"+")
  qxt2=exp(temp2)/(exp(temp2)+1)
  rownames(qxt2)=ages
  mxt2=-log(1-qxt2)
  Out=list(label1=label1,label2=label2,ages=ages,years=years,delta=delta,
           mxt1=mxt1,mxt2=mxt2,qxt1=qxt1,qxt2=qxt2)
}


#' Validation SSE for each bivariate-population CBD-1 with Delta model
#'
#' @param Obj
#' @param Series
#' @param Dataset_Validate
#'
ValidationSSE_CBD1withDelta=function(Obj,Series, Dataset_Validate){
  Seires=Series
  Deltat=Obj$delta
  target1=Obj$label1
  target2=Obj$label2
  data.validate1=Dataset_Validate[[which(names(Dataset_Validate)==target1)]]
  data.real1.central=StMoMoData(data.validate1,series=Series,type="central")
  data.real1.initial=central2initial(data.real1.central)
  qxt1.real=data.real1.initial$Dxt/data.real1.initial$Ext

  data.validate2=Dataset_Validate[[which(names(Dataset_Validate)==target2)]]
  data.real2.central=StMoMoData(data.validate2,series=Series,type="central")
  data.real2.initial=central2initial(data.real2.central)
  qxt2.real=data.real2.initial$Dxt/data.real2.initial$Ext

  SSE_V.1=sum((log(qxt1.real/(1-qxt1.real))-log(Obj$qxt1/(1-Obj$qxt1)))^2)
  SSE_V.2=sum((log(qxt2.real/(1-qxt2.real))-log(Obj$qxt2/(1-Obj$qxt2)))^2)
  SSE_V=SSE_V.1+SSE_V.2

  Out=data.frame(label1=target1,label2=target2,Deltat=Deltat,SSE_V=SSE_V)
  return(Out)
}

#' relevantscale between K(t) and k_j(t) for models
#' relevantscale is the ratio bewteen mean_delta_K(t) and mean_delta_k_j(t)
#'
#' @param Obj
#'
ObtainRelevantScale=function(Obj){
  mean_delta_Kt=mean(Obj$Kt[-1]-Obj$Kt[-length(Obj$Kt)])
  mean_delta_k1t=mean(Obj$kt_1[-1]-Obj$kt_1[-length(Obj$kt_1)])
  mean_delta_k2t=mean(Obj$kt_2[-1]-Obj$kt_2[-length(Obj$kt_2)])
  RS1=abs(mean_delta_k1t/mean_delta_Kt)
  RS2=abs(mean_delta_k1t/mean_delta_Kt)
  Out=data.frame(label1=Obj$label1,label2=Obj$label2,Deltat=Obj$Deltat,RS1=RS1,RS2=RS2,RS=mean(c(RS1,RS2)))
  return(Out)
}



## RankAvg Validation SSE ####
#' PerformanceValidationCBD1withDelta
#'
#' @param Obj bivarite
#' @param Series female / male
#' @param Dataset_Validate validate set
#' @import StMoMo
#' @return performance with delta
#'
PerformanceValidationCBD1withDelta=function(Obj,Series,Dataset_Validate){
  Seires=Series
  Deltat=Obj$delta
  target1=Obj$label1
  target2=Obj$label2

  data.validate1=Dataset_Validate[[which(names(Dataset_Validate)==target1)]]
  data.real1.central=StMoMoData(data.validate1,series=Series,type="central")
  data.real1.initial=central2initial(data.real1.central)
  qxt1.real=data.real1.initial$Dxt/data.real1.initial$Ext

  data.validate2=Dataset_Validate[[which(names(Dataset_Validate)==target2)]]
  data.real2.central=StMoMoData(data.validate2,series=Series,type="central")
  data.real2.initial=central2initial(data.real2.central)
  qxt2.real=data.real2.initial$Dxt/data.real2.initial$Ext

  SSE_V.1=sum((log(qxt1.real/(1-qxt1.real))-log(Obj$qxt1/(1-Obj$qxt1)))^2)
  SSE_V.2=sum((log(qxt2.real/(1-qxt2.real))-log(Obj$qxt2/(1-Obj$qxt2)))^2)
  SSE_V=SSE_V.1+SSE_V.2
  Out=data.frame(label1=target1,label2=target2,Deltat=Deltat,SSE_V=SSE_V,SSE_V.1=SSE_V.1,SSE_V.2=SSE_V.2)
  return(Out)
}



#' Choose best group size
#'
#' @param Obj obj
#' @param Series female / male
#' @param Dataset_Validate Valid
#' @import StMoMo
#' @return the best size
#'
ChooseGroupSize=function(Obj,Series, Dataset_Validate, ModelInfo,ValidateBivarateCBD1withDelta){
  Series=Series
  target=as.character(Obj$target)[1]
  reference=as.character(Obj$reference)
  data.validate1=Dataset_Validate[[which(names(Dataset_Validate)==target)]]
  data.real1.central=StMoMoData(data.validate1,series=Series,type="central")
  data.real1.initial=central2initial(data.real1.central)
  qxt1.real=data.real1.initial$Dxt/data.real1.initial$Ext
  SSE_log=vector()
  for(Groupsize in seq_len(length(reference))){
    RefenceChoice=reference[1:Groupsize]
    Delta=Obj$delta[1:Groupsize]
    logitqxt_fcast=0
    for(i in 1:length(RefenceChoice)){
      Index1=which(ModelInfo$target1==target & ModelInfo$target2==RefenceChoice[i] & ModelInfo$delta==Delta[i])
      if(length(Index1)>0){
        qxt=ValidateBivarateCBD1withDelta[[Index1]]$qxt1
        Forecast=log(qxt/(1-qxt))
        logitqxt_fcast=logitqxt_fcast+Forecast
      }
      Index2=which(ModelInfo$target2==target & ModelInfo$target1==RefenceChoice[i] & ModelInfo$delta==-Delta[i])
      if(length(Index2)>0){
        qxt=ValidateBivarateCBD1withDelta[[Index2]]$qxt2
        Forecast=log(qxt/(1-qxt))
        logitqxt_fcast=logitqxt_fcast+Forecast
      }
    }
    logitqxt_fcast=logitqxt_fcast/(length(RefenceChoice))
    SSE_log[Groupsize]=sum((log(qxt1.real/(1-qxt1.real))-logitqxt_fcast)^2)
  }
  GroupSizeChoice=which(SSE_log==min(SSE_log))
  Out=list(target=target,RefGroupSizeChoice=GroupSizeChoice,OptimSSE=min(SSE_log),SSE_V=SSE_log,
           RefreenceGroup=reference[1:GroupSizeChoice])
  return(Out)
}



#' PerformanceMeasure_RankAvg_SSE
#'
#' @param Obj obj
#' @param Series male / female
#' @param Dataset_Test test set
#' @param PredictionBivarateCBD1withDelta prediction
#'
#' @return Rand_Avg SSE
PerformanceMeasure_CBD_RankAvg_SSE=function(Obj,Series,Dataset_Test,PredictionBivarateCBD1withDelta){
  Series=Series
  target=Obj$target
  OptimSize=Obj$RefGroupSizeChoice
  RefreenceGroup=Obj$RefreenceGroup
  data.test=Dataset_Test[[which(names(Dataset_Test)==target)]]
  data.real.central=StMoMoData(data.test,series=Series,type="central")
  data.real.initial=central2initial(data.real.central)
  qxt.real=data.real.initial$Dxt/data.real.initial$Ext
  logit_qxt.real=log(qxt.real/(1-qxt.real))
  logitqxt_fcast=0
  S=0
  for(i in 1:length(PredictionBivarateCBD1withDelta)){
    if(PredictionBivarateCBD1withDelta[[i]]$label1==target & PredictionBivarateCBD1withDelta[[i]]$label2 %in% RefreenceGroup){
      qxt=PredictionBivarateCBD1withDelta[[i]]$qxt1
      Forecast=log(qxt/(1-qxt))
      logitqxt_fcast=logitqxt_fcast+Forecast
      S=S+1
    }
    if(PredictionBivarateCBD1withDelta[[i]]$label2==target & PredictionBivarateCBD1withDelta[[i]]$label1 %in% RefreenceGroup){
      qxt=PredictionBivarateCBD1withDelta[[i]]$qxt2
      Forecast=log(qxt/(1-qxt))
      logitqxt_fcast=logitqxt_fcast+Forecast
      S=S+1
    }
  }
  logitqxt_fcast=logitqxt_fcast/(length(RefreenceGroup))
  pred = exp(logitqxt_fcast) / (1+exp(logitqxt_fcast))
  #SSE=sum((logit_qxt.real-logitqxt_fcast)^2)
  Out=list(target=target,Series=Series,RefGroupSizeChoice=length(RefreenceGroup),
           forecast = pred)
  # Out=list(target=target,Series=Series,RefGroupSizeChoice=length(RefreenceGroup),SSE=SSE,
  #          forecast = pred)
  return(Out)
}


