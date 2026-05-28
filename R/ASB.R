#' Define AgeBand
#'
#' @param TargetAge target age
#' @param BandWidth bandwidth
#'
#' @return the age band
AgeBand=function(TargetAge,BandWidth){
  MinAge=max(TargetAge-BandWidth,0)
  MaxAge=min(TargetAge+BandWidth,100)
  Out=list(TargetAge=TargetAge,BandWidth=BandWidth,MinAge=MinAge,MaxAge=MaxAge,AgeBand=c(MinAge:MaxAge))
  return(Out)
}

.BCV_4_ageband_generic <- function(TargetPopulation, TargetAge, Series, K, RealMortality) {
  AgeBandSelection <- AgeBand(TargetAge, K)$AgeBand
  n_years <- ncol(RealMortality)
  n_folds <- min(4, n_years)
  fold_id <- cut(seq_len(n_years), breaks = n_folds, labels = FALSE)
  sse <- vapply(seq_len(n_folds), function(fold) {
    validate_idx <- which(fold_id == fold)
    model_idx <- setdiff(seq_len(n_years), validate_idx)
    selected <- which(rownames(RealMortality) %in% AgeBandSelection)
    target <- which(rownames(RealMortality) %in% TargetAge)

    log_model <- log(RealMortality[selected, model_idx, drop = FALSE])
    log_validate <- log(RealMortality[target, validate_idx, drop = FALSE])
    h <- length(validate_idx)

    if (nrow(log_model) == 1) {
      ax <- mean(log_model)
      kt <- as.numeric(log_model - ax)
      kt_f <- rwf(kt, h = h, drift = TRUE)$mean
      pred <- matrix(ax + kt_f, nrow = 1)
    } else {
      ax <- rowMeans(log_model)
      residual <- sweep(log_model, 1, ax, "-")
      svd.mx <- svd(residual)
      bx <- svd.mx$u[, 1] / sum(svd.mx$u[, 1])
      kt <- svd.mx$d[1] * svd.mx$v[, 1] * sum(svd.mx$u[, 1])
      kt_f <- rwf(kt, h = h, drift = TRUE)$mean
      pred_all <- ax + outer(bx, kt_f)
      rownames(pred_all) <- rownames(log_model)
      pred <- pred_all[which(rownames(pred_all) %in% TargetAge), , drop = FALSE]
    }
    sum((log_validate - pred)^2)
  }, numeric(1))

  data.frame(
    TargetPopulation = TargetPopulation,
    TargetAge = TargetAge,
    Series = Series,
    K = K,
    SSE.BCV = mean(sse)
  )
}



#' Choose Optimal BandWidth via blocked-cross-validation (BCV)
#'
#' @param TargetPopulation target country
#' @param TargetAge target age
#' @param Series female / male
#' @param K the bandwidth
#' @param Dataset_training traing set
#'
#' @return 4 fold cross validation
BCV_4=function(TargetPopulation,TargetAge,Series,K, Dataset_training){
  Series=Series
  TargetPopulation=TargetPopulation
  TargetAge=TargetAge
  RealMortality=Dataset_training[[which(names(Dataset_training)==TargetPopulation)]]$rate[[which(names(Dataset_training[[which(names(Dataset_training)==TargetPopulation)]]$rate)==Series)]]
  if (nrow(RealMortality) >= 101) rownames(RealMortality)[101]=100

  AgeBandSelection=AgeBand(TargetAge,K)$AgeBand # The AgeBand for targetage with a given bandwidth
  if (ncol(RealMortality) != 33) {
    return(.BCV_4_ageband_generic(TargetPopulation, TargetAge, Series, K, RealMortality))
  }

  # validation 1 ####
  data.model=log(RealMortality[which(rownames(RealMortality) %in% AgeBandSelection),c(14:ncol(RealMortality))])
  data.validate=log(RealMortality[which(rownames(RealMortality) %in% AgeBandSelection),c(1:8)])
  if(length(AgeBandSelection)==1)
    ## k=0, no borrowing information
  {
    ax=mean(data.model)
    Residual_c=data.model-ax
    svd.mx=svd(t(Residual_c))
    bx=svd.mx$u[,1]/sum(svd.mx$u[,1])
    kt=svd.mx$d[1]*svd.mx$v[,1]*sum(svd.mx$u[,1])
    ts.model=Arima(kt,order=c(0,1,0),include.drift=TRUE)
    kt.validate=kt[1]-ts.model$coef*c(13:6)
    logrates.validate=ax+outer(bx,kt.validate)
    SSE.V1=sum((data.validate-logrates.validate)^2)
  }
  else
  {
    ax=rowMeans(data.model)
    Residual_c=sweep(data.model,1,ax,"-")
    svd.mx=svd(Residual_c)
    bx=svd.mx$u[,1]/sum(svd.mx$u[,1])
    kt=svd.mx$d[1]*svd.mx$v[,1]*sum(svd.mx$u[,1])
    ts.model=Arima(kt,order=c(0,1,0),include.drift=TRUE)
    kt.validate=kt[1]-ts.model$coef*c(13:6)
    logrates.validate=ax+outer(bx,kt.validate)
    row.names(logrates.validate)=AgeBandSelection
    SSE.V1=sum((data.validate[which(rownames(data.validate) %in% TargetAge),]-logrates.validate[which(rownames(logrates.validate) %in% TargetAge),])^2)
  }
  # validation 2 ####
  data.model=log(RealMortality[which(rownames(RealMortality) %in% AgeBandSelection),c(c(1:3),c(22:ncol(RealMortality)))])
  data.validate=log(RealMortality[which(rownames(RealMortality) %in% AgeBandSelection),c(9:16)])
  if(length(AgeBandSelection)==1)
    ## k=0, no borrowing information
  {
    ax=mean(data.model)
    Residual_c=data.model-ax
    svd.mx=svd(t(Residual_c))
    bx=svd.mx$u[,1]/sum(svd.mx$u[,1])
    kt=svd.mx$d[1]*svd.mx$v[,1]*sum(svd.mx$u[,1])
    kt.adjust=c(kt[1:3],rep(NA,18),kt[4:length(kt)])
    ts.model=Arima(kt.adjust,order=c(0,1,0),include.drift=TRUE)
    kt.validate=kt[3]+ts.model$coef*c(6:13)
    logrates.validate=ax+outer(bx,kt.validate)
    SSE.V2=sum((data.validate-logrates.validate)^2)
  }
  else
  {
    ax=rowMeans(data.model)
    Residual_c=sweep(data.model,1,ax,"-")
    svd.mx=svd(Residual_c)
    bx=svd.mx$u[,1]/sum(svd.mx$u[,1])
    kt=svd.mx$d[1]*svd.mx$v[,1]*sum(svd.mx$u[,1])
    kt.adjust=c(kt[1:3],rep(NA,18),kt[4:length(kt)])
    ts.model=Arima(kt.adjust,order=c(0,1,0),include.drift=TRUE)
    kt.validate=kt[3]+ts.model$coef*c(6:13)
    logrates.validate=ax+outer(bx,kt.validate)
    row.names(logrates.validate)=AgeBandSelection
    SSE.V2=sum((data.validate[which(rownames(data.validate) %in% TargetAge),]-logrates.validate[which(rownames(logrates.validate) %in% TargetAge),])^2)
  }
  # validation 3 ####
  data.model=log(RealMortality[which(rownames(RealMortality) %in% AgeBandSelection),c(c(1:11),c(30:ncol(RealMortality)))])
  data.validate=log(RealMortality[which(rownames(RealMortality) %in% AgeBandSelection),c(17:24)])
  if(length(AgeBandSelection)==1)
    ## k=0, no borrowing information
  {
    ax=mean(data.model)
    Residual_c=data.model-ax
    svd.mx=svd(t(Residual_c))
    bx=svd.mx$u[,1]/sum(svd.mx$u[,1])
    kt=svd.mx$d[1]*svd.mx$v[,1]*sum(svd.mx$u[,1])
    kt.adjust=c(kt[1:11],rep(NA,18),kt[12:length(kt)])
    ts.model=Arima(kt.adjust,order=c(0,1,0),include.drift=TRUE)
    kt.validate=kt[11]+ts.model$coef*c(6:13)
    logrates.validate=ax+outer(bx,kt.validate)
    SSE.V3=sum((data.validate-logrates.validate)^2)
  }
  else
  {
    ax=rowMeans(data.model)
    Residual_c=sweep(data.model,1,ax,"-")
    svd.mx=svd(Residual_c)
    bx=svd.mx$u[,1]/sum(svd.mx$u[,1])
    kt=svd.mx$d[1]*svd.mx$v[,1]*sum(svd.mx$u[,1])
    kt.adjust=c(kt[1:11],rep(NA,18),kt[12:length(kt)])
    ts.model=Arima(kt.adjust,order=c(0,1,0),include.drift=TRUE)
    kt.validate=kt[11]+ts.model$coef*c(6:13)
    logrates.validate=ax+outer(bx,kt.validate)
    row.names(logrates.validate)=AgeBandSelection
    SSE.V3=sum((data.validate[which(rownames(data.validate) %in% TargetAge),]-logrates.validate[which(rownames(logrates.validate) %in% TargetAge),])^2)
  }
  # validation 4 ####
  data.model=log(RealMortality[which(rownames(RealMortality) %in% AgeBandSelection),c(1:19)])
  data.validate=log(RealMortality[which(rownames(RealMortality) %in% AgeBandSelection),c(25:32)])
  if(length(AgeBandSelection)==1)
    ## k=0, no borrowing information
  {
    ax=mean(data.model)
    Residual_c=data.model-ax
    svd.mx=svd(t(Residual_c))
    bx=svd.mx$u[,1]/sum(svd.mx$u[,1])
    kt=svd.mx$d[1]*svd.mx$v[,1]*sum(svd.mx$u[,1])
    ts.model=Arima(kt,order=c(0,1,0),include.drift=TRUE)
    kt.validate=kt[length(kt)]+ts.model$coef*c(6:13)
    logrates.validate=ax+outer(bx,kt.validate)
    SSE.V4=sum((data.validate-logrates.validate)^2)
  }
  else
  {
    ax=rowMeans(data.model)
    Residual_c=sweep(data.model,1,ax,"-")
    svd.mx=svd(Residual_c)
    bx=svd.mx$u[,1]/sum(svd.mx$u[,1])
    kt=svd.mx$d[1]*svd.mx$v[,1]*sum(svd.mx$u[,1])
    ts.model=Arima(kt,order=c(0,1,0),include.drift=TRUE)
    kt.validate=kt[length(kt)]+ts.model$coef*c(6:13)
    logrates.validate=ax+outer(bx,kt.validate)
    row.names(logrates.validate)=AgeBandSelection
    SSE.V4=sum((data.validate[which(rownames(data.validate) %in% TargetAge),]-logrates.validate[which(rownames(logrates.validate) %in% TargetAge),])^2)
  }
  ## averaging ####
  SSE.BCV=mean(c(SSE.V1,SSE.V2,SSE.V3,SSE.V4))
  Out=data.frame(TargetPopulation=TargetPopulation,TargetAge=TargetAge,Series=Series,K=K,SSE.BCV=SSE.BCV)
  return(Out)
}




#' UpperLimit regulate the largest allowed Bandwidth
#'
#' @param TargetPopulation target country
#' @param TargetAge target age
#' @param Series female / male
#' @param UpperLimit upperlimit bandwidth
#' @param Dataset_training traing set
#'
#' @return the best bandwidth
ChoseK_BCV=function(TargetPopulation,TargetAge,Series,UpperLimit, Dataset_training){
  Series=Series
  TargetPopulation=TargetPopulation
  TargetAge=TargetAge
  UpperLimit=UpperLimit

  target=paste(TargetPopulation,TargetAge,sep="")
  message(target)
  ValidationSSE=data.frame()

  for(k in 0:UpperLimit){
    temp=BCV_4(TargetPopulation,TargetAge,Series,k, Dataset_training)
    ValidationSSE=rbind(ValidationSSE,temp)
  }
  ChosenK=ValidationSSE$K[which(ValidationSSE$SSE.BCV==min(ValidationSSE$SSE.BCV))[1]]
  message("Chosen Optim K is ",ChosenK)
  ChosenAgeBand=AgeBand(TargetAge,ChosenK)
  Out=data.frame(TargetPopulation=TargetPopulation,TargetAge=TargetAge,ChosenK=ChosenK,optimSSE_V=min(ValidationSSE$SSE.BCV),
                 minAge=ChosenAgeBand$MinAge,maxAge=ChosenAgeBand$MaxAge)
  return(Out)
}


#' Smoothing.mean
#'
#' @param Obj Obj
#' @param h h
#' @import caTools
#'
#' @return smoothed mean
Smoothing.mean=function(Obj,h){
  temp=Obj$ChosenK
  smooth.K=round(runmean(temp,h,endrule = "mean"))
  return(smooth.K)
}



#' The SSE
#'
#' @param TargetPopulation target country
#' @param TargetAge the target age
#' @param Series female / male
#' @param Dataset_Test test set
#' @param Dataset_training training set.
#' @param Smooth.ChosenK The best K
#'
#' @return The SSE
performance_test=function(TargetPopulation,TargetAge,Series,Dataset_Test, Dataset_training, Smooth.ChosenK, h = 8){
  Series=Series
  TargetPopulation=TargetPopulation
  TargetAge=TargetAge

  data.test=Dataset_Test[[which(names(Dataset_Test)==TargetPopulation)]]$rate[[which(names(Dataset_Test[[which(names(Dataset_Test)==TargetPopulation)]]$rate)==Series)]]
  if (nrow(data.test) >= 101) rownames(data.test)[101]=100
  RealMortality=data.test[which(rownames(data.test) %in% TargetAge),]

  # lc sse
  data.train=Dataset_training[[which(names(Dataset_training)==TargetPopulation)]]
  if (nrow(data.train$rate[[Series]]) >= 101) {
    rownames(data.train$rate$total)[101]=rownames(data.train$rate$female)[101]=rownames(data.train$rate$male)[101]=100
  }
  ax=lca(data.train,Series)$ax
  log_mx=log(data.train$rate[[which(names(data.train$rate)==Series)]])
  Residual_c=sweep(log_mx,1,ax,"-")
  svd.mx=svd(Residual_c)
  bx=svd.mx$u[,1]/sum(svd.mx$u[,1])
  kt=svd.mx$d[1]*svd.mx$v[,1]*sum(svd.mx$u[,1])

  kt_f=rwf(kt,h=h,drift=TRUE)
  logrates_fcast=ax+outer(bx,kt_f$mean)
  rownames(logrates_fcast)=rownames(data.train$rate$total)
  RealMortality=data.test[which(rownames(data.test) %in% TargetAge),]
  #sse_lc=sum((log(RealMortality)-logrates_fcast[which(rownames(logrates_fcast) %in% TargetAge),])^2)


  ## optim chosen smmothed k for the target
  ChosenKvalue=Smooth.ChosenK$smooth.K[which(Smooth.ChosenK$TargetPopulation==TargetPopulation & Smooth.ChosenK$TargetAge ==TargetAge )]
  ## find the optim chosen ages
  if(ChosenKvalue==0)
    ## k=0, no borrowing information
  {
    population_temp=Dataset_training[[which(names(Dataset_training)==TargetPopulation)]]
    data_temp=population_temp$rate[[which(names(population_temp$rate)==Series)]]
    Mortalitydata_temp=data_temp[which(rownames(data_temp) %in% TargetAge),]
    data.combined=Mortalitydata_temp
    log_mx=log(data.combined)
    ax=mean(log_mx)
    Residual_c=log_mx-ax
    svd.mx=svd(t(Residual_c))
    bx=svd.mx$u[,1]/sum(svd.mx$u[,1])
    kt=svd.mx$d[1]*svd.mx$v[,1]*sum(svd.mx$u[,1])
    kt_f=rwf(kt,h=h,drift=TRUE)
    logrates_fcast=ax+outer(bx,kt_f$mean)
    rownames(logrates_fcast)=rownames(data.combined)
    #sse_AgeBand=sum((log(RealMortality)-logrates_fcast)^2)
    fcast.logrates=logrates_fcast
  }else{
    ## k>0
    #message(k)
    AgeBandSelection=AgeBand(TargetAge,ChosenKvalue)$AgeBand # The AgeBand for targetage with a given bandwidth
    #merge training mortality data
    data.combined=data.frame()

    population_temp=Dataset_training[[which(names(Dataset_training)==TargetPopulation)]]
    data_temp=extract.ages(population_temp,AgeBandSelection,combine.upper = FALSE)
    data.combined=data_temp$rate[[which(names(population_temp$rate)==Series)]]
    #colnames(data.combined)=population_temp$year

    log_mx=log(data.combined)
    ax=rowMeans(log_mx)
    Residual_c=sweep(log_mx,1,ax,"-")
    svd.mx=svd(Residual_c)
    bx=svd.mx$u[,1]/sum(svd.mx$u[,1])
    kt=svd.mx$d[1]*svd.mx$v[,1]*sum(svd.mx$u[,1])
    kt_f=rwf(kt,h=h,drift=TRUE)
    logrates_fcast=ax+outer(bx,kt_f$mean)
    rownames(logrates_fcast)=rownames(data.combined)
    #sse_AgeBand=sum((log(RealMortality)-logrates_fcast[which(rownames(logrates_fcast) %in% TargetAge),])^2)
    fcast.logrates=logrates_fcast[which(rownames(logrates_fcast) %in% TargetAge),]
  }
  Out=list(TargetPopulation=TargetPopulation,TargetAge=TargetAge,Series=Series,ChosenKvalue=ChosenKvalue,
           logrates_fcast= exp(fcast.logrates))
  # Out=list(TargetPopulation=TargetPopulation,TargetAge=TargetAge,Series=Series,sse_lc=sse_lc,sse_AgeBand=sse_AgeBand,ChosenKvalue=ChosenKvalue,
  #          logrates_fcast=fcast.logrates)
  return(Out)
}











