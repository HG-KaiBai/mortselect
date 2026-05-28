#' Random sampling data
#'
#' @param Target target country
#' @param Size size
#' @param Series female / male
#' @param Dictionary country names
#' @param Training.rate training rate
#'
#' @return the random sampled data
RandomSampleData.training=function(Target,Size,Series, Dictionary,Training.rate){
  Series=Series
  target.index=which(Dictionary$population==Target$population & Dictionary$age==Target$age)
  Dictionary.temp=c(1:length(Dictionary$population))[-target.index]
  sample.index=sample(Dictionary.temp, Size, replace = FALSE, prob = NULL)
  Index=c(target.index,sample.index)
  SampleData=Training.rate[Index,]
  Out=list(Target=Target,Series=Series,Size=Size,target.index=target.index,SampleIndex=Index,SampleData=SampleData)
  return(Out)
}


#' Calculate Distance
#'
#' @param Data data
#' @param Target target country
#'
#' @return The distance
CalculateDistance=function(Data,Target){
  Distance=vector()
  s=1
  for(index1 in 1:nrow(Data)){
    if(rownames(Data)[index1] != Target){
      # find data for Age1
      RealMortality.Target=(Data[which(rownames(Data)==Target),])
      RealMortality.index1=(Data[index1,])

      #Calculate distance
      dist_temp=var(RealMortality.Target-RealMortality.index1)
      Distance[s]=dist_temp
      names(Distance)[s]=rownames(Data)[index1]
      s=s+1
    }
  }
  return(Distance)
}


#' Select age
#'
#' @param Data data
#' @param Target target country
#' @param K hyperpar K
#'
#' @return the age and distances
AgeSelectionK=function(Data,Target,K){
  if(K==0){
    AgeSelection=Target
    DistanceMeasure=NULL
  }else{
    DistanceMeasure=CalculateDistance(Data,Target)
    AgeSelection=c(Target,names(sort(DistanceMeasure))[1:K])
  }
  Out=list(AgeSelection=AgeSelection,DistanceMeasure=DistanceMeasure)
  return(Out)
}

.BCV_4_multi_generic <- function(TargetPopulation, TargetAge, Series, K, RealMortality) {
  Target <- paste(TargetPopulation, TargetAge, sep = ".")
  n_years <- ncol(RealMortality)
  n_folds <- min(4, n_years)
  fold_id <- cut(seq_len(n_years), breaks = n_folds, labels = FALSE)
  sse <- vapply(seq_len(n_folds), function(fold) {
    validate_idx <- which(fold_id == fold)
    model_idx <- setdiff(seq_len(n_years), validate_idx)
    data.model <- log(RealMortality[, model_idx, drop = FALSE])
    data.validate <- log(RealMortality[, validate_idx, drop = FALSE])
    AgeBandSelection <- AgeSelectionK(data.model, Target, K)$AgeSelection
    data.model.knn <- data.model[which(rownames(data.model) %in% AgeBandSelection), , drop = FALSE]
    h <- length(validate_idx)

    if (nrow(data.model.knn) == 1) {
      ax <- mean(data.model.knn)
      kt <- as.numeric(data.model.knn - ax)
      kt_f <- rwf(kt, h = h, drift = TRUE)$mean
      pred <- matrix(ax + kt_f, nrow = 1)
    } else {
      ax <- rowMeans(data.model.knn)
      residual <- sweep(data.model.knn, 1, ax, "-")
      svd.mx <- svd(residual)
      bx <- svd.mx$u[, 1] / sum(svd.mx$u[, 1])
      kt <- svd.mx$d[1] * svd.mx$v[, 1] * sum(svd.mx$u[, 1])
      kt_f <- rwf(kt, h = h, drift = TRUE)$mean
      pred_all <- ax + outer(bx, kt_f)
      row.names(pred_all) <- rownames(data.model.knn)
      pred <- pred_all[which(rownames(pred_all) %in% Target), , drop = FALSE]
    }
    sum((data.validate[which(rownames(data.validate) %in% Target), , drop = FALSE] - pred)^2)
  }, numeric(1))

  data.frame(
    TargetPopulation = TargetPopulation,
    TargetAge = TargetAge,
    Target = Target,
    Series = Series,
    K = K,
    SSE.BCV = mean(sse),
    AllData = I(list(RealMortality))
  )
}


#' Choose Optimal hyperparameter K via blocked-cross-validation (BCV)
#'
#' @param TargetPopulation  TargetPopulation
#' @param TargetAge  TargetAge
#' @param Series female / male
#' @param K hyper K
#' @param Training.rate training rate
#'
#' @return sse of corss validation
BCV_4withAllData=function(TargetPopulation,TargetAge,Series,K, Training.rate){
  Series=Series
  # get information
  TargetPopulation=TargetPopulation
  TargetAge=TargetAge
  Target=paste(TargetPopulation,TargetAge,sep=".")

  AllData=Training.rate

  RealMortality=AllData
  if (ncol(RealMortality) != 33) {
    return(.BCV_4_multi_generic(TargetPopulation, TargetAge, Series, K, RealMortality))
  }

  # validation 1 ####
  data.model=log(RealMortality[,c(14:ncol(RealMortality))])
  data.validate=log(RealMortality[,c(1:8)])
  AgeBandSelection=AgeSelectionK(data.model,Target,K)$AgeSelection
  data.model.knn=data.model[which(rownames(data.model) %in% AgeBandSelection),]
  if(length(AgeBandSelection)==1)
    ## k=0, no borrowing information
  {
    ax=mean(data.model.knn)
    Residual_c=data.model.knn-ax
    svd.mx=svd(t(Residual_c))
    bx=svd.mx$u[,1]/sum(svd.mx$u[,1])
    kt=svd.mx$d[1]*svd.mx$v[,1]*sum(svd.mx$u[,1])
    ts.model=Arima(kt,order=c(0,1,0),include.drift=TRUE)
    kt.validate=kt[1]-ts.model$coef*c(13:6)
    logrates.validate=ax+outer(bx,kt.validate)
    SSE.V1=sum((data.validate[which(rownames(data.validate) %in% Target),]-logrates.validate)^2)
  }else{
    ax=rowMeans(data.model.knn)
    Residual_c=sweep(data.model.knn,1,ax,"-")
    svd.mx=svd(Residual_c)
    bx=svd.mx$u[,1]/sum(svd.mx$u[,1])
    kt=svd.mx$d[1]*svd.mx$v[,1]*sum(svd.mx$u[,1])
    ts.model=Arima(kt,order=c(0,1,0),include.drift=TRUE)
    kt.validate=kt[1]-ts.model$coef*c(13:6)
    logrates.validate=ax+outer(bx,kt.validate)
    row.names(logrates.validate)=rownames(data.model.knn)
    SSE.V1=sum((data.validate[which(rownames(data.validate) %in% Target),]-logrates.validate[which(rownames(logrates.validate) %in% Target),])^2)
  }
  # validation 2 ####
  data.model=log(RealMortality[,c(c(1:3),c(22:ncol(RealMortality)))])
  data.validate=log(RealMortality[,c(9:16)])
  AgeBandSelection=AgeSelectionK(data.model,Target,K)$AgeSelection
  data.model.knn=data.model[which(rownames(data.model) %in% AgeBandSelection),]
  if(length(AgeBandSelection)==1)
    ## k=0, no borrowing information
  {
    ax=mean(data.model.knn)
    Residual_c=data.model.knn-ax
    svd.mx=svd(t(Residual_c))
    bx=svd.mx$u[,1]/sum(svd.mx$u[,1])
    kt=svd.mx$d[1]*svd.mx$v[,1]*sum(svd.mx$u[,1])
    kt.adjust=c(kt[1:3],rep(NA,18),kt[4:length(kt)])
    ts.model=Arima(kt.adjust,order=c(0,1,0),include.drift=TRUE)
    kt.validate=kt[3]+ts.model$coef*c(6:13)
    logrates.validate=ax+outer(bx,kt.validate)
    SSE.V2=sum((data.validate[which(rownames(data.validate) %in% Target),]-logrates.validate)^2)
  }else{
    ax=rowMeans(data.model.knn)
    Residual_c=sweep(data.model.knn,1,ax,"-")
    svd.mx=svd(Residual_c)
    bx=svd.mx$u[,1]/sum(svd.mx$u[,1])
    kt=svd.mx$d[1]*svd.mx$v[,1]*sum(svd.mx$u[,1])
    kt.adjust=c(kt[1:3],rep(NA,18),kt[4:length(kt)])
    ts.model=Arima(kt.adjust,order=c(0,1,0),include.drift=TRUE)
    kt.validate=kt[3]+ts.model$coef*c(6:13)
    logrates.validate=ax+outer(bx,kt.validate)
    row.names(logrates.validate)=rownames(data.model.knn)
    SSE.V2=sum((data.validate[which(rownames(data.validate) %in% Target),]-logrates.validate[which(rownames(logrates.validate) %in% Target),])^2)
  }

  # validation 3 ####
  data.model=log(RealMortality[,c(c(1:11),c(30:ncol(RealMortality)))])
  data.validate=log(RealMortality[,c(17:24)])
  AgeBandSelection=AgeSelectionK(data.model,Target,K)$AgeSelection
  data.model.knn=data.model[which(rownames(data.model) %in% AgeBandSelection),]
  if(length(AgeBandSelection)==1)
    ## k=0, no borrowing information
  {
    ax=mean(data.model.knn)
    Residual_c=data.model.knn-ax
    svd.mx=svd(t(Residual_c))
    bx=svd.mx$u[,1]/sum(svd.mx$u[,1])
    kt=svd.mx$d[1]*svd.mx$v[,1]*sum(svd.mx$u[,1])
    kt.adjust=c(kt[1:11],rep(NA,18),kt[12:length(kt)])
    ts.model=Arima(kt.adjust,order=c(0,1,0),include.drift=TRUE)
    kt.validate=kt[11]+ts.model$coef*c(6:13)
    logrates.validate=ax+outer(bx,kt.validate)
    SSE.V3=sum((data.validate[which(rownames(data.validate) %in% Target),]-logrates.validate)^2)
  }else{
    ax=rowMeans(data.model.knn)
    Residual_c=sweep(data.model.knn,1,ax,"-")
    svd.mx=svd(Residual_c)
    bx=svd.mx$u[,1]/sum(svd.mx$u[,1])
    kt=svd.mx$d[1]*svd.mx$v[,1]*sum(svd.mx$u[,1])
    kt.adjust=c(kt[1:11],rep(NA,18),kt[12:length(kt)])
    ts.model=Arima(kt.adjust,order=c(0,1,0),include.drift=TRUE)
    kt.validate=kt[11]+ts.model$coef*c(6:13)
    logrates.validate=ax+outer(bx,kt.validate)
    row.names(logrates.validate)=rownames(data.model.knn)
    SSE.V3=sum((data.validate[which(rownames(data.validate) %in% Target),]-logrates.validate[which(rownames(logrates.validate) %in% Target),])^2)
  }

  # validation 4 ####
  data.model=log(RealMortality[,c(1:19)])
  data.validate=log(RealMortality[,c(25:32)])
  AgeBandSelection=AgeSelectionK(data.model,Target,K)$AgeSelection
  data.model.knn=data.model[which(rownames(data.model) %in% AgeBandSelection),]

  if(length(AgeBandSelection)==1)
    ## k=0, no borrowing information
  {
    ax=mean(data.model.knn)
    Residual_c=data.model.knn-ax
    svd.mx=svd(t(Residual_c))
    bx=svd.mx$u[,1]/sum(svd.mx$u[,1])
    kt=svd.mx$d[1]*svd.mx$v[,1]*sum(svd.mx$u[,1])
    ts.model=Arima(kt,order=c(0,1,0),include.drift=TRUE)
    kt.validate=kt[length(kt)]+ts.model$coef*c(6:13)
    logrates.validate=ax+outer(bx,kt.validate)
    SSE.V4=sum((data.validate[which(rownames(data.validate) %in% Target),]-logrates.validate)^2)
  }else{
    ax=rowMeans(data.model.knn)
    Residual_c=sweep(data.model.knn,1,ax,"-")
    svd.mx=svd(Residual_c)
    bx=svd.mx$u[,1]/sum(svd.mx$u[,1])
    kt=svd.mx$d[1]*svd.mx$v[,1]*sum(svd.mx$u[,1])
    ts.model=Arima(kt,order=c(0,1,0),include.drift=TRUE)
    kt.validate=kt[length(kt)]+ts.model$coef*c(6:13)
    logrates.validate=ax+outer(bx,kt.validate)
    row.names(logrates.validate)=rownames(data.model.knn)
    SSE.V4=sum((data.validate[which(rownames(data.validate) %in% Target),]-logrates.validate[which(rownames(logrates.validate) %in% Target),])^2)
  }
  ## averaging ####
  SSE.BCV=mean(c(SSE.V1,SSE.V2,SSE.V3,SSE.V4))
  Out=data.frame(TargetPopulation=TargetPopulation,TargetAge=TargetAge,Target=Target,Series=Series,K=K,SSE.BCV=SSE.BCV,AllData=AllData)
  return(Out)
}


#' the best K
#'
#' @param TargetPopulation TargetPopulation
#' @param TargetAge TargetAge
#' @param Series female / male
#' @param UpperLimit upper limite of bandwidth
#' @param Training.rate trainging rate
#'
#' @return The best K
ChoseK_BCV_multi=function(TargetPopulation,TargetAge,Series,UpperLimit, Training.rate)
  #UpperLimit regulate the largest allowed Bandwidth
{
  Series=Series
  TargetPopulation=TargetPopulation
  TargetAge=TargetAge
  UpperLimit=UpperLimit

  target=paste(TargetPopulation,TargetAge,sep=".")
  message(target)
  ValidationSSE=data.frame()

  for(k in 0:UpperLimit)
  {
    temp=BCV_4withAllData(TargetPopulation,TargetAge,Series,k, Training.rate)
    ValidationSSE=rbind(ValidationSSE,temp)
  }
  ChosenK=ValidationSSE$K[which(ValidationSSE$SSE.BCV==min(ValidationSSE$SSE.BCV))[1]]
  message("Chosen Optim K is ",ChosenK)
  Out=data.frame(TargetPopulation=TargetPopulation,TargetAge=TargetAge,ChosenK=ChosenK,optimSSE_V=min(ValidationSSE$SSE.BCV))
  return(Out)
}




#' performance_test_multi
#'
#' @param TargetPopulation TargetPopulation
#' @param TargetAge TargetAge
#' @param Series  female / male
#' @param Dataset_Test test set
#' @param Dataset_training train set
#' @param Smooth.ChosenK the best K
#'
#' @return
#' @export
#'
#' @examples
performance_test_multi=function(TargetPopulation,TargetAge,Series,Training.rate,Dataset_Test,Dataset_training, Smooth.ChosenK,Dictionary,h = 8){
  Series=Series
  TargetPopulation=TargetPopulation
  TargetAge=TargetAge
  Target=paste(TargetPopulation,TargetAge,sep=".")

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


  ## optim chosen k for the target
  ChosenKvalue=Smooth.ChosenK$smooth.K[which(Smooth.ChosenK$TargetPopulation==TargetPopulation & Smooth.ChosenK$TargetAge ==TargetAge )]
  ## find the optim chosen ages
  if(ChosenKvalue==0)
    ## k=0, no borrowing information
  {
    DomesticK=0;ForeignK=0

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
    sse_AgeBand=sum((log(RealMortality)-logrates_fcast)^2)
    fcast.logrates=logrates_fcast
  }else
  {
    ## k>0
    #message(k)
    population_temp=Dataset_training[[which(names(Dataset_training)==TargetPopulation)]]
    data_temp=population_temp$rate[[which(names(population_temp$rate)==Series)]]

    AgeBandSelection=AgeSelectionK(log(Training.rate),Target,ChosenKvalue)$AgeSelection # The AgeBand for targetage with a given bandwidth
    DomesticK=sum(Dictionary$population[which(Dictionary$target %in% AgeBandSelection[-1])]==TargetPopulation)
    ForeignK=ChosenKvalue-DomesticK
    #merge training mortality data

    data.combined=Training.rate[which(rownames(Training.rate) %in% AgeBandSelection),]
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
    #sse_AgeBand=sum((log(RealMortality)-logrates_fcast[which(rownames(logrates_fcast) %in% Target),])^2)
    fcast.logrates=logrates_fcast[which(rownames(logrates_fcast) %in% Target),]
  }
  #SE=(log(RealMortality)-fcast.logrates)^2
  Out=list(TargetPopulation=TargetPopulation,TargetAge=TargetAge,Series=Series,ChosenKvalue=ChosenKvalue,DomesticK=DomesticK,ForeignK=ForeignK,
           logrates_fcast=exp(fcast.logrates))
  # Out=list(TargetPopulation=TargetPopulation,TargetAge=TargetAge,Series=Series,sse_lc=sse_lc,sse_AgeBand=sse_AgeBand,ChosenKvalue=ChosenKvalue,DomesticK=DomesticK,ForeignK=ForeignK,
  #          logrates_fcast=fcast.logrates,SE=SE)
  return(Out)
}


