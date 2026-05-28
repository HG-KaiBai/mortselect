
#' Forecast with the deletion-substitution-addition algorithm
#'
#' Implements the DSA population-selection workflow of Diao, Meng, and Weng
#' (2020) with the augmented common factor fitting routine used in the
#' accompanying research code. Call [prepare_mortselect_data()] first to set
#' age ranges, year windows, exclusions, and missing-rate handling.
#'
#' @param prepared A prepared mortality dataset returned by
#'   [prepare_mortselect_data()]. For DSA, prepare with `relabel_sex = TRUE`.
#' @param COG Maximum selected group size.
#' @param N Maximum number of country-specific components.
#' @param LOSS Selection loss. One of `"SSE"`, `"MAE"`, `"MSPE"`, or `"MAPE"`.
#' @param t Number of years to forecast.
#' @param save_file Optional `.RData` path. If `NULL`, no file is written.
#'
#' @return A named list of forecasts. Each element contains `label` and
#'   `Prediction`, where `Prediction` is a mortality-rate matrix.
#' @export
#'
forecast_dsa = function(prepared,
                        COG = 30,
                        N = 5,
                        LOSS = "SSE",
                        t = 15,
                        save_file = NULL)
{
  validate_prepared_mortselect_data(
    prepared,
    series = NULL,
    min_populations = 2,
    require_relabel_sex = TRUE
  )

  Training_relabel = prepared$full
  Candidate_relabel = prepared$learning
  Dataset_Validate_relabel = prepared$validation
  Dataset_Test_relabel = prepared$test
  if (length(Candidate_relabel) < 2) {
    stop("At least two sex-specific populations are required for DSA selection.", call. = FALSE)
  }
  COG <- min(COG, length(Candidate_relabel))
  ########################################### model Fitting
  Labels = names(Dataset_Test_relabel)
  Ans = list()
  for (country in Labels){
    message("Current predicting ", country)
    CV_example = DSA_V5(Dataset_Validate_relabel[[country]], Candidate_relabel, COG, N, LOSS)
    best_idx = which.min(CV_example$LOSS)
    Best_group = CV_example$Best[[best_idx]]
    group = Best_group[[1]]
    if(length(Best_group) > 1){
      for (i in 2:length(Best_group)){
        group = Groupdemogdata_V2(group,Best_group[[i]],
                                  Best_group[[i]]$year[1],Best_group[[i]]$year[length(Best_group[[i]]$year)])
      }
    }
    ### prediction
    predict = prefit_v2(group,Dataset_Test_relabel[[country]], (best_idx-1)%%N+1, Training_relabel, t)
    #Err = loss(predict)[[LOSS]]
    Predict = predict[[1]]
    rownames(Predict) = rownames(predict[[2]])
    last_year <- max(Training_relabel[[country]]$year)
    colnames(Predict) = c((last_year + 1):(last_year + t))
    Rate_predict = exp(Predict)
    Ans[[country]] = list(label = country, Group = Best_group, Prediction = Rate_predict)
    #Ans[[country]] = list(label = country, Error = Err, Prediction = Rate_predict)
  }
  if (!is.null(save_file)) {
    save(Ans, file = save_file)
  }
  Ans
}



#' Bivariate CBD model-based ensemble forecast
#'
#' Applies the CBD-type bivariate model-based ensemble described in the BMBE
#' paper. The default age range, 55--90, follows the senior-age CBD study.
#'
#' @param prepared A prepared mortality dataset returned by
#'   [prepare_mortselect_data()].
#' @param Series Mortality series to forecast.
#' @param t Number of years to forecast.
#' @inheritParams forecast_dsa
#' @return A named list of CBD ensemble forecast objects.
#' @export
#'
forecast_bmbe_cbd = function(prepared,
                             Series = "female",
                             t = 15,
                             save_file = NULL){
  validate_prepared_mortselect_data(prepared, series = Series, min_populations = 2)
  Training = prepared$full
  Candidate = prepared$learning
  Dataset_Validate = prepared$validation
  Dataset_Test = prepared$test
  validation_h <- ncol(Dataset_Validate[[1]]$rate[[Series]])
  if (length(Candidate) < 2) {
    stop("At least two populations are required for the bivariate CBD ensemble.", call. = FALSE)
  }

  ### Fit all bivariate models and obtain optim Delta ####
  BivarateCBD1withDelta_modeling=list() # record fitting results of all models
  ValidateBivarateCBD1withDelta=list() # record extrapolation on validation set for all models
  ValidationSSESummary=data.frame() # record validation SSE
  OptimDelta=data.frame() # record deltat choice with smallest SSE_V

  s=1
  for(Index1 in 1:(length(Candidate)-1)){
    for(Index2 in (Index1+1):length(Candidate)){
      ## indication current chosen populations
      message(names(Candidate)[c(Index1,Index2)])
      Test=data.frame()
      for(Deltat in -10:10){
        F=FitBivariateCBD1withDelta_modeling(Index1,Index2,Series,Deltat, Candidate)
        BivarateCBD1withDelta_modeling[[s]]=F
        E=ExtrapolateBivarateCBD1withDelta(F, h = validation_h)
        ValidateBivarateCBD1withDelta[[s]]=E
        temp=ValidationSSE_CBD1withDelta(E,Series, Dataset_Validate)
        ValidationSSESummary=rbind(ValidationSSESummary,temp)
        Test=rbind(Test,temp)
        s=s+1
      }
      optim=Test[which(Test$SSE_V==min(Test$SSE_V)),]
      OptimDelta=rbind(OptimDelta,optim)
    }
  }
  ####################################################
  BivarateCBD1withDelta=list() # record fitting results of all models
  PredictionBivarateCBD1withDelta=list() # record extrapolation on validation set for all models
  RelevantScale=data.frame() # record relevantscale between K(t) and k_j(t) for all models

  s=1
  for(Index1 in 1:(length(Candidate)-1)){
    for(Index2 in (Index1+1):length(Candidate)){
      ## indication current chosen populations
      message(names(Candidate)[c(Index1,Index2)])
      Index=which(OptimDelta$label1== names(Candidate)[Index1] & OptimDelta$label2== names(Candidate)[Index2])
      DeltatChoice=OptimDelta$Deltat[Index]
      F=FitBivariateCBD1withDelta_modeling(Index1,Index2,Series,DeltatChoice, Training)
      message(DeltatChoice)
      temp_RS=ObtainRelevantScale(F)
      RelevantScale=rbind(RelevantScale,temp_RS)
      E=ExtrapolateBivarateCBD1withDelta(F, h = t)
      BivarateCBD1withDelta[[s]]=F
      PredictionBivarateCBD1withDelta[[s]]=E
      s=s+1
    }
  }
  ####################################################
  # Rand_Avg
  ModelInfo=data.frame()
  for(i in 1:length(BivarateCBD1withDelta_modeling)){
    temp=data.frame(Index=i,target1=BivarateCBD1withDelta_modeling[[i]]$label1,target2=BivarateCBD1withDelta_modeling[[i]]$label2,
                    delta=BivarateCBD1withDelta_modeling[[i]]$Deltat)
    ModelInfo=rbind(ModelInfo,temp)
  }

  ValidationSSE_Rank=data.frame()
  for(i in 1:length(OptimDelta$label1)){
    Index=which(ModelInfo$target1==OptimDelta$label1[i] & ModelInfo$target2==OptimDelta$label2[i] & ModelInfo$delta==OptimDelta$Deltat[i])
    temp=PerformanceValidationCBD1withDelta(ValidateBivarateCBD1withDelta[[Index]],Series,Dataset_Validate)
    ValidationSSE_Rank=rbind(ValidationSSE_Rank,temp)
  }

  ValidationSSE_Rank_byPopulation=list()
  for(i in 1:length(Candidate)){
    target=names(Candidate)[i]
    Index1=which(ValidationSSE_Rank$label1==target)
    Index2=which(ValidationSSE_Rank$label2==target)
    temp1=data.frame(target=ValidationSSE_Rank$label1[Index1],reference=ValidationSSE_Rank$label2[Index1],
                     delta=ValidationSSE_Rank$Deltat[Index1],SSE=ValidationSSE_Rank$SSE_V.1[Index1])
    temp2=data.frame(target=ValidationSSE_Rank$label2[Index2],reference=ValidationSSE_Rank$label1[Index2],
                     delta=-ValidationSSE_Rank$Deltat[Index2],SSE=ValidationSSE_Rank$SSE_V.2[Index2])
    temp=rbind(temp1,temp2)
    temp.sort=temp[order(temp$SSE),]
    ValidationSSE_Rank_byPopulation[[i]]=temp.sort
    names(ValidationSSE_Rank_byPopulation)[i]=target
  }

  PredictionRule_RankAvgSSE=list()
  for(i in 1:length(Candidate)){
    PredictionRule_RankAvgSSE[[i]]=ChooseGroupSize(ValidationSSE_Rank_byPopulation[[i]],Series, Dataset_Validate,
                                                   ModelInfo,ValidateBivarateCBD1withDelta)
    names(PredictionRule_RankAvgSSE)[i]=names(ValidationSSE_Rank_byPopulation)[i]
  }

  CBD=list()
  for(i in 1:length(Candidate)){
    temp = PerformanceMeasure_CBD_RankAvg_SSE(PredictionRule_RankAvgSSE[[i]],Series,Dataset_Test,PredictionBivarateCBD1withDelta)
    CBD[[i]]= temp
  }
  if (!is.null(save_file)) {
    save(CBD, file = save_file)
  }
  CBD
}



#' Bivariate ACF model-based ensemble forecast
#'
#' Applies the ACF-type bivariate model-based ensemble described in the BMBE
#' paper.
#'
#' @param prepared A prepared mortality dataset returned by
#'   [prepare_mortselect_data()].
#' @param Series Mortality series to forecast.
#' @param t Number of years to forecast.
#' @inheritParams forecast_dsa
#' @return A named list of ACF ensemble forecast objects.
#' @export
#'
forecast_bmbe_acf = function(prepared,
                             Series = 'female',
                             t = 15,
                             save_file = NULL){
  validate_prepared_mortselect_data(prepared, series = Series, min_populations = 2)
  Training = prepared$full
  Candidate = prepared$learning
  Dataset_Validate = prepared$validation
  Dataset_Test = prepared$test
  validation_h <- ncol(Dataset_Validate[[1]]$rate[[Series]])
  if (length(Candidate) < 2) {
    stop("At least two populations are required for the bivariate ACF ensemble.", call. = FALSE)
  }

  Result_bivariateACF=data.frame()
  for(i in 1:(length(Candidate)-1))
  {
    for(j in (i+1):length(Candidate))
    {
      message(names(Candidate)[i])
      message(names(Candidate)[j])
      temp=JointModelACF(Candidate[[i]],Candidate[[j]],Series, Candidate, Dataset_Validate,
                         h = validation_h)
      Result_bivariateACF=rbind(Result_bivariateACF,temp)
    }
  }

  ## Obtain Optim Delta choice based on sum of validation SSE ####
  Result_DeltaChoice_SumSSE_V=data.frame()
  for(i in 1:(length(Candidate)-1)){
    for(j in (i+1):length(Candidate)){
      Index=which(Result_bivariateACF$target==names(Candidate)[i] & Result_bivariateACF$reference==names(Candidate)[j])
      temp=Result_bivariateACF[Index,]
      DeltaChoice=which(temp$SSE_V1+temp$SSE_V2==min(temp$SSE_V1+temp$SSE_V2))
      Out.raw=temp[DeltaChoice,]
      Out1=data.frame(target=Out.raw$target,reference=Out.raw$reference,delta=Out.raw$delta,SSE=Out.raw$SSE_V1)
      Out2=data.frame(target=Out1$reference,reference=Out1$target,delta=-Out1$delta,SSE=Out.raw$SSE_V2)
      Result_DeltaChoice_SumSSE_V=rbind(Result_DeltaChoice_SumSSE_V,Out1,Out2)
    }
  }

  LogRatePrediction=list()
  for(i in 1:length(Result_DeltaChoice_SumSSE_V[,1])){
    LogRatePrediction[[i]]=Predict(Result_DeltaChoice_SumSSE_V[i,],Series, Training, Dataset_Test, t)
  }


  Validate_Forecast_Result=list()
  s=1
  for(i in 1:length(Result_DeltaChoice_SumSSE_V[,1])){
    target_label=Result_DeltaChoice_SumSSE_V[i,]$target
    reference_label=Result_DeltaChoice_SumSSE_V[i,]$reference
    message(target_label)
    message(reference_label)
    target=Candidate[[which(names(Candidate)==target_label)]]
    reference=Candidate[[which(names(Candidate)==reference_label)]]
    delta=Result_DeltaChoice_SumSSE_V[i,]$delta
    G=GroupdemogdataWithDrift(target,reference,target$year[1],target$year[length(target$year)],delta)
    message(G$year)
    V=Validate_Forecast(target,G,Series, Dataset_Validate, h = validation_h)
    Validate_Forecast_Result[[s]]=V
    s=s+1
  }

  PredictionRule_RankAvg=list()
  for(i in 1:length(Candidate)){
    PredictionRule_RankAvg[[i]]=RankAndAvg(Candidate[[i]], Series,
                                           Result_DeltaChoice_SumSSE_V,
                                           Dataset_Validate,
                                           Validate_Forecast_Result)
    names(PredictionRule_RankAvg)[i]=names(Candidate)[i]
  }

  ACF=list()
  for(i in 1:length(Candidate)){
    temp=PerformanceMeasure_RankAvg(Candidate[[i]],Series,
                                    Dataset_Test,PredictionRule_RankAvg,
                                    Result_DeltaChoice_SumSSE_V,LogRatePrediction)
    ACF[[i]]=temp
  }
  if (!is.null(save_file)) {
    save(ACF, file = save_file)
  }
  ACF
}


#' Single-population age-specific band forecast
#'
#' Implements the age-specific band-selection method for each population
#' separately.
#'
#' @param prepared A prepared mortality dataset returned by
#'   [prepare_mortselect_data()].
#' @param Series Mortality series to forecast.
#' @param t Number of years to forecast.
#' @param smoothing_window Window size for smoothing selected bandwidths.
#' @inheritParams forecast_dsa
#' @return A named list of log-rate forecast matrices.
#' @export
#'
forecast_ass_single = function(prepared,
                               Series = 'female',
                               t = 15,
                               smoothing_window = 5,
                               save_file = NULL){
  validate_prepared_mortselect_data(prepared, series = Series, min_populations = 1)
  Training = prepared$full
  Candidate = prepared$learning
  Dataset_training = prepared$training
  Dataset_Test = prepared$test
  ages <- prepared$metadata$ages %||% Candidate[[1]]$age

  ChosenK=data.frame()
  for(i in seq_along(Candidate)){
    TargetPopulation=names(Candidate)[i]
    for(Age in ages){
      TargetAge=Age
      target=paste(TargetPopulation,TargetAge,sep="")
      message(message("the target is ",target))
      UpperLimit_temp=max(abs(TargetAge - range(ages)))
      temp=ChoseK_BCV(TargetPopulation,TargetAge,Series,UpperLimit=UpperLimit_temp, Dataset_training)
      ChosenK=rbind(ChosenK,temp)
    }
  }

  Smooth.ChosenK=data.frame()
  for (Population in names(Candidate)){
    ChosenK.temp=subset(ChosenK,TargetPopulation==Population)
    ChosenK.temp$smooth.K=Smoothing.mean(ChosenK.temp,smoothing_window)
    Smooth.ChosenK=rbind(Smooth.ChosenK,ChosenK.temp)
  }

  Result_allTargetPopulation=list()
  #testSSE_AgeBand=data.frame()
  LogRate.Fcast.matrix=list()
  for(i in seq_along(Candidate)){
    TargetPopulation=names(Candidate)[i]
    message(TargetPopulation)
    matrix.temp=matrix(nrow=length(ages),ncol=t)
    #matrix.temp=matrix(nrow=101,ncol=8)
    s=1
    for (TargetAge in ages){
      message(TargetAge)
      out=performance_test(TargetPopulation,TargetAge,Series,Dataset_Test, Training, Smooth.ChosenK, t)
      Result_allTargetPopulation[[i]]=out
      names(Result_allTargetPopulation)[i]=TargetPopulation
      #temp=data.frame(TargetPopulation=TargetPopulation,TargetAge=TargetAge,sse_lc=out$sse_lc,sse_AgeBand=out$sse_AgeBand,BandWidth=out$ChosenKvalue)
      #testSSE_AgeBand=rbind(testSSE_AgeBand,temp)
      matrix.temp[s,]=as.vector(out$logrates_fcast)
      s=s+1
    }
    colnames(matrix.temp) = c((Training[[i]]$year[length(Training[[i]]$year)]+1):
                                (Training[[i]]$year[length(Training[[i]]$year)]+t))
    LogRate.Fcast.matrix[[i]]=matrix.temp
    names(LogRate.Fcast.matrix)[i]=paste0(TargetPopulation,"_",Series)
  }
  # Summary.Result=as.data.frame(testSSE_AgeBand %>%
  #                                group_by(TargetPopulation) %>%
  #                                summarise(sse_lc = sum(sse_lc),sse_AgeBand=sum(sse_AgeBand)))
  if (!is.null(save_file)) {
    save(LogRate.Fcast.matrix, file = save_file)
  }
  LogRate.Fcast.matrix
  #resave(Summary.Result,file="Single_ASS_Prediction.RData")
}



#' Multi-population age-specific band forecast
#'
#' Extends age-specific band selection by selecting information across both
#' ages and populations.
#'
#' @param prepared A prepared mortality dataset returned by
#'   [prepare_mortselect_data()].
#' @param Series Mortality series to forecast.
#' @param t Number of years to forecast.
#' @param smoothing_window Window size for smoothing selected bandwidths.
#' @inheritParams forecast_dsa
#' @return A named list of log-rate forecast matrices.
#' @export
#'
forecast_ass_multi = function(prepared,
                              Series = 'female',
                              t = 15,
                              smoothing_window = 9,
                              save_file = NULL){
  validate_prepared_mortselect_data(prepared, series = Series, min_populations = 1)
  Training = prepared$full
  Candidate = prepared$learning
  Dataset_training = prepared$training
  Dataset_Test = prepared$test
  ages <- prepared$metadata$ages %||% Candidate[[1]]$age

  ## Population-Age dictionary ####
  Dictionary=create_population_age_dictionary(Candidate, ages = ages)

  ## All Traning Data
  Training.rate=bind_training_rates(Dataset_training, series = Series, ages = ages)

  ChosenK=data.frame()
  for(i in seq_along(Candidate)){
    TargetPopulation=names(Candidate)[i]
    for(Age in ages){
      TargetAge=Age
      target=paste(TargetPopulation,TargetAge,sep="")
      message(message("the target is ",target))
      UpperLimit_temp=length(ages) * length(Candidate) - 1
      temp=ChoseK_BCV_multi(TargetPopulation,TargetAge,Series,UpperLimit=UpperLimit_temp, Training.rate)
      ChosenK=rbind(ChosenK,temp)
    }
  }
  Smooth.ChosenK=data.frame()
  for (Population in names(Candidate)){
    ChosenK.temp=subset(ChosenK,TargetPopulation==Population)
    ChosenK.temp$smooth.K=Smoothing.mean(ChosenK.temp,smoothing_window)
    Smooth.ChosenK=rbind(Smooth.ChosenK,ChosenK.temp)
  }

  Result_allTargetPopulation=list()
  #testSSE_AgeBand=data.frame()
  LogRate.Fcast.matrix=list()
  #Error.matrix=list()
  for(i in seq_along(Candidate)){
    TargetPopulation=names(Candidate)[i]
    message(TargetPopulation)
    matrix.temp=matrix(nrow=length(ages),ncol=t)
    #matrix.temp=matrix(nrow=101,ncol=8)
    #error.matrix.temp=matrix(nrow=101,ncol=8)
    s=1
    for (TargetAge in ages){
      message(TargetAge)
      out=performance_test_multi(TargetPopulation,TargetAge,Series,Training.rate,Dataset_Test,Training, Smooth.ChosenK,Dictionary, t)
      Result_allTargetPopulation[[i]]=out
      names(Result_allTargetPopulation)[i]=TargetPopulation
      #temp=data.frame(TargetPopulation=TargetPopulation,TargetAge=TargetAge,sse_lc=out$sse_lc,sse_AgeBand=out$sse_AgeBand,BandWidth=out$ChosenKvalue,DomesticK=out$DomesticK,ForeignK=out$ForeignK)
      #testSSE_AgeBand=rbind(testSSE_AgeBand,temp)
      matrix.temp[s,]=as.vector(out$logrates_fcast)
      #error.matrix.temp[s,]=as.vector(out$SE)
      s=s+1
    }
    colnames(matrix.temp) = c((Training[[i]]$year[length(Training[[i]]$year)]+1):
                                (Training[[i]]$year[length(Training[[i]]$year)]+t))
    LogRate.Fcast.matrix[[i]]=matrix.temp
    #Error.matrix[[i]]=error.matrix.temp
    names(LogRate.Fcast.matrix)[i]=paste0(TargetPopulation,"_",Series)
    #names(LogRate.Fcast.matrix)[i]=names(Error.matrix)[i]=TargetPopulation
  }
  # Summary.Result.smooth=as.data.frame(testSSE_AgeBand %>%
  #                                       group_by(TargetPopulation) %>%
  #                                       summarise(sse_lc = sum(sse_lc),sse_AgeBand=sum(sse_AgeBand)))
  if (!is.null(save_file)) {
    save(LogRate.Fcast.matrix, file = save_file)
  }
  LogRate.Fcast.matrix
  #resave(Summary.Result.smooth,file="Multi_ASS_Prediction.RData")
}





