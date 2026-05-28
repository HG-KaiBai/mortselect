
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
#' @param parallel Logical; if `TRUE`, run independent target or model-pair
#'   tasks in parallel.
#' @param workers Optional number of parallel workers. Defaults to one fewer
#'   than the number of physical CPU cores, capped at the number of tasks.
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
                        save_file = NULL,
                        parallel = FALSE,
                        workers = NULL)
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
  Ans <- .mort_lapply(Labels, function(country) {
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
    list(label = country, Group = Best_group, Prediction = Rate_predict)
    #Ans[[country]] = list(label = country, Error = Err, Prediction = Rate_predict)
  }, parallel = parallel, workers = workers)
  names(Ans) <- Labels
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
                             save_file = NULL,
                             parallel = FALSE,
                             workers = NULL){
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
  PairIndex <- t(utils::combn(seq_along(Candidate), 2))
  pair_modeling <- .mort_lapply(seq_len(nrow(PairIndex)), function(pair_row) {
      Index1 <- PairIndex[pair_row, 1]
      Index2 <- PairIndex[pair_row, 2]
      message(names(Candidate)[c(Index1, Index2)])
      Test=data.frame()
      models <- list()
      validations <- list()
      for(Deltat in -10:10){
        F=FitBivariateCBD1withDelta_modeling(Index1,Index2,Series,Deltat, Candidate)
        E=ExtrapolateBivarateCBD1withDelta(F, h = validation_h)
        temp=ValidationSSE_CBD1withDelta(E,Series, Dataset_Validate)
        models[[length(models) + 1]] <- F
        validations[[length(validations) + 1]] <- E
        Test=rbind(Test,temp)
      }
      optim=Test[which(Test$SSE_V==min(Test$SSE_V)),]
      list(models = models, validations = validations, validation = Test, optim = optim)
  }, parallel = parallel, workers = workers)
  BivarateCBD1withDelta_modeling <- unlist(lapply(pair_modeling, `[[`, "models"), recursive = FALSE)
  ValidateBivarateCBD1withDelta <- unlist(lapply(pair_modeling, `[[`, "validations"), recursive = FALSE)
  ValidationSSESummary <- do.call(rbind, lapply(pair_modeling, `[[`, "validation"))
  OptimDelta <- do.call(rbind, lapply(pair_modeling, `[[`, "optim"))
  ####################################################
  final_pair_models <- .mort_lapply(seq_len(nrow(PairIndex)), function(pair_row) {
      Index1 <- PairIndex[pair_row, 1]
      Index2 <- PairIndex[pair_row, 2]
      message(names(Candidate)[c(Index1, Index2)])
      Index=which(OptimDelta$label1== names(Candidate)[Index1] & OptimDelta$label2== names(Candidate)[Index2])
      DeltatChoice=OptimDelta$Deltat[Index]
      F=FitBivariateCBD1withDelta_modeling(Index1,Index2,Series,DeltatChoice, Training)
      message(DeltatChoice)
      temp_RS=ObtainRelevantScale(F)
      E=ExtrapolateBivarateCBD1withDelta(F, h = t)
      list(model = F, prediction = E, relevant_scale = temp_RS)
  }, parallel = parallel, workers = workers)
  BivarateCBD1withDelta <- lapply(final_pair_models, `[[`, "model")
  PredictionBivarateCBD1withDelta <- lapply(final_pair_models, `[[`, "prediction")
  RelevantScale <- do.call(rbind, lapply(final_pair_models, `[[`, "relevant_scale"))
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

  PredictionRule_RankAvgSSE <- .mort_lapply(seq_along(Candidate), function(i) {
    ChooseGroupSize(ValidationSSE_Rank_byPopulation[[i]],Series, Dataset_Validate,
                    ModelInfo,ValidateBivarateCBD1withDelta)
  }, parallel = parallel, workers = workers)
  names(PredictionRule_RankAvgSSE)=names(ValidationSSE_Rank_byPopulation)

  CBD <- .mort_lapply(seq_along(Candidate), function(i) {
    PerformanceMeasure_CBD_RankAvg_SSE(PredictionRule_RankAvgSSE[[i]],Series,Dataset_Test,PredictionBivarateCBD1withDelta)
  }, parallel = parallel, workers = workers)
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
                             save_file = NULL,
                             parallel = FALSE,
                             workers = NULL){
  validate_prepared_mortselect_data(prepared, series = Series, min_populations = 2)
  Training = prepared$full
  Candidate = prepared$learning
  Dataset_Validate = prepared$validation
  Dataset_Test = prepared$test
  validation_h <- ncol(Dataset_Validate[[1]]$rate[[Series]])
  if (length(Candidate) < 2) {
    stop("At least two populations are required for the bivariate ACF ensemble.", call. = FALSE)
  }

  PairIndex <- t(utils::combn(seq_along(Candidate), 2))
  Result_bivariateACF <- do.call(rbind, .mort_lapply(seq_len(nrow(PairIndex)), function(pair_row) {
      i <- PairIndex[pair_row, 1]
      j <- PairIndex[pair_row, 2]
      message(names(Candidate)[i])
      message(names(Candidate)[j])
      JointModelACF(Candidate[[i]],Candidate[[j]],Series, Candidate, Dataset_Validate,
                    h = validation_h)
  }, parallel = parallel, workers = workers))

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

  LogRatePrediction <- .mort_lapply(seq_len(nrow(Result_DeltaChoice_SumSSE_V)), function(i) {
    Predict(Result_DeltaChoice_SumSSE_V[i,],Series, Training, Dataset_Test, t)
  }, parallel = parallel, workers = workers)


  Validate_Forecast_Result <- .mort_lapply(seq_len(nrow(Result_DeltaChoice_SumSSE_V)), function(i) {
    target_label=Result_DeltaChoice_SumSSE_V[i,]$target
    reference_label=Result_DeltaChoice_SumSSE_V[i,]$reference
    message(target_label)
    message(reference_label)
    target=Candidate[[which(names(Candidate)==target_label)]]
    reference=Candidate[[which(names(Candidate)==reference_label)]]
    delta=Result_DeltaChoice_SumSSE_V[i,]$delta
    G=GroupdemogdataWithDrift(target,reference,target$year[1],target$year[length(target$year)],delta)
    message(G$year)
    Validate_Forecast(target,G,Series, Dataset_Validate, h = validation_h)
  }, parallel = parallel, workers = workers)

  PredictionRule_RankAvg <- .mort_lapply(seq_along(Candidate), function(i) {
    RankAndAvg(Candidate[[i]], Series,
               Result_DeltaChoice_SumSSE_V,
               Dataset_Validate,
               Validate_Forecast_Result)
  }, parallel = parallel, workers = workers)
  names(PredictionRule_RankAvg)=names(Candidate)

  ACF <- .mort_lapply(seq_along(Candidate), function(i) {
    PerformanceMeasure_RankAvg(Candidate[[i]],Series,
                               Dataset_Test,PredictionRule_RankAvg,
                               Result_DeltaChoice_SumSSE_V,LogRatePrediction)
  }, parallel = parallel, workers = workers)
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
forecast_asb_single = function(prepared,
                               Series = 'female',
                               t = 15,
                               smoothing_window = 5,
                               save_file = NULL,
                               parallel = FALSE,
                               workers = NULL){
  validate_prepared_mortselect_data(prepared, series = Series, min_populations = 1)
  Training = prepared$full
  Candidate = prepared$learning
  Dataset_training = prepared$training
  Dataset_Test = prepared$test
  ages <- prepared$metadata$ages %||% Candidate[[1]]$age

  ChosenK <- do.call(rbind, .mort_lapply(names(Candidate), function(TargetPopulation) {
    out <- data.frame()
    for(Age in ages){
      TargetAge=Age
      target=paste(TargetPopulation,TargetAge,sep="")
      message(message("the target is ",target))
      UpperLimit_temp=max(abs(TargetAge - range(ages)))
      temp=ChoseK_BCV(TargetPopulation,TargetAge,Series,UpperLimit=UpperLimit_temp, Dataset_training)
      out=rbind(out,temp)
    }
    out
  }, parallel = parallel, workers = workers))

  Smooth.ChosenK=data.frame()
  for (Population in names(Candidate)){
    ChosenK.temp=ChosenK[ChosenK$TargetPopulation == Population, , drop = FALSE]
    ChosenK.temp$smooth.K=Smoothing.mean(ChosenK.temp,smoothing_window)
    Smooth.ChosenK=rbind(Smooth.ChosenK,ChosenK.temp)
  }

  LogRate.Fcast.matrix <- .mort_lapply(seq_along(Candidate), function(i) {
    TargetPopulation=names(Candidate)[i]
    message(TargetPopulation)
    matrix.temp=matrix(nrow=length(ages),ncol=t)
    #matrix.temp=matrix(nrow=101,ncol=8)
    s=1
    for (TargetAge in ages){
      message(TargetAge)
      out=performance_test(TargetPopulation,TargetAge,Series,Dataset_Test, Training, Smooth.ChosenK, t)
      #temp=data.frame(TargetPopulation=TargetPopulation,TargetAge=TargetAge,sse_lc=out$sse_lc,sse_AgeBand=out$sse_AgeBand,BandWidth=out$ChosenKvalue)
      #testSSE_AgeBand=rbind(testSSE_AgeBand,temp)
      matrix.temp[s,]=as.vector(out$logrates_fcast)
      s=s+1
    }
    colnames(matrix.temp) = c((Training[[i]]$year[length(Training[[i]]$year)]+1):
                                (Training[[i]]$year[length(Training[[i]]$year)]+t))
    matrix.temp
  }, parallel = parallel, workers = workers)
  names(LogRate.Fcast.matrix)=paste0(names(Candidate),"_",Series)
  # Summary.Result=as.data.frame(testSSE_AgeBand %>%
  #                                group_by(TargetPopulation) %>%
  #                                summarise(sse_lc = sum(sse_lc),sse_AgeBand=sum(sse_AgeBand)))
  if (!is.null(save_file)) {
    save(LogRate.Fcast.matrix, file = save_file)
  }
  LogRate.Fcast.matrix
  #resave(Summary.Result,file="Single_ASB_Prediction.RData")
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
forecast_asb_multi = function(prepared,
                              Series = 'female',
                              t = 15,
                              smoothing_window = 9,
                              save_file = NULL,
                              parallel = FALSE,
                              workers = NULL){
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

  ChosenK <- do.call(rbind, .mort_lapply(names(Candidate), function(TargetPopulation) {
    out <- data.frame()
    for(Age in ages){
      TargetAge=Age
      target=paste(TargetPopulation,TargetAge,sep="")
      message(message("the target is ",target))
      UpperLimit_temp=length(ages) * length(Candidate) - 1
      temp=ChoseK_BCV_multi(TargetPopulation,TargetAge,Series,UpperLimit=UpperLimit_temp, Training.rate)
      out=rbind(out,temp)
    }
    out
  }, parallel = parallel, workers = workers))
  Smooth.ChosenK=data.frame()
  for (Population in names(Candidate)){
    ChosenK.temp=ChosenK[ChosenK$TargetPopulation == Population, , drop = FALSE]
    ChosenK.temp$smooth.K=Smoothing.mean(ChosenK.temp,smoothing_window)
    Smooth.ChosenK=rbind(Smooth.ChosenK,ChosenK.temp)
  }

  LogRate.Fcast.matrix <- .mort_lapply(seq_along(Candidate), function(i) {
    TargetPopulation=names(Candidate)[i]
    message(TargetPopulation)
    matrix.temp=matrix(nrow=length(ages),ncol=t)
    #matrix.temp=matrix(nrow=101,ncol=8)
    #error.matrix.temp=matrix(nrow=101,ncol=8)
    s=1
    for (TargetAge in ages){
      message(TargetAge)
      out=performance_test_multi(TargetPopulation,TargetAge,Series,Training.rate,Dataset_Test,Training, Smooth.ChosenK,Dictionary, t)
      #temp=data.frame(TargetPopulation=TargetPopulation,TargetAge=TargetAge,sse_lc=out$sse_lc,sse_AgeBand=out$sse_AgeBand,BandWidth=out$ChosenKvalue,DomesticK=out$DomesticK,ForeignK=out$ForeignK)
      #testSSE_AgeBand=rbind(testSSE_AgeBand,temp)
      matrix.temp[s,]=as.vector(out$logrates_fcast)
      #error.matrix.temp[s,]=as.vector(out$SE)
      s=s+1
    }
    colnames(matrix.temp) = c((Training[[i]]$year[length(Training[[i]]$year)]+1):
                                (Training[[i]]$year[length(Training[[i]]$year)]+t))
    matrix.temp
  }, parallel = parallel, workers = workers)
  names(LogRate.Fcast.matrix)=paste0(names(Candidate),"_",Series)
  # Summary.Result.smooth=as.data.frame(testSSE_AgeBand %>%
  #                                       group_by(TargetPopulation) %>%
  #                                       summarise(sse_lc = sum(sse_lc),sse_AgeBand=sum(sse_AgeBand)))
  if (!is.null(save_file)) {
    save(LogRate.Fcast.matrix, file = save_file)
  }
  LogRate.Fcast.matrix
  #resave(Summary.Result.smooth,file="Multi_ASB_Prediction.RData")
}





