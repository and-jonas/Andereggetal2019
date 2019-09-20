#====================================================================================== -

#HEADER ----

# Author: Jonas Anderegg, ETH Zürich
# Copyright (C) 2019  ETH Zürich, Jonas Anderegg (jonas.anderegg@usys.ethz.ch)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

#====================================================================================== -

#Helper functions ----
######################################################################################## -

#calculate RMSE
get_rmse = function(actual, predicted) {
  sqrt(mean((actual - predicted) ^ 2))
}

#perform up- or downsampling
perform_sampling <- function(data, Trait, method = "upsampling"){
  data[, Trait] <- as.factor(data[, Trait])
  if(method == "upsampling"){
    data <- upSample(x = data,
                     y = data[, Trait],
                     yname = Trait)
  } else if (method == "downsampling"){
    data <- downSample(x = data,
                       y = data[, Trait],
                       yname = Trait)
  } else print("specify sampling procedure for test data")
  data[, Trait] <- as.numeric(as.character(data[, Trait]))
  data <- data[, -ncol(data)]
  return(data)
}

#calculate RMSE on upsampled validation dataset
get_RMSE_sample <- function(sample, fit, Trait){
  ID_sub <- sample %>% dplyr::select(-contains("rflt"), -SnsCnp, -SnsFl0)
  data_sub <- sample %>% dplyr::select(contains("rflt"), SnsCnp, -SnsFl0)
  #predict new data using trained model
  pred <- caret::predict.train(fit, newdata = data_sub, type = "raw", na.action = na.omit)
  predobs <- cbind(ID_sub, sample[, Trait], pred) %>% as_tibble()
  names(predobs)[length(predobs)-1] <- "obs"
  RMSE <- round(get_rmse(predobs$obs, predobs$pred), 2)
  return(RMSE)
}

#generate predictions and extract corresponding real observations
extract_predobs <- function(fit, data, Trait){
  ID <- data %>% dplyr::select(-contains("rflt"), -SnsCnp, -SnsFl0)
  dat <- data %>% dplyr::select(contains("rflt"), SnsCnp, -SnsFl0)
  #create predictions
  pred <- caret::predict.train(fit, newdata = dat, type = "raw", na.action = na.omit)
  #create output tibble
  predobs <- cbind(ID, dat[, Trait], pred) %>% as_tibble()
  names(predobs)[length(predobs)-1] <- "obs"
  return(predobs)
}

#Main function ----
######################################################################################## -

#evaluate model performance within and across years
eval_full_spc_cross <- function(Trait, method,
                                    maxcomp = switch(method == "pls", 20), 
                                    tuneGrid = switch(method == "cubist", Grid), 
                                    data. = data, data_type = "rflt", 
                                    train. = train, test. = test, 
                                    trainsample = "fullsample", 
                                    testsample = "upsample",
                                    preProc = "center_scale"
                                    ){
  
  #select type of pre-processsing
  DATA = data.[[data_type]]
  
  #for each of the possible training datasets
  #i.e. for each year and all combinations of years
  #train and evaluate full-spectrum models
  out_train <- out_all <- list()
  for (i in 1:length(train.)) {
    
    #Verbose
    print(paste("> train on", paste(train.[[i]], collapse = "_")))
    
    #====================================================================================== -
    #TRAIN ALGORITHM, VALIDATE ON THE SAME YEAR#
    
    dat_train0 <- DATA %>% filter(Exp %in% train.[[i]]) %>% as.data.frame()
    
    if(trainsample == "downsample"){
      dat_train <- perform_sampling(dat_train0, Trait, method = "downsampling")
    } else if (trainsample == "upsample"){
      dat_train <- perform_sampling(dat_train0, Trait, method = "upsampling")
    } else if(trainsample == "fullsample") {
      dat_train <- dat_train0
    } else print("specify sampling procedure for train data")
    
    nsamp <- nrow(dat_train)
    ID_train <- dat_train %>% dplyr::select(-contains("rflt"), -SnsCnp, -SnsFl0)
    dat <- dat_train %>% dplyr::select(Trait, contains("rflt"))
    
    preProcess <- switch(!is.null(preProc), strsplit(preProc, "_") %>% unlist()) #translate pre-processing procedure for trainControl function
    train_index <- caret::createMultiFolds(y = dat[, Trait], k = 10, times = 1)
    ctrl <- caret::trainControl(method = "repeatedcv", 
                                index = train_index,
                                savePredictions = TRUE, 
                                selectionFunction = "oneSE",
                                allowParallel = TRUE,
                                verboseIter = FALSE)
    fit <- caret::train(as.formula(paste(Trait, "~.", sep = "")),
                        data = dat, 
                        preProcess = preProcess,
                        method = method,
                        tuneLength = switch(method == "pls", maxcomp), 
                        tuneGrid = switch(method == "cubist", Grid), 
                        trControl = ctrl,
                        importance = TRUE)
    #get model performance for training experiment
    predobs <- plyr::match_df(fit$pred, fit$bestTune, on = names(fit$bestTune)) %>% 
      group_by(rowIndex) %>%
      summarize(obs = mean(obs), mean_pred = mean(pred))
    ID_train$rowIndex <- 1:nrow(ID_train)
    predobs <- dplyr::full_join(ID_train, predobs, by = "rowIndex") %>% 
      dplyr::select(-rowIndex) %>% as_tibble()
    RMSE_train <- caret::getTrainPerf(fit) %>% pull(TrainRMSE)
    clusterExport(cluster, c("dat_train0", "get_rmse", "get_RMSE_sample", "perform_sampling", "Trait", "fit"), envir=environment())
      RMSE_train_up <- parSapply(cluster, 1:10, function(i, ...) {x <- perform_sampling(dat_train0, Trait, method = "upsampling") %>% 
        get_RMSE_sample(., fit, Trait)}) %>% 
        mean() %>% round(., 2)
    if(trainsample == "downsample"){
      predobs_all <- extract_predobs(fit, dat_train0, Trait)
    } else {
      predobs_all <- predobs
    }
    out_train[[i]] <-  tibble(Trait = Trait, 
                              train = paste(train.[[i]], collapse = "_"),
                              trainsample = trainsample,
                              data_type = data_type,
                              method = method, 
                              tunepar = list(fit$bestTune[1, 1]), 
                              preProc = preProc,
                              RMSE_train = RMSE_train,
                              RMSE_train_up = RMSE_train_up,
                              nsamp_train = nsamp,
                              predobs_train_all = list(predobs_all))
    
    #====================================================================================== -
    
    #for each validation dataset (i.e. different years), evaluate model performance
    out_test <- list()
    for(j in 1:length(test.[[i]])){
      print(paste0("validating on ", test.[[i]][[j]], collapse = "_"))
      dat_test0 <- DATA %>% 
        filter(Exp %in% test.[[i]][[j]]) %>% as.data.frame()
      nsamp_all <- nrow(dat_test0)
      predobs <- extract_predobs(fit, dat_test0, Trait)
      clusterExport(cluster, c("dat_test0", "fit"), envir=environment())
        RMSE_test_up <- parSapply(cluster, 1:10, function(i, ...) {x <- perform_sampling(dat_test0, Trait, method = "upsampling") %>%
          get_RMSE_sample(., fit, Trait)}) %>%
          mean() %>% round(., 2)
      predobs_all <- extract_predobs(fit, dat_test0, Trait)
      RMSE_test <- round(get_rmse(predobs_all$obs, predobs_all$pred), 2)
      out_test[[j]] <-  tibble(test = paste(test.[[i]][[j]], collapse = "_"),
                               RMSE_test = RMSE_test,
                               RMSE_test_up = RMSE_test_up,
                               nsamp_test = nsamp_all,
                               predobs_test_all = list(predobs_all))
    }#end of validation loop
    out_test_done <- bind_rows(out_test)
    out_all[[i]] <- tibble(train = list(out_train[[i]]), test = list(out_test_done)) %>% 
      unnest(train) %>% unnest(test, .preserve = c(predobs_train_all, tunepar))
  }#end of training set loop
  out <- bind_rows(out_all)
  return(out)
}#end of function
