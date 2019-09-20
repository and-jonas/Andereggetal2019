#====================================================================================== -

#HEADER ----

# Author: Jonas Anderegg, ETH Zürich
# Copyright (C) 2019  ETH Zürich, Jonas Anderegg (jonas.anderegg@usys.ethz.ch)

# Extract error and senescence dynamics parameters from model predictions and scorings
# Extract correlations, mean errors and model bias


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

.libPaths("") #Specify path to R Libraries

library(tidyverse)
library(future)
library(furrr)

path_to_data <- "" #Specify path to data
path_to_utils <- "" #Specify path to functions

source(paste0(path_to_utils, "001_Spectra_utils.R"))

#====================================================================================== -

# Prepare data ----
 
#load model output 
data <- as.list(path = paste0(path_to_data, "OUT/model_output/"), list.files(pattern = "result_")) %>% 
  lapply(., readRDS) %>% lapply(., unnest) %>% bind_rows() %>% dplyr::select(-contains("1")) %>% 
  dplyr::select(1:7, 11, 13, 17)

# rearrange output
## results from within-year validation
inyear <- data %>% dplyr::select(1:8) %>% 
  mutate(test = train) %>% dplyr::select(1:7, 9, 8) %>% 
  unique() %>% 
  dplyr::rename("predobs" = predobs_train_all)
## results from across-year validation
acyear <- data %>% dplyr::select(1:7, 9:10) %>% 
  dplyr::rename("predobs" = predobs_test_all) %>% 
  #combine
  bind_rows(inyear, .) %>% 
  #scale predictions and scorings
  mutate(rel = purrr::map(predobs, `[`, c(2, 8, 18, 19)) %>% 
           purrr::map(., f_scale_si)) 

# replace the unscaled scorings and predictions by scaled
# rename as required by function
for(i in 1:nrow(acyear)){
  acyear$predobs[[i]][18] <- acyear$rel[[i]][3]
  acyear$predobs[[i]][19] <- acyear$rel[[i]][4]
  names(acyear$predobs[[i]])[19] <- "value"
  names(acyear$predobs[[i]])[18] <- "scoring"
}

#====================================================================================== -

# model performance ----
# > Error, bias ----

# extract senescence dynamics parameters from model predictions and scorings;
# from linear interpolations or parametric models;
# calculate area between the interpolation curves;
# calculate bias as difference;
plan("multiprocess")
perf <- acyear %>% dplyr::select(-rel) %>% 
  unnest() %>% group_by(Trait, method, data_type, trainsample, testsample, 
                        preProc, train, test, Exp, Plot_ID) %>% 
  nest() %>% 
  mutate(eval = furrr::future_map(.x = data, method = "lin",
                                  .f = possibly(get_errors_and_dynpars,
                                                otherwise = NA_real_)))
plan("sequential")

#====================================================================================== -

#get index performance metrics per experiment and trait
metrics_exp <- perf %>% dplyr::select(-data) %>%
  #helper, to check successful evaluation
  mutate(ctrl = purrr::map_dbl(eval, length)) %>% filter(ctrl > 1) %>% dplyr::select(-ctrl) %>% 
  unnest() %>% dplyr::select(1:10, contains("d_"), Error) %>%  
  gather(., metric, value, d_onsen:Error) %>% 
  group_by(Trait, method, data_type, trainsample, testsample, preProc, train, test, Exp, metric) %>% 
  dplyr::summarize(mean = mean(value),
            sd = sd(value))

#get errors per experiment
err_exp <- metrics_exp %>% filter(metric == "Error") %>% 
  arrange(Exp, Trait, mean)

#====================================================================================== -

# > Correlations ----

#rearrange output from function
perf_long <- perf %>% dplyr::select(-data) %>% 
  #helper, to check successful evaluation
  mutate(ctrl = purrr::map_dbl(eval, length)) %>% filter(ctrl > 1) %>% dplyr::select(-ctrl) %>% 
  unnest() %>% dplyr::select(-contains("d_"), -Error) %>%  
  gather(., Dynpar, value, onsen_SI:tsen_Trait) %>% 
  mutate(level = lapply(strsplit(Dynpar, "_"), "[[", 2) %>% unlist(),
         Dynpar = lapply(strsplit(Dynpar, "_"), "[[", 1) %>% unlist()) %>% 
  mutate(level = ifelse(level == "Trait", "scoring", level)) %>% 
  spread(., level, value) %>% 
  arrange(Dynpar)

#get correlations per experiment and trait
corr_exp <- perf_long %>% 
  group_by(Trait, method, data_type, trainsample, testsample, preProc, train, test, Exp, Dynpar) %>% 
  nest() %>% 
  mutate(pearson_r = purrr::map_dbl(.x = data, x = "scoring", y = "SI", 
                                    .f = do_cor_test),
         p_val = purrr::map_dbl(.x = data, x = "scoring", y = "SI", return = "p.value",
                                .f = do_cor_test)) %>% 
  dplyr::select(-data) %>% 
  arrange(Trait, Exp, Dynpar, desc(abs(pearson_r)))

#====================================================================================== -

## THESE ARE FINAL RESULTS 

#====================================================================================== -
