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
source(paste0(path_to_utils, "000_helper_funcs.R"))

#====================================================================================== -

# Prepare data ----

#scaled visual scorings for SI and model performance
scr_sca <- data.table::fread(paste0(path_to_data, "scr_data/scr_sca_corr.csv")) %>% 
  dplyr::mutate(Plot_ID = as.factor(Plot_ID),
                Gen_Name = as.factor(Gen_Name),
                grading_date = as.Date(grading_date),
                heading_date = as.Date(heading_date)) %>% 
  as_tibble() %>% 
  rename("scoring" = SnsCnp)

#GDDAH for measurements and scorings (plot-specific)
gddah <- read.csv(paste0(path_to_data, "helper_files/gddah_data.csv")) %>% 
  mutate(meas_date = meas_date %>% as.Date()) %>% 
  as_tibble()

#template of matching dates to join spectral and scoring data
match_dates <- read.csv(paste0(path_to_data, "helper_files/match_dates.csv"), sep = ";") #load matching dates template

#load model output 
data <- as.list(list.files(path = paste0(path_to_data, "OUT/model_output/"), pattern = "result_4", full.names = TRUE)) %>% 
  lapply(., readRDS) %>% lapply(., unnest) %>% bind_rows() %>% dplyr::select(-contains("1"))

data_perf <- data %>% dplyr::select(Trait, method, data_type, trainsample, testsample, preProc, train, test, contains("RMSE")) %>% 
  mutate(RMSE_train = round(RMSE_train, 2))

data <- data %>%  
  dplyr::select(1:7, 11, 13, 17)

# rearrange output
## results from within-year validation
inyear <- data %>% dplyr::select(1:8) %>% 
  mutate(test = train) %>% dplyr::select(1:7, 9, 8) %>% 
  unique() %>% 
  dplyr::rename("predobs" = predobs_train_all)
## results from across-year validation
all <- data %>% dplyr::select(1:7, 9:10) %>% 
  dplyr::rename("predobs" = predobs_test_all) %>% 
  #combine
  bind_rows(inyear, .) %>% 
  #scale predictions
  mutate(rel = purrr::map(predobs, `[`, c(2, 8, 19)) %>% 
           purrr::map(., f_scale_si, sub = "post_heading"))

# replace the unscaled scorings and predictions by scaled
# rename as required by function
for(i in 1:nrow(all)){
  all$predobs[[i]][19] <- all$rel[[i]][3]
  names(all$predobs[[i]])[19] <- "value"
}

final <- all %>% mutate(predobs = purrr::map(predobs, dplyr::select, -Exp, -obs))

data <- final %>% dplyr::select(-rel, -trainsample, -testsample, -preProc) %>% 
  #filter results
  filter(data_type == "rflt_bin3") %>%
  filter(!grepl("_", test)) %>% 
  filter(!train == "FPWW012_FPWW018_FPWW022") %>% 
  mutate(predobs = purrr::map(predobs, dplyr::select, Plot_ID, meas_date, value)) %>% 
  #add scaled scorings: needs to be the full sequence! Only this allows direct comparison to SI performance
  mutate(predobs = purrr::map(predobs, f_match_join, scr_sca, gddah, match_dates)) %>% 
  mutate(predobs = purrr::map(predobs, filter, !is.na(Exp))) %>% 
  #remove Experiments added by match-joining with scorings
  mutate(predobs = purrr::map(predobs, remove_exps, value = "value", missing = 5000)) %>% 
  dplyr::select(-Trait, -data_type) %>% 
  arrange(method, test) %>% 
  #remove plots without a complete series of measurements (as with SI evaluation)
  mutate(predobs = purrr::map(predobs, remove_incomplete_series)) 

#====================================================================================== -

# model performance ----
# > Error, bias ----

# extract senescence dynamics parameters from model predictions and scorings;
# from linear interpolations or parametric models;
# calculate area between the interpolation curves;
# calculate bias as difference;
plan("multiprocess")
perf <- data %>%  
  unnest() %>% group_by(method, train, test, Plot_ID) %>% 
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
  group_by(method,train, test, metric) %>% 
  dplyr::summarize(mean = mean(value),
            sd = sd(value))

#get errors per experiment
err_exp <- metrics_exp %>% filter(metric == "Error") %>% 
  arrange(mean)

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
  group_by(method, train, test, Dynpar) %>% 
  nest() %>% 
  mutate(pearson_r = purrr::map_dbl(.x = data, x = "scoring", y = "SI", 
                                    .f = do_cor_test),
         p_val = purrr::map_dbl(.x = data, x = "scoring", y = "SI", return = "p.value",
                                .f = do_cor_test)) %>% 
  dplyr::select(-data) %>% 
  arrange(Dynpar, desc(abs(pearson_r)))

#====================================================================================== -

## THESE ARE FINAL RESULTS 

#====================================================================================== -
