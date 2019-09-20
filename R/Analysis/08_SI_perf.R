#====================================================================================== -

#HEADER ----

# Author: Jonas Anderegg, ETH Zürich
# Copyright (C) 2019  ETH Zürich, Jonas Anderegg (jonas.anderegg@usys.ethz.ch)

# Evaluate performance of spectral indices in tracking visually observed senescence dynamics
# extract senescence dynamics parameters from spectral indices and scorings,
# from linear interpolations or parametric models;
# calculate area between the interpolation curves;
# calculate bias as difference.


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

list.of.packages <- c("doParallel", "Rmpi", "future", "furrr", "tidyverse", "tictoc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE, repos='https://stat.ethz.ch/CRAN/')

library(doParallel)
library(tictoc)
library(future)
library(furrr)
library(tidyverse)

path_to_data <- "" #Specify path to data
path_to_utils <- "" #Specify path to functions

#create sink directory
if(!dir.exists(paste0(path_to_data, "OUT"))){
  dir.create(paste0(path_to_data, "OUT"))
} else NULL
sinkdir <- paste0(path_to_data, "OUT/") #Specify output directory

source(paste0(path_to_utils, "001_spectra_utils.R"))

#====================================================================================== -

#Performance evaluation ----

#> Selected SI ----
#Evaluate performance of selected subset of SI
data <- readRDS(paste0(path_to_data, "OUT/SI/SI_forperf.rds"))
data <- data[!data$Lot == 6, ] %>% droplevels() #remove incomplete time series
data <- data[!(data$Lot == 3 & data$Exp == "FPWW022"), ] #remove incomplete time series
ind <- readRDS(paste0(path_to_data, "helper_files/SI_decorr.rds")) #Define subset of SI for which to carry out the analysis
ind <- c(ind, "SI_NDVI_nb_ASD") #add narrow-band NDVI to the dataset as benchmark

plan("multiprocess")
perf <- data %>%
  gather(SVI, value, contains("SI_")) %>%
  gather(Trait, scoring, contains("Sns")) %>%
  #use the selected subset of spectral indices
  filter(SVI %in% ind) %>%
  group_by(Exp, Trait, SVI, Plot_ID) %>%
  nest() %>%
  mutate(eval = furrr::future_map(.x = data, method = "lin",
                                  .f = possibly(get_errors_and_dynpars, otherwise = NA_real_)))
plan("sequential")

saveRDS(perf, pasteo(sinkdir, "SI_perf.rds"))

#====================================================================================== -

#> PSRI ----
##Evaluate performance of PSRI with different constituting wavelengths 
data <- readRDS(paste0(path_to_data, "OUT/SI/SI_PSRI.rds"))
data <- data[!data$Lot == 6, ] %>% droplevels()
data <- data[!(data$Lot == 3 & data$Exp == "FPWW022"), ]

plan("multiprocess")
perf <- data %>% 
  gather(SVI, value, starts_with("SI")) %>% 
  gather(Trait, scoring, contains("Sns")) %>%
  group_by(Exp, Trait, SVI, Plot_ID) %>% 
  nest() %>% 
  mutate(eval = furrr::future_map(.x = data, method = "lin",
                                  .f = possibly(get_errors_and_dynpars, otherwise = NA_real_)))

plan("sequential")

saveRDS(perf, paste0(sinkdir, "SI_PSRI_perf.rds"))

# ====================================================================================== -

#Summarize results ----
#> selected SI ----

#rearrange output
perf_long <- readRDS(paste0(sinkdir, "SI_perf.rds")) %>%
  dplyr::select(-data) %>%
  #drop plots with missing results (NA in eval)
  mutate(remove = purrr::map_chr(eval, length)) %>% filter(remove > 1) %>% dplyr::select(-remove) %>%
  unnest() %>% dplyr::select(1:4, contains("d_"), Error) %>%
  gather(., metric, value, d_onsen:Error)

#get index performance metrics per experiment and trait
metrics_exp <- perf_long %>%
  group_by(Exp, Trait, SVI, metric) %>%
  summarize(mean = mean(value),
            sd = sd(value))

#get errors per experiment
err_exp <- metrics_exp %>% filter(metric == "Error") %>%
  arrange(Exp, Trait, mean)

#get index performance metrics per trait (across Experiments)
metrics_overall <- perf_long %>%
  group_by(Trait, SVI, metric) %>%
  summarize(mean = mean(value),
            sd = sd(value))

#get errors overall
err_overall <- metrics_overall %>% filter(metric == "Error") %>%
  arrange(Trait, mean)

#====================================================================================== -

#rearrange output
perf_long <- readRDS(paste0(sinkdir, "SI_perf.rds")) %>%
  dplyr::select(-data) %>%
  #drop plots with missing results (NA in eval)
  mutate(remove = purrr::map_chr(eval, length)) %>% filter(remove > 1) %>% dplyr::select(-remove) %>%
  unnest() %>% dplyr::select(-contains("d_"), -Error) %>%
  gather(., Dynpar, value, onsen_SI:tsen_Trait) %>%
  mutate(level = lapply(strsplit(Dynpar, "_"), "[[", 2) %>% unlist(),
         Dynpar = lapply(strsplit(Dynpar, "_"), "[[", 1) %>% unlist()) %>%
  mutate(level = ifelse(level == "Trait", "scoring", level)) %>%
  spread(., level, value) %>%
  arrange(Dynpar)

#get correlations per experiment and trait
corr_exp <- perf_long %>%
  group_by(Exp, Trait, SVI, Dynpar) %>% nest() %>%
  mutate(pearson_r = purrr::map_dbl(.x = data, x = "scoring", y = "SI",
                                    .f = do_cor_test),
         p_val = purrr::map_dbl(.x = data, x = "scoring", y = "SI", return = "p.value",
                                .f = do_cor_test)) %>%
  dplyr::select(-data) %>%
  arrange(Trait, Exp, Dynpar, desc(abs(pearson_r)))

#get correlations overall
corr_overall <- perf_long %>%
  group_by(Trait, SVI, Dynpar) %>% nest() %>%
  mutate(pearson_r = purrr::map_dbl(.x = data, x = "scoring", y = "SI",
                                    .f = do_cor_test),
         p_val = purrr::map_dbl(.x = data, x = "scoring", y = "SI", return = "p.value",
                                .f = do_cor_test)) %>%
  dplyr::select(-data) %>%
  arrange(Trait, Dynpar, desc(abs(pearson_r)))

#====================================================================================== -

#> PSRI' ----

#rearrange output
perf_long <- readRDS(paste0(sinkdir, "SI_PSRI_perf.rds")) %>%
  dplyr::select(-data) %>%
  #drop plots with missing results (NA in eval)
  mutate(remove = purrr::map_chr(eval, length)) %>% filter(remove > 1) %>% dplyr::select(-remove) %>%
  unnest() %>% dplyr::select(-contains("d_"), -Error) %>%
  gather(., Dynpar, value, onsen_SI:tsen_Trait) %>%
  mutate(level = lapply(strsplit(Dynpar, "_"), "[[", 2) %>% unlist(),
         Dynpar = lapply(strsplit(Dynpar, "_"), "[[", 1) %>% unlist()) %>%
  mutate(level = ifelse(level == "Trait", "scoring", level)) %>%
  spread(., level, value) %>%
  arrange(Dynpar)

#get correlations per experiment and trait
corr_exp <- perf_long %>%
  group_by(Exp, Trait, SVI, Dynpar) %>% nest() %>%
  mutate(pearson_r = purrr::map_dbl(.x = data, x = "scoring", y = "SI",
                                    .f = do_cor_test),
         p_val = purrr::map_dbl(.x = data, x = "scoring", y = "SI", return = "p.value",
                                .f = do_cor_test)) %>%
  dplyr::select(-data) %>%
  arrange(Trait, Exp, Dynpar, desc(abs(pearson_r)))

#====================================================================================== -

## THESE ARE FINAL RESULTS 

#====================================================================================== -
