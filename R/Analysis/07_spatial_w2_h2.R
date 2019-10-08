#====================================================================================== -

#HEADER ----

# Author: Jonas Anderegg, ETH Zürich
# Copyright (C) 2019  ETH Zürich, Jonas Anderegg (jonas.anderegg@usys.ethz.ch)

# Calculate within-year repeatbility (w2)
# Calculate BLUEs per year
# Calculate spatially correct plot values
# Calculate across-year heritability (h2)


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
library(SpATS)
library(asreml)

#set working directory
path_to_data <- "" #Specify path to data
path_to_utils <- "" #Specify path to functions

#create sink directory
if(!dir.exists(paste0(path_to_data, "OUT"))){
  dir.create(paste0(path_to_data, "OUT"))
}
sinkdir <- paste0(path_to_data, "OUT/") #Specify output directory

#load required function
source(paste0(path_to_utils, "004_spat_aov_utils.R"))
source(paste0(path_to_utils, "000_helper_funcs.R"))

#====================================================================================== -

#load data
Data <- readRDS(paste0(path_to_data, "other_data/data_aov.rds"))

#w2, BLUE, spatially corrected plot value ----

#calculate within-year repeatability and extract BLUEs for all senesence dynamics traits,
#for GY, GPC and heading date
Exp_wise <- Data %>% 
  dplyr::select(-heading_date) %>% 
  gather(Trait, value, 16:length(.)) %>%
  filter(grepl("PSRI|GY|GPC|lin_Cnp", Trait)) %>%
  distinct() %>% 
  #remove lots with missing data
  group_by(Exp, Trait, Rep) %>% nest() %>% 
  mutate(remove = purrr::map_chr(data, remove_reps, missing = 200)) %>% filter(remove == FALSE) %>% dplyr::select(-remove) %>% unnest() %>% ungroup() %>% 
  group_by(Exp, Trait) %>% nest() %>%  
  #calculate within-year repeatablity
  mutate(w2 = purrr::map(.x = data, response = "value", random = "~ Xf + Yf", fixed = "~ check", genotype.as.random = TRUE, genotype = "Gen_Name",
                         .f = possibly(f_spats, otherwise = NA_real_)) %>%
           purrr::map_dbl(.x = .,
                          .f = possibly(get_h2, otherwise = NA_real_))) %>%
  #extract BLUEs and spatially corrected plot values
  mutate(obj = purrr::map(.x = data, response = "value", random = "~ Xf + Yf", fixed = "~ NULL", genotype.as.random = FALSE, genotype = "Gen_Name", 
                           .f = possibly(f_spats, otherwise = NA_real_))) %>%  
  mutate(BLUE =  purrr::map(.x = obj,
                      .f = possibly(get_BLUE_spats, otherwise = NA_real_))) %>% 
  mutate(spat_corr = purrr::map(.x = obj, 
                                .f = possibly(spat_corr_spats, otherwise = NA_real_)))

BLUEs <- Exp_wise %>% dplyr::select(Exp, Trait, BLUE) %>% 
  #remove year-trait combinations where no data exists
  mutate(remove = purrr::map_chr(BLUE, length)) %>% filter(remove > 1) %>% dplyr::select(-remove) %>% 
  unnest() %>% 
  spread(., Trait, BLUE)

#export for feature selection
saveRDS(BLUEs, paste0(sinkdir, "BLUEs.rds"))

#====================================================================================== -

#TWO STAGE ----
#heritability from best linear unbiased estimators
h2_BLUE <- Exp_wise %>%
  dplyr::select(Trait, Exp, BLUE) %>% 
  #no data for GPC in 2018, exclude
  filter(Trait!="GPC" | Exp != "FPWW022") %>% unnest() %>% 
  mutate(Gen_Name = as.factor(Gen_Name)) %>% 
  group_by(Trait) %>% nest() %>% 
  mutate(h2_BLUE = purrr::map_dbl(data, possibly(get_h2_asreml2, otherwise = NA_real_), #after bug-fix, replace by map_dbl
                                         fixed = "BLUE ~ Exp",
                                         random = "~Gen_Name",
                                         residual = "~NULL",
                                         cullis = FALSE)) %>% 
  mutate(h2_BLUE_c = purrr::map_dbl(data, possibly(get_h2_asreml2, otherwise = NA_real_), #after bug-fix, replace by map_dbl
                                  fixed = "BLUE ~ Exp",
                                  random = "~Gen_Name",
                                  residual = "~NULL",
                                  cullis = TRUE)) 

#====================================================================================== -

#ONE STAGE ----
#heritability from spatially corrected plot values
#scorings and agronomic traits
h2_spatcorr <- Exp_wise %>%
  dplyr::select(Trait, Exp, spat_corr) %>% 
  #no data for GPC in 2018, exclude
  filter(Trait!="GPC" | Exp != "FPWW022") %>% unnest() %>%
  filter(grepl("GY|lin_Cnp", Trait)) %>% 
  mutate(Gen_Name = as.factor(Gen_Name)) %>% 
  group_by(Trait) %>% nest() %>% 
  mutate(h2_spatcorr = purrr::map_dbl(data, possibly(get_h2_asreml, otherwise = NA_real_),
                                      fixed = "spat_corr ~ Exp ",
                                      random  = "~ Gen_Name + Gen_Name:Exp + Rep:at(Exp)",
                                      residual = "~dsum(~id(units) | Exp)",
                                      cullis = TRUE))

#for spectral features
h2_spatcorr <- Exp_wise %>%
  dplyr::select(Trait, Exp, spat_corr) %>% 
  #no data for GPC in 2018, exclude
  filter(Trait!="GPC" | Exp != "FPWW022") %>% unnest() %>%
  filter(grepl("PSRI", Trait)) %>% 
  mutate(Gen_Name = as.factor(Gen_Name)) %>% 
  group_by(Trait) %>% nest() %>% 
  mutate(h2_spatcorr = purrr::map_dbl(data, possibly(get_h2_asreml, otherwise = NA_real_),
                                      fixed = "spat_corr ~ Exp",
                                      random  = "~ Gen_Name + Gen_Name:Exp + Rep:at(Exp, 1)",
                                      residual = "~dsum(~id(units) | Exp)",
                                      cullis = TRUE))

#====================================================================================== -

## THESE ARE FINAL RESULTS 

#====================================================================================== -