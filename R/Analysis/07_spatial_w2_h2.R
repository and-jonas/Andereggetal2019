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
} else NULL
sinkdir <- paste0(path_to_data, "OUT/") #Specify output directory

#load required function
source(paste0(path_to_utils, "004_spat_aov_utils.R"))

#====================================================================================== -

#load data
Data <- readRDS(paste0(path_to_data, "other_data/data_aov.rds"))

#w2, BLUE, spatially corrected plot value ----

#calculate within-year repeatability and extract BLUEs for all senesence dynamics traits,
#for GY, GPC and heading date
Exp_wise <- Data %>% 
  dplyr::select(-heading_date) %>% 
  gather(Trait, value, 15:length(.)) %>%
  filter(Trait %in% c("heading_GDDAS", "GY", "GPC")) %>%
  distinct() %>% 
  group_by(Exp, Trait) %>% 
  nest() %>%  
  #calculate within-year repeatablity
  mutate(w2 = purrr::map(.x = data, genotype.as.random = TRUE,
                         .f = possibly(f_spats, otherwise = NA_real_)) %>%
           purrr::map_dbl(.x = .,
                          .f = possibly(getHeritability, otherwise = NA_real_))) %>%
  #extract BLUEs and spatially corrected plot values
  mutate(obj = purrr::map(.x = data, genotype.as.random = FALSE,
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

saveRDS(BLUEs, paste0(sinkdir, "BLUEs.rds"))

#====================================================================================== -
  
#h2 ----

across <- Exp_wise %>%
  dplyr::select(Trait, Exp, spat_corr) %>% 
  #no data for GPC in 2018, exclude
  filter(Trait!="GPC" | Exp != "FPWW022") %>% unnest() %>% 
  mutate(Gen_Name = as.factor(Gen_Name)) %>% 
  group_by(Trait, Exp, Lot) %>% nest() %>% 
  #check whether there is data for both lots
  mutate(remove = purrr::map_chr(data, remove_reps, value = "spat_corr", missing = 200)) %>% filter(remove == FALSE) %>% dplyr::select(-remove) %>% 
  ungroup() %>% unnest() %>% 
  group_by(Trait) %>% nest() %>% 
  mutate(h2_spat_corr_c = purrr::map_dbl(data, possibly(get_h2_asreml, otherwise = NA_real_), #after bug-fix, replace by map_dbl
                                         fixed = "spat_corr ~ Exp",
                                         random = "Gen_Name + Gen_Name:Exp",
                                         residual = "NULL",
                                         cullis = TRUE))

#====================================================================================== -

## THESE ARE FINAL RESULTS 

#====================================================================================== -