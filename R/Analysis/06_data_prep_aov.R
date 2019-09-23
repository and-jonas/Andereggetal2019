#====================================================================================== -

#HEADER ----

# Author: Jonas Anderegg, ETH Zürich
# Copyright (C) 2019  ETH Zürich, Jonas Anderegg (jonas.anderegg@usys.ethz.ch)

# Preparation agronomic data and senescence dynamics data 
# for calculation of heritability and BLUEs


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
library(gapminder)

path_to_data <- ""#Specify path to data

#====================================================================================== -

grnyld <- read.csv(paste0(path_to_data, "other_data/gy.csv")) #grain yield
grnprt <- read.csv(paste0(path_to_data, "other_data/gpc.csv")) #grain protein concentration
head <- read.csv(paste0(path_to_data, "other_data/heading.csv")) %>% # heading date
  mutate(heading_date = as.Date(heading_date),
         heading_DAS = as.numeric(heading_DAS))

#rearrange output of parameter extraction (SI)
excl_plots <- read.csv(paste0(path_to_data, "other_data/lodging.csv"), sep = ";") %>% filter(Lodging %in% c("l", "ll")) %>% pull(Plot_ID)
pars_SI <- readRDS(paste0(path_to_data, "OUT/dynpars/dynpars_SI.rds")) %>% #read output of parameter extrcation
  dplyr::select(-data, -Exp) %>% 
  mutate(remove = purrr::map_chr(dynpars, length)) %>% filter(remove > 1) %>% dplyr::select(-remove) %>% #remove plots with missing data
  unnest() %>% 
  mutate_if(is.integer, as.numeric) %>% 
  gather(., parameter, value, onsen_SI:tsen_SI, factor_key = TRUE) %>% 
  #replace with NA, where affected by lodging  
  mutate(value = ifelse(Plot_ID %in% excl_plots, NA_real_, value),
         parameter = gsub("_SI", "_lin", parameter),
         Trait = paste(parameter, SVI, sep = "_")) %>% 
  dplyr::select(-parameter, -SVI) %>% 
  spread(., Trait, value) %>% 
  arrange(Plot_ID)

#rearrange output of parameter extraction (scr)
pars_scr <- readRDS(paste0(path_to_data, "OUT/dynpars/dynpars_scr.rds")) %>% 
  dplyr::select(-convInfo) %>% 
  #transform to wide
  #in order to analyze scoring at each date seperately
  #i.e. as a separate trait
  ##first to FULL long format
  gather(., parameter, value, onsen_gom:tsen_lin, factor_key = TRUE) %>% 
  mutate(level = gsub("Sns", "", Trait)) %>% 
  mutate(Trait = paste(parameter, level, sep = "_")) %>% 
  dplyr::select(-level, -parameter) %>% 
  ##then to wide
  spread(., Trait, value)

# Merge all
Data <- list(head, grnyld, grnprt, pars_scr, pars_SI) %>% 
  Reduce(function(dtf1,dtf2) full_join(dtf1, dtf2, by = "Plot_ID"), .)

# Add design
exp <- read.csv(paste0(path_to_data, "helper_files/exp_design.csv")) %>% 
  as_tibble() %>% 
  mutate_at(vars(Lot, RangeLot, RowLot, Range, Row, RowBL, RangeBL), funs(as.numeric)) %>% 
  mutate_at(vars(Gen_ID, Rep, Xf, Yf), funs(as.factor))

Data <- full_join(exp_design, Data, by = "Plot_ID") %>%
  as_tibble() %>% 
  arrange(Exp, Lot, RangeLot, RowLot)

saveRDS(Data, paste0(path_to_data, "other_data/data_aov.rds"))

#====================================================================================== -
