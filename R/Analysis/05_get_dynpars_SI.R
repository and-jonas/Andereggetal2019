#====================================================================================== -

#HEADER ----

# Author: Jonas Anderegg, ETH Zürich
# Copyright (C) 2019  ETH Zürich, Jonas Anderegg (jonas.anderegg@usys.ethz.ch)

# Extract dynamics parameters from spectral indices
# using linear interpolation.


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

list.of.packages <- c("doParallel", "Rmpi", "future", "furrr", "tidyverse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE, repos = "https://stat.ethz.ch/CRAN/")

library(doParallel)
library(future)
library(furrr)
library(tidyverse)

path_to_data <- "" #Specify path to data
path_to_utils <- "" #Specify path to functions

#create sink directory
if(!dir.exists(paste0(path_to_data, "OUT/dynpars"))){
  dir.create(paste0(path_to_data, "OUT/dynpars"))
} else NULL
sinkdir <- paste0(path_to_data, "OUT/dynpars/") #Specify output directory

source(paste0(path_to_utils, "001_spectra_utils.R"))

#====================================================================================== -

data <- readRDS(paste0(path_to_data, "OUT/SI/SI_forpars.rds"))

#remove incomplete time series
#no meaningful parameters can be extracted from these plots
data <- data[!data$Lot == 6, ] %>% droplevels()
data <- data[!(data$Lot == 3 & data$Exp == "FPWW022"), ]

plan("multiprocess")
pars <- data %>% gather(SVI, value, contains("SI_")) %>% 
  group_by(Exp, SVI, Plot_ID) %>% 
  nest() %>% 
  mutate(eval = furrr::future_map(.x = data, method = "lin",
                                  .f = possibly(get_dynpars_SI, otherwise = NA_real_)))
plan("sequential")

#tidy up output
pars_out <- pars %>% dplyr::select(-data) %>% 
  mutate(remove = purrr::map_dbl(.x = dynpars, .f = length)) %>% 
  #remove plots where no data is available
  dplyr::filter(remove == 4) %>% dplyr::select(-remove) %>% 
  unnest()

#save to directory
saveRDS(pars_out, paste0(sinkdir, "dynpars_SI.rds"))

#====================================================================================== -
