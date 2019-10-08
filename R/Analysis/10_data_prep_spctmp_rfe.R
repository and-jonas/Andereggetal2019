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
library(readxl)

path_to_data <- "" #Specify path to data

#create sink directory
if(!dir.exists(paste0(path_to_data, "other_data/spctempfeat"))){
  dir.create(paste0(path_to_data, "other_data/spctempfeat"))
}
sinkdir <- paste0(path_to_data, "other_data/spctempfeat/") #Specify output directory

#====================================================================================== -

# read data
BLUEs <- readRDS(paste0(path_to_data, "OUT/BLUEs.rds"))
dynpars <- read.csv(paste0(path_to_data, "helper_files/dynpars_decorr.csv")) %>% pull(x) %>% as.character()
dynpars <- c(dynpars, "onsen_lin_SI_NDVI_nb_ASD",
             "midsen_lin_SI_NDVI_nb_ASD",
             "endsen_lin_SI_NDVI_nb_ASD",
             "tsen_lin_SI_NDVI_nb_ASD")

# GY FPWW012
data <- BLUEs %>%
  filter(Exp == "FPWW012") %>%
  dplyr::select(one_of(dynpars), contains("SnsCnp"), GY) %>%
  filter(complete.cases(.))
data <- data[!data$GY < 3, ]
saveRDS(data, paste0(sinkdir, "spctempgy2016.rds"))

# GY FPWW018
data <- BLUEs %>%
  filter(Exp == "FPWW018") %>%
  dplyr::select(one_of(dynpars), contains("SnsCnp"), GY) %>%
  filter(complete.cases(.))
saveRDS(data, paste0(sinkdir, "spctempgy2017.rds"))

# GY FPWW022
data <- BLUEs %>%
  filter(Exp == "FPWW022") %>%
  dplyr::select(one_of(dynpars), contains("SnsCnp"), GY) %>%
  filter(complete.cases(.))
saveRDS(data, paste0(sinkdir, "spctempgy2018.rds"))

# GPC FPWW012
data <- BLUEs %>%
  filter(Exp == "FPWW012") %>%
  dplyr::select(one_of(dynpars), contains("SnsCnp"), GPC) %>%
  filter(complete.cases(.))
saveRDS(data, paste0(sinkdir, "spctempgpc2016.rds"))

# GPC FPWW018
data <- BLUEs %>%
  filter(Exp == "FPWW018") %>%
  dplyr::select(one_of(dynpars), contains("SnsCnp"), GPC) %>%
  filter(complete.cases(.))
saveRDS(data, paste0(sinkdir, "spctempgpc2017.rds"))

#====================================================================================== -