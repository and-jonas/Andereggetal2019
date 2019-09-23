#====================================================================================== -

#HEADER ----

# AUTHOR: Jonas Anderegg

# Preparation agronomic data and senescence dynamics data 
# for calculation of heritability and BLUEs

#====================================================================================== -

.libPaths("") #Specify path to R Libraries

library(tidyverse)
library(readxl)

path_to_data <- "" #Specify path to data

#create sink directory
if(!dir.exists(paste0(path_to_data, "other_data/spctempfeat"))){
  dir.create(paste0(path_to_data, "other_data/spctempfeat"))
} else NULL
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
  filter(Year == "2017") %>%
  dplyr::select(one_of(dynpars), contains("SnsCnp"), GY) %>%
  filter(complete.cases(.))
saveRDS(data, paste0(sinkdir, "spctempgy2017.rds"))

# GY FPWW022
data <- BLUEs %>%
  filter(Year == "2018") %>%
  dplyr::select(one_of(dynpars), contains("SnsCnp"), GY) %>%
  filter(complete.cases(.))
saveRDS(data, paste0(sinkdir, "spctempgy2018.rds"))

# GPC FPWW012
data <- BLUEs %>%
  filter(Year == "2016") %>%
  dplyr::select(one_of(dynpars), contains("SnsCnp"), GPC) %>%
  filter(complete.cases(.))
saveRDS(data, paste0(sinkdir, "spctempgpc2016.rds"))

# GPC FPWW018
data <- BLUEs %>%
  filter(Year == "2017") %>%
  dplyr::select(one_of(dynpars), contains("SnsCnp"), GPC) %>%
  filter(complete.cases(.))
saveRDS(data, paste0(sinkdir, "spctempgpc2017.rds"))

#====================================================================================== -