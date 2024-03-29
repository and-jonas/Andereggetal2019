#====================================================================================== -

#HEADER ----

# Author: Jonas Anderegg, ETH Z�rich
# Copyright (C) 2019  ETH Z�rich, Jonas Anderegg (jonas.anderegg@usys.ethz.ch)

# Prepare spectral datasets for subsequent analyses


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
library(prospectr)
library(stringr)
library(data.table)

path_to_utils <- "" #Specify path to utils folder
path_to_data <- "" #Specify path to data folder

#create sink directory
if(!dir.exists(paste0(path_to_data, "OUT"))){
  dir.create(paste0(path_to_data, "OUT"))
}
sinkdir <- paste0(path_to_data, "OUT/") #Specify output directory

source(paste0(path_to_utils, "001_spectra_Utils.R"))

#====================================================================================== -

#Load data ----

#spectral reflectance data
spc <- data.table::fread(paste0(path_to_data, "spc_data/spc_raw.csv")) %>% 
  mutate(Plot_ID = as.factor(Plot_ID),
         meas_date = as.Date(meas_date),
         replicate = as.factor(replicate)) %>% 
  as_tibble()

#scaled visual scorings for SI and model performance
scr_sca <- data.table::fread(paste0(path_to_data, "scr_data/scr_sca.csv")) %>% 
  dplyr::mutate(Plot_ID = as.factor(Plot_ID),
         Gen_Name = as.factor(Gen_Name),
         grading_date = as.Date(grading_date),
         heading_date = as.Date(heading_date)) %>% 
  as_tibble()

#unscaled visual scorings for training and validating spectral models
scr_usca <- data.table::fread(paste0(path_to_data, "scr_data/scr_usca.csv")) %>% 
  as_tibble() %>% 
  mutate(Plot_ID = as.factor(Plot_ID),
         Gen_Name = as.factor(Gen_Name),
         grading_date = as.Date(grading_date),
         heading_date = as.Date(heading_date),
         SnsFl0 = as.double(SnsFl0),
         SnsCnp = as.double(SnsCnp)) %>% 
  as_tibble()

#GDDAH for measurements and scorings (plot-specific)
gddah <- read.csv(paste0(path_to_data, "helper_files/gddah_data.csv")) %>% 
  mutate(meas_date = meas_date %>% as.Date()) %>% 
  as_tibble()

#template of matching dates to join spectral and scoring data
match_dates <- read.csv(paste0(path_to_data, "helper_files/match_dates.csv"), sep = ";") #load matching dates template

#lodging plots, to be excluded from the analysis
lodg_plots <- read.csv(paste0(path_to_data, "other_data/lodging.csv"), sep = ";") %>% 
  filter(Lodging %in% c("l", "ll")) %>% pull(Plot_ID)

#====================================================================================== -

#Spectral Index datasets ----

#create sink directory
if(!dir.exists(paste0(sinkdir, "SI"))){
  dir.create(paste0(sinkdir, "SI"))
}

#remove lodging plots
spc <- spc %>% 
  filter(!Plot_ID %in% lodg_plots)

#create dataset to extract dynamics parameters
##This dataset comprises all measurement dates available,
##except for those carried out prior to heading
spc %>%
  #smooth spectra using the Savitzky-Golay smoothing filter
  f_spc_smooth(3, 11, 0) %>%
  #average spectra for each plot
  f_spc_avg() %>% 
  #calculate spectral vegetation indices
  f_calc_si() %>%
  #scale values of spectral vegetation indices;
  #use only measurements carried out AFTER heading
  f_scale_si(sub = "post_heading") %>% 
  #add gddah data
  left_join(., gddah, by = c("Plot_ID", "meas_date")) %>% 
  #rearrange
  dplyr::select(-contains("SI"), everything()) %>% 
  saveRDS(paste0(sinkdir, "SI/SI_forpars.rds"))

#create dataset to quantify index performance
#in approximating visually observed senescence dynamics
##this dataset is limited to the common temporal range 
##of scorings and spectral measurements in each year
spc %>%
  #smooth spectra using the Savitzky-Golay smoothing filter
  f_spc_smooth(3, 11, 0) %>%
  #average spectra for each plot
  f_spc_avg() %>% 
  #calculate spectral vegetation indices
  f_calc_si() %>%
  #scale values of spectral vegetation indices
  f_scale_si(sub = "post_heading") %>% 
  #add scorings and gddah data
  f_match_join(., scr_sca, gddah, match_dates) %>% 
  #rearrange
  dplyr::select(-contains("SI_"), everything()) %>% 
  # #drop rows with missing data
  filter(!is.na(Exp)) %>%
  as_tibble() %>% 
  saveRDS(paste0(sinkdir, "SI/SI_forperf.rds"))

#create dataset of PSRI' with different constituting wavelengths
##this dataset is limited to the common temporal range 
##of scorings and spectral measurements in each year
spc %>%
  #smooth spectra using the Savitzky-Golay smoothing filter
  f_spc_smooth(3, 11, 0) %>%
  #average spectra for each plot
  f_spc_avg() %>% 
  #calculate spectral vegetation indices
  f_calc_sen_si(wvlt_d = seq(400, 850, 3), wvlt_n = seq(400, 850, 3)) %>%
  #scale values of spectral vegetation indices
  f_scale_si(sub = "post_heading") %>% 
  #add scorings and gddah data
  f_match_join(., scr_sca, gddah, match_dates) %>% 
  #rearrange
  dplyr::select(-starts_with("SI"), everything()) %>% 
  as_tibble() %>% 
  #some wvlt combinations produce highly unstable values
  #these can be recognized by their non-reversion
  #these can be excluded
  dplyr::select(Exp:ref_date, ends_with("_r")) %>% 
  saveRDS(paste0(sinkdir, "SI/SI_PSRI.rds"))
  
#====================================================================================== -

#Datasets for full-spc models ----

#create sink directory
if(!dir.exists(paste0(sinkdir, "MM"))){
  dir.create(paste0(sinkdir, "MM"))
}

#average smoothed reflectance spectra,
#rflt values at adjacent wavelengths are binned
spc %>%
  f_spc_smooth(3, 11, 0) %>%
  f_spc_avg() %>%
  f_spc_bin(bin_size = 3) %>%
  f_spc_trim() %>%
  f_match_join(., scr_usca, gddah, match_dates) %>%
  filter(complete.cases(.)) %>% 
  as_tibble() %>% 
  saveRDS(paste0(sinkdir, "MM/smth_avg_rflt_bin3.rds")) 

#average smoothed first derivative spectra,
#rflt values at adjacent wavelengths are binned
spc %>%
  f_spc_smooth(3, 11, 1) %>%
  f_spc_avg() %>%
  f_spc_bin(bin_size = 3) %>%
  f_spc_trim() %>%
  f_match_join(.,  scr_usca, gddah, match_dates) %>%
  filter(complete.cases(.)) %>% 
  as_tibble() %>% 
  saveRDS(paste0(sinkdir, "MM/avg_bin3_smth_der1.rds"))

#average smoothed continuum-removed reflectance spectra,
#rflt values at adjacent wavelengths are binned
spc %>%
  f_spc_avg() %>%
  f_spc_bin(bin_size = 3) %>%
  f_spc_trim() %>%
  f_cont_rem() %>%
  f_spc_smooth(3, 11, 0) %>%
  f_match_join(.,  scr_usca, gddah, match_dates) %>%
  filter(complete.cases(.)) %>%   
  as_tibble() %>% 
  saveRDS(paste0(sinkdir, "MM/avg_rflt_bin3_cr_smth.rds"))

#average smoothed reflectance spectra,
#rflt values in the spectral range between 500nm and 700nm 
spc %>%
  f_spc_avg() %>%
  f_spc_trim() %>%
  f_match_join(.,  scr_usca, gddah, match_dates) %>%
  filter(complete.cases(.)) %>%
  dplyr::select(-contains("rflt_"), rflt_500:rflt_700) %>% 
  as_tibble() %>% 
  saveRDS(paste0(sinkdir, "MM/avg_rflt_restr.rds"))

#====================================================================================== -
