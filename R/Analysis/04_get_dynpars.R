#====================================================================================== -

#HEADER ----

# Author: Jonas Anderegg, ETH Zürich
# Copyright (C) 2019  ETH Zürich, Jonas Anderegg (jonas.anderegg@usys.ethz.ch)

# Extract senescence dynamics parameters from scorings
# using linear interpolation and gompertz models.


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
library(nls.multstart)
library(doParallel)

path_to_data <- "" #Specify path to data
path_to_utils <- "" #Specify path to functions

#create sink directory
if(!dir.exists(paste0(path_to_data, "OUT/dynpars"))){
  dir.create(paste0(path_to_data, "OUT/dynpars"))
} else NULL
sinkdir <- paste0(path_to_data, "OUT/dynpars/") #Specify output directory

source(paste0(path_to_utils, "001_spectra_utils.R"))

#====================================================================================== -

#load data
data <- readRDS(paste0(path_to_data, "scr_sca.rds"))

# extract DynPars: Scorings ----
data <- data %>% 
  tidyr::gather(Trait, Score, SnsFl0:SnsCnp, factor_key = TRUE) %>% 
  filter(complete.cases(.)) %>% 
  data.frame()

# do linear interpolation of scorings
# and fit nls using multstart
# Performs 250 times repeated NLS fitting (Levenberg-Marquardt algorithm)
# with random-search start parameter sets randomly sampled from a uniform
# distribution between upper and lower starting parameter bounds
data_fits <- data %>%
  nest(-c(Plot_ID, Trait)) %>%
  group_by(Plot_ID, Trait) %>%
  mutate(fit_lin = purrr::map(data, lin_approx)) %>% 
  mutate(fit_cgom = purrr::map(data, ~ nls_multstart(Score ~ Gompertz_constrained(b, M, tempsum = grading_GDDAH),
                                                     data = .x,
                                                     iter = 250,
                                                     start_lower = c(b = -0.1, M = 550),
                                                     start_upper = c(b = 0, M = 750),
                                                     convergence_count = 100)))

# new data frame of predictions
new_preds <- data %>%
  do(., data.frame(grading_GDDAH = seq(min(.$grading_GDDAH), max(.$grading_GDDAH), length.out = 1000), stringsAsFactors = FALSE))

# max and min for each curve
max_min <- group_by(data, Plot_ID) %>%
  summarise(., min_gGDDAH = min(grading_GDDAH), max_gGDDAH = max(grading_GDDAH)) %>%
  ungroup()

# create new predictions
preds2 <- data_fits %>%
  tidyr::unnest(fit_cgom %>% purrr::map(broom::augment, newdata = new_preds)) %>%
  merge(., max_min, by = "Plot_ID") %>%
  group_by(., Plot_ID) %>%
  filter(., grading_GDDAH > unique(min_gGDDAH) & grading_GDDAH < unique(max_gGDDAH)) %>%
  arrange(., Plot_ID, Trait, grading_GDDAH) %>%
  rename(., Score = .fitted) %>%
  ungroup()

# check whether model converged
convInfo <- data_fits %>%
  transmute(convInfo = purrr::map(purrr::map(fit_cgom, "convInfo"), "isConv"))

#extract parameters from nls fits
preds_gom <- preds2 %>% 
  group_by(Plot_ID, Trait) %>%
  nest() %>%
  mutate(pars_gom = purrr::map(data, extract_pars)) %>%
  select(-data)

#extract parameters from linear interpolation
preds_lin <- data_fits %>%
  unnest(fit_lin) %>% 
  mutate(data_lin = purrr::map(fit_lin, bind_cols)) %>%
  transmute(pars_lin = purrr::map(data_lin, extract_pars))

#combine predictions and convInfo
pred <- list(preds_gom, preds_lin, convInfo) %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by=c("Plot_ID", "Trait")), .)

#combine predictions and convInfo
df <- pred %>% unnest(pars_gom, pars_lin, convInfo)
names(df)[4:length(df)] <- c("onsen_gom", "midsen_gom", "endsen_gom", "tsen_gom",
                             "onsen_lin", "midsen_lin", "endsen_lin", "tsen_lin")

#replace parameter value with NA if model did not converge
df[df$convInfo == FALSE, c(4:7)] <- NA

saveRDS(df, file = paste0(sinkdir, "dynpars_scr.rds"))

#====================================================================================== -
