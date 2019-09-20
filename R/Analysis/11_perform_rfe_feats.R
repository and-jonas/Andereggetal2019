#====================================================================================== -

#HEADER ----

# Author: Jonas Anderegg, ETH Zürich
# Copyright (C) 2019  ETH Zürich, Jonas Anderegg (jonas.anderegg@usys.ethz.ch)

# Perform recursive feature elimination using rf as base learner
# Get robust feature ranks
# create performance profile plots


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

list.of.packages <- c("doParallel", "Rmpi", "caret", "Cubist", "tidyverse", "ranger")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE, repos='https://stat.ethz.ch/CRAN/')

library(doParallel)
library(caret)
library(Cubist)
library(ranger)
library(tidyverse)

path_to_data <- "" #Specify path to data
path_to_utils <- "" #Specify path to functions

#create sink directory
if(!dir.exists(paste0(path_to_data, "OUT/out_rfe_feats"))){
  dir.create(paste0(path_to_data, "OUT/out_rfe_feats"))
} else NULL
sinkdir <- paste0(path_to_data, "OUT/out_rfe_feats/") #Specify output directory

#====================================================================================== -

data <- readRDS(paste0(path_to_data, "other_data/spctempfeat/spctempgy2016.rds"))

# The candidate set of the number of predictors to evaluate                      
subsets <- c(length(data)-1, seq(90, 20, -5), 17, 14, 12:1)

# Set up a cluster
p = strtoi(Sys.getenv('LSB_DJOB_NUMPROC'))
print(paste0("Detected cores: ", p)) #check cluster set up
cluster <- makeCluster(p - 1, type = "MPI", outfile = "")
registerDoParallel(cluster)
clusterEvalQ(cluster, {
  library(MASS)
  library(Cubist)
  library(caret)
  library(tidyverse)
})

# Perform recursive feature elimination
rfe <- perform_rfe(response = "GY", base_learner = "ranger",
                   p = 0.8, times = 30, groups = 5, 
                   subsets = subsets, data = data,
                   importance = "permutation",
                   num.trees = 10000)

stopCluster(cluster)
registerDoSEQ()

saveRDS(rfe, paste0(sinkdir, "rfe_out_feats.rds"))

#====================================================================================== -
