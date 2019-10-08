#====================================================================================== -

#HEADER ----

# Author: Jonas Anderegg, ETH Zürich
# Copyright (C) 2019  ETH Zürich, Jonas Anderegg (jonas.anderegg@usys.ethz.ch)

# Fit full-spectrum models to predict senescence scorings from reflectance spectra
# Evaluate across-year applicability of models through leave-year(s)-out cross-validation

# This analysis was performed on a high performance computing cluster of ETH Zürich;
# Script is designed to run on a cluster with >= 11 cores;
# Serialize if such a cluster is not available to you, but expect long computation times. 


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

list.of.packages <- c("doParallel", "Rmpi", "caret", "Cubist", "pls", "plyr", "tidyverse", "purrrlyr", "tictoc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE, repos='https://stat.ethz.ch/CRAN/')

library(tidyverse)
library(caret)
library(parallel)
library(doParallel)
library(Cubist)
library(pls)
library(purrrlyr)
library(tictoc)

path_to_utils <- "" #Specify path to functions 
path_to_data <- "" #Specify path to data

#create sink directory
if(!dir.exists(paste0(path_to_data, "OUT/model_output"))){
  dir.create(paste0(path_to_data, "OUT/model_output"))
}
sinkdir <- paste0(path_to_data, "OUT/model_output/") #Specify output directory

source(paste0(path_to_utils, "002_full_spc_utils.R"))

#====================================================================================== -

#Load data ----

data_rflt_bin3 <- readRDS(paste0(path_to_data, "OUT/MM/smth_avg_rflt_bin3.rds"))
data <- list("rflt_bin3" = data_rflt_bin3)

#====================================================================================== -

#Specify function arguments ----

# Experiments
Exp <- c("FPWW012", "FPWW018", "FPWW022")

# create train list of length 7
train <- list(Exp[1], Exp[2], Exp[3], Exp[1:2], Exp[2:3], Exp[c(1, 3)], Exp[1:3])

# create validation list of length 7
test <- list(c(list(Exp[2], Exp[3], Exp[2:3])),
             c(list(Exp[1], Exp[3], Exp[c(1, 3)])),
             c(list(Exp[1], Exp[2], Exp[1:2])),
             c(list(Exp[3])),
             c(list(Exp[1])),
             c(list(Exp[2])),
             c(list(Exp[1], Exp[2], Exp[3])))

Grid <- expand.grid(committees = c(1, 2, 5, 10, 20),
                        neighbors = 0)

#====================================================================================== -

#Run ----
# Create a grid with all possible argument combinations
args <- expand.grid(Trait = "SnsCnp",
                    method = c("pls", "cubist"),
                    data_type = "rflt_bin3", #other datatypes can be added here, to apply function sequentially
                    trainsample = "fullsample",
                    testsample = "upsample",
                    preProc = "center_scale",
                    stringsAsFactors = FALSE)

# Set up a cluster
p = strtoi(Sys.getenv('LSB_DJOB_NUMPROC'))
print(paste0("Detected cores: ", p)) #check cluster set up
cluster <- makeCluster(p - 1, type = "MPI", outfile = "")
registerDoParallel(cluster)
clusterEvalQ(cluster, {
  library(MASS)
  library(caret)
  library(tidyverse)
})

# sequentially, run function for all argument combinations
result <- args %>% purrrlyr::invoke_rows(.f = eval_full_spc_cross)

# close cluster
stopCluster(cluster)
registerDoSEQ()

# save result
saveRDS(result, paste0(sinkdir, "result_1.rds")) #save results for other trait combinations in the same folder, name output files result_*.rds

#====================================================================================== -