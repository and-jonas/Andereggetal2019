#====================================================================================== -

#HEADER ----

# Author: Jonas Anderegg, ETH Zürich
# Copyright (C) 2019  ETH Zürich, Jonas Anderegg (jonas.anderegg@usys.ethz.ch)

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

# main functions ---- 

#Wrapper for SpATS
f_spats <- function(data, response = response, random = random, fixed = NULL, 
                    genotype.as.random = genotype.as.random, genotype = genotype) {
  SpATS(response = response, random = as.formula(random), fixed = as.formula(fixed), 
        spatial = ~PSANOVA(RangeBL, RowBL, nseg = c(20,20), nest.div=c(2,2)), 
        genotype = genotype, genotype.as.random = genotype.as.random, data = data,
        control = list(maxit = 100, tolerance = 1e-03, monitoring = 0))
}

get_h2 <- function(obj){
  ifelse(length(unique(obj$data$Rep))> 1, SpATS::getHeritability(obj), NA)
}

#extract BLUE from SpATS
get_BLUE_spats <- function(object) {
  fitted <- object$fitted
  intercept <- object$coeff['Intercept']
  gen_mod_mat <- construct_genotype_prediction_matrix(object, object$data)
  gen_coeff <- object$coeff[1:ncol(gen_mod_mat)]
  geno_pred <- as.vector(gen_mod_mat %*% gen_coeff)
  BLUE <- as.data.frame(intercept + gen_coeff) %>% 
    rownames_to_column() %>% 
    as_tibble() %>% 
    dplyr::rename(Gen_Name = 1, BLUE = 2)
}

#extract spatially corrected plot values from SpATS
spat_corr_spats <- function(object) {
  fitted <- object$fitted
  intercept <- object$coeff['Intercept']
  gen_mod_mat <- construct_genotype_prediction_matrix(object, object$data)
  gen_coeff <- object$coeff[1:ncol(gen_mod_mat)]
  geno_pred <- as.vector(gen_mod_mat %*% gen_coeff)
  residual <- object$residuals
  plot_corr <- as.data.frame(intercept + geno_pred + residual) %>% rename(spat_corr = 1) %>% 
    bind_cols(object$data, .) %>% as_tibble() %>% dplyr::select(-value)
}

#calculate heritability, using spatially corrected plot values
get_h2_asreml <- function(data, fixed, random, residual, cullis = TRUE){
  if(!is.factor(data$Gen_Name)){
    data$Gen_Name <- as.factor(data$Gen_Name)
  }
  as <- asreml(fixed = as.formula(fixed),
               random = as.formula(random),
               residual = as.formula(residual),
               data = data,
               trace = FALSE
               ,na.action = na.method(x = "include", y = "include")
               )
  
  if(cullis){
    GenVar <- summary(as)$varcomp["Gen_Name","component"]
    #Predict genotypic mean values to estimate avsed
    pv.G <- predict(as, classify = "Gen_Name")
    avsed.G<-pv.G$avsed
    #The Cullis version of heritability based on avsed^2 and V.G
    h2<- 1-(avsed.G^2/(2*GenVar))
  } else {
    #number of Years
    n1 <- data %>% group_by(Exp) %>% nest() %>% nrow
    #number of Year:Replicate combinations
    n2 <- data %>% group_by(Exp, Lot) %>% nest() %>% nrow
    GenVar <- summary(as)$varcomp["Gen_Name","component"]
    ErrVar <- summary(as)$varcomp["units!R","component"]
    GenExpVar <- summary(as)$varcomp["Gen_Name:Exp","component"]
    h2 <- GenVar/(GenVar +  GenExpVar/n1 + ErrVar/n2)
  }
}

#calculate heritability, using year-wise BLUEs
get_h2_asreml2 <- function(data, fixed, random, residual, cullis = TRUE){
  if(!is.factor(data$Gen_Name)){
    data$Gen_Name <- as.factor(data$Gen_Name)
  }
  as <- asreml(fixed = as.formula(fixed),
               random = as.formula(random),
               residual = as.formula(residual),
               data = data,
               trace = FALSE)
  if(cullis){
    GenVar <- summary(as)$varcomp["Gen_Name","component"]
    #Predict genotypic mean values to estimate avsed
    pv.G <- predict(as, classify = "Gen_Name")
    avsed.G<-pv.G$avsed
    #The Cullis version of heritability based on avsed^2 and V.G
    h2<- 1-(avsed.G^2/(2*GenVar))
  } else {
    #number of Years
    n <- data %>% group_by(Exp) %>% nest() %>% nrow
    GenVar <- summary(as)$varcomp["Gen_Name","component"]
    ErrVar <- summary(as)$varcomp["units!R","component"]
    h2 <- GenVar/(GenVar + ErrVar)
  }
}



# helper functions ----

#remove reps not measured
remove_reps <- function(data, value, missing){
  remove <- ifelse(sum(is.na(data %>% pull(value)), na.rm = TRUE) > missing, TRUE, FALSE)
}

#helper function
make_design_matrix <- function(geno, names) {
  Nrow <- length(geno)
  Ncol <- length(names)
  col <- match(geno, names)
  frame <- data.frame(i = c(1:Nrow), j = col, v = rep(1,Nrow))
  frame <- subset(frame, is.na(col) == FALSE)
  L <- as.list(frame)
  X <- spam::spam(L, nrow = Nrow, ncol = Ncol)
  return(X)
}

#helper function
construct_genotype_prediction_matrix <- function(object, newdata) {
  Z_geno = make_design_matrix(newdata[,object$model$geno$genotype], object$terms$geno$geno_names)
  if(object$model$geno$as.random)
    Z_geno <- Z_geno[,object$terms$geno$ndx]
  else
    Z_geno <- Z_geno[, object$terms$geno$ndx[2:length(object$terms$geno$ndx)]]
  Z_geno	
}
