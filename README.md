# Phenotyping Senescence Dynamics

This repository contains functions and analysis scripts required to reproduce results reported by Anderegg *et al*. (2019). 

## Author


> Jonas Anderegg  
> Crop Science Group  
> ETH Zürich  


## Dependencies

Analysis of variance is performed using the R-package `asreml`. Computationally demanding steps in the analysis are designed to run on a cluster, requiring the R-packages `doParallel`, `future`, and `furrr`. Scripts build primarily on functions from the R packages `prospectr`, `caret` and `tidyverse`.

## Content  

Folder `Utils` contains functions for the pre-processing of spectra, calculation and evaluation of spectral indices, training and evaluation of full-spectrum models, recursive feature elimination and analysis of variance. 

Folder `Analysis`contains scripts to implement the analysis and obtain results contained in the study. 

### Utils

#### 000_helper_funcs.R

1.-3. Helper functions to remove incomplete datasets (Lots, Years, measurement series). 

#### 001_spectra_utils.R

1. `f_spc_smooth` applies the Savitzky-Golay smoothing filter to raw spectra (wrapper for `prospectr::SavitzkyGolay`)
2. `f_spc_avg` averages replicate measurements
3. `f_calc_si`calculates spectral indices
4. `f_calc_sen_si` calculates all spectral indices of the form of the PSRI, to enable the wavelength sensitivity analysis
5. `f_scale_si` scales spectral indices and predictions obtained from full-spectrum model to range from 0 to 10
6. `f_spc_bin` computes average values of a signal in pre-determined bins (wrapper for `prospectr::binning`)
7. `f_spc_trim` removes noisy parts of the signal in pre-determined ranges
8. `f_cont_rem` performs continuum removal (wrapper for `prospectr::ContinuumRemoval`)
9. `f_match_join` joins spectral and scoring datasets by assessment/measurement dates, using a supplied template of matching dates. This function also adds growing degree days for scoring and measurement time points.
10. `get_dynpars_SI` extracts dynamics parameters from spectral indices. 
11. `get_errors_and_dynpars`Performs linear interpolation of scorings and index values (data in growing degree days after heading),  extracts dynamics parameters and calculates overall error between the two fitted curves. Parametric models instead of linear interpolation are also supported.
12. Several data wrangling helper functions. 


#### 002_full_spc_utils.R

1. `get_rmse` calculates the root mean square error (RMSE) 
2. `perform_sampling` performs up- or down-sampling to create balanced evaluation datasets (wrapper for `caret::upSample` or `caret::downSample`)
3. `get_RMSE_sample` calculates the root mean square error (RMSE) of up- or downsamples
4. `extract_predobs` creates predictions using a fitted model and extracts the corresponding real observations
5. `eval_full_spc_cross` fits full-spectrum models, extracts performance metrics, performs sampling (up- or downsampling of both training and test data) and extracts predictions and real observations for these samples, performs all possible variants of leave-year(s)-out cross-validation. This is essentially a wrapper function for `pls::plsr` and `Cubist::cubist` via the `caret::train` interface. This function is designed to run in parallel on the high performance computing cluster of ETH Zürich. Must be serialized if no access to such a facility is available. Very long computation times should be expected. 
6. Several data wrangling helper functions. 

#### 003_rfe_utils.R

1. `perform_rfe` performs recursive feature elimination using random forest or cubist regression as base learners. This is essentially a wrapper function for `ranger::ranger` and `Cubist::cubist` using the`caret::train` interface.
2. `tidy_rfe_output` gathers results of resamples and creates summary statistics of performance and feature ranks.
3. `plot_perf_profile` creates simple performance profile plots.

#### 004_spat_aov_utils.R

1. `f_spats` fits the SpATS model to raw plot-based data (wrapper for `SpATS::SpATS`)
2. `get_h2` reports repeatability, or `NA`if no replicate measurements are available. This is a wrapper for `SpATS::getHeritability`. 
3. `spat_corr_spats` computes spatially corrected values for each experimental plot from the spats object (`object$intercept` + `object$coeff`+ `object$residuals`).  
4. `get_BLUE_spats` computes the best linear unbiased estimator (BLUE) for each genotype. 
5. `get_h2_asreml` calculates heritability across years using spatially corrected plot values and the formula proposed by Cullis. 
6. `get_h2_asreml2` calculates heritability across years using best linear unbiased estimators. 
7. Data wrangling helper functions. 

### Scripts

1. `01_prep_spcdat.R` Prepare spectral datasets for evaluation of spectral indices and full-spectrum models.
2. `02_full_spc_cv.R` Train and evaluate full-spectrum models
3. `03_perform_rfe.R` An example script that performs rfe using cubist as base learner for the selection of the most important wavelengths to predict visual senescence scorings for the 2016 experiment. This script is designed to run on the high performance computing cluster of ETH Zürich. It must be serialized if no access to such a facility is available. Expect very long computation times. 
4. `04_get_dynpars.R` Extract dynamics parameters from visual senescence scorings. 
5. `05_get_dynpars_SI.R` Extract dynamics parameters from spectral indices.
6. `06_data_prep_aov.R` Assemble all data for correction of spatial heterogeneity. 
7. `07_spatial_w2_h2.R` Correct for spatial heterogeneity, calculate repeatability (w2) and across-year heritability (h2). Calculate best linear unbiased estimators. 
8. `08_SI_perf.R` Evaluate performance of spectral indices in tracking visually observed senescence dynamics. 
9. `09_full_spc_perf.R` Evaluate performance of full-spectrum models in tracking visually observed senescence dynamics. 
10. `10_data_prep_spctemp_rfe.R` Assemble data sets for feature selection
11. `11_perform_rfe_feats.R` Perform feature selection. 
