# Phenotyping Senescence Dynamics

## Author


> Jonas Anderegg  
> Crop Science Group  
> ETH Zürich  

## Content  

Folder `Utils` contains functions for the pre-processing of spectra, calculation and evaluation of spectral indices, training and evaluation of full-spectrum models, recursive feature elimination and analysis of variance. 

Folder `Analysis`contains scripts to implement the analysis and obtain results contained in the study. 

## Function and Script description

`f_spc_smooth` applies the Savitzky-Golay smoothing filter to raw spectra
`f_spc_avg` averages replicate measurements 
`f_calc_si`calculates spectral indices
`f_calc_sen_si` calculates all spectral indices of the form of the PSRI, to enable the wavelength sensitivity analysis
`f_scale_si` scales spectral indices and predictions obtained from full-spectrum model to range from 0 to 10
`f_spc_bin`
