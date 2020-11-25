# ACT_PAS_OHC_VAR

This repository provides the necessary source code, run files and diagnostic scripts to recreate the experiments in our manuscript
"The active and passive roles of the ocean in generating basin-scale heat content variability"
Also required are model outputs from the IPSL-CM5a run used to create our stochastic representation.

There are three directories:
## `1_CREATE_STOCHASTIC_REPRESENTATION`
This contains three python scripts for taking IPSL-CM5a model output, creating a covariance matrix of fluxes and an array of local flux decorrelation times, and combining the two to create the stochastic representation (Sigma in our manuscript)

## `2_RUN_ADJOINT_MODEL`
This contains two directories: 
- `MY_SRC`, which contains modified source code for the adjoint model and instructions for installation
- `EXPERIMENTS`, which contains namelist files to generate a trajectory and run the adjoint model, as well as a python script to create the adjoint input (heat content cost functions for different basins and depth levels). Instructions are contained in the directory.

## `3_DIAGNOSE_RESPONSE_VARIANCE`
This contains the python script `calculate_variance.py` which takes the adjoint output produced in step 2 and the stochastic representation produced in step 1 and produces a netCDF file containing the distributions (due to buoyancy and wind forcing) of variance accumulated during each of the 60 years of the adjoint run.

Instructions and README files can be found in the relevant locations within the repository, and I am also happy to be contacted at the corresponding author address to provide any further assistance.
