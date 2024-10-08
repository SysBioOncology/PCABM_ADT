## PCABM_ADT

Model used in the paper 'Agent-based modeling of the prostate tumor microenvironment uncovers spatial tumor growth constraints and immunomodulatory properties' by MNG van Genderen et. al. https://doi.org/10.1038/s41540-024-00344-6

## Main Folder  
The main folder contains examples of how to run and use PCABM. 

***ABM_run***
Example of how to run a simple simulation of PCABM. Function getHyperParameters determines which set of parameter is used (with or without R1881 or castration resistance). 

***ABM_CRPC_variation***
Run CRPC growth in the TME with different defined ratios of cells.

***ABM_CRPC_spontaneous***
Run 'spontaneous' development of CRPC in the TME. 
This simulation contains a 'full' cycle of PCa growth. Starting with growth in R1881, moving to DMSO (castration) and eventually growth
of resistant cells in DMSO. Resistant tumor cells develop spontaneously upon proliferation of non-resistant cells with probability TUpres.

## Subroutines_ABM
Contains all files to run the prostate cancer agent based model (PCABM). Model was based on ABM for CRC by Jakob Nikolas Kather et al. (2017)
PCABM was adapted to be more prostate cancer specific and contains tumor cells, fibroblast, M1 and M2 type macrophages. 
For detailed explenation on interactions between cell types, please see project report. 

***growTumor*** is the main function that bundles all other functions and runs all parts. 

***getHyperParameters***
Get parameters for condition with R1881, DMSO, CRPC, or LNCaP-abl only in DMSO

***getSystemParamters***
Get all general parameters
Although not used in project report - tumor cells can posibly be stem cells (set parameter TUps to non-zero in getSystemParams.m).

***getAdjacent***
Get cell neighborhood

***shuffleCells***
Randomly order vector of cells

***TU_go_grow_die***
Function for modeling tumor cell migration, growth and death.

***F_go_grow_die***
Function for modeling firoblast growth and migration

***mCellRound***
Function for modeling macrophage round/behavior. Combines all seperate macrophages functions.

***M_go***
Function for modeling macrophage migration

***M_go_die***
Function for modeling macropahge migration & death

***M_kill***
Function for modeling macropahge killing of tumor cells

***M_promote***
NOT USED: Model macrophages promoting tumor cells, but only for tumor cells close to the macrophage. Could be added/replaced in mCellRound in future use.

***updateSystem***
Function to update sytem in each time step (iteration) after all cells performed their round of actions.

***visualizeSystem***
Function to visualize simulation

***writeMyVideo***
Function to write visulized simulation to video file

## Functions
Contains functions for optimizing parameters and visualizing outputs of PCABM. 

***ABM_PSO_function***
Function that runs particle swarm optimization (PSO). PSO fits parameters by comparing best Mean Squared Error (MSE) between Relative Tumor Cell Numbers from Incucyte data to model output for each fit.

***MSE***
Objective function for PSO: the MSE.  

***shadedErrorBar***
Function by Rob Campbell - November 2009, for plotting shaded errobar lineplots. 

## Acknowledgements
Part of the code is adapted from a model by Kather et. al. (2017) http://dx.doi.org/10.5281/zenodo.853342
