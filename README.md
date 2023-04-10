# Description

This repository contains the Matlab prototype for the VEXAT muscle model described in the paper Millard et al. The model is named after the viscoelastic (VE) crossbridge (X) active-titin (AT) elements contained in the model

Matthew Millard, David W. Franklin, Walter Herzog. A three filament mechanistic model of musculotendon %force and impedance. bioRxiv 2023.03.27.534347; doi: https://doi.org/10.1101/2023.03.27.534347 

# Quick start guide:

Execute 'main_OuterLoop.m' from Matlab to run everything. Roughly 6 hours (Intel i7-3630QM @ 2.40 GHz, Ubuntu 22 8 GB ram, SSD harddrive) is needed to run all of the experiments. The experiments require roughly 500 simulations, some of which are numerically stiff. At the end if this round of simulation, you can find all of the figures that appear in Millard et al., and many more besides, in the folders

output/
- HerzogLeonard2002
- InitializationBenchmark
- KirschBoskovRymer1994
- LeonardJoumaaHerzog2010
- MuscleCurves
- NettiDamoreRoncaAmbrosioNicolais1996
- StandardTests

# Repostiory Overview

1. The muscle models used for the experiments consist of a cat soleus, a rabbit psoas, and a human soleus. These models are created and fitted in the function main_CreateModels_OuterLoop.m.

2. Each experiment is accompanied by a script (e.g. main_HerzogLeonard2002.m) and an script-harness to run the script with multiple configurations (e.g. main_HerzogLeonard2002_OuterLoop.m). Please refer to these files to see exactly how these experiments are simulated.

3. These scripts are not run as a part of main_OuterLoop.m, but are important nonetheless:
  - This function compares the response of the viscoelastic tendon model to the data of Netti et al.: main_NettiDamoreRoncaAmbrosioNicolais1996_TendonDamping.m
  - This function numerically evaluates the quality of the initialization routine: main_InitializationBenchmark.m

4. The model is implemented in 2 files: models/calcMillard2023VexatMuscleInfo.m and models/updateMillard2023VexatCache.m. The first file (models/calcMillard2023VexatMuscleInfo.m) acts as a wrapper file which takes in the inputs from the user, solves for the solution of the models state derivative, and populates the output structures. All of the numerical number crunching needed to evaluate the various computational stages of the model appear in  (models/updateMillard2023VexatCache.m).

5. The Hill model that is used for comparison purposes appears in model/calcMillard2012DampedEquilibriumMuscleInfo.m, and this is a Matlab implementation of the model described in 

    Millard, M., Uchida, T., Seth, A., & Delp, S. L. (2013). 
    Flexing computational muscle: modeling and simulation of 
    musculotendon dynamics. Journal of biomechanical engineering, 
    135(2), 021005.

6. Here is a quick overview of all of the folders that appear
  - curves: Contains a library of functions to create and evaluate the various Bezier curves used for the models
  - experiments: Contains a series of subfolders that contain experimental data (either raw or digitized from the paper) needed both to simulate experiments and generate comparison plots
  - models: Contains the code needed to evaluate the state derivative of the various models that are simulated.
  - output:	All of the plots, structs, and tables created during the process of simulation are written to this folder	
  - parameters: Contains the scripts needed to set all of the parameters needed to simulate the cat soleus used in this work.	
  - postprocessing:	Contains the scripts needed to generate the various custom plots and tables that are generated during the process of simulation
  - simulation: Contains the scripts that are needed to run the various simulations that are applied to each model.
  - LICENSE: A folder that contains the licenses that apply to the files in this project. This project's licensing will be compliant with the license auditing tool provided by https://api.reuse.software/

#Readme Details

- author: Matthew Millard
- date: 10 April 2023
- version: 0.0
