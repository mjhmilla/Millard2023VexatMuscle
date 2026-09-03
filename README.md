# Description

This repository contains the Matlab prototype for the VEXAT muscle model described in Millard et al. and published in eLife and the bioRxiv. The contents of both versions are identical though the eLife version includes the eLife assessment text and links to the public reviews of the paper. The model is named after the viscoelastic (VE) crossbridge (X) active-titin (AT) elements contained in the model.

Millard Matthew, Franklin David W., Herzog Walter (2023) A three filament mechanistic model of musculotendon force and impedance eLife 12:RP88344 https://doi.org/10.7554/eLife.88344.3

Matthew Millard, David W. Franklin, Walter Herzog. A three filament mechanistic model of musculotendon force and impedance. bioRxiv 2023.03.27.534347; doi: https://doi.org/10.1101/2023.03.27.534347 


# Quick Start: Analysis and simulation of active muscle fibers during lengthening across the force-length relation

This part of the README.md refers specifically to code that accompanies a paper titled <em>Analysis and simulation of active muscle fibers during lengthening across the force-length relation</em> in which an analysis and simulation of Tomalka et al. is presented.

Tomalka A, Rode C, Schumacher J, Siebert T. The active force–length relationship is invisible during extensive eccentric contractions in skinned skeletal muscle fibres. Proceedings of the Royal Society B: Biological Sciences. 2017 May 17;284(1854):20162497. https://doi.org/10.1098/rspb.2016.2497

## Running the simulations

To run the simulations that are presented in <em>Analysis and simulation of active muscle fibers during lengthening across the force-length relation</em> please:

1. Run main_TomalkaRodeSchumacherSiebert2017_OuterLoop.m in Matlab. Please note that the simulations presented in <em>Analysis and simulation ...</em> were run on the following systems
  - Ubuntu 20.04.6 LTS + Matlab version 9.11.0.1873467 (R2021b) Update 3 
  - Windows 11 Pro + Matlab version 9.13.0.2105380 (R2022b) Update 2
2. The function main_TomalkaRodeSchumacherSiebert2017_OuterLoop.m requires about 10 minutes to terminate on a 2.80 GHz Intel machine with 32 GB ram, and a SSD.
3. These simulations generate the following plots found in the paper:
  - Figures 1,6,7A,7B from the paper in <em>output/plots/TomalkaRodeSchumacherSiebert2017/fig_Sim_TRSS2017_123_Fl_Fv_QToF.pdf</em>
  - Figures 7C, and 7D from the paper in <em>output/plots/TomalkaRodeSchumacherSiebert2017/fig_Sim_TRSS2017_123_Fl_Fv_QToF_i.pdf</em>
  - Figures 4 and 9 from the paper in <em>output/plots/MuscleCurves/fig_Pub_RatMuscleCurves_TRSS2017_0_Pub.pdf</em>
  - Figure 2 in <em>output/plots/MuscleCurves/fig_Pub_RatMuscleCurves_TRSS2017_0.pdf</em>
4. Initial model parameters (before fitting) can be found: 
  - <em>output/structs/FittedModels/ratTRSS2017EDLFibrilActiveTitin_0.mat</em>. 
  - Note that <em>ratTRSS2017EDLFibrilActiveTitin_1.mat</em> and <em>ratTRSS2017EDLFibrilActiveTitin_2.mat</em> are identical to <em>ratTRSS2017EDLFibrilActiveTitin_0.mat</em>
5. Outputs related to fitting with one Q value: 
  - Fitted model parameters: 
    - <em>output/structs/FittedModels/ratTRSS2017EDLFibrilActiveTitinFitted_123_Fl_Fv_QToF.mat</em> 
  - Time-series data of the fitted simulations: 
    - <em>output/structs/TomalkaRodeSchumacherSiebert2017/benchRecordVexat_TRSS2017_fitted_123_Fl_Fv_QToF.mat</em>
  - Fitted parameters and error information
    - <em>output/tables/TomalkaRodeSchumacherSiebert2017/fittingInfo_ratTRSS2017EDLFibrilActiveTitinFitted_123_Fl_Fv_QToF.tex</em>   
    - <em>output/structs/TomalkaRodeSchumacherSiebert2017/fittingInfo_ratTRSS2017EDLFibrilActiveTitinFitted_123_Fl_Fv_QToF.mat</em>  
    - <em>output/structs/fittingLog__123_Fl_Fv_QToF.txt</em>
6. Outputs related to fitting with a Q value for each trial (Q_1, Q_2, and Q_3). Note that the file names have <em>_i</em> appended
  - Fitted model parameters: 
    - <em>output/structs/FittedModels/ratTRSS2017EDLFibrilActiveTitinFitted_123_Fl_Fv_QToF_i.mat</em> 
  - Time-series data of the fitted simulations: 
    - <em>output/structs/TomalkaRodeSchumacherSiebert2017/benchRecordVexat_TRSS2017_fitted_123_Fl_Fv_QToF_i.mat</em>
  - Fitted parameters and error information
    - <em>output/tables/TomalkaRodeSchumacherSiebert2017/fittingInfo_ratTRSS2017EDLFibrilActiveTitinFitted_123_Fl_Fv_QToF_i.tex</em>         
    - <em>output/structs/TomalkaRodeSchumacherSiebert2017/fittingInfo_ratTRSS2017EDLFibrilActiveTitinFitted_123_Fl_Fv_QToF_i.mat</em>  
    - <em>output/structs/fittingLog__123_Fl_Fv_QToF_i.txt</em>    
7. To generate Figure 8:
  - Open main_TomalkaRodeSchumacherSiebert2017_OuterLoop.m
  - Set all occurances of simConfigInput.manuallySetTitinParameters to 1. Its a little annoying to have 3 copies of this varible, but this is necessary because everything is cleared between calls to main_TomalkaRodeSchumacherSiebert2017.m 
      - line 19: simConfigInput.manuallySetTitinParameters = 1;
      - line 27: simConfigInput.manuallySetTitinParameters = 1;
      - line 43: simConfigInput.manuallySetTitinParameters = 1;
  - Run main_TomalkaRodeSchumacherSiebert2017_OuterLoop.m.
  - Note 1: This will over-write all plots, structs, and generated tex files using the manually set titin parameters. Note that these parameters have been identified using optimisation. Due to the structure of the software written for this project it is easiest to update the model and simulations by manually copying these parameters over and re-running the data.
  - Note 2: To anyone working on the software, both main_CreateRatMuscleModel.m and main_TomalkaRodeSchumacherSiebert2017.m have structs that will manually set the force-length relation of titin that must be manually set to be identical.
  - Note 3: If you would like to solve for the passive force-length parameters of titin that were used in Figure 8,
      - Open main_TomalkaRodeSchumacherSiebert2018.m
      - line 7 update: flag_OuterLoopMode=0;      
      - line 89 update: fittingConfig.fitQToF=0;
      - line 94 update: fittingConfig.fitFpeQToF=1; 
      - Re-run main_TomalkaRodeSchumacherSiebert2018.m
      - As before, this may overwrite plots/tables/structs produced by this script.
      
## Code overview

Here is a list of files to look at if the implementation interests you:

- Model generation
  - <em>main_CreateRatMuscleModel.m</em>
  - Generate the reference human soleus titin model model. Used in the methods described in Appendix B 
    - <em>parameters/createHumanSoleusModel2025.m</em> 
  - Generate the rat EDL model 
    - <em>parameters/createRatSkeletalMuscleModel.m</em> 
  - Get Sarcomere related parameters
    - <em>parameters/getMammalianSkeletalMuscleNormalizedSarcomereProperties2025.m</em>
  - Titin curve scaling from Appendix B.
    - <em>parameters/scaleTitinElongationFunction2025.m</em>      
  - Get Musculotendon parameters
    - <em>parameters/getRatMusculotendonProperties.m</em>
  - Bezier curve generation
    - <em>curves/createFittedMuscleCurvesTRSS2017.m</em>
  - Active force length relation curves
    - <em>curves/MusculoTendonCurves/createFiberActiveForceLengthCurve2025.m</em>
  - Titin curves
    - <em>curves/MusculoTendonCurves/createTitinCurves2025.m</em>

- Fitting to Tomalka et al. 2017
  - <em>main_TomalkaRodeSchumacherSiebert2017.m</em>
  - Fit one set of parameters
    - <em>parameters/fitRatFibrilTRSS2017.m</em>
  - Evaluate the error of the force-length and force-velocity curves
    - <em>parameters/calcErrorTRSS2017ForceLengthRelationAscendingLimb.m</em>
  - Evalaute the error of the Q parameters (among others)
    - <em>parameters/calcErrorTRSS2017RampFraction.m</em>


# Quick start guide: Millard, Franklin, Herzog (2024)

The instructions here are related to the simulations presented in Millard et al. (2024)

Millard M, Franklin DW, Herzog W. A three filament mechanistic model of musculotendon force and impedance. Elife. 2024 Sep 10;12:RP88344.(https://doi.org/10.7554/eLife.88344.4)

Execute 'main_OuterLoop.m' from Matlab to run everything. Roughly 5 hours and 45 minutes is needed to run all of the experiments on my old 2013 Samsung laptop (Intel i7-3630QM @ 2.40 GHz, Ubuntu 22 8 GB ram, SSD harddrive), while the newer Lenovo machine (Windows 10) requires nearly an hour less time. The experiments require roughly 500 simulations, some of which are numerically stiff. At the end if this round of simulation, you can find all of the figures that appear in Millard et al., and many more besides, in the folders:

output/
- HerzogLeonard2002
- InitializationBenchmark
- KirschBoskovRymer1994
- LeonardJoumaaHerzog2010
- MuscleCurves
- NettiDamoreRoncaAmbrosioNicolais1996
- StandardTests

Please be aware that there might be slight differences between the results contained in Millard et al. and the simulation. These small differences can arise for two reasons:

1. The properties of titin are fitted by optimization, the solution of which may differ from one run to the next.
2. The simulations of Kirsch et al. require a pseudo-random perturbation waveform. This waveform is generated from scratch each time the script is run.

## Repository Overview

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


