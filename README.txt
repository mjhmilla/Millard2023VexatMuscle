@author Matthew Millard
@date 6/9/2021

This repository contains the Matlab prototype for the muscle model described in the paper:

Milard M, Franklin D, and Herzog W. A three filament mechanistic model of musculotendon force and impedance (2021-2022). 


Quick start guide:

1. Execute 'main_OuterLoop.m' from Matlab to run everything. Several hours are needed to run all of the experiments due to the number and numerical stiffness of the simulations performed.

2. The muscle model used for all of the experiments that of a cat soleus. The script that creates all of the parameters needed for the model is 'main_createDefaultFelineSoleusModel.m', while all of the sub-scripts used to create the cat soleus parameters appears in the 'parameters' folder.

3. Each experiment is accompanied by a script (e.g. main_HerzogLeonard2002.m) and an script-harness to run the script with multiple configurations (e.g. main_HerzogLeonard2002_OuterLoop.m). Please refer to these files to see exactly how these experiments are simulated.

4. The model is implemented in 2 files: models/calcMillard2019MuscleInfoOpus31.m and models/calcEquilibriumErrorOpus31.m. The first file (models/calcMillard2019MuscleInfoOpus31.m) acts as a wrapper file which takes in the inputs from the user, solves for the solution of the models state derivative, and populates the output structures. The state derivative requires numerical iteration, and so all of the models code that might require iteration is implemented in the second file (models/calcEquilibriumErrorOpus31.m).

5. The Hill model that is used for comparison purposes appears in model/calcMillard2012DampedEquilibriumMuscleInfo.

6. Here is a quick overview of all of the folders that appear

curves: 
  Contains a library of functions to create the various Bezier curves used for the models

experiments:
  Contains a series of subfolders that contain experimental data (either raw or digitzed from the paper) needed both to simulate experiments and generate comparison plots

models:	
  Contains the code needed to evaluate the state dervative of the various models that are simulated.

output:	
  All of the plots, structs, and tables created during the process of simulation are written to this folder	

parameters:
  Contains the scripts needed to set all of the parameters needed to simulate the cat soleus used in this work.	

postprocessing:	
  Contains the scripts needed to generate the various custom plots and tables that are generated during the process of simulation

simulation:
  Contains the scripts that are needed to run the various simulations that are applied to each model.

