%%
% SPDX-FileCopyrightText: 2023 Matthew Millard <millard.matthew@gmail.com>
%
% SPDX-License-Identifier: MIT
%
% If you use this code in your work please cite the pre-print of this paper
% or the most recent peer-reviewed version of this paper:
%
%    Matthew Millard, David W. Franklin, Walter Herzog. 
%    A three filament mechanistic model of musculotendon force and impedance. 
%    bioRxiv 2023.03.27.534347; doi: https://doi.org/10.1101/2023.03.27.534347 
%
%%

%This flag allows us to avoid the memory clearing functions so that
%this can be timed using tic and tock from within main_OuterLoop
flag_OuterOuterLoopMode =1;
if(flag_OuterOuterLoopMode ==0)
    clc;
    close all;
    clear all;
end

disp('----------------------------------------');
disp(' running main_KirschBoskovRymer1994_OuterLoop');
disp('----------------------------------------');
disp('   :run-time: 3h 45min (226 min) ');
disp('            :Why so long? Nearly 400 simulations are performed');
disp('            *Intel i7-3630QM @ 2.40 GHz, Ubuntu 22');
disp '             8 GB ram, SSD harddrive');
disp('----------------------------------------');

rootDir         = getRootProjectDirectory();
projectFolders  = getProjectFolders(rootDir);

%%
%Run all of the simulations
%%

flag_runSimulations     = 1;
flag_frequencyAnalysis  = 1;
flag_generatePlots      = 1;
flag_generateTables     = 1;

flag_useSameRandomPerturbationAsPublication   = 1;

%Parameters that are tuned
%updSlidingTimeConstant              = 0.0005;
%updForceVelocityCalibrationFactor    = 0.95;

flag_useCalibratedVexatCurves = 1;
flag_useTwoSidedTitinCurves   = 0;

if(flag_runSimulations == 1)
  %Hill and the proposed model
  %Elastic tendon
  flag_simulateHillModel                        = 1; 
  flag_simulateVexatModel                      = 1;
  flag_fitToFig3KirchBoskovRymer1994            = 0;
  flag_useElasticTendon                         = 1;
  flag_useFiberDamping                          = 1;
  flag_frequencyAnalysisMuscleModels            = 0;
  flag_plotAccelerationEquationFactors          = 0;  
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 0;
  flag_pubTabulateStiffnessDampingVariation     = 0;

  flag_generateRandomInput    = 1;
  flag_processInputFunctions  = 1;

  main_KirschBoskovRymer1994;

  close all;

  %Hill and the proposed model
  %Rigid tendon
  flag_simulateHillModel                        = 1; 
  flag_simulateVexatModel                      = 1;
  flag_fitToFig3KirchBoskovRymer1994            = 0;
  flag_useElasticTendon                         = 0;
  flag_useFiberDamping                          = 1;
  flag_frequencyAnalysisMuscleModels            = 0;
  flag_plotAccelerationEquationFactors          = 0;  
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 0;
  flag_pubTabulateStiffnessDampingVariation     = 0;

  flag_generateRandomInput    = 0;
  flag_processInputFunctions  = 0;

  main_KirschBoskovRymer1994;

  close all;  

  %Proposed model only
  %Elastic tendon
  flag_simulateHillModel                        = 1; 
  flag_simulateVexatModel                      = 1;
  flag_fitToFig3KirchBoskovRymer1994            = 1;
  flag_useElasticTendon                         = 1;
  flag_useFiberDamping                          = 1;
  flag_frequencyAnalysisMuscleModels            = 0;
  flag_plotAccelerationEquationFactors          = 0;  
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 0;
  flag_pubTabulateStiffnessDampingVariation     = 0;

  flag_generateRandomInput    = 0;
  flag_processInputFunctions  = 0;

  main_KirschBoskovRymer1994;

  close all;  

  %Proposed model only
  %Rigid tendon  
  flag_simulateHillModel                        = 1; 
  flag_simulateVexatModel                      = 1;
  flag_fitToFig3KirchBoskovRymer1994            = 1;
  flag_useElasticTendon                         = 0;
  flag_useFiberDamping                          = 1;
  flag_frequencyAnalysisMuscleModels            = 0;
  flag_plotAccelerationEquationFactors          = 0;  
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 0;
  flag_pubTabulateStiffnessDampingVariation     = 0;

  flag_generateRandomInput    = 0;
  flag_processInputFunctions  = 0;

  main_KirschBoskovRymer1994;

  close all;





end

%%
%Frequency Analysis
%%

if(flag_frequencyAnalysis==1)
  flag_generateRandomInput                      = 0;
  flag_processInputFunctions                    = 0;
  %Rigid tendon - do not fit Kx and betaX to Kirsch et al. Fig 3
  flag_simulateHillModel                        = 0; 
  flag_simulateVexatModel                      = 0;
  flag_fitToFig3KirchBoskovRymer1994            = 0; 
  flag_useElasticTendon                         = 0;
  flag_useFiberDamping                          = 1;
  flag_frequencyAnalysisMuscleModels            = 1;
  flag_plotAccelerationEquationFactors          = 0;  
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 0;
  flag_pubTabulateStiffnessDampingVariation     = 0;

  main_KirschBoskovRymer1994;

  close all;

  %Elastic tendon - do not fit Kx and betaX to Kirsch et al. Fig 3
  flag_simulateHillModel                        = 0; 
  flag_simulateVexatModel                      = 0;
  flag_fitToFig3KirchBoskovRymer1994            = 0; 
  flag_useElasticTendon                         = 1;
  flag_useFiberDamping                          = 1;
  flag_frequencyAnalysisMuscleModels            = 1;
  flag_plotAccelerationEquationFactors          = 0;  
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 0;
  flag_pubTabulateStiffnessDampingVariation     = 0;

  main_KirschBoskovRymer1994;

  close all;
  
  %Rigid tendon - fit Kx and betaX to Kirsch et al. Fig 3  
  flag_simulateHillModel                        = 0; 
  flag_simulateVexatModel                      = 0;
  flag_fitToFig3KirchBoskovRymer1994            = 1; 
  flag_useElasticTendon                         = 0;
  flag_useFiberDamping                          = 1;
  flag_frequencyAnalysisMuscleModels            = 1;
  flag_plotAccelerationEquationFactors          = 0;  
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 0;
  flag_pubTabulateStiffnessDampingVariation     = 0;

  main_KirschBoskovRymer1994;

  close all;

%  Elastic tendon - fit Kx and betaX to Kirsch et al. Fig 3    
  flag_simulateHillModel                        = 0; 
  flag_simulateVexatModel                      = 0;
  flag_fitToFig3KirchBoskovRymer1994            = 1; 
  flag_useElasticTendon                         = 1;
  flag_useFiberDamping                          = 1;
  flag_frequencyAnalysisMuscleModels            = 1;
  flag_plotAccelerationEquationFactors          = 0;  
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 0;
  flag_pubTabulateStiffnessDampingVariation     = 0;

  main_KirschBoskovRymer1994;

  close all;  
end

%%
%Plots
%%

if(flag_generatePlots==1)
  flag_generateRandomInput    = 0;
  flag_processInputFunctions  = 0;



  %Elastic tendon
  flag_simulateHillModel                        = 0; 
  flag_simulateVexatModel                       = 0;
  flag_fitToFig3KirchBoskovRymer1994            = 0; 
  flag_useElasticTendon                         = 1;
  flag_useFiberDamping                          = 1;
  flag_frequencyAnalysisMuscleModels            = 0;
  flag_plotAccelerationEquationFactors          = 1;  
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 1;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 1;
  flag_pubTabulateStiffnessDampingVariation     = 0;

  main_KirschBoskovRymer1994;

  %Rigid tendon  
  flag_simulateHillModel                        = 0; 
  flag_simulateVexatModel                       = 0;
  flag_fitToFig3KirchBoskovRymer1994            = 0; 
  flag_useElasticTendon                         = 0;
  flag_useFiberDamping                          = 1;
  flag_frequencyAnalysisMuscleModels            = 0;
  flag_plotAccelerationEquationFactors          = 1;  
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 1;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 1;
  flag_pubTabulateStiffnessDampingVariation     = 0;

  main_KirschBoskovRymer1994;

  close all;

  %When flag_pubPlotFrequencyResponseKBR1994Fig3 is 1 both 
  %elastic tendon and rigid tendon plots are made simultaneously
  flag_simulateHillModel                        = 0; 
  flag_simulateVexatModel                       = 0;
  flag_fitToFig3KirchBoskovRymer1994            = 0; 
  flag_useElasticTendon                         = 0;
  flag_useFiberDamping                          = 1;
  flag_frequencyAnalysisMuscleModels            = 0;
  flag_plotAccelerationEquationFactors          = 1;
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 1;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 0;
  flag_pubTabulateStiffnessDampingVariation     = 0;

  main_KirschBoskovRymer1994;

  close all;


end

if(flag_generateTables)

  flag_generateRandomInput    = 0;
  flag_processInputFunctions  = 0;

  %Elastic tendon
  flag_simulateHillModel                        = 0; 
  flag_simulateVexatModel                       = 0;
  flag_fitToFig3KirchBoskovRymer1994            = 0; 
  flag_useElasticTendon                         = 1;
  flag_useFiberDamping                          = 1;
  flag_frequencyAnalysisMuscleModels            = 0;
  flag_plotAccelerationEquationFactors          = 0;  
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 0;
  flag_pubTabulateStiffnessDampingVariation     = 1;

  main_KirschBoskovRymer1994;

  %Rigid tendon  
  flag_simulateHillModel                        = 0; 
  flag_simulateVexatModel                       = 0;
  flag_fitToFig3KirchBoskovRymer1994            = 0; 
  flag_useElasticTendon                         = 0;
  flag_useFiberDamping                          = 1;
  flag_frequencyAnalysisMuscleModels            = 0;
  flag_plotAccelerationEquationFactors          = 0;  
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 0;
  flag_pubTabulateStiffnessDampingVariation     = 1;

  main_KirschBoskovRymer1994;

  close all;
  

end
