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
disp(' running main_HerzogLeonard2002_OuterLoop');
disp('----------------------------------------');
disp('   :run-time: ~52 minutes ');
disp('            *Intel i7-3630QM @ 2.40 GHz, Ubuntu 22');
disp '             8 GB ram, SSD harddrive');
disp('----------------------------------------');

%%
% Generic setup
%%
rootDir         = getRootProjectDirectory();
projectFolders  = getProjectFolders(rootDir);


flag_useFig3KirchBoskovRymer1994              = 0; 
flag_useFiberDamping                          = 1;

figHerzogLeonard2002Fig7Comparision = figure;

figDescendingCombined = figure;

flag_testVexatDerivativeFunction = 0;
flag_useOctave                    = 0;


%%
% Runs
%%
flag_simulateActiveStretch  = 1;
flag_simulatePassiveStretch = 1;
flag_simulateStatic         = 1;

flag_simulateHillModel       = 1; 
flag_simulateVexatModel      = 1;
flag_useCalibratedVexatCurves= 1;
flag_useTitinCurvesWithRigidIgDSegment   = 0;
flag_useTwoSidedTitinCurves   = 0;
flag_plotDataComparison       = 1;
flag_useHorizontalPlotLayout  = 0;

simParamsHL2002File = fullfile( projectFolders.experiments_HL2002,...
                                'simulationParametersHerzogLeonard2002.csv');
simParamsHL2002 = csvread(simParamsHL2002File);
nominalNormalizedFiberLength = simParamsHL2002(1,1);


scaleMaximumIsometricForce     = [];
scaleOptimalFiberLength        = [];%[]: take the default. Prev. 1.0;
scaleSlidingTimeConstant       = [];%[]: take the default. Prev. 0.2;
normActiveTitinToActinDamping  = [];%[]: take the default. Now 100
normPassiveTitinToActinDamping = [];%[]: take the default. Prev. 1.5; 


if( flag_simulateHillModel == 1 ... 
    || flag_simulateVexatModel == 1)
  
  for indexTendonType=1:1:2

    flag_useElasticTendon       = 0;  
    if(indexTendonType == 1)
      flag_useElasticTendon     = 1;  
    end
    for indexSubFigureNumber=1:1:3

      flag_savePlotsToFile      = 0;
      
      for indexTrialNumber=1:1:3
        
        flag_savePlotsToFile      = 0;
        flag_plotData             = 0;

        figureNumber       = 7;
        subFigureNumber    = indexSubFigureNumber;
        trialNumber        = indexTrialNumber;  

        flag_buildCombinedPlot = 1;
        main_HerzogLeonard2002;

      end
   end
  end
  
end

if(flag_plotDataComparison==1)

  flag_simulateHillModel      = 0; 
  flag_simulateVexatModel    = 0;

  flag_plotData=1;
  for indexSubFigureNumber=1:1:3
    clf(figHerzogLeonard2002Fig7Comparision);   
    clf(figDescendingCombined)
    for indexTendonType=1:1:2

      flag_useElasticTendon       = 0;  
      if(indexTendonType == 2)
        flag_useElasticTendon     = 1;  
      end



      indexTrialNumber=3;

      figureNumber       = 7;
      subFigureNumber    = indexSubFigureNumber;
      trialNumber        = indexTrialNumber;  

      flag_buildCombinedPlot = 1;      
      flag_savePlotsToFile   = 0;
      
      if(indexTendonType==2)
        flag_buildCombinedPlot = 2;      
        flag_savePlotsToFile   = 1;        
      end

      main_HerzogLeonard2002;
      here=1;
    end

  end
end

