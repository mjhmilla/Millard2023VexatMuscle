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
flag_OuterOuterLoopMode =0;
if(flag_OuterOuterLoopMode ==0)
    clc;
    close all;
    clear all;
end


disp('----------------------------------------');
disp(' running main_LeonardJoumaaHerzog2010_OuterLoop');
disp('----------------------------------------');
disp('   :run-time: 6 minutes ');
disp('            *Intel i7-3630QM @ 2.40 GHz, Ubuntu 22');
disp '             8 GB ram, SSD harddrive');
disp('----------------------------------------');


rootDir         = getRootProjectDirectory();
projectFolders  = getProjectFolders(rootDir);

disp('Set these three flags to 1 if running from scratch.');
flag_defaultSimulation       =1;
flag_tunedSimulation         =1;
flag_defaultHillSimulation   =1;


%Generic config.
flag_simulateActiveStretch      = 1;
flag_simulatePassiveStretch     = 1;

flag_useCalibratedVexatCurves  = 1;
flag_useTwoSidedTitinCurves     = 0;

nominalNormalizedFiberLength = 1.0;

flag_useOctave                    = 0;
%tunedVexatResults = 'benchRecordVexat_RigidTendon_K39p47D0p35Tau_LJH2010__TiAD1000p00_TiPD0p25_NomLen1p00_90Hz_TiAdj';

flag_useFig3KirchBoskovRymer1994  = 0; 
flag_useElasticTendon             = 0; 
flag_useFiberDamping              = 1;
fiberDampingCoefficient           = 0.01;

%Simulation specific config: adj. models
if(flag_tunedSimulation==1)
    flag_simulateHillModel          = 0; 
    flag_simulateVexatModel        = 1;
    flag_plotData                   = 0;  
    flag_savePlotsToFile            = 0;    
    flag_useTunedRabbitPsoasModel   = 1;
    main_LeonardJoumaaHerzog2010;
end

%Simulation specific config: default models
if(flag_defaultSimulation==1)
    flag_simulateHillModel      = 0; 
    flag_simulateVexatModel    = 1;
    flag_plotData               = 0;  
    flag_savePlotsToFile        = 0;
    flag_useTunedRabbitPsoasModel = 0;
    

    main_LeonardJoumaaHerzog2010;
end

if(flag_defaultHillSimulation==1)
    flag_simulateHillModel      = 1; 
    flag_simulateVexatModel    = 0;
    flag_plotData               = 0;  
    flag_savePlotsToFile        = 0;
    flag_useTunedRabbitPsoasModel = 0;
    

    main_LeonardJoumaaHerzog2010;
end





%Plot.
flag_simulateHillModel          = 0; 
flag_simulateVexatModel        = 0;
flag_plotData                   = 1;  
flag_savePlotsToFile            = 1;
flag_useTunedRabbitPsoasModel   = 0;

figLeonardJoumaaHerzog2010Fig2Comparision = figure;
main_LeonardJoumaaHerzog2010;
