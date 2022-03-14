clc;
close all;
clear all;


%Generic config.
flag_simulateActiveStretch  = 1;
flag_simulatePassiveStretch = 1;

flag_useCalibratedOpus31Curves= 1;

nominalNormalizedFiberLength = 1.0;

normActiveTitinToActinDamping   = [];% take the default. Prev. 6; 
normPassiveTitinToActinDamping  = []; %[]: take the default. Prev. 1.5; 

flag_useOctave                    = 0;
tunedOpus31Results = 'experiments/LeonardJoumaaHerzog2010/benchRecordOpus31_RigidTendon_K44p31D0p50Tau_LJH2010__TiAD1000p00_TiPD1p50_NomLen1p00_90Hz_TiAdj.mat';


flag_useFig3KirchBoskovRymer1994              = 0; 
flag_useElasticTendon                         = 0; 
flag_useFiberDamping                          = 1;
fiberDampingCoefficient = 0.01;

%Simulation specific config: default models

flag_simulateHillModel      = 1; 
flag_simulateOpus31Model    = 1;
flag_plotData         = 0;  
flag_savePlotsToFile  = 0;
flag_fitTitin                   = 0;  

main_LeonardJoumaaHerzog2010;

%Simulation specific config: adj. models
flag_simulateHillModel      = 0; 
flag_simulateOpus31Model    = 1;
flag_fitTitin                   = 1;
normActiveTitinToActinDamping   = 1000;
main_LeonardJoumaaHerzog2010;

%Plot.
flag_simulateHillModel      = 0; 
flag_simulateOpus31Model    = 0;
flag_plotData         = 1;  
flag_savePlotsToFile  = 1;
flag_fitTitin                   = 0;  
normActiveTitinToActinDamping = [];

figLeonardJoumaaHerzog2010Fig2Comparision = figure;
main_LeonardJoumaaHerzog2010;
