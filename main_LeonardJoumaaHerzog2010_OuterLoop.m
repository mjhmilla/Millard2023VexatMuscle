clc;
close all;
clear all;


disp('Set these three flags to 1 if running from scratch.');
flag_defaultSimulation       =1;
flag_tunedSimulation         =1;
flag_defaultHillSimulation   =1;


%Generic config.
flag_simulateActiveStretch      = 1;
flag_simulatePassiveStretch     = 1;

flag_useCalibratedOpus31Curves  = 1;
flag_useTwoSidedTitinCurves     = 0;

nominalNormalizedFiberLength = 1.0;

flag_useOctave                    = 0;
%tunedOpus31Results = 'benchRecordOpus31_RigidTendon_K39p47D0p35Tau_LJH2010__TiAD1000p00_TiPD0p25_NomLen1p00_90Hz_TiAdj';

flag_useFig3KirchBoskovRymer1994  = 0; 
flag_useElasticTendon             = 0; 
flag_useFiberDamping              = 1;
fiberDampingCoefficient           = 0.01;

%Simulation specific config: adj. models
if(flag_tunedSimulation==1)
    flag_simulateHillModel          = 0; 
    flag_simulateOpus31Model        = 1;
    flag_plotData                   = 0;  
    flag_savePlotsToFile            = 0;    
    flag_useTunedRabbitPsoasModel   = 1;
    main_LeonardJoumaaHerzog2010;
end

%Simulation specific config: default models
if(flag_defaultSimulation==1)
    flag_simulateHillModel      = 0; 
    flag_simulateOpus31Model    = 1;
    flag_plotData               = 0;  
    flag_savePlotsToFile        = 0;
    flag_useTunedRabbitPsoasModel = 0;
    

    main_LeonardJoumaaHerzog2010;
end

if(flag_defaultHillSimulation==1)
    flag_simulateHillModel      = 1; 
    flag_simulateOpus31Model    = 0;
    flag_plotData               = 0;  
    flag_savePlotsToFile        = 0;
    flag_useTunedRabbitPsoasModel = 0;
    

    main_LeonardJoumaaHerzog2010;
end





%Plot.
flag_simulateHillModel          = 0; 
flag_simulateOpus31Model        = 0;
flag_plotData                   = 1;  
flag_savePlotsToFile            = 1;
flag_useTunedRabbitPsoasModel   = 0;

figLeonardJoumaaHerzog2010Fig2Comparision = figure;
main_LeonardJoumaaHerzog2010;
