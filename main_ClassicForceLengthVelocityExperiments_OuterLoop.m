clc;
close all;
clear all;

%Parameters that do not change
flag_simulateHillModel            = 1; 
flag_simulateOpus31Model          = 1;


flag_plotDataComparison           = 1;  


flag_useFiberDamping              = 1;
fiberDampingCoefficient           = 0.01;

flag_useFig3KirchBoskovRymer1994  = 0; 
flag_useOctave                    = 0;

flag_activeForceLengthSimulations  = 1;
flag_passiveForceLengthSimulations = 1;  
flag_forceVelocitySimulations      = 1;

normFiberLengthAtForceVelocitySample  = 1.;

maxShorteningVelocity = 3.0;
forceVelocityNormFiberHalfLength = 0.05;

numberOfLengthSteps   = 20;
numberOfVelocitySteps = 10;

%Common parameters


if(   flag_simulateHillModel   == 1 ...
   || flag_simulateOpus31Model == 1)

  flag_plotData = 0;
  flag_buildCombinedPlot = 0;  
  flag_useElasticTendon             = 1; 
  main_ClassicForceLengthVelocityExperiments;

  flag_useElasticTendon             = 0; 
  main_ClassicForceLengthVelocityExperiments;  
  
end

if(flag_plotDataComparison==1)
  flag_plotData         = 1;
  flag_savePlotsToFile  = 1;
  
  figClassicFl = figure;
  figClassicFv = figure;  
  
  flag_buildCombinedPlot  = 1;  
  flag_useElasticTendon   = 0;   
  main_ClassicForceLengthVelocityExperiments;

  flag_buildCombinedPlot  = 2;  
  flag_useElasticTendon   = 1;   
  main_ClassicForceLengthVelocityExperiments;
  
  
end

% 1 Feb 2021
% The model gives descent results, with a slower time constant (which
% improves the results of Kirch et al., and speeds up the simulations)
% if the sarcomere properties are massaged a bit:
%
% sarcomerePropertiesOpus31.fvNVelocityScaling = 1;
% sarcomerePropertiesOpus31.slidingTimeConstantSlow = 0.002;
% sarcomerePropertiesOpus31.slidingTimeConstantFast = 0.002;
% 
% sarcomerePropertiesOpus31.normCrossBridgeCyclingDampingSlow = 10;%120/10;
% sarcomerePropertiesOpus31.normCrossBridgeCyclingDampingFast = 1;%20/10;
