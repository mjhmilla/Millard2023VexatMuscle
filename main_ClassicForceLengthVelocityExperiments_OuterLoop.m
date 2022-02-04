clc;
close all;
clear all;

%Parameters that do not change
flag_simulateHillModel            = 0;  
flag_simulateOpus31Model          = 1;
flag_simulateRigidTendon          = 1;
flag_simulateElasticTendon        = 0;



flag_plotDataComparison           = 1;  


flag_useFiberDamping              = 1;
fiberDampingCoefficient           = 0.01;

flag_useFig3KirchBoskovRymer1994  = 0; 
flag_useOctave                    = 0;


flag_activeForceLengthSimulations  = 1;
flag_passiveForceLengthSimulations = 1;  
flag_forceVelocitySimulations      = 1;

normFiberLengthAtForceVelocitySample  = 1.;

maxShorteningVelocity = 4.0;
forceVelocityNormFiberHalfLength = 0.1;

numberOfLengthSteps   = 3;
numberOfVelocitySteps = 3;

%Common parameters



flag_plotData           = 0;
flag_buildCombinedPlot  = 0;  
figClassicFl = figure;
figClassicFv = figure; 

if(flag_plotDataComparison==1)
  flag_plotData         = 1;
  flag_savePlotsToFile  = 1;
end

if(flag_simulateRigidTendon==1)
  if(flag_plotDataComparison==1)
      flag_buildCombinedPlot  = 1;  
  end

  flag_useElasticTendon   = 0;  
  flag_useCrossbridgeStiffnessCalibratedCurves = 0; 
  main_ClassicForceLengthVelocityExperiments;


  if(flag_simulateOpus31Model==1)
    flag_useCrossbridgeStiffnessCalibratedCurves = 1;  
    flag_simulateHillModelStashed = flag_simulateHillModel;
    
    flag_simulateHillModel=0;   
    main_ClassicForceLengthVelocityExperiments;
    flag_simulateHillModel=flag_simulateHillModelStashed
  end
end
if(flag_simulateElasticTendon==1)
  if(flag_plotDataComparison==1)
      flag_buildCombinedPlot  = 2; 
  end   
  flag_useElasticTendon   = 1;   

  flag_useCrossbridgeStiffnessCalibratedCurves = 0; 
  main_ClassicForceLengthVelocityExperiments;

  if(flag_simulateOpus31Model==1)
    flag_useCrossbridgeStiffnessCalibratedCurves = 1;  
    flag_simulateHillModelStashed = flag_simulateHillModel;
    
    flag_simulateHillModel=0;   
    main_ClassicForceLengthVelocityExperiments;
    flag_simulateHillModel=flag_simulateHillModelStashed
  end
end
  
