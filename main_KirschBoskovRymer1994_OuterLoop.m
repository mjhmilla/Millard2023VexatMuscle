clc;
close all;
clear all;

%%
%Run all of the simulations
%%

flag_runSimulations     = 1;
flag_frequencyAnalysis  = 1;
flag_generatePlots      = 1;




if(flag_runSimulations == 1)
  
  flag_simulateHillModel      = 1; 
  flag_simulateOpus31Model    = 1;
  flag_fitToFig3KirchBoskovRymer1994            = 0;
  flag_useElasticTendon       = 1;
  flag_useFiberDamping       = 1;
  flag_frequencyAnalysisMuscleModels            = 0;
  flag_plotAccelerationEquationFactors          = 0;  
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 0;
  flag_pubTabulateStiffnessDampingVariation     = 0;

  main_KirschBoskovRymer1994;

  close all;

  flag_simulateHillModel      = 1; 
  flag_simulateOpus31Model    = 1;
  flag_fitToFig3KirchBoskovRymer1994            = 0;
  flag_useElasticTendon       = 0;
  flag_useFiberDamping       = 1;
  flag_frequencyAnalysisMuscleModels            = 0;
  flag_plotAccelerationEquationFactors          = 0;  
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 0;
  flag_pubTabulateStiffnessDampingVariation     = 0;

  main_KirschBoskovRymer1994;

  close all;  

  flag_simulateHillModel      = 1; 
  flag_simulateOpus31Model    = 1;
  flag_fitToFig3KirchBoskovRymer1994            = 1;
  flag_useElasticTendon       = 1;
  flag_useFiberDamping       = 1;
  flag_frequencyAnalysisMuscleModels            = 0;
  flag_plotAccelerationEquationFactors          = 0;  
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 0;
  flag_pubTabulateStiffnessDampingVariation     = 0;

  main_KirschBoskovRymer1994;

  close all;  
  
  flag_simulateHillModel      = 1; 
  flag_simulateOpus31Model    = 1;
  flag_fitToFig3KirchBoskovRymer1994            = 1;
  flag_useElasticTendon       = 0;
  flag_useFiberDamping       = 1;
  flag_frequencyAnalysisMuscleModels            = 0;
  flag_plotAccelerationEquationFactors          = 0;  
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 0;
  flag_pubTabulateStiffnessDampingVariation     = 0;

  main_KirschBoskovRymer1994;

  close all;





end

%%
%Frequency Analysis
%%

if(flag_frequencyAnalysis==1)
  flag_simulateHillModel      = 0; 
  flag_simulateOpus31Model    = 0;
  flag_fitToFig3KirchBoskovRymer1994            = 0; 
  flag_useElasticTendon       = 0;
  flag_useFiberDamping       = 1;
  flag_frequencyAnalysisMuscleModels            = 1;
  flag_plotAccelerationEquationFactors          = 0;  
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 0;
  flag_pubTabulateStiffnessDampingVariation     = 0;

  main_KirschBoskovRymer1994;

  close all;

  flag_simulateHillModel      = 0; 
  flag_simulateOpus31Model    = 0;
  flag_fitToFig3KirchBoskovRymer1994            = 0; 
  flag_useElasticTendon       = 1;
  flag_useFiberDamping       = 1;
  flag_frequencyAnalysisMuscleModels            = 1;
  flag_plotAccelerationEquationFactors          = 0;  
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 0;
  flag_pubTabulateStiffnessDampingVariation     = 0;

  main_KirschBoskovRymer1994;

  close all;
  
  flag_simulateHillModel      = 0; 
  flag_simulateOpus31Model    = 0;
  flag_fitToFig3KirchBoskovRymer1994            = 1; 
  flag_useElasticTendon       = 0;
  flag_useFiberDamping       = 1;
  flag_frequencyAnalysisMuscleModels            = 1;
  flag_plotAccelerationEquationFactors          = 0;  
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 0;
  flag_pubTabulateStiffnessDampingVariation     = 0;

  main_KirschBoskovRymer1994;

  close all;

  flag_simulateHillModel      = 0; 
  flag_simulateOpus31Model    = 0;
  flag_fitToFig3KirchBoskovRymer1994            = 1; 
  flag_useElasticTendon       = 1;
  flag_useFiberDamping       = 1;
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
  
  flag_simulateHillModel      = 0; 
  flag_simulateOpus31Model    = 0;
  flag_fitToFig3KirchBoskovRymer1994               = 0; 
  flag_useElasticTendon       = 1;
  flag_useFiberDamping       = 1;
  flag_frequencyAnalysisMuscleModels            = 0;
  flag_plotAccelerationEquationFactors          = 1;  
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 1;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 1;
  flag_pubTabulateStiffnessDampingVariation     = 1;

  main_KirschBoskovRymer1994;
  
  flag_simulateHillModel      = 0; 
  flag_simulateOpus31Model    = 0;
  flag_fitToFig3KirchBoskovRymer1994               = 0; 
  flag_useElasticTendon       = 0;
  flag_useFiberDamping       = 1;
  flag_frequencyAnalysisMuscleModels            = 0;
  flag_plotAccelerationEquationFactors          = 1;  
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 1;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 1;
  flag_pubTabulateStiffnessDampingVariation     = 1;

  main_KirschBoskovRymer1994;

  close all;



  
  
  close all;
  flag_simulateHillModel      = 0; 
  flag_simulateOpus31Model    = 0;
  flag_fitToFig3KirchBoskovRymer1994               = 0; 
  flag_useElasticTendon       = 0;
  flag_useFiberDamping       = 1;
  flag_frequencyAnalysisMuscleModels            = 0;
  flag_plotAccelerationEquationFactors          = 1;
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 1;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 0;
  flag_pubTabulateStiffnessDampingVariation     = 0;

  main_KirschBoskovRymer1994;

  close all;

  flag_simulateHillModel      = 0; 
  flag_simulateOpus31Model    = 0;
  flag_fitToFig3KirchBoskovRymer1994               = 0; 
  flag_useElasticTendon       = 1;
  flag_useFiberDamping       = 1;
  flag_frequencyAnalysisMuscleModels            = 0;
  flag_plotAccelerationEquationFactors          = 1;  
  flag_pubPlotFrequencyResponseKBR1994Fig3      = 1;
  flag_pubPlotStiffnessDampingKBR1994Fig9Fig10  = 0;
  flag_pubPlotStiffnessDampingKBR1994Fig12      = 0;
  flag_pubTabulateStiffnessDampingVariation     = 0;

  main_KirschBoskovRymer1994;

  close all;



  close all;
end

