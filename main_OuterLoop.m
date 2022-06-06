clc;
close all;
clear all;


fitCrossBridgeStiffnessDampingToKirch199490Hz = 1;
flag_useFixedLambdaECM      = 0; %Deprecated: flag_useFixedLambdaECM
flag_makeAndSavePubPlots    = 1;
flag_fitActiveTitinProperties=0; 
%0: Use previously computed values
%1: Solve for the titin-actin attachment point and active damping
%coefficent (this is time consuming ~ 30 min)

flag_fitToFig3KirchBoskovRymer1994 = 0;
%0: will fit to data from Fig. 12
%1: will fit to data from Fig. 3

main_createDefaultFelineSoleusModel;

clc;
close all;
clear all;

main_KirschBoskovRymer1994_OuterLoop;
 
clc;
close all;
clear all;

main_HerzogLeonard2002_OuterLoop;

clc;
close all;
clear all;

main_LeonardJoumaaHerzog2010_OuterLoop;

clc;
close all;
clear all;

main_ClassicForceLengthVelocityExperiments_OuterLoop;