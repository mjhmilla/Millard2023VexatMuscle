clc;
close all;
clear all;


fitCrossBridgeStiffnessDampingToKirch199490Hz = 1;
flag_useFixedLambdaECM = 0;
flag_makeAndSavePubPlots = 1;
disp('Deprecated: flag_useFixedLambdaECM');

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