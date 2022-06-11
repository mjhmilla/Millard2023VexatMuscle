flag_fitActiveTitinProperties              =0; 


disp(['Note: flag_fitActiveTitinProperties must be set to 1 the first']);
disp([' time you run this script.']);
%0: Use previously computed values
%1: Solve for the titin-actin attachment point and active damping
%coefficent (this is time consuming ~ 30 min)

fitCrossBridgeStiffnessDampingToKirch199490Hz = 1;
flag_useFixedLambdaECM      = 0; %Deprecated: leave as 0


rigidTendonReferenceModel = [];%...
%    'output/structs/felineSoleusRigidTendonKBR1994Fig12.mat';
elasticTendonReferenceModel=[];%...
%    'output/structs/felineSoleusElasticTendonKBR1994Fig12.mat';

flag_fitToFig3KirchBoskovRymer1994 = 0;
flag_makeAndSavePubPlots = 1;
main_createDefaultFelineSoleusModel;

% Generate the model with Kx and Bx fitted to Figure 12 of 
% Kirsch Boskov and Rymer. Models fitted to both Figure 12 and Figure 3 of 
% Kirsch, Boskov, and Rymer are need for later simulations in 
% main_KirschBoskovRymer1994_OuterLoop
%
%0. Will generate models with Kx and Bx fitted to Figure 12 of Kirch,
%Boskov, and Rymer
%1. Will generate models with Kx and Bx fitted to Figure 3 of Kirsch,
%Boskov, and Rymer

close all;

flag_fitToFig3KirchBoskovRymer1994 = 1;
% Generate the model with Kx and Bx fitted to Figure 3 of 
% Kirsch Boskov and Rymer.

flag_fitActiveTitinProperties=0; %Takes previously generated versions 
flag_makeAndSavePubPlots    = 0; 
main_createDefaultFelineSoleusModel;