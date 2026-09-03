clc; 
close all;
clear all;

%%
%
% Before running this script, make sure that
%
% 1. Open main_CreateRatMuscleModel.m and set
%    flag_OuterLoopMode=1;
%    on line  8
%
% 2. Open main_TomalkaRodeSchumacherSiebert2017.m and set
%    flag_OuterLoopMode=1;
%    on line  7
%%

main_CreateRatMuscleModel;

simConfigInput.runFitting              = 1; 
simConfigInput.generatePlots           = 1;
simConfigInput.fitToIndividualTrials   = 1; %Individually fitted Q's
simConfigInput.manuallySetTimeConstant = 0;

main_TomalkaRodeSchumacherSiebert2017;

pause(1.0)

clc;
close all;
clear all;

simConfigInput.runFitting              = 1; 
simConfigInput.generatePlots           = 1;
simConfigInput.fitToIndividualTrials   = 0;  %One Q for all trials
simConfigInput.manuallySetTimeConstant = 0;

main_TomalkaRodeSchumacherSiebert2017;

