clc;
close all;
clear all;

main_createFelineSoleusModels_OuterLoop.m;

main_createRabbitPsoasModels_OuterLoop.m;

main_createHumanSoleusModels_OuterLoop.m;


%%
% Generate publication quality plots
%%
if(flag_makeAndSavePubPlots==1)
  plotMuscleCurves( felineSoleusNormMuscleCurvesUpd_ET,...
                      activeForceLengthCurveAnnotationPoints,...
                      musculotendonPropertiesOpus31_ET,...
                      sarcomerePropertiesOpus31_ET,...
                      felineSoleusActiveForceLengthDataDefault,...
                      felineSoleusPassiveForceLengthDataDefault,...
                      pubOutputFolder);
end