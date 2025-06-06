%%
% SPDX-FileCopyrightText: 2023 Matthew Millard <millard.matthew@gmail.com>
%
% SPDX-License-Identifier: MIT
%
% If you use this code in your work please cite the pre-print of this paper
% or the most recent peer-reviewed version of this paper:
%
%    Matthew Millard, David W. Franklin, Walter Herzog. 
%    A three filament mechanistic model of musculotendon force and impedance. 
%    bioRxiv 2023.03.27.534347; doi: https://doi.org/10.1101/2023.03.27.534347 
%
%%

clc;
close all;
clear all;

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

disp(['Note: flag_fitActiveTitinProperties must be set to 1 the first']);
disp([' time you run this script.']);
flag_fitActiveTitinProperties              = 0; 
normPevkToActinAttachmentPoint             = 0.5;
%0: Use previously computed values
%1: Solve for the titin-actin attachment point and active damping
%   coefficent (this is time consuming ~ 30 min)

normFiberLengthAtOneNormPassiveForce=1.367732948060934e+00;

useCalibratedCurves = 1;

flag_fitToFig3KirchBoskovRymer1994            = 0;
fitCrossBridgeStiffnessDampingToKirch199490Hz = 1;

rigidTendonReferenceModel = [];%...
%    'output/structs/felineSoleusRigidTendonKBR1994Fig12.mat';
elasticTendonReferenceModel=[];%...
%    'output/structs/felineSoleusElasticTendonKBR1994Fig12.mat';
flag_makeAndSavePubPlots = 1;

%%
% Generate rigid-tendon and elastic-tendon cat soleus model with
%   1. The titin-actin attachment point fitted to match Herzog & Leonard
%      as closely as possible.
%   2. The lumped cross-bridge stiffness and damping tuned to best match
%      the frequency response of Figure 3 of Kirsch, Boskov, & Rymer 
%      as closely as possible
%%
%main_createDefaultFelineSoleusModel;

[defaultFelineSoleusModel,...
felineSoleusPassiveForceLengthCurveSettings]...
   = createFelineSoleusModel(normPevkToActinAttachmentPoint,...
                            normFiberLengthAtOneNormPassiveForce,... 
                            useCalibratedCurves,...
                            smallNumericallyNonZeroNumber,...
                            flag_enableNumericallyNonZeroGradients,...
                            rigidTendonReferenceModel,...
                            elasticTendonReferenceModel,...
                            flag_plotAllCurves,...
                            flag_useOctave);


save('/output/structs/defaultFelineSoleus.mat',...
     'defaultFelineSoleus');  

close all;

flag_fitToFig3KirchBoskovRymer1994              = 1;
fitCrossBridgeStiffnessDampingToKirch199490Hz   = 0;

flag_fitActiveTitinProperties  = 0; %Takes previously generated versions 
flag_makeAndSavePubPlots       = 0; 

%%
% Same as before but now:
%   1. ... same as before
%   2. ... best match the frequency response of Figure 12 of Kirsch, 
%      Boskov, & Rymer as closely as possible
%
% Generate rigid-tendon and elastic-tendon cat soleus model with
%   1. The titin-actin attachment point fitted to match Herzog & Leonard
%      as closely as possible.
%   2. The lumped cross-bridge stiffness and damping tuned to best match
%      the frequency response of Figure 12 of Kirsch, Boskov, & Rymer 
%      as closely as possible.
%%
main_createDefaultFelineSoleusModel;