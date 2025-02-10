%%
% SPDX-FileCopyrightText: 2023 Matthew Millard <millard.matthew@gmail.com>
%
% SPDX-License-Identifier: MIT
%
%%

rootDir         = getRootProjectDirectory();
projectFolders  = getProjectFolders(rootDir);


%%
% Rat soleus fibril Model
%%

fprintf('\n\nCreating: default rat soleus fibril model\n');
fprintf('  used to simulate Tomalka, Weider, Hahn, Seiberl, Siebert 2020.\n\n');
%From Fig. 2 Leonard, Joumaa, Herzog, 2010




normCrossBridgeStiffness    = fittedFelineSoleusKBR1994Fig12_RT.sarcomere.normCrossBridgeStiffness;
normCrossBridgeDamping      = fittedFelineSoleusKBR1994Fig12_RT.sarcomere.normCrossBridgeDamping;

titinMolecularWeightInkDDefault =[];


%For now this is being set to a length at which the passive force of titin
%reaches 0.5 fiso. To compensate, the ecm fraction is set to 0.5 thereby
%achieving the desired effect: titin reaches a force of 0.5 iso at
% normFiberLengthAtOneNormPassiveForce.
%
%This work around has to be done because, by default, I assume that the
% toe force of the passive curve has a value of 1. Since rabbit titin is so 
% stiff and short and its contour length is shorter than the length at 
% which the fibrils reach 1 fiso ... and this causes the functions
% charged with making the titin curves to fail.
passiveForceKeyPoints = [1.0,    0;...
                         2.86,1.31];  

normFiberLengthAtOneNormPassiveForceRabbitFibril = ...
    interp1(passiveForceKeyPoints(:,2),passiveForceKeyPoints(:,1),0.5);

% Now we adjust how much of the passive force is comprised of the ECM in
% order to make the titin force-length curve more compliant.
ecmForceFractionRabbitPsoasFitted = 0.675;% 

normPevkToActinAttachmentPointRabbitPsoasFitted=0.675;


rabbitPsoasFibrilWLC = createRabbitPsoasFibrilModel(...
                              normCrossBridgeStiffness,...
                              normCrossBridgeDamping,...
                              normPevkToActinAttachmentPointRabbitPsoasFitted,...
                              normMaxActiveTitinToActinDamping,...
                              normFiberLengthAtOneNormPassiveForceRabbitFibril,...
                              ecmForceFractionRabbitPsoasFitted,...
                              titinMolecularWeightInkDDefault,...
                              wlcTitinModel,...
                              useCalibratedCurves,...
                              useTwoSidedTitinCurves,...
                              smallNumericallyNonZeroNumber,...
                              flag_enableNumericallyNonZeroGradients,...
                              scaleOptimalFiberLengthRabbitPsoas,...
                              scaleMaximumIsometricTensionRabbitPsoas, ...
                              projectFolders,...
                              flag_useOctave);

filePathRabbitPsoas = fullfile(projectFolders.output_structs_FittedModels,...
                                'rabbitPsoasFibrilWLC.mat');

save(filePathRabbitPsoas,'rabbitPsoasFibrilWLC');