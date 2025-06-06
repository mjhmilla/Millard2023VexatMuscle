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

function [defaultFelineSoleusModel,...
          activeForceLengthCurveAnnotationPoints,...
          felineSoleusActiveForceLengthDataDefault,...
          felineSoleusPassiveForceLengthDataDefault,...
          felineSoleusPassiveForceLengthCurveSettings]...
           = createFelineSoleusModel(normPevkToActinAttachmentPoint,...
                                    normMaxActiveTitinToActinDamping,...
                                    normFiberLengthAtOneNormPassiveForce,... 
                                    ecmForceFraction,...
                                    useWLCTitinModel,...
                                    useCalibratedCurves,...
                                    useTwoSidedTitinCurves,...
                                    smallNumericallyNonZeroNumber,...
                                    flag_enableNumericallyNonZeroGradients,... 
                                    scaleOptimalFiberLength,... 
                                    scaleMaximumIsometricTension,...
                                    projectFolders,...
                                    flag_useOctave)

rigidTendonReferenceModel   = [];
elasticTendonReferenceModel = [];

%  flag_makeAndSavePubPlots                      = 1;



% 1. Active force length curve vs. data
% Solution: There were some initial descrepencies between the experimental force
%length data and a theoretical curve. These errors almost completely go
%away if it is assumed that the experimental recordings are of total 
%path length, rather than fiber length. In this case, when the elastiticy
%of the tendon is taken into account the theoretical active-force-length 
%curve and the transformed data nicely align.

%Failed attempt:
%This creates a cat soleus with an optimal fiber length of 58 mm: this
%is simply way too big to be realistic (given the data I'm working from)
flag_solveForOptimalFiberLengthOfBestFit  = 0; 

%Failed attempt:
shiftLengthActiveForceLengthCurveDescendingCurve = 0.;%...
%  (1/3)*( (1.154-1.087) + (1.23-1.162) + (1.077-1.039) );


%%
% Add the directories needed to run this script
%%

%if(exist('fitCrossBridgeStiffnessDampingToKirch199490Hz','var')==0)
%  fitCrossBridgeStiffnessDampingToKirch199490Hz = 1;
%end

[felineSoleusMusculotendonProperties, ...
 felineSoleusSarcomereProperties,...
 felineSoleusActiveForceLengthData,...
 felineSoleusPassiveForceLengthData] = ...
 createFelineSoleusParameters(    scaleOptimalFiberLength,...
                        scaleMaximumIsometricTension,...
                        normFiberLengthAtOneNormPassiveForce,...
                        normPevkToActinAttachmentPoint,...
                        normMaxActiveTitinToActinDamping,...
                        ecmForceFraction,...
                        projectFolders,...
                        flag_useOctave);

createMusculoTendonFcn = ...
  @(argScaleFiberLength,argScaleFiso)createFelineSoleusParameters(...
                                        argScaleFiberLength,...
                                        argScaleFiso,...
                                        normFiberLengthAtOneNormPassiveForce,...
                                        normPevkToActinAttachmentPoint,...
                                        normMaxActiveTitinToActinDamping,...
                                        ecmForceFraction,...
                                        flag_useOctave); 
                                        
[felineSoleusNormMuscleCurvesDefault,...
 felineSoleusMusculotendonPropertiesDefault,...
 felineSoleusSarcomerePropertiesDefault,... 
 activeForceLengthCurveAnnotationPoints,...
 felineSoleusActiveForceLengthDataDefault,...
 felineSoleusPassiveForceLengthDataDefault,...
 felineSoleusPassiveForceLengthCurveSettings]= ...
    createFittedMuscleCurves( ...
      felineSoleusMusculotendonProperties,...
      felineSoleusSarcomereProperties,...
      useWLCTitinModel,...
      useCalibratedCurves,...
      useTwoSidedTitinCurves,...
      felineSoleusActiveForceLengthData,...
      felineSoleusPassiveForceLengthData,...
      shiftLengthActiveForceLengthCurveDescendingCurve,...
      flag_enableNumericallyNonZeroGradients,...
      smallNumericallyNonZeroNumber,...
      flag_solveForOptimalFiberLengthOfBestFit,...
      createMusculoTendonFcn,...
      flag_useOctave);


% Some muscles appear to have a minimum length for developing linear
% eccentric force profiles, others not. Here we set the default to 
% be the start of the passive-force-length curve and adjust as needed
%
%Tomalka A. Eccentric muscle contractions: from single muscle fibre to 
%whole muscle mechanics. Pflügers Archiv-European Journal of Physiology. 
%2023 Apr;475(4):421-35.

felineSoleusSarcomerePropertiesDefault.normLengthTitinActinBondMinimum = ...
  felineSoleusNormMuscleCurvesDefault.fiberForceLengthCurve.xEnd(1,1);


%%
%Check to make sure that
% The normFiberLengthAtOneNormPassiveForce used to create keypoints for the
% titin model actually matches the norm fiber length at which the passive
% curve develops 1 normalized force.
%%
fpe2 = felineSoleusNormMuscleCurvesDefault.fiberForceLengthCurve.yEnd(1,2);
lpe2 = felineSoleusNormMuscleCurvesDefault.fiberForceLengthCurve.xEnd(1,2);
assert(abs(fpe2-1)<1e-3);
errLp2 = 1e-3;
if(isempty(felineSoleusPassiveForceLengthData)==0)
    errLp2 = 0.05;
end
assert(abs(lpe2-felineSoleusSarcomerePropertiesDefault.normFiberLengthAtOneNormPassiveForce)<errLp2);



defaultFelineSoleusModel = struct('musculotendon',...
                            felineSoleusMusculotendonPropertiesDefault,...
                            'sarcomere',...
                            felineSoleusSarcomerePropertiesDefault,...
                            'falData',...
                            felineSoleusActiveForceLengthDataDefault,...
                            'fpeData',...
                            felineSoleusPassiveForceLengthDataDefault,...
                            'curves',...
                            felineSoleusNormMuscleCurvesDefault,...
                            'fitting',...
                            []);
                      



