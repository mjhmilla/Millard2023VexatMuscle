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
                                    useElasticTendon,...
                                    elasticTendonReferenceModel,...
                                    projectFolders,...
                                    flag_useOctave)



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
                                        projectFolders,...
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
      useElasticTendon,...
      elasticTendonReferenceModel,...
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


%%
% Populate the output structure
%%
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


%%
% Check to make sure the the passive curves of the rigid and elastic 
% tendon models match
%%
flag_examineQualityOfRigidToElasticForceLengthCurveMatch=1;
if(flag_examineQualityOfRigidToElasticForceLengthCurveMatch == 1 && ...
        useElasticTendon == 0 && ...
        isempty(elasticTendonReferenceModel)==0)

    lopt    = elasticTendonReferenceModel.musculotendon.optimalFiberLength;
    penOpt  = elasticTendonReferenceModel.musculotendon.pennationAngle;
    ltSlk   = elasticTendonReferenceModel.musculotendon.tendonSlackLength;
    
    
    lceN0 = elasticTendonReferenceModel.curves.fiberForceLengthCurve.xEnd(1,1);
    lceN1 = elasticTendonReferenceModel.curves.fiberForceLengthCurve.xEnd(1,2);

    lceNV = zeros(100,1);
    lceNVRT=zeros(100,1);

    lpV = zeros(100,1);
    fpeNATV=zeros(100,1);
    fpeNRTATV=zeros(100,1);
    for(i=1:1:100)
       n = (i-1)/99;
       nMin = 0.05;    
       lceN_ET =lceN0 + (n*(1-nMin) + nMin)*(lceN1-lceN0);
       fpeN_ET = calcBezierYFcnXDerivative(lceN_ET, ...
           elasticTendonReferenceModel.curves.fiberForceLengthCurve,0);
    
       fiberKin = calcFixedWidthPennatedFiberKinematicsAlongTendon(lceN_ET*lopt,0,lopt,penOpt);
       alpha = fiberKin.pennationAngle;

       fpeNATV(i,1)=fpeN_ET*cos(alpha);

    
       ftN_ET = fpeN_ET*cos(alpha);
    
       ltN_ET = calcBezierYFcnXDerivative(ftN_ET, ...
                elasticTendonReferenceModel.curves.tendonForceLengthInverseCurve, 0);
       lp = (lopt * lceN_ET)*cos(alpha) + ltN_ET*ltSlk;
       lpV(i,1)=lp;
    
       lceAT_RT = lp - ltSlk;
       
       fiberKinRT = calcFixedWidthPennatedFiberKinematics(lceAT_RT,0,lopt,penOpt);
       lce_RT = fiberKinRT.fiberLength;
       alpha_RT = fiberKinRT.pennationAngle;
    
       lceN_RT = lce_RT/lopt;
    
       fpeN_RT = calcBezierYFcnXDerivative(lceN_RT, ...
                    defaultFelineSoleusModel.curves.fiberForceLengthCurve,0);
       fpeNRTATV(i,1)=fpeN_RT*cos(alpha_RT);
       
    end        
    figFpeETvsRT=figure;
    subplot(2,2,1);
        plot(lpV,fpeNATV,'b');
        hold on;
        plot(lpV,fpeNRTATV,'r');
        box off;
        xlabel('Path Length (m)');
        ylabel('Norm. Force');
        title('Passive force length: ET vs RT');
    subplot(2,2,2);
        plot(lpV,fpeNRTATV-fpeNATV,'m');
        box off;
        xlabel('Path Length (m)');
        ylabel('Norm. Force');
        title('Passive force error: ET vs RT');
    here=1;


    lTitinFixedHN = elasticTendonReferenceModel.sarcomere.ZLineToT12NormLengthAtOptimalFiberLength ...
                  + elasticTendonReferenceModel.sarcomere.IGDFixedNormLengthAtOptimalFiberLength;    
    ecmFrac = elasticTendonReferenceModel.sarcomere.extraCellularMatrixPassiveForceFraction;

    lTitinFixedHN_RT = defaultFelineSoleusModel.sarcomere.ZLineToT12NormLengthAtOptimalFiberLength ...
                     + defaultFelineSoleusModel.sarcomere.IGDFixedNormLengthAtOptimalFiberLength;    

    ecmFrac_RT = defaultFelineSoleusModel.sarcomere.extraCellularMatrixPassiveForceFraction;

    fecm12NTiV=zeros(100,1);
    fecm12NTiV_RT=zeros(100,1);
    lpTiV = zeros(100,1);

    lceN0 = elasticTendonReferenceModel.curves.fiberForceLengthCurve.xEnd(1,1);
    lceN1 = elasticTendonReferenceModel.curves.fiberForceLengthCurve.xEnd(1,2);    

    for(i=1:1:100)
       n = (i-1)/99;
       nMin = 0.1;    
       lceN_ET =lceN0 + (n*(1-nMin) + nMin)*(lceN1-lceN0);
        
        %Evaluate the ECM and titin, assuming the prox and distal 
        %segments of titin are in a force equilibrium
        fecmN = calcBezierYFcnXDerivative(lceN_ET*0.5, ...
                elasticTendonReferenceModel.curves.forceLengthECMHalfCurve, 0);
        l12N = lceN_ET*0.5 - lTitinFixedHN;
        if(i==99)
            here=1;
        end
        [x1, x2, y12] = calcSeriesSpringStretch(l12N, ...
            elasticTendonReferenceModel.curves.forceLengthProximalTitinCurve,...
            elasticTendonReferenceModel.curves.forceLengthProximalTitinInverseCurve,...
            elasticTendonReferenceModel.curves.forceLengthDistalTitinCurve,...
            elasticTendonReferenceModel.curves.forceLengthDistalTitinInverseCurve);
        l1N=x1;
        l2N=x2;
        f1N=y12;
        f2N=y12;        
        %Set the ce strain of the elastic tendon model
        fiberKin = calcFixedWidthPennatedFiberKinematicsAlongTendon(lceN_ET*lopt,0,lopt,penOpt);
        alpha = fiberKin.pennationAngle;    

        if(i==90)
            here=1;
        end

        fecm12NTiV(i,1)=(fecmN+f2N)*cos(alpha);
        %Evaluate the tendon strain
        ltN  = calcBezierYFcnXDerivative(fecm12NTiV(i,1), ...
                elasticTendonReferenceModel.curves.tendonForceLengthInverseCurve, 0);
        lpTiV(i,1) = (lceN_ET*lopt)*cos(alpha) + ltN*ltSlk;

        %Evaluate the rigid-tendon model fiber kinematics
        lceAT_RT = lpTiV(i,1) - ltSlk;
       
        fiberKinRT = calcFixedWidthPennatedFiberKinematics(lceAT_RT,0,lopt,penOpt);
        lce_RT = fiberKinRT.fiberLength;
        alpha_RT= fiberKinRT.pennationAngle;
    
        lceN_RT = lce_RT/lopt;
        fecmN_RT = calcBezierYFcnXDerivative(lceN_RT*0.5, ...
                defaultFelineSoleusModel.curves.forceLengthECMHalfCurve, 0);
        l12N = lceN_RT*0.5 - lTitinFixedHN_RT;

        [x1, x2, y12] = calcSeriesSpringStretch(l12N, ...
            defaultFelineSoleusModel.curves.forceLengthProximalTitinCurve,...
            defaultFelineSoleusModel.curves.forceLengthProximalTitinInverseCurve,...
            defaultFelineSoleusModel.curves.forceLengthDistalTitinCurve,...
            defaultFelineSoleusModel.curves.forceLengthDistalTitinInverseCurve);
        l1N_RT=x1;
        l2N_RT=x2;

        f1N_RT=y12;
        f2N_RT=y12;

        fecm12NTiV_RT(i,1)= (fecmN_RT + f2N_RT)*cos(alpha_RT);
    
    end
    subplot(2,2,3);
        plot(lpTiV,fecm12NTiV,'b');
        hold on;
        plot(lpTiV,fecm12NTiV_RT,'r');
        box off;
        xlabel('Path Length (m)');
        ylabel('Norm. Force');
        title('ECM + Titin force-length: ET vs RT');
    subplot(2,2,4);
        plot(lpTiV,fecm12NTiV-fecm12NTiV_RT,'m');
        box off;
        xlabel('Path Length (m)');
        ylabel('Norm. Force');
        title('ECM + Titin force error: ET vs RT');
    here=1;


end


                      




