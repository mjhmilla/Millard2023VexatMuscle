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

function [muscleArchitecture,...
          crossBridgeParameters,...
          sarcomereParameters,...
          muscleCurveData,...
          excitationFcn,...
          activationFcn, ...
          pathLengthFcn, ...
          expForceData,...
          expLengthData,...
          initialState] = config_HerzogLeonard2002(...
                              figureNumber, subFigureNumber, trialNumber,...
                              maximumNormalizedFiberVelocity,...
                              forceVelocityMultiplierAtHalfOmegaMax,...
                              flagElasticTendon, ...
                              tendonStrainAtOneNormForce,...
                              scaleOptimalFiberLength,...
                              scaleSrsStiffness,...
                              scaleSrsMaxForce, ...
                              scaleTitinStiffness,...
                              scaleTitinDamping,...
                              fisoScale, ...
                              scaleFastForceVelocityLimit, ...
                              scaleSlowForceVelocityLimit,...
                              scaleRatioSlowForceVelocityToShortRangeStiffness)

assert(figureNumber == 7, ['Only Fig. 7 is allowed: Fig. 5 & 6 are ', ...
                           'from different cats']);                          
                          
muscleCurveData = struct('fpeData',[],'falData',[]);
expForceData = [];



%%
%
% Muscle Architecture & CrossBridge Properties: Extract from Fig. 7A.
%
%%
fitData          = csvread('dataHerzogLeonard2002Figure6and7A.csv',1,0);
idxFitF43=2;   % Length 0mm, Activation 0-1-0 
idxFitL43=3;
idxFitF44=4;   % Length 9mm, Activation 0-1-0
idxFitL44=5;
idxFitF45=6;   % Length 6-9mm, Activation 0-1-0
idxFitL45=7;
idxFitF47=8;   % Length 3-9mm, Activation 0-1-0
idxFitL47=9;
idxFitF49=10;  % Length 0-9mm, Activation 0 (Passive)
idxFitL49=11;
idxFitF50=12;  % Length 0-9mm, Activation 0-1-0
idxFitL50=13;

idxFitMax = 1500; %All columns except 50 series stop at 1500

%extracted from Herzog & Leonard 2002
lceOptExp = (27/0.63)/1000; % From the text: pg 1277 paragraph 2;
lceOpt   = lceOptExp*scaleOptimalFiberLength;

alphaOpt = 7.5*(pi/180);   % from Scott et al. 1996  
ltSlk    = 65/1000;        % from Scott et al. 1996



[ftMax idxFtMax] = max(fitData(:,idxFitF43));
lp      = fitData(idxFtMax, idxFitL43);
fp      = 0;
idxBest = 0;
fpBest  = 0;
lerr    = Inf;
for k=1:1:size(fitData,1)
    if( abs(lp - fitData(k,idxFitL50)) < lerr)
        idxBest = k;
        lerr    = abs(lp - fitData(k,idxFitL50));
        fpBest  = fitData(k,idxFitF50);
    end
end
fp = fpBest;


fisoExp  = (ftMax-fp)/cos(alphaOpt); 
fiso     = fisoScale*fisoExp;
% From Fig.7A from the 9mm stretch just before the stretch begins

etIso = 0.0;
normTendonLength =1.0;
if(flagElasticTendon==1)
  etIso = tendonStrainAtOneNormForce;
  normTendonLength = 1+etIso;
end


%Set later when the curves are created
%        'minimumNormHillDamping', NaN,...
%        'maximumHillDamping'    , NaN, ...

%titinNormStiffness = 1*scaleTitinStiffness; %As in 1 fiso/lopt.
%titinNormDamping   = 0.1*scaleTitinDamping;
%titinCoilingVelocityThreshold = 0.05;

minimumFiberLengthAlongTendon = sqrt(eps);

lcePerp    = lceOpt*sin(alphaOpt);
lceAT      = minimumFiberLengthAlongTendon;

alphaMax        = atan2(lcePerp, ...
                       lceAT);
minimumFiberLength = sqrt(lcePerp*lcePerp + lceAT*lceAT);

dlceMaxN=maximumNormalizedFiberVelocity;

muscleArchitecture = struct(...
        'name'                        , 'catSoleus',  ... 
        'abbr'                        , 'catSol',     ...
        'fiso'                        , fiso,         ...  
        'fisoExp'                     , fisoExp,      ...
        'optimalFiberLength'          , lceOpt,       ... 
        'optimalFiberLengthExp'       , lceOptExp,    ...
        'pennationAngle'              , alphaOpt,     ... 
        'tendonSlackLength'           , ltSlk,        ...  
        'tendonStrainAtOneNormForce'  , etIso,        ...
        'normTendonDampingLinear'     , 0.05,          ...
        'normTendonDampingConstant'   , 0.05,           ...
        'normNumericalDamping'        , 0.001,         ...
        'normNumericalStiffness'      , 0.001,        ...
        'maximumNormalizedFiberVelocity', maximumNormalizedFiberVelocity,...
        'minimumFiberLength'          ,       minimumFiberLength,    ...
        'minimumFiberLengthAlongTendon',      minimumFiberLengthAlongTendon, ...
        'pennationAngleAtMinimumFiberLength', alphaMax,     ...
        'forceVelocityMultiplierAtHalfMaximumFiberVelocity' , ...
            forceVelocityMultiplierAtHalfOmegaMax,...
        'normMass'              , 0.001,...
        'normPassiveForceLengthCurveStiffnessAtSlack', 0.001, ...
        'normFpeDamping' , NaN, ...
        'normECMDamping' , 0.0, ...
        'scalePEVK', 1.0,...
        'scaleECM', 1.0,...
        'scaleIG1',1.0);


%%
%
% X-Bridge Properties
%
%%


%From Fig. 7A. These properties are taken from the graph without regard
%for the changes in pennation angle, and tendon stretch. 
%Just to get in ball park

idx0    = 229; %Index beginning of short-range-stiffness
idx1    = 253; %Index end of short-range, beginning of long-range
idx2    = 517; %Index at the end of the long range.

timeFit = [fitData(idx0,1);...
           fitData(idx1,1);...
           fitData(idx2,1)];

ftATFit = [fitData(idx0, idxFitF50);...
           fitData(idx1, idxFitF50);...
           fitData(idx2, idxFitF50)];
fpeATFit = zeros(3,1);
lpATFitErr=zeros(3,1);
       
fceFit = zeros(3,1);                                  

lpDeltaFit = [fitData(idx0,idxFitL50);...
            fitData(idx1,idxFitL50);...
            fitData(idx2,idxFitL50)];

[fpeMax, idxMaxFpe] = max( fitData(:,idxFitL49));
        
for k=1:1:3
   lerr = Inf;
   fpeBest=0;
   for z=1:1:idxMaxFpe
       err = abs(lpDeltaFit(k)-fitData(z,idxFitL49));
       if(err < lerr)
           fpeATFit(k) = fitData(z,idxFitF49);
           lerr = err;
           lpATFitErr(k) = err;
       end
   end
end
        
lceATFit = lpDeltaFit./1000 + (lceOptExp*cos(alphaOpt)).*ones(3,1);

lceFit = zeros(3,1);
alphaFit =zeros(3,1);

flagPlotCrossBridgeFittingCurves =0;
if(flagPlotCrossBridgeFittingCurves == 1)
   figCrossBridgeFit = figure;
   subplot(2,1,1)    
    plot( fitData(1:idxFitMax,1), fitData(1:idxFitMax,idxFitF50),'r');    
    hold on;
    plot( timeFit, ftATFit,'xb');
    hold on;
    xlabel('Time (s)');
    ylabel('Force (N)');
       
   subplot(2,1,2)
    plot( fitData(1:idxFitMax,1), fitData(1:idxFitMax,idxFitL50),'k');
    hold on;
    plot( timeFit, lpDeltaFit,'xb');
    xlabel('Time (s)');
    ylabel('Length (mm)');    
    hold on;
end

%The change in load is so small that we can ignore the change in
%tendon length
fceAT = zeros(3,1);
for k=1:1:3
    fiberKinematics = ...
     calcFixedWidthPennatedFiberKinematics(lceATFit(k),...
                                          NaN,...
                                          lceOptExp,...
                                          alphaOpt);
    lceFit(k)   = fiberKinematics.fiberLength;
    alphaFit(k) = fiberKinematics.pennationAngle;
    fceAT(k)    = ftATFit(k)-fpeATFit(k);
    fceFit(k)   = fceAT(k)/(cos(alphaFit(k)));
end

lceNormFit = lceFit./lceOptExp;
fceNormFit = fceFit./fisoExp;
                                         
normCrossBridgeStiffnessShortRange = ...
    scaleSrsStiffness*(fceNormFit(2)-fceNormFit(1)) ...
                      /(lceNormFit(2)-lceNormFit(1));

normCrossBridgeStiffnessLongRange = ...
    (fceNormFit(3)-fceNormFit(2)) ...
   /(lceNormFit(3)-lceNormFit(2));

maxNormShortRangeCrossBridgeForce = 1 + scaleSrsMaxForce*(fceNormFit(2)-1);

%maxNormShortRangeCrossBridgeForce = 1 + 0.5*(maxNormShortRangeCrossBridgeForce-1);
%disp('**max Norm Short Range Force region scaled**');

normFastForceVelocityForceLimit = ...
  scaleFastForceVelocityLimit*(maxNormShortRangeCrossBridgeForce-1);

normSlowForceVelocityForceLimit = ...
  scaleSlowForceVelocityLimit*(1.6-1-normFastForceVelocityForceLimit);


maxNormLongRangeCrossBridgeForce = fceNormFit(3);

maxNormCrossBridgeLengthChange = ...
    ((maxNormLongRangeCrossBridgeForce-maxNormShortRangeCrossBridgeForce) ...
    /normCrossBridgeStiffnessLongRange);



maxNormCrossBridgeLengtheningRate = maxNormCrossBridgeLengthChange/3.0;

%normCrossBridgeLengthAtFiso = ...
%  maxNormShortRangeCrossBridgeForce/normCrossBridgeStiffnessShortRange;

normFiberForceDelta   = (fceNormFit(3) -fceNormFit(2));
normFiberLengthDelta  = (lceNormFit(3) -lceNormFit(2));
deltaT                = (timeFit(3) - timeFit(2));
normFiberVelocityDelta= (normFiberLengthDelta/deltaT)...
                         /(maximumNormalizedFiberVelocity);

normCrossBridgeVelocity    = ...
  ((normFiberForceDelta/normCrossBridgeStiffnessShortRange)/deltaT) ...
  /maximumNormalizedFiberVelocity;

%normSlidingElementVelocity = normFiberVelocityDelta - normCrossBridgeVelocity; 


normEccentricForceVelocitySlope = (normFiberForceDelta ...
                                 / normFiberVelocityDelta);

normEccentricForceVelocitySlopeMax = 2*normEccentricForceVelocitySlope;
normEccentricForceVelocitySlopeMin = 0.05*normEccentricForceVelocitySlopeMax;

ratioShortRangeStiffnessToSlowForceVelocity = ...
  scaleRatioSlowForceVelocityToShortRangeStiffness;


normTitinDamping = (maxNormShortRangeCrossBridgeForce-1)...
                  /normFiberVelocityDelta;

normTitinActinCompressiveLengthAtZeroForce = 0.05;  %As in 0.05*lceOpt            
normTitinActinCompressiveLengthAtFiso      = 0;

normOptimalTitinLength = 0.5;
%normTitinLengthAtFiso = normTitinLengthAtFiso/normOptimalTitinLength;             
normTitinStiffnessLongRange =  normCrossBridgeStiffnessLongRange ...
                               /normOptimalTitinLength;
%%
%If you ever pick new names, stick to this convention
%
% physicalQuantityPartName
%
% e.g. normStiffnessShortRangeCrossBridge
%      maxForceShortRangeCrossbridge
%
% Or maybe not. I'm not sure this is any clearer than what I have.
%%

optSarcomereLength = 0.5*(2.64+2.81); %Rassier et al. 1999
normTitinSlackLength = 0.5; %As in 0.5*lceOpt
normTitinLengthAtFiso= 0.5*1.6; 
normHalfMyosinLength = 0.5*(4.24 - optSarcomereLength)/optSarcomereLength; %Rassier et al.
normHalfMyosinBareLength = (0.2/optSarcomereLength)*0.5; %Rassier et al. 
normActinLength      = 0.5*(4.24 - 2*normHalfMyosinLength*optSarcomereLength)...
                      /optSarcomereLength;
normZLineLength      = 0.05/optSarcomereLength; %Rassier et al.                     
%1.27/optSarcomereLength; %Rassier et al.

normLengthSteepToShallowAscendingCurve = 1.7/(optSarcomereLength);
%  'normTitinActinOffsetSteepAscendingCurve'     , normTitinActinOffsetSteepAscendingCurve,...

%Note All IG1, IG2, PEVK, and T12 lengths come from Trombitas et al. 1997
%     by assuming lopt is 2.725 um, and that fiso is hit at 1.6 lopt.

loptT12  = 0.1/optSarcomereLength;
loptIG1  = 0.19/optSarcomereLength; %Should be renamed IGP: proximal Ig
loptPEVK = 0.19/optSarcomereLength;
loptIG2  = 1*0.5-(loptT12+loptIG1+loptPEVK); %0.5: applies to a half sarcomere
loptIGD  = loptIG2 - normHalfMyosinLength; %IGD: the distal Ig section that
% is between the PEVK section and myosin. This is the small section of Ig that
% is still free to stretch. The rest is bound to myosin and does not flex at 
% all.


lfisoIG1  = 0.45812/optSarcomereLength;
lfisoPEVK = 0.6454/optSarcomereLength;
lfisoIG2  = (4.36*0.5/optSarcomereLength-(loptT12 + lfisoIG1+lfisoPEVK));
lfisoIGD  = lfisoIG2 - normHalfMyosinLength;

activationTimeConstant = 0.1;
deactivationTimeConstant= 0.1;

ig1NormDamping = 10*(4/fiso)/(0.63*0.5);
ig2NormDamping = 10*(4/fiso)/(0.63*0.5);
pevkNormDamping =10*(4/fiso)/(0.63*0.5);
fpeNormDamping = 1/( (1/ig1NormDamping)...
                    +(1/ig2NormDamping)...
                    +(1/pevkNormDamping));

%Numbers from Prado                  
extraCellularMatrixPassiveForceFraction = 1-0.5*(0.24+0.57);

%Extracted from Fig 7A first trial of Herzog & Leonard 2002
%

normActiveFiberForce  = [1.222,1.506];
normFiberLength       = [1.022,1.175];

kNormPevkIGd = 2*diff(normActiveFiberForce)/diff(normFiberLength);
%The 2* (roughly) accomodates for the fact that this measurement is made on the 
%descending limb

sarcomereParameters = ...
    struct( 'normTitinSlackLength'              , normTitinSlackLength,...
            'normActinLength'                   , normActinLength,...  
            'normMyosinHalfWidth'               , normHalfMyosinLength,...
            'normMyosinBareHalfLength'          , normHalfMyosinBareLength,...
            'normZLineLength'                   , normZLineLength,...
            'normTitinNormLengthAtFiso'         , normTitinLengthAtFiso,...
            'normMyosinActinCompressiveDamping' , 1,...
            'normTitinDamping'                  , NaN,...
            'normTransitionForce'               , 1.1,...
            'normTransitionVelocity'            , NaN,...
            'MlineToMyosinTitinAttachmentNormLength'  , normHalfMyosinLength*0.05,...
            'ZLineToT12NormLengthAtOptimalFiberLength', loptT12,...
            'IG1NormLengthAtOptimalFiberLength'       , loptIG1,...
            'IG2NormLengthAtOptimalFiberLength'       , loptIG2,...            
            'PEVKNormLengthAtOptimalFiberLength'      , loptPEVK,...
            'IGDNormLengthAtOptimalFiberLength'       , loptIGD, ...
            'IG1NormLengthAtFiso'                     , lfisoIG1,...
            'IG2NormLengthAtFiso'                     , lfisoIG2,...
            'PEVKNormLengthAtFiso'                    , lfisoPEVK,...
            'IGDNormLengthAtFiso'                     , lfisoIGD,...      
            'PEVKIGDNormStiffness'                    , kNormPevkIGd,...      
            'IG1NormDamping'                          , ig1NormDamping , ...
            'IG2NormDamping'                          , ig2NormDamping , ...
            'PEVKNormDamping'                         , pevkNormDamping, ...
            'fpeNormDamping'                          , fpeNormDamping,...
            'activationTimeConstant'                  , 0.1, ...
            'deactivationTimeConstant'                , 0.1,...
            'normMyosinTitinDampingShortening'        , 2,...
            'normMyosinTitinDampingLengthening'       , 20,...
            'activationSwitchScaling'                 , 0.001,...
            'lengtheningSwitchScaling'                , 0.05,...
            'activeTitinFraction'                     , 1.0 , ...
            'extraCellularMatrixPassiveForceFraction' , ...
                extraCellularMatrixPassiveForceFraction, ...
            'normCrossBridgeStiffness', 8*2,...%21.0*2, ...
            'normCrossBridgeDamping', 0.126/(0.5*maximumNormalizedFiberVelocity),...
            'normCrossBridgeCyclingDamping', 75,...
            'normMaxActiveTitinToActinDamping' , 60,... %60.0, ...
            'normPassiveTitinToActinDamping'   , 0.01/(0.5*maximumNormalizedFiberVelocity), ...
            'slidingTimeConstant', 0.005 ,...
            'activeLengtheningStiffening', 0.0, ...
            'titinActiveRelaxationTimeConstant',100, ...
            'minimumCrossBridgeStiffness',0.1);

% 'normCrossBridgeStiffness'  : 21.0   fiso/lopt 
% 'normCrossBridgeDamping'    :  0.126/vmax fiso/(vmax)
%
% These parameters come from Fig. 12 of Kirch, Boskov, & Rymer assuming
% an fiso of 12.5 N (from the same paper, estimate from 1st paragraph of)
% the results, and a 42 mm fiber length (estimated from Herzog & Leonard 1997, 
% soleus, 1 year old cats, min mass 3.5 kg compared to Kirsch et al. cats
% with a mass of 2-5 kg - same average mass)
%
%'normMaxActiveTitinToActinDamping' , 20.0, ...
%'normPassiveTitinToActinDamping'   , 0.10/(0.5*maximumNormalizedFiberVelocity), ...
%  note: maximumNormalizedFiberVelocity is in units of fiber lengths/sec
%
% Scaled for a half sarcomere: 
%   normMaxActiveTitinToActinDamping can sustain 10 fiso at a dlce of 1 lopt/sec
%   normPassiveTitinToActinDamping can sustain 0.1 fiso at dlce of vmax
%  note: maximumNormalizedFiberVelocity is in units of fiber lengths/sec

crossBridgeParameters = struct(...
  'normCrossBridgeStiffnessShortRange'          , normCrossBridgeStiffnessShortRange      ,...
  'normOptimalTitinLength'                      , normOptimalTitinLength,... 
  'normHalfMyosinLength'                        , normHalfMyosinLength,...
  'normActinLength'                             , normActinLength,...
  'normLengthSteepToShallowAscendingCurve'      , normLengthSteepToShallowAscendingCurve,...
  'normTitinLengthAtFiso'                       , normTitinLengthAtFiso,...
  'normTitinStiffnessLongRange'                 , scaleTitinStiffness*normTitinStiffnessLongRange,...  
  'normTitinDamping'                            , scaleTitinDamping*normTitinDamping, ...
  'normTitinActinCompressiveLengthAtZeroForce'  , normTitinActinCompressiveLengthAtZeroForce, ...
  'normTitinActinCompressiveLengthAtFiso'       , normTitinActinCompressiveLengthAtFiso, ...
  'normTitinActinCompressiveCurveNormLowStiffnessWidth', 0.5,...
  'maxNormShortRangeCrossBridgeForce'           , maxNormShortRangeCrossBridgeForce       ,...  
  'normVelocityDuringShortRangeCrossBridgeForce', normFiberVelocityDelta,...      
  'normShortRangeForcePreload'                  , 1,...  
  'crossBridgeNaturalFrequency'                 , 15,... 
  'crossBridgeDamping'                          , 4,...
  'fiberDamping'                                , 0.01,...
  'ratioOmegaSlowToOmegaFast'                   , 0.01,...
  'fastForceVelocityFastTimeConstant'           , 1/15, ...
  'fastForceVelocityRatioConcentricToEccentric' , 2.0, ...
  'fastForceVelocitySlowTimeConstant'           , 0.1, ...
  'fastForceVelocityNormFiberVelocityThreshold' , 0.01, ...
  'normFastForceVelocityForceLimit'             , normFastForceVelocityForceLimit     ,...
  'normSlowForceVelocityForceLimit'             , normSlowForceVelocityForceLimit     ,...  
  'ratioShortRangeStiffnessToSlowForceVelocity' , ratioShortRangeStiffnessToSlowForceVelocity ,...  
  'slowForceVelocityRatioConcentricToEccentric' , 0.1, ...
  'slowForceVelocityTimeConstantSlow'           , 4, ...
  'slowForceVelocityTimeConstantFast'           , 0.040, ...
  'normCrossBridgeStiffnessShortRangeScaling'   , scaleSrsStiffness      ,...    
  'maxNormShortRangeCrossBridgeForceScaling'    , scaleSrsMaxForce);     


%%
%
% Passive & Active Element Fitting Data
%
%%



idx49LStretch = find(fitData(:,idxFitL49) > 1.0);
idx49LStart   = idx49LStretch(1);     

[maxF49 idxMaxF49] = max(fitData(:,idxFitF49));



lmOpt = lceOptExp*cos(alphaOpt) + ltSlk*normTendonLength;
idx   = 1;
idxDelta = floor((idxMaxF49-idx49LStart)/10); %So that we get ~10 pts

muscleCurveData.fpeData = zeros(floor((idxMaxF49-idx49LStart)/idxDelta),2);
muscleCurveData.falData = zeros(4, 2);

for(i=idx49LStart:idxDelta:idxMaxF49)
  muscleCurveData.fpeData(idx,1) = lmOpt + fitData(i,idxFitL49)/1000;
  muscleCurveData.fpeData(idx,2) = fitData(i,idxFitF49);    
  idx = idx+1;
end          


colFalF =[idxFitF43; idxFitF44; idxFitF45; idxFitF47];
colFalL =[idxFitL43; idxFitL44; idxFitL45; idxFitL47];
rowFal  =zeros(length(colFalF),1);
[tmp rowFal(1)] = max(fitData(:,idxFitF43));
[tmp rowFal(2)] = max(fitData(:,idxFitF44));
rowFal(3) = 422;
rowFal(4) = 326;


falLpDelta    = zeros(length(colFalL),1);
falFt         = zeros(length(colFalL),1);
falFpe        = zeros(length(colFalL),1);
falFpErr      = zeros(length(colFalL),1);
for(k=1:1:length(colFalF))
    falLpDelta(k) = fitData(rowFal(k), colFalL(k));
    falFt(k)      = fitData(rowFal(k), colFalF(k));
    
    lerr = Inf;
    for(z=1:1:idxMaxFpe)
        err = abs(falLpDelta(k)-fitData(z,idxFitL49));
        if(err < lerr)
           lerr = err;
           falFpe(k) = fitData(z,idxFitF49);
           falFpErr(k) = err;
        end        
    end
end

flagPlotFalFitCurves =0;
if(flagPlotFalFitCurves == 1)
  figCurveFit = figure;
  color0 = [0;0;0];
  color1 = [1;0;0];   
  subplot(2,1,1)    
  for(i=1:1:length(colFalF))
    n = (i-1)/length(colFalF);
    color01 = color0.*(1-n) + n.*color1;
    plot( fitData(1:idxFitMax,1), fitData(1:idxFitMax,colFalF(i)),...
              'Color',color01);    
    hold on;        
  end
  plot(fitData(rowFal,1), falFt,'xb');

  xlabel('Time (s)');
  ylabel('Force (N)');        
  
  subplot(2,1,2)      
  for(i=1:1:length(colFalF))
    n = (i-1)/length(colFalF);
    color01 = color0.*(1-n) + n.*color1;
    plot( fitData(1:idxFitMax,1), fitData(1:idxFitMax,colFalL(i)),...
              'Color',color01);    
    hold on;        
  end
  plot(fitData(rowFal,1), falLpDelta,'xb');
  
  xlabel('Time (s)');
  ylabel('Length (mm)');    
  hold on;
end

for(k=1:1:length(colFalF))  
  muscleCurveData.falData(k,1) = falLpDelta(k)./1000 + lmOpt;    
  muscleCurveData.falData(k,2) = falFt(k)-falFpe(k);
end

%%
%
% Simulation Data
%
%%

rampSlope      = 0;
rampStartTime  = 0;
rampStretch    = 0;

actT0          = 0;
actT1          = NaN;

colLength =0;
colForce = 0;
colNoStretch=0;

    


switch(figureNumber)

  case 7
    switch(subFigureNumber)
      case 1
        expDataFile = 'dataHerzogLeonard2002Figure6and7A.csv';     
        tmp       = csvread(expDataFile,1,0);
        idxMax = 0;
        for(z=1:1:size(tmp,1))
           if(tmp(z,1) > 0)
              idxMax = z; 
           end
        end
        expData = tmp(1:1:idxMax,:);        
        rampSlope      = (3/1000); 
        colNoStretch = 2;
        switch trialNumber
            case 1
                colLength = 7;
                colForce  = colLength-1; 
            case 2
                colLength = 9;
                colForce  = colLength-1;        
            case 3
                colLength = 13;
                colForce  = colLength-1;        
            otherwise
                assert(0,'trialNumber must be 1,2, or 3');
        end
      case 2
        expDataFile = 'dataHerzogLeonard2002Figure7B.csv';        
        expData       = csvread(expDataFile,1,0);        
        rampSlope      = (9/1000);  
        idxFitMax = size(expData,1);
        colNoStretch = 2;
        switch trialNumber
            case 1
                colLength = 7;
                colForce  = colLength-1; 
            case 2
                colLength = 9;
                colForce  = colLength-1;        
            case 3
                colLength = 11;
                colForce  = colLength-1;        
            otherwise
                assert(0,'trialNumber must be 1,2, or 3');
        end
      case 3
        expDataFile = 'dataHerzogLeonard2002Figure7C.csv';        
        expData       = csvread(expDataFile,1,0);        
        rampSlope      = (27/1000);  
        idxFitMax = size(expData,1);
        colNoStretch = 2;
        switch trialNumber
            case 1
                colLength = 7;
                colForce  = colLength-1; 
            case 2
                colLength = 9;
                colForce  = colLength-1;        
            case 3
                colLength = 11;
                colForce  = colLength-1;        
            otherwise
                assert(0,'trialNumber must be 1,2, or 3');
        end
      otherwise 
        msg = sprintf('Fig %i (%i) : Does not exist',...
                      figureNumber,subFigureNumber);
        assert(0,msg);    
    end
  otherwise
    msg = sprintf('Fig %i (%i) : Does not exist',...
                figureNumber,subFigureNumber);
    assert(0,msg);    
end

switch(trialNumber)
    case 1
        rampStretch=3/1000;
    case 2
        rampStretch=6/1000;
    case 3
        rampStretch=9/1000;
    otherwise
        assert(0,'trialNumber must be 1, 2, or 3');
end


expLengthData      = [  expData(1:idxFitMax,1),...
                        expData(1:idxFitMax,colLength)];
expDeltaLengthData = expLengthData;                    
expLengthData(:,2) = expLengthData(:,2)./1000 ...
                   + lmOpt.*(ones(size(expLengthData,1),1)); 
expForceData       = [  expData(1:idxFitMax,1),...
                        expData(1:idxFitMax,colForce )];

expForceStaticData = [  expData(1:idxFitMax,1),...
                        expData(1:idxFitMax,colNoStretch )];
                    
%%
% Extract the times the activation signal turns on/off
% Extract the time the ramps starts and ends.
%%
diffTime = diff(expData(:,1));
meanSampleTime = mean(diffTime);
stdSampleTime  = std(diffTime);
assert( abs(stdSampleTime/meanSampleTime) < 1e-3, ...
        'Sample time is not regular');

freq = 1/meanSampleTime;
filtFreq = freq*0.125;
[b,a] = butter(2,filtFreq/freq,'low');

expLengthDataFilt = zeros(size(expLengthData,1),3);
expForceDataFilt  = zeros(size(expLengthData,1),3);
expForceStaticDataFilt  = zeros(size(expLengthData,1),3);


expLengthDataFilt(:,1)  = expLengthData(:,1);
expForceDataFilt(:,1)   = expForceData(:,1);
expForceStaticDataFilt(:,1)   = expForceData(:,1);

expLengthDataFilt(:,2)  = filtfilt(b,a,expLengthData(:,2));
expForceDataFilt(:,2)   = filtfilt(b,a,expForceData(:,2));
expForceStaticDataFilt(:,2)   = filtfilt(b,a,expForceStaticData(:,2));

dt = meanSampleTime;
for k=2:1:(size(expLengthData,1)-1)
    expLengthDataFilt(k,3) = ...
        ((expLengthDataFilt(k,2)-expLengthDataFilt(k-1,2)) ...
      + (expLengthDataFilt(k+1,2)-expLengthDataFilt(k,2))) / (2*dt);

    expForceDataFilt(k,3) = ...
        ((expForceDataFilt(k,2) - expForceDataFilt(k-1,2)) ...
      + (expForceDataFilt(k+1,2)- expForceDataFilt(k,2))) / (2*dt);  
  
    expForceStaticDataFilt(k,3) =...
        ((expForceStaticDataFilt(k,2) - expForceStaticDataFilt(k-1,2)) ...
      + (expForceStaticDataFilt(k+1,2)- expForceStaticDataFilt(k,2)) )/(2*dt);  
end


%%
% Get the activation signal start and end times
%%


[dFMaxVal, idxDfMax] = max(expForceStaticDataFilt(:,3));
[dFMaxVal, idxDfMin] = min(expForceStaticDataFilt(:,3));
[maxFVal  , idxFMax]  = max(expForceData(:,2));

actT0 = NaN;
actT1 = NaN;
idxActT0 = NaN;
idxActT1 = NaN;
flagSet = 0;
k = idxDfMax;
while(flagSet == 0 && k > 2)
   if(expForceData(k,2) <= expForceData(k-1,2) && flagSet == 0)
       actT0 = expForceData(k+1,1);
       idxActT0 = k+1;
       flagSet =1;
   end
   k = k-1;
end
flagSet=0;
k = idxDfMin;
while(flagSet==0 && k > idxFMax)
   if(expForceData(k,2) >= expForceData(k-1,2) && flagSet==0)
       actT1 = expForceData(k+1,1);
       idxActT1 = k+1;
       flagSet =1;
   end
   k = k-1;
end

%%
% Get the ramp start time
%%
[maxDlVal, idxMaxDL] = max(expLengthDataFilt(:,3));

k=idxMaxDL;
flagSet=0;
lenT0 = 0;
idxRampT0 = 0;
while(flagSet==0 && k > 2)
   if(expLengthDataFilt(k,3) <= 0 && flagSet == 0)
    lenT0 = expLengthDataFilt(k+5,1);
    idxRampT0 = k+5;
    flagSet=1;
   end
   k=k-1;
end


flagPlotFilteredSignals =0;
if(flagPlotFilteredSignals==1)
   figFilt = figure;
   subplot(2,1,1);
       plot(expForceData(:,1),expForceData(:,2),...
            'Color',[0.25,0.25,0.25],'LineWidth',2);
       hold on;
       plot(expForceDataFilt(:,1),expForceDataFilt(:,2),'r');
       hold on;
       plot(expForceDataFilt(:,1),expForceDataFilt(:,3),'b');
       hold on;
       plot(expForceData(idxActT0,1),expForceData(idxActT0,2),'x');
       hold on;
       plot(expForceData(idxActT1,1),expForceData(idxActT1,2),'x');
       hold on;

       xlabel('Time (s)');
       ylabel('Force (N)');
   subplot(2,1,2);
       plot(expLengthData(:,1),expLengthData(:,2),...
            'Color',[0.25,0.25,0.25],'LineWidth',2);
       hold on;
       plot(expLengthDataFilt(:,1),expLengthDataFilt(:,2),'r');
       hold on;
       plot(expLengthDataFilt(:,1),expLengthDataFilt(:,3),'b');
       hold on;      
       plot(expLengthData(idxRampT0,1),expLengthData(idxRampT0,2),'xk');
       hold on;
       xlabel('Time (s)');
       ylabel('Length (mm)');          
end




here=1;


%excitationFcn = @(argt)calcStepFunction(argt,actT0,actT1,1);
excitationFcn = @(argt)calcStepFunction(argt,actT0,actT1,1);
activationFcn = @(argU,argA)calcFirstOrderActivationDerivative(...
  argU,argA,activationTimeConstant,deactivationTimeConstant,0);



lenT1   = lenT0 + rampStretch/rampSlope;
length0 = expLengthData(idxRampT0,2);
pathLengthFcn = @(argT)calcRampStateSharp(...
                         argT,lenT0,lenT1,length0,rampSlope);

%%
%
%%
tmp = pathLengthFcn(0);

lerr  = Inf;
fpe   = 0;
idxFpe= 0;
for k=idx49LStart:idxDelta:idxMaxF49
  err = abs(expDeltaLengthData(idxRampT0,2)-fitData(k,idxFitL49));
  if(err < lerr)
    lerr = err;
    fpe = fitData(k,idxFitF49);
    idxFpe=k;
  end
end


initialState.pathLength       = tmp(2);
initialState.pathVelocity     = tmp(1);
initialState.forceAlongTendon = expForceData(idxRampT0,2);
initialState.activeFiberForceAlongTendon = expForceData(idxRampT0,2)-fpe;
initialState.activeFiberForce = NaN;
initialState.fiberLengthAlongTendon   = NaN;
initialState.fiberLength              = NaN;
initialState.fiberForce               = NaN;
initialState.pennationAngle           = NaN;
initialState.fiberVelocityAlongTendon = NaN;
initialState.fiberVelocity            = NaN;
initialState.pennationAngularVelocity = NaN;
initialState.activation               = 0;



