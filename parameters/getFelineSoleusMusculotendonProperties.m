function [musculotendonProperties,...
          musculotendonPropertiesExp] = getFelineSoleusMusculotendonProperties(...
                                      maximumNormalizedFiberVelocity,...
                                      forceVelocityMultiplierAtHalfMaximumFiberVelocity,...
                                      tendonStrainAtOneNormForce,...
                                      scaleOptimalFiberLength,...                                          
                                      fisoScale,...
                                      normPlateauOffset,...
                                      useElasticTendon,...
                                      flag_useOctave)                        
%%
% This function uses data from Herzog & Leonard to build a model of the
% cat soleus that is used in Figure 7 of the paper. Missing details
% (pennation angle, tendon slack length) are filled in by Sacks et al.
%
% The following parameters are numerical in nature and as such are chosen
% to be small and non-zero.
%
%        'normNumericalDamping'        
%        'minimumFiberLength'         
%        'minimumFiberLengthAlongTendon'
%        'normPassiveForceLengthCurveStiffnessAtSlack'
%
% References
%   Herzog W, Leonard TR. Force enhancement following stretching of skeletal 
%   muscle: a new mechanism. Journal of Experimental Biology. 2002 May 
%   1;205(9):1275-83.
%
%   Sacks RD, Roy RR. Architecture of the hind limb muscles of cats: functional 
%   significance. Journal of Morphology. 1982 Aug;173(2):185-95.
%
%
%%

%Load in the passive and active force-length trials that to fit the musculotendon
%parameters to
fitData =csvread('experiments/HerzogLeonard2002/data/dataHerzogLeonard2002Figure6and7A.csv',1,0);

%Note: 'F43', 'L43', etc are the column headings in the file.
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
mm2m = 1/1000;

%extracted from Herzog & Leonard 2002
lceOptExp = (27/0.63)/1000; % From the text: pg 1277 paragraph 2;
lceOpt   = lceOptExp*scaleOptimalFiberLength;

% 9 April 2020
% Previously I was using Sacks & Roy to fill in all additional information 
% about the cat soleus muscle architecture that I was missing. It has bencome
% clear that the soleus musculotendon length (84.7mm) and fiber length (41.7mm) 
% they report results in a tendon that is 43 mm long. All of their measurements
% were made on slack tissue (I believe).
%
% Scott & Loeb made detailed measurements on the length changes of tendon,
% apneurosis, and on muscle fibers in response to loading using piezoelectric
% crystals. The cats used in the study are of comparable size to those used
% by Sacks & Roy. The basic architectural properties differ in a few respects:
% soleus musculotendon length (38 + 27 = 65 mm), fiber length (38 mm +/- 0.6),
% and tendon (27mm +/- 0.3). These measurements were made between loads of zero
% to one maximum isometric force allowing them to establish the elasticity of
% the tendon and the apneurosis.
%
% At least for the length of the tendon I will be using Scott & Loeb's estimates
% which should be more accurate because they made measurements across a range
% of physiologically relevant loads.
%
% Sacks RD, Roy RR. Architecture of the hind limb muscles of cats: functional 
% significance. Journal of Morphology. 1982 Aug;173(2):185-95.
%
% Scott SH, Loeb GE. Mechanical properties of aponeurosis and tendon of the 
% cat soleus muscle during whole‐muscle isometric contractions. Journal of 
% Morphology. 1995 Apr;224(1):73-86.
%

%Additional parameters from Sacks et al.
% From Table 1
alphaOpt = 7.0*(pi/180); 

% Taking the mean tendon length from Table 1 of Scott et al. and scaling it
% by the ratio of the optimal fiber length and the mean fiber length we can
% get an estimate of the soleus tendon length for the cat used in Herzog 
% and Leonard's experiment.
% 
% Note: here we are treating the apneurosis as rigid.              
mtLengthSacksRoy      = 84.7/1000;
fiberLengthSacksRoy   = 41.7/1000;
lengthScalingSacksRoy = lceOpt/(fiberLengthSacksRoy);
ltSlkSacksRoy         = (mtLengthSacksRoy-fiberLengthSacksRoy)*lengthScalingSacksRoy;  

fiberLengthScottLoeb   = 38.0/1000;
lengthScalingScottLoeb = lceOpt / (fiberLengthScottLoeb);
ltSlkScottLoeb         = (27.0/1000)*lengthScalingScottLoeb;

ltSlk = ltSlkScottLoeb;

%Set the elasticity of the tendon as desired
etIso = tendonStrainAtOneNormForce;
normTendonLength = 1+etIso;
if(useElasticTendon==0)
  etIso = 0;
  normTendonLength = 1;
end


%%
%Indentify the maximum-isometric-force of the soleus from the 
%data of Herzog & Leonard. This has a few steps:
%
% 1. Create a tendon model to estimate the fiber length given the
%    path length and the load.
% 2. Use the passive-force-length data to estimate the passive-force-length
%    curve of the CE
% 3. Use the tendon and fiber passive-force-length models to get the
%    length of the CE and the active force developed by the CE
%
%%

% 1. Create a tendon model to estimate the fiber length given the
%    path length and the load.
tendonForceLengthCurve = [];

if(useElasticTendon==1)
  kIso            = 1.375/etIso;
  fToe            = 2.0/3.0;
  curviness       = 0.5;
  computeIntegral = 0;
  minimumSlope = sqrt(eps);
  flag_enableNumericallyNonZeroGradients = 0;
  smallNumericallyNonZeroNumber = sqrt(sqrt(eps));
  tendonForceLengthCurve = ...
    createTendonForceLengthCurve2019( ...
        etIso, kIso, fToe, curviness, computeIntegral, ...
        flag_enableNumericallyNonZeroGradients,...
        smallNumericallyNonZeroNumber,...
        'temporary',flag_useOctave);

end 

    
    
%Back to using Herzog & Leonard
%Extract the maximum tension recorded at the optimal fiber length
[ftMax idxFtMax] = max(fitData(:,idxFitF43));

lceOptLeft = lceOpt*(1-normPlateauOffset);
fibKin = calcFixedWidthPennatedFiberKinematicsAlongTendon(...
                                        lceOptLeft,...
                                        0,...
                                        lceOpt,...
                                        alphaOpt);

alphaOptLeft = fibKin.pennationAngle;  

ftN = 1*cos(alphaOptLeft);
ltN = 1;
if(useElasticTendon==1)
  ltN = calcBezierFcnXGivenY(ftN,tendonForceLengthCurve); 
end

lmtOptFiso = fitData(idxFtMax, idxFitL43)*mm2m ...
        + (lceOptLeft)*cos(alphaOptLeft)...
        + (ltN*ltSlk);
lmtOptZero = lmtOptFiso;



% 2. Use the passive-force-length data to estimate the passive force developed
%    by the CE at its optimal fiber length

fp      = 0;
idxBest = 0;
fpBest  = 0;
lerr    = Inf;
lceBest = 0;

idxStart = 2440;
idxEnd   = 5454;

for k=idxStart:1:idxEnd
    ftN = ( fitData(       k, idxFitF50)*cos(alphaOptLeft) ...
           /fitData(idxFtMax, idxFitF43));
    ltN=1;
    if(useElasticTendon)
      ltN = calcBezierFcnXGivenY(ftN,tendonForceLengthCurve); 
    end    
    
    lceAT = (lmtOptZero+fitData(k,idxFitL50)*mm2m)-ltN*ltSlk;
    fibKin  = calcFixedWidthPennatedFiberKinematics(lceAT,0,lceOpt,alphaOpt);     
    lce = fibKin.fiberLength;  
    
    if( abs(lce - lceOptLeft) < abs(lerr))
        idxBest = k;
        lerr    = lce - lceOptLeft;
        fpBest  = fitData(k,idxFitF50);
        if(fpBest > 1)
          here=1;
        end
        lceBest = lce;
    end
end
fp = fpBest;

% 3. Use the tendon and fiber passive-force-length models to get the
%    length of the CE and the active force developed by the CE

fisoExp  = (ftMax-fp)/cos(alphaOptLeft); 

%Scale it as desired
fiso     = fisoScale*fisoExp;



minimumFiberLengthAlongTendon = sqrt(eps);

lcePerp    = lceOpt*sin(alphaOpt);
lceAT      = minimumFiberLengthAlongTendon;

alphaMax        = atan2(lcePerp, ...
                       lceAT);
minimumFiberLength = sqrt(lcePerp*lcePerp + lceAT*lceAT);

dlceMaxN=maximumNormalizedFiberVelocity;

musculotendonProperties = struct(...
        'name'                        , 'Feline Soleus',  ... 
        'abbr'                        , 'fSol',     ...
        'fiso'                        , fiso,         ...  
        'optimalFiberLength'          , lceOpt,       ... 
        'pennationAngle'              , alphaOpt,     ... 
        'pennationAngleAtMinimumFiberLength',alphaMax,...
        'tendonSlackLength'           , ltSlk,        ...  
        'tendonStrainAtOneNormForce'  , etIso,        ...
        'normTendonDampingLinear'     , 0.0565,       ...
        'normTendonDampingConstant'   , 0, ...
        'normNumericalDamping'        , 1e-4,         ...
        'minimumFiberLength'          ,       minimumFiberLength,    ...
        'minimumFiberLengthAlongTendon',      minimumFiberLengthAlongTendon, ...
        'normPassiveForceLengthCurveStiffnessAtSlack', 0.001,...
        'maximumNormalizedFiberVelocity',     maximumNormalizedFiberVelocity,...
        'forceVelocityMultiplierAtHalfMaximumFiberVelocity',...
          forceVelocityMultiplierAtHalfMaximumFiberVelocity,...
        'forceVelocityMultiplierAtLowEccentricFiberVelocity',...
        1.30,... %1.15
        'forceVelocityMultiplierAtMaximumEccentricFiberVelocity',...
        1.45,... %1.30
        'appliedFiberLengthScaling', scaleOptimalFiberLength,...
        'appliedMaxIsometricForceScaling', fisoScale);

musculotendonPropertiesExp = struct(...
        'fiso'                     , fisoExp,      ...
        'optimalFiberLength'       , lceOptExp,    ...
        'name'                        , 'HerzogLeonard2002');




