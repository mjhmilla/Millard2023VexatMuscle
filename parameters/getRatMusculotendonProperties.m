%%
% SPDX-FileCopyrightText: 2023 Matthew Millard <millard.matthew@gmail.com>
%
% SPDX-License-Identifier: MIT
%
%
%%
function [musculotendonProperties] = ...
            getRatMusculotendonProperties(...
              ratSoleusSarcomereProperties,...
              maximumNormalizedFiberVelocity,...
              forceVelocityMultiplierAtHalfMaximumFiberVelocity,...
              forceVelocityMultiplierAtLowEccentricFiberVelocity,...
              forceVelocityMultiplierAtMaximumEccentricFiberVelocity,...
              tendonStrainAtOneNormForce,...
              scaleOptimalFiberLength,...                                          
              scaleMaximumIsometricTension,...
              normPlateauOffset,...
              useElasticTendon,...
              useFibrilModel,...
              muscleName,...
              flag_useOctave)

assert(strcmp(muscleName,'SOL') || strcmp(muscleName,'EDL'),...
       'Error: muscleName is not SOL or EDL');

%disp('getRatSoleusMusculotendonProperties:');
mm2m    = 0.001;
g2N     = 9.81/1000;
deg2rad = pi/180;

lceOpt   = nan;
alphaOpt = nan; 
fiso     = nan;
ltSlk    = nan; 
etIso    = nan;
normTendonLength = nan;

if(useFibrilModel==1)
  assert(useElasticTendon==0,'Error: a fibril has no tendon');

  lceOpt   = ratSoleusSarcomereProperties.optimalSarcomereLength;
  alphaOpt = 0; 
  fiso     = 1;
  ltSlk    = 0; 
  etIso    = 0;
  normTendonLength = 1;

end

if(useFibrilModel==0)



  % Lemaire KK, Baan GC, Jaspers RT, van Soest AK. Comparison of the validity of 
  % Hill and Huxley muscle–tendon complex models using experimental data  
  % obtained from rat m. soleus in situ. Journal of experimental biology. 
  % 2016 Apr 1;219(7):977-87.

  table5Lemaire2016.fceMax = [1.27; 1.25; 1.26];          % N
  table5Lemaire2016.lceOpt = [17.1; 21.6; 18.6].*(mm2m); % mm/1000
  table5Lemaire2016.cSee   = [1.30; 0.85;  1.5];          % kN/m^2
  table5Lemaire2016.lSee0  = [19.5; 13.3; 18.6].*(mm2m); % mm/1000
  table5Lemaire2016.cPee   = [0.87; 1.21; 0.68].*(mm2m); % mm/1000
  table5Lemaire2016.lPee0  = [ 9.8; 12.4;  9.7].*(mm2m); % mm/1000

  idxRat = 1;

  lceOpt   = table5Lemaire2016.lceOpt(idxRat,1);
  alphaOpt = 0; 
  fiso     = table5Lemaire2016.fceMax(idxRat,1);
  ltSlk    = table5Lemaire2016.lSee0(idxRat,1); 

  %Taking the inverse of Eqn. 3 of Lemaire et al. to solve for the tendon
  %length under fiso
  lSee1 = (sqrt( table5Lemaire2016.fceMax(idxRat,1)...
                  /table5Lemaire2016.cSee(idxRat,1)) ...
           +  table5Lemaire2016.lSee0(idxRat,1)      );

  etIso = (lSee1-lSee0)/lSee0;
  
  if(isempty(tendonStrainAtOneNormForce)==0 ...
    || isnan(tendonStrainAtOneNormForce)==1)
    etIso = tendonStrainAtOneNormForce;
  end

  normTendonLength = 1+etIso;

  if(strcmp(muscleName,'EDL')==1)
      %From: Table II of
      %
      %W. L. Johnson, D. L. Jindrich, H. Zhong, R. R. Roy and V. R. Edgerton, 
      % "Application of a Rat Hindlimb Model: A Prediction of Force Spaces 
      % Reachable Through Stimulation of Nerve Fascicles," in IEEE 
      % Transactions on Biomedical Engineering, vol. 58, no. 12, pp. 
      % 3328-3338, Dec. 2011, doi: 10.1109/TBME.2011.2106784. 
    
      fiso     = 225*g2N; %225g 
      lceOpt   = 13.7*mm2m; 
      alphaOpt = 10*deg2rad;
      ltSlk    = 9*mm2m;
      etIso    = 0.033; % Johnson et al. took the default value from Zajac

      if(isempty(tendonStrainAtOneNormForce)==0 ...
        || isnan(tendonStrainAtOneNormForce)==1)
        etIso = tendonStrainAtOneNormForce;
      end      


  end


end

lceOpt = lceOpt*scaleOptimalFiberLength;
fiso   = fiso*scaleMaximumIsometricTension;

minimumFiberLengthAlongTendon = sqrt(eps);

lcePerp    = lceOpt*sin(alphaOpt);
lceAT      = minimumFiberLengthAlongTendon;

alphaMax        = atan2(lcePerp, ...
                       lceAT);
minimumFiberLength = sqrt(lcePerp*lcePerp + lceAT*lceAT);

dlceMaxN=maximumNormalizedFiberVelocity;

nameMTU = ['Rat ',muscleName];
abbrMTU = ['rat',muscleName];


musculotendonProperties = struct(...
        'name'                              , nameMTU,  ... 
        'abbr'                              , abbrMTU,     ...
        'fiso'                              , fiso,         ...  
        'optimalFiberLength'                , lceOpt,       ... 
        'pennationAngle'                    , alphaOpt,     ... 
        'pennationAngleAtMinimumFiberLength', alphaMax,...
        'tendonSlackLength'                 , ltSlk,        ...  
        'tendonStrainAtOneNormForce'        , etIso,        ...
        'normTendonDampingLinear'           , 0.0565,          ...
        'normTendonDampingConstant'         , 0.0,           ...
        'normNumericalDamping'              , 1e-4,         ...
        'minimumFiberLength'                , minimumFiberLength,    ...
        'minimumFiberLengthAlongTendon'     , minimumFiberLengthAlongTendon, ...
        'normPassiveForceLengthCurveStiffnessAtSlack', 0.001,...
        'maximumNormalizedFiberVelocity',     maximumNormalizedFiberVelocity,...
        'forceVelocityMultiplierAtHalfMaximumFiberVelocity',...
          forceVelocityMultiplierAtHalfMaximumFiberVelocity,...
        'forceVelocityMultiplierAtLowEccentricFiberVelocity',...
         forceVelocityMultiplierAtLowEccentricFiberVelocity,... %1.15
        'forceVelocityMultiplierAtMaximumEccentricFiberVelocity',...
         forceVelocityMultiplierAtMaximumEccentricFiberVelocity,... %1.30
        'appliedFiberLengthScaling', ...
         scaleOptimalFiberLength,...
        'appliedMaxIsometricForceScaling', ...
         scaleMaximumIsometricTension);


