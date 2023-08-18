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

function [musculotendonProperties,...
          musculotendonPropertiesExp] = ...
            getRabbitTAMusculotendonProperties(...
              maximumNormalizedFiberVelocity,...
              forceVelocityMultiplierAtHalfMaximumFiberVelocity,...
              tendonStrainAtOneNormForce,...
              scaleOptimalFiberLength,...                                          
              fisoScale,...
              normPlateauOffset,...
              useElasticTendon,...
              projectFolders,...
              flag_useOctave)                        
%%
% This function uses data from Siebert et al. and Winters et al. to build 
% the table of musculotendon properties for a Rabbit TA
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
%   Siebert T, Leichsenring K, Rode C, Wick C, Stutzig N, Schubert H, 
%   Blickhan R, Böl M. Three-dimensional muscle architecture and 
%   comprehensive dynamic properties of rabbit gastrocnemius, plantaris and 
%   soleus: input for simulation studies. PLoSlceOptExp* one. 2015 Jun 26;10(6):e0130985.
%
%   Winters TM, Takahashi M, Lieber RL, Ward SR. Whole muscle length-tension 
%   relationships are accurately modeled as scaled sarcomeres in rabbit 
%   hindlimb muscles. Journal of biomechanics. 2011 Jan 4;44(1):109-15.%
%
%%


%extracted from Siebert et al.
lceOptExpSiebert    = 36.7/1000;

%lopt from Siebert et al. 2015 is defined as the shortest length of the
%plateau. All of the curves in this code base assume that lopt is in the
%middle of the plateau, and so we must add an offset.
lceOptExp           = lceOptExpSiebert*(1+0.5*(1.16-1.0));% 
lceOpt              = lceOptExp*scaleOptimalFiberLength;

alphaOpt    = 0; %No data       
ltSlk       = 56.2/1000;
      
%Scale it as desired
%From Winters et al.
fisoExp     = 18.5;
fiso        = fisoScale*fisoExp;

minimumFiberLengthAlongTendon = sqrt(eps);

lcePerp     = lceOpt*sin(alphaOpt);
lceAT       = minimumFiberLengthAlongTendon;

alphaMax        = atan2(lcePerp, ...
                       lceAT);
minimumFiberLength = sqrt(lcePerp*lcePerp + lceAT*lceAT);

if(isempty(forceVelocityMultiplierAtHalfMaximumFiberVelocity))

    
    %Use the force-velocity model of Siebert et al to evaluate the normalized
    %force at half the maximum contraction velocity
    vmax = -maximumNormalizedFiberVelocity;
    curv = 0.33;
    vcc =  vmax*0.5;
    
    fvHalf = (vmax-vcc)/(vmax + vcc/curv);
    
    forceVelocityMultiplierAtHalfMaximumFiberVelocity = fvHalf;  

end

etIso = tendonStrainAtOneNormForce;
if(isempty(tendonStrainAtOneNormForce))
    %From Siebert et al S1_Table.
    lsec0 = 0.0562;
    lsec1 = lsec0*0.027;
    f1N   = 0.2;
    k     = 9/(1/1000); %N/mm -> N/m
    kN    = k/fisoExp;
    dlsec = (1-f1N)/kN;
    lsecAtFiso = dlsec + lsec1;
    
    etIso      = lsecAtFiso /lsec0;
end

if(useElasticTendon==0)
    etIso=1;
end


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
        'optimalFiberLengthSiebert'   , lceOptExpSiebert,...
        'name'                        , 'Siebert2015');




