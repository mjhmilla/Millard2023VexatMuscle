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

function [ratSoleusMusculotendonProperties, ...
          ratSoleusSarcomereProperties,...
          activeForceLengthDataNormalized,...
          passiveForceLengthDataNormalized] = ...
          createRatSoleusParameters(...
            scaleOptimalFiberLength,...
            scaleMaximumIsometricTension,...
            activeForceLengthData,...
            passiveForceLengthData,...            
            normFiberLengthAtZeroPassiveForce,...
            normFiberLengthAtOneNormPassiveForce,...
            normFiberStiffnessAtOneNormPassiveForce,...
            normPevkToActinAttachmentPoint,...
            normMaxActiveTitinToActinDamping,...
            ecmForceFraction,...
            titinMolecularWeightInkD,...
            specimenTemperature,...
            makeFibrilModel,...
            useElasticTendon,...
            mapToEDLModel,...
            projectFolders,...
            flag_useOctave)
%%
% This function uses data from the literature to return a series of structs that
% contain the necessary musculotendon properties, sarcomere properties, and
% processed data to construct a rat psoas fibril model for Opus 31 
% which is a rather detailed muscle model.
% This is a minimal extension of createFelineSoleus: actin, myosin, and 
% titin geometry do reflect a rat soleus. 
%
% @param scaleOptimalFiberLength: scales the optimal fiber length
%
% @param scaleMaximumIsometricTension: scales the maximum isometric tension
%
% @param flag_useOctave
%   Setting this to 1 will ensure that no parts of the code that are 
%   incompatible with octave are called.
%
% @return Four structs:
%   ratSoleusMusculotendonProperties
%     architectural properties and some gross mechanical properites
%   ratSoleusSarcomereProperties
%     lengths of all of the various filaments of a rat soleus along with 
%     detailed information regarding the segment lengths of titin
%   ratSoleusActiveForceLengthData
%     an n by 2 matrix that contains in the first column the normalized fiber
%     length and in the second column the normalized active force developed by
%     the muscle. (presently empty)
%   ratSoleusPassiveForceLengthData
%     an n by 2 matrix that contains in the first column the normalized fiber
%     length and in the second column the normalized passive force developed by
%     the muscle. (presently empty)
%%
%scaleOptimalFiberLength               = 1.0; %user-settable parameter
%scaleMaximumIsometricTension          = 1.0; %user-settable parameter


%Get the default sarcomere properties for a rat soleus          
             
nameSarcomerePropertiers = 'ratSOL';
if(mapToEDLModel==1)
    nameSarcomerePropertiers = 'ratEDL';
end

[ratSoleusSarcomereProperties] =...
  getMammalianSkeletalMuscleNormalizedSarcomereProperties(...
    nameSarcomerePropertiers,...
    normFiberLengthAtOneNormPassiveForce,...
    normPevkToActinAttachmentPoint,...
    normMaxActiveTitinToActinDamping,...
    ecmForceFraction,...
    titinMolecularWeightInkD,...
    projectFolders);
%
% Note: ratSol & ratEDL have the same 
%

%
% This is a new parameter for sarcomere parameters
%
ratSoleusSarcomereProperties.normFiberLengthAtZeroPassiveForce=...
    normFiberLengthAtZeroPassiveForce;

ratSoleusSarcomereProperties.normFiberStiffnessAtOneNormPassiveForce =...
    normFiberStiffnessAtOneNormPassiveForce;
% 
% From Tomalka et al. (skinned fibril at 12'C)
%
% vceMax      = 0.46
% fvCurvature = 0.07 
%
% From Degens et al. (skinned fibril at 12'C)
%
% YM 1.10 +/- 0.44
% YF 1.02 +/- 0.35
% OM 0.62 +/- 0.27
% OF 0.65 +/- 0.53
%
% From Ranatunga et al. (whole muscle), the curvature of the soleus (a/Po)
% goes from
%
% 35'C    30'C    25'C    20'C
% 0.212   0.18    0.157   0.137
%
% Fitting a linear (lsq) model to this data and extrapolating to 12'C we would 
% have:
%
% c = 0.0946
%
% Weidner suggests that I use Degens et al. as his measurements on subsequent
% experiments recorded maximum shortening velocities of ~1 lo/s.
%
% Tomalka et al. used 3 month old female rats, so I'm using the YF data 
% from Degens et al. 
%
% Tomalka A, Weidner S, Hahn D, Seiberl W, Siebert T. Power amplification 
% increases with contraction velocity during stretch-shortening cycles of 
% skinned muscle fibers. Frontiers in physiology. 2021 Mar 31;12:644981.
%
% Degens H, Yu F, Li X, Larsson L. Effects of age and gender on shortening 
% velocity and myosin isoforms in single rat muscle fibres. Acta physiologica 
% scandinavica. 1998 May;163(1):33-40.
%
% Ranatunga KW. The force‐velocity relation of rat fast‐and slow‐twitch muscles 
% examined at different temperatures. The Journal of physiology. 1984 
% Jun 1;351(1):517-29.
%
% Ranatunga KW. Temperature‐dependence of shortening velocity and rate of 
% isometric tension development in rat skeletal muscle. The Journal of 
% Physiology. 1982 Aug 1;329(1):465-83.
%


maximumNormalizedFiberVelocity = 1.02; % in units of norm fiber lengths/second

if(mapToEDLModel==1)
    maximumNormalizedFiberVelocity = 2.25;    
end
% 2.25 Lo/s from pg 4, column 1, paragraph 3
%
% Tomalka A, Rode C, Schumacher J, Siebert T. The active force–length 
% relationship is invisible during extensive eccentric contractions in 
% skinned skeletal muscle fibres. Proceedings of the Royal Society B: 
% Biological Sciences. 2017 May 17;284(1854):20162497.

%
% To map curvature to fv-at-half-vceMax
%
% given Hill's f-v hyperbola 
%
%   (P+a)(V+b) = (Po+a)b
%
% with a vceMax of 1.02 and a curvature of 0.0946 we end up with a value of
% fvN of 0.0796 at 0.5*vceMax.
%
% as described in
%
% Alcazar J, Csapo R, Ara I, Alegre LM. On the Shape of the Force-Velocity 
% Relationship in Skeletal Muscles: The Linear, the Hyperbolic, and the 
% Double-Hyperbolic. Front Physiol. 2019 Jun 19;10:769. 
% doi: 10.3389/fphys.2019.00769. PMID: 31275173; PMCID: PMC6593051.
%

% Fitting a linear model to the termperature-curvature data in
% Ranatunga 1982 (Table 1) to get the soleus curvature at 12 C. Note
% that this data comes from whole muscle.
%
% Ranatunga KW. Temperature‐dependence of shortening velocity and rate of 
% isometric tension development in rat skeletal muscle. The Journal of 
% Physiology. 1982 Aug 1;329(1):465-83


t=[35, 30, 25, 20]';
c=[0.212, 0.18, 0.157, 0.137]';

if(mapToEDLModel==1)
    t=[35, 30, 25, 20]';
    c=[0.415, 0.429, 0.385, 0.283]';
end

% c0 + c1*t = c
% A*x = b 
% [35, 1](c1) = 0.212 
% [30, 1](c0) = 0.18 
% [25, 1]     = 0.157 
% [20, 1]     = 0.137

A = [t , ones(size(t))];
b = c;

% Ax = b
% A'Ax = A'b
% x = (A'A)\(A'b)

x = (A'*A)\(A'*b);

% Expected curvature value at 12 C in whole soleus muscle
c12 = [specimenTemperature,1]*x;

% Solving for fv at vceMax*0.5


% vmax for YF
%
% Degens H, Yu F, Li X, Larsson L. Effects of age and gender on shortening 
% velocity and myosin isoforms in single rat muscle fibres. Acta physiologica 
% scandinavica. 1998 May;163(1):33-40.

halfMaximumNormalizedFiberVelocity = maximumNormalizedFiberVelocity*0.5;

Po = 1;
c = c12;
b  = c*maximumNormalizedFiberVelocity;

forceVelocityMultiplierAtHalfMaximumFiberVelocity = ...
  ((1+c)*b - c.*(halfMaximumNormalizedFiberVelocity+b)) ...
  ./ (halfMaximumNormalizedFiberVelocity+b);  

flag_debugFv = 0;
if(flag_debugFv==1)
  figDebugFv = figure;
  v = [0.01:0.01:1].*maximumNormalizedFiberVelocity;
  fv = ((1+c)*b - c.*(v+b)) ./ (v+b);
  plot(v,fv,'-k');
  hold on;
  plot(halfMaximumNormalizedFiberVelocity,...
       forceVelocityMultiplierAtHalfMaximumFiberVelocity,...
       'xb');
  xlabel('Norm. Velocity ($$\ell / \ell^M_o$$)');
  ylabel('Norm. Force ($$f / f^M_o$$)'); 
  title('Force-velocity curve and curvature'); 
end

forceVelocityMultiplierAtLowEccentricFiberVelocity     = 1.35;
forceVelocityMultiplierAtMaximumEccentricFiberVelocity = 1.40;

kisoScott                       = nan;
tendonStrainAtOneNormForce      = nan;



normPlateauOffset = ...
  ratSoleusSarcomereProperties.normMyosinBareHalfLength;



%Get the default musculotendon properties for the feline soleus
[ratSoleusMusculotendonProperties] = ...
  getRatSoleusMusculotendonProperties(...
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
            makeFibrilModel,...
            mapToEDLModel,...
            flag_useOctave);


ratSoleusMusculotendonProperties.temperature=specimenTemperature;

activeForceLengthDataNormalized = [];
passiveForceLengthDataNormalized = [];

if(isempty(activeForceLengthData)==0)
    activeForceLengthDataNormalized=activeForceLengthData;
    
    activeForceLengthDataNormalized= ...
           activeForceLengthDataNormalized...
           ./[ratSoleusMusculotendonProperties.optimalFiberLength,...
              ratSoleusMusculotendonProperties.fiso];
end

if(isempty(passiveForceLengthData)==0)
    passiveForceLengthDataNormalized=passiveForceLengthData;          
    
    passiveForceLengthDataNormalized= ...
           passiveForceLengthDataNormalized...
           ./[ratSoleusMusculotendonProperties.optimalFiberLength,...
              ratSoleusMusculotendonProperties.fiso];
end

  

