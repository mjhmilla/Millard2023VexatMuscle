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
          ratSoleusSarcomereProperties] = ...
          createRatSoleusFibrilParameters(...
            scaleOptimalFiberLength,...
            scaleMaximumIsometricTension,...
            normFiberLengthAtOneNormPassiveForce,...
            normPevkToActinAttachmentPoint,...
            normMaxActiveTitinToActinDamping,...
            ecmForceFraction,...
            titinMolecularWeightInkD,...
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
                          
[ratSoleusSarcomereProperties] =...
  getMammalianSkeletalMuscleNormalizedSarcomereProperties(...
    'ratSOL',...
    normFiberLengthAtOneNormPassiveForce,...
    normPevkToActinAttachmentPoint,...
    normMaxActiveTitinToActinDamping,...
    ecmForceFraction,...
    titinMolecularWeightInkD,...
    projectFolders);

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

forceVelocityMultiplierAtHalfMaximumFiberVelocity = 0.0796;  



kisoScott                       = nan;
tendonStrainAtOneNormForce      = nan;

useElasticTendonExp = 0;

normPlateauOffset = ...
  ratSoleusSarcomereProperties.normMyosinBareHalfLength;

%Get the default musculotendon properties for the feline soleus
[ratSoleusMusculotendonProperties] = ...
  getRatSoleusMusculotendonProperties(...
            maximumNormalizedFiberVelocity,...
            forceVelocityMultiplierAtHalfMaximumFiberVelocity,...
            tendonStrainAtOneNormForce,...
            scaleOptimalFiberLength,...                              
            scaleMaximumIsometricTension,...
            normPlateauOffset,...
            useElasticTendonExp,...
            flag_useOctave);
          



  

