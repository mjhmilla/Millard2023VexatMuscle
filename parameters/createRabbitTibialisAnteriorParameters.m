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

function [rabbitTAMusculotendonProperties, ...
          rabbitTASarcomereProperties,...
          rabbitTAActiveForceLengthData,...
          rabbitTAPassiveForceLengthData] = ...
          createRabbitTibialisAnteriorParameters(...
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
% This function uses data from Siebert et al. and Prado et al. to return a 
% series of structs that contain the necessary musculotendon properties, 
% sarcomere properties, and processed data to construct Opus 31 which is 
% a rather detailed muscle model. Please review the script for details and 
% references to the literature. Default values that have been determined 
% by simulation have been noted as such in the script
%
% Prado LG, Makarenko I, Andresen C, Krüger M, Opitz CA, Linke WA. Isoform 
% diversity of giant proteins in relation to passive and active contractile 
% properties of rabbit skeletal muscles. The Journal of general physiology. 
% 2005 Nov 1;126(5):461-80.
%
% Siebert T, Leichsenring K, Rode C, Wick C, Stutzig N, Schubert H, 
% Blickhan R, Böl M. Three-dimensional muscle architecture and 
% comprehensive dynamic properties of rabbit gastrocnemius, plantaris and 
% soleus: input for simulation studies. PLoS one. 2015 Jun 26;10(6):e0130985.
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
%   rabbitTAMusculotendonProperties
%     architectural properties and some gross mechanical properites
%   rabbitTASarcomereProperties
%     lengths of all of the various filaments of a feline soleus along with 
%     detailed information regarding the segment lengths of titin
%   rabbitTAActiveForceLengthData
%     an n by 2 matrix that contains in the first column the normalized fiber
%     length and in the second column the normalized active force developed by
%     the muscle.
%   rabbitTAPassiveForceLengthData
%     an n by 2 matrix that contains in the first column the normalized fiber
%     length and in the second column the normalized passive force developed by
%     the muscle.
%%
%scaleOptimalFiberLength               = 1.0; %user-settable parameter
%scaleMaximumIsometricTension          = 1.0; %user-settable parameter

% This is the normalized fiber length at which the extrapolated 
% passive-force-length curve is expected to develop 1 normalized force.
% Here this is the passive force length curve that is fitted to the data
% of 

%Get the default sarcomere properties for a feline soles                               
[rabbitTASarcomereProperties] =...
  getMammalianSkeletalMuscleNormalizedSarcomereProperties(...
    'rabbitTA',...
    normFiberLengthAtOneNormPassiveForce,...
    normPevkToActinAttachmentPoint,...
    normMaxActiveTitinToActinDamping,...
    ecmForceFraction,...
    titinMolecularWeightInkD,...
    projectFolders);


%Siebert T, Leichsenring K, Rode C, Wick C, Stutzig N, Schubert H, 
% Blickhan R, Böl M. Three-dimensional muscle architecture and 
% comprehensive dynamic properties of rabbit gastrocnemius, plantaris and 
% soleus: input for simulation studies. PLoS one. 2015 Jun 26;10(6):e0130985.



%Take the default
forceVelocityMultiplierAtHalfMaximumFiberVelocity = [];

%From Siebert et al S1_Table.
maximumNormalizedFiberVelocity = 16.4; % in units of norm fiber lengths/second



tendonStrainAtOneNormForce = [];

%Get the (formatted) experimental data on the active/passive
%force-length curves
useElasticTendonExp = 1;


normPlateauOffset = ...
  rabbitTASarcomereProperties.normMyosinBareHalfLength;


%Get the default musculotendon properties for the feline soleus
[rabbitTAMusculotendonProperties,...
 rabbitTAMusculotendonPropertiesExp ] = ...
  getRabbitTAMusculotendonProperties(...
            maximumNormalizedFiberVelocity,...
            forceVelocityMultiplierAtHalfMaximumFiberVelocity,...
            tendonStrainAtOneNormForce,...
            scaleOptimalFiberLength,...                              
            scaleMaximumIsometricTension,...
            normPlateauOffset,...
            useElasticTendonExp,...
            projectFolders,...
            flag_useOctave);
          

[rabbitTAActiveForceLengthData,...
rabbitTAPassiveForceLengthData] = getRabbitTAMusculotendonData(...
                              rabbitTAMusculotendonPropertiesExp,...
                              rabbitTAMusculotendonProperties,...
                              normPlateauOffset,...
                              useElasticTendonExp,...
                              projectFolders,...
                              flag_useOctave);


  

