function [rabbitPsoasMusculotendonProperties, ...
          rabbitPsoasSarcomereProperties] = ...
          createRabbitPsoasFibrilParameters(...
            scaleOptimalFiberLength,...
            scaleMaximumIsometricTension,...
            normFiberLengthAtOneNormPassiveForce,...
            normPevkToActinAttachmentPoint,...
            normMaxActiveTitinToActinDamping,...
            ecmForceFraction,...
            titinMolecularWeightInkD,...
            flag_useOctave)
%%
% This function uses data from the literature to return a series of structs that
% contain the necessary musculotendon properties, sarcomere properties, and
% processed data to construct a rabbit psoas fibril model for Opus 31 
% which is a rather detailed muscle model.
% This is a minimal extension of createFelineSoleus: actin, myosin, and 
% titin geometry do reflect a rabbit psoas. However, the architectural
% properties are just stand-in-values since I cannot find any data in the 
% literature on the architectural properties of a rabbit psoas. This does
% not matter too much, since this model is being used to replicate fibril
% experiments.
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
%   rabbitPsoasMusculotendonProperties
%     architectural properties and some gross mechanical properites
%   rabbitPsoasSarcomereProperties
%     lengths of all of the various filaments of a feline soleus along with 
%     detailed information regarding the segment lengths of titin
%   rabbitPsoasActiveForceLengthData
%     an n by 2 matrix that contains in the first column the normalized fiber
%     length and in the second column the normalized active force developed by
%     the muscle. (presently empty)
%   rabbitPsoasPassiveForceLengthData
%     an n by 2 matrix that contains in the first column the normalized fiber
%     length and in the second column the normalized passive force developed by
%     the muscle. (presently empty)
%%
%scaleOptimalFiberLength               = 1.0; %user-settable parameter
%scaleMaximumIsometricTension          = 1.0; %user-settable parameter


%Get the default sarcomere properties for a feline soles          
                          
[rabbitPsoasSarcomereProperties] =...
  getMammalianSkeletalMuscleNormalizedSarcomereProperties(...
    'rabbit',...
    normFiberLengthAtOneNormPassiveForce,...
    normPevkToActinAttachmentPoint,...
    normMaxActiveTitinToActinDamping,...
    ecmForceFraction,...
    titinMolecularWeightInkD);

% I have no data on the force-velocity characteristics of a rabbit psoas 
% specifically, so I'm taking the value of 4.5 fiber lengths/s from Scott
% et al. for a cat soleus.
%
% Scott SH, Brown IE, Loeb GE. Mechanics of feline soleus: I. Effect of 
% fascicle length and velocity on force output. Journal of Muscle Research & 
% Cell Motility. 1996 Apr 1;17(2):207-19.
maximumNormalizedFiberVelocity = 4.5; % in units of norm fiber lengths/second


% I have no data on the force-velocity characteristics of a rabbit psoas 
% specifically, so I'm taking the value of 0.1 fiber lengths/s from 
% slow-twitch fibers plotted in Fig. 3 of Ranatunga 1984 develop a 
% normalized force of 0.1 at half the maximum contraction velocity.
%
% Ranatunga KW. The force‐velocity relation of rat fast‐and slow‐twitch muscles 
% examined at different temperatures. The Journal of physiology. 1984 Jun 1;
% 351(1):517-29.
forceVelocityMultiplierAtHalfMaximumFiberVelocity = 0.1;  


% From Scott et al. pg 211 column 2, 2nd last paragraph
forceVelocityNormalizedFittingData =[-0.5/maximumNormalizedFiberVelocity,0.5,...
                                       -1/maximumNormalizedFiberVelocity,0.3];

% I have no data on the force-velocity characteristics of a rabbit psoas 
% specifically, so I'm using the values from Scott. These are pretty
% similar even to the values typically used for a human Achilles tendon, so
% I'm reasonably comfortable that the rabbit psoas tendon strains will not
% be too different.
%
% Of course, this parameter has no impact on the skinned fibril experiments
% I'm replicating: there is no tendon.
%
% Scott SH, Loeb GE. Mechanical properties of aponeurosis and tendon of the 
% cat soleus muscle during whole‐muscle isometric contractions. Journal of 
% Morphology. 1995 Apr;224(1):73-86.
%
kisoScott = 30; %Scott & Loeb 1995: pg 80 paragraph 1
tendonStrainAtOneNormForce      = 1.375/kisoScott; 

%Get the (formatted) experimental data on the active/passive
%force-length curves
useElasticTendonExp = 1;

normPlateauOffset = ...
  rabbitPsoasSarcomereProperties.normMyosinBareHalfLength;

%Get the default musculotendon properties for the feline soleus
[rabbitPsoasMusculotendonProperties] = ...
  getRabbitPsoasMusculotendonProperties(...
            maximumNormalizedFiberVelocity,...
            forceVelocityMultiplierAtHalfMaximumFiberVelocity,...
            tendonStrainAtOneNormForce,...
            scaleOptimalFiberLength,...                              
            scaleMaximumIsometricTension,...
            normPlateauOffset,...
            useElasticTendonExp,...
            flag_useOctave);
          



  

