function [felineSoleusMusculotendonProperties, ...
          felineSoleusSarcomereProperties,...
          felineSoleusActiveForceLengthData,...
          felineSoleusPassiveForceLengthData] = ...
          createFelineSoleus(scaleOptimalFiberLength,...
                             scaleMaximumIsometricTension,...
                             fitCrossBridgeStiffnessDampingToKirch199490Hz,...
                             flag_useOctave)
%%
% This function uses data from the literature to return a series of structs that
% contain the necessary musculotendon properties, sarcomere properties, and
% processed data to construct Opus 31 which is a rather detailed muscle model.
% Please review the script for details and references to the literature.
% Default values that have been determined by simulation have been noted
% as such in the script
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
%   felineSoleusMusculotendonProperties
%     architectural properties and some gross mechanical properites
%   felineSoleusSarcomereProperties
%     lengths of all of the various filaments of a feline soleus along with 
%     detailed information regarding the segment lengths of titin
%   felineSoleusActiveForceLengthData
%     an n by 2 matrix that contains in the first column the normalized fiber
%     length and in the second column the normalized active force developed by
%     the muscle.
%   felineSoleusPassiveForceLengthData
%     an n by 2 matrix that contains in the first column the normalized fiber
%     length and in the second column the normalized passive force developed by
%     the muscle.
%%
%scaleOptimalFiberLength               = 1.0; %user-settable parameter
%scaleMaximumIsometricTension          = 1.0; %user-settable parameter


%Get the default sarcomere properties for a feline soles          
flag_Cat1_Human2 = 1;                           
[felineSoleusSarcomereProperties] =...
  getMammalianSkeletalMuscleNormalizedSarcomereProperties(...
    scaleOptimalFiberLength,...
    flag_Cat1_Human2,...
    fitCrossBridgeStiffnessDampingToKirch199490Hz);

% From Scott et al. pg 211, column 2, 2nd last paragraph. Note that Scott et al.
% did their experiment on a cat soleus and that the experiments were done at
% normal body temperature.
%
% Scott SH, Brown IE, Loeb GE. Mechanics of feline soleus: I. Effect of 
% fascicle length and velocity on force output. Journal of Muscle Research & 
% Cell Motility. 1996 Apr 1;17(2):207-19.
maximumNormalizedFiberVelocity = 4.5; % in units of norm fiber lengths/second

% The slow-twitch fibers plotted in Fig. 3 of Ranatunga 1984 develop a 
% normalized force of 0.1 at half the maximum contraction velocity. I am
% assuming that the rats slow-twich normalized force-velocity curve will also
% fit that of a feline. Given that both a rat and a cat are mammals they should 
% be similar.
%
% Ranatunga KW. The force‐velocity relation of rat fast‐and slow‐twitch muscles 
% examined at different temperatures. The Journal of physiology. 1984 Jun 1;
% 351(1):517-29.
forceVelocityMultiplierAtHalfMaximumFiberVelocity = 0.1;  


% From Scott et al. pg 211 column 2, 2nd last paragraph
forceVelocityNormalizedFittingData =[-0.5/maximumNormalizedFiberVelocity,0.5,...
                                       -1/maximumNormalizedFiberVelocity,0.3];

% Fit to in-vivo data of the human Achilles tendon from Magnusson et al.
% I'd like to use something that is fitted to a feline soleus, but I haven't
% been able to find this parameter in the literature.
%
% Magnusson, S. P., Aagaard, P., Rosager, S., Dyhre-Poulsen, P., and Kjaer, M.,
% 2001, “Load–Displacement Properties of the Human Triceps Surae Aponeurosis
% In Vivo,” J. Physiol., 531(1), pp. 277–288.
%
% A better references for a cat soleus tendon:
%
% Scott SH, Loeb GE. Mechanical properties of aponeurosis and tendon of the 
% cat soleus muscle during whole‐muscle isometric contractions. Journal of 
% Morphology. 1995 Apr;224(1):73-86.
%
kisoScott = 30; %Scott & Loeb 1995: pg 80 paragraph 1
tendonStrainAtOneNormForce      = 1.375/kisoScott; %0.049; 

%Get the (formatted) experimental data on the active/passive
%force-length curves
useElasticTendonExp = 1;

normPlateauOffset = ...
  felineSoleusSarcomereProperties.normMyosinBareHalfLength*2 ...
  +0.004877-0.000022;



%Get the default musculotendon properties for the feline soleus
[felineSoleusMusculotendonProperties,...
 felineSoleusMusculotendonPropertiesExp ] = ...
  getFelineSoleusMusculotendonProperties(...
            maximumNormalizedFiberVelocity,...
            forceVelocityMultiplierAtHalfMaximumFiberVelocity,...
            tendonStrainAtOneNormForce,...
            scaleOptimalFiberLength,...                              
            scaleMaximumIsometricTension,...
            normPlateauOffset,...
            useElasticTendonExp,...
            flag_useOctave);
          

[felineSoleusActiveForceLengthData,...
felineSoleusPassiveForceLengthData] = getFelineSoleusMusculotendonData(...
                              felineSoleusMusculotendonPropertiesExp,...
                              felineSoleusMusculotendonProperties,...
                              normPlateauOffset,...
                              useElasticTendonExp,...
                              flag_useOctave);


  

