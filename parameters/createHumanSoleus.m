function [humanSoleusMusculotendonProperties, ...
          humanSoleusSarcomereProperties,...
          humanSoleusActiveForceLengthData,...
          humanSoleusPassiveForceLengthData] = ...
          createHumanSoleus(scaleOptimalFiberLength,...
                             scaleMaximumIsometricTension,...
                             fitCrossBridgeStiffnessDampingToKirch199490Hz,...
                             flag_useOctave)
%%
% This function will generate the structs and parameter values that are
% consistent with the geometry of a human. Specifically the optimal fiber length
% active-force-length curve, and sarcomere filament geometry (actin, myosin, and
% titin) have been taken from published human soleus data in the literature.
% Where necessary, normalized data has been taken from a cat soleus and scaled
% to be appropriate for a human soleus fiber. At the moment I have used scaled
% feline soleus data for the passive-force-length curve, and the force velocity
% curve.
%
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
%   humanSoleusMusculotendonProperties
%     architectural properties and some gross mechanical properites
%   humanSoleusSarcomereProperties
%     lengths of all of the various filaments of a feline soleus along with 
%     detailed information regarding the segment lengths of titin
%   humanSoleusActiveForceLengthData
%     an n by 2 matrix that contains in the first column the normalized fiber
%     length and in the second column the normalized active force developed by
%     the muscle.
%   humanSoleusPassiveForceLengthData
%     an n by 2 matrix that contains in the first column the normalized fiber
%     length and in the second column the normalized passive force developed by
%     the muscle.
%%
%scaleOptimalFiberLength               = 1.0; %user-settable parameter
%scaleMaximumIsometricTension          = 1.0; %user-settable parameter


%Get the default sarcomere properties for a feline soles          
flag_Cat1_Human2 = 2;                           
[humanSoleusSarcomereProperties] =...
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


tendonStrainAtOneNormForce      = 0.049; 
%Average of the data from Maganaris 2002
% Maganaris, C. N., and Paul, J. P., 2002, “Tensile Properties of the In Vivo
% Human Gastrocnemius Tendon,” J. Biomech., 35(12), pp. 1639–1646.
%

%Get the (formatted) experimental data on the active/passive
%force-length curves
useElasticTendonExp = 1;

normPlateauOffset = ...
  humanSoleusSarcomereProperties.normMyosinBareHalfLength*2 ...
  +0.004877-0.000022;



%Get the default musculotendon properties for the feline soleus
[humanSoleusMusculotendonProperties,...
 humanSoleusMusculotendonPropertiesExp ] = ...
  getFelineSoleusMusculotendonProperties(...
            maximumNormalizedFiberVelocity,...
            forceVelocityMultiplierAtHalfMaximumFiberVelocity,...
            tendonStrainAtOneNormForce,...
            scaleOptimalFiberLength,...                              
            scaleMaximumIsometricTension,...
            normPlateauOffset,...
            useElasticTendonExp,...
            flag_useOctave);
          

[humanSoleusActiveForceLengthData,...
humanSoleusPassiveForceLengthData] = getFelineSoleusMusculotendonData(...
                              felineSoleusMusculotendonPropertiesExp,...
                              felineSoleusMusculotendonProperties,...
                              normPlateauOffset,...
                              useElasticTendonExp,...
                              flag_useOctave);


  

