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

% A simulation of versions of the experiments of Gordon et al., and Hill.
%
% Gordon, A. M., Huxley, A. F., & Julian, F. J. (1966). The variation 
% in isometric tension with sarcomere length in vertebrate muscle 
% fibres. The Journal of physiology, 184(1), 170-192.
%
% Hill AV. The heat of shortening and the dynamic constants of muscle. 
% Proceedings of the Royal Society of London. Series B-Biological Sciences. 
% 1938 Oct 10;126(843):136-95.

%This flag allows us to avoid the memory clearing functions so that
%this can be timed using tic and tock from within main_OuterLoop
flag_OuterOuterLoopMode =1;
if(flag_OuterOuterLoopMode ==0)
    clc;
    close all;
    clear all;
end

disp('----------------------------------------');
disp(' running main_ClassicForceLengthVelocityExperiments_OuterLoop');
disp('----------------------------------------');
disp('   :run-time: ~26 minutes ');
disp('            *Intel i7-3630QM @ 2.40 GHz, Ubuntu 22');
disp '             8 GB ram, SSD harddrive');
disp('----------------------------------------');

%Parameters that are tuned
% updSlidingTimeConstant              = 0.0025;
% updForceVelocityCalibrationFactor   = 0.95;

%Parameters that do not change
flag_simulateHillModel            = 1;  
flag_simulateVexatModel          = 1;

flag_useSimulatePlotRigidTendon   = 1;
flag_useSimulatePlotElasticTendon = 1;

flag_removeActiveTitinForces      = 0;
flag_useTwoSidedTitinCurves       = 0;

flag_plotDataComparison           = 1;  
flag_plotHillModel                = 0;

flag_useFiberDamping              = 1;
fiberDampingCoefficient           = 0.01;

flag_useFig3KirchBoskovRymer1994  = 0; 
flag_useOctave                    = 0;


flag_activeForceLengthSimulations  = 1;
flag_passiveForceLengthSimulations = 1;  
flag_forceVelocitySimulations      = 1;

normFiberLengthAtForceVelocitySample  = 1.;

maxShorteningVelocity = 4.0;
forceVelocityNormFiberHalfLength = 0.1;

numberOfLengthSteps   = 20;
numberOfVelocitySteps = 10;

%Common parameters



flag_plotData           = 0;
flag_buildCombinedPlot  = 0;  
figClassicFl = figure;
figClassicFv = figure; 

if(flag_plotDataComparison==1)
  flag_plotData         = 1;
  flag_savePlotsToFile  = 1;
end

if(flag_useSimulatePlotRigidTendon==1)
    if(flag_plotDataComparison==1)
      flag_buildCombinedPlot  = 1;  
    end
    flag_useElasticTendon=0;
    main_ClassicForceLengthVelocityExperiments;
end

if(flag_useSimulatePlotElasticTendon==1)
    if(flag_plotDataComparison==1)
      flag_buildCombinedPlot  = 2; 
    end   
    flag_useElasticTendon   = 1;   
    main_ClassicForceLengthVelocityExperiments;
end