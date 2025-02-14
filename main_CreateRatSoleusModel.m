%%
% SPDX-FileCopyrightText: 2023 Matthew Millard <millard.matthew@gmail.com>
%
% SPDX-License-Identifier: MIT
%
%%

rootDir         = getRootProjectDirectory();
projectFolders  = getProjectFolders(rootDir);

addpath( genpath(projectFolders.parameters)     );
addpath( genpath(projectFolders.curves)         );
addpath( genpath(projectFolders.experiments)    );
addpath( genpath(projectFolders.simulation)     );
addpath( genpath(projectFolders.models)         );
addpath( genpath(projectFolders.postprocessing) );

%%
% Parameters
%%

makeFibrilModel         = 1;
useElasticTendon        = 1 && ~makeFibrilModel;
useWlcTitinModel        = 0;
useCalibratedCurves     = 0;
useTwoSidedTitinCurves  = 0;

flag_enableNumericallyNonZeroGradients  = 1;
scaleOptimalFiberLengthRatSoleus        = 1;
scaleMaximumIsometricTensionRatSoleus   = 1;
flag_useOctave                          = 0;

%%
% Rat soleus fibril Model
%%

fprintf('\n\nCreating: default rat soleus fibril model\n');
fprintf('  used to simulate Tomalka, Weider, Hahn, Seiberl, Siebert 2020.\n\n');

fprintf('\n\nTo do:');
fprintf('\n1. Update createFiberActiveForceLengthCurve:');
fprintf('\na. Option 1: Add optional parameters to set the coordinates of the');
fprintf('\n             transition from the steep-to-shallow ascending limb.');
fprintf('\nb. Option 2: Add a new function that takes the saromere data + exp');
fprintf('\n             measurements and makes a fit');
fprintf('\nc. Option 3: Add a new function that takes the saromere data + exp');
fprintf('\n             uses Guenter and Rockenfellers model and makes a fit');
fprintf('\nMotivation: Stephenson and Williams present a mean +/- 1sd force-');
fprintf('\n            length curve that is far below the theoretical curve on');
fprintf('\n            the shallow ascending limb.');
fprintf('\n\n2. Look at Prado again: there are a lot of references related to');
fprintf(  '\n   rat muscle.')
% Stephenson DG, Williams DA. Effects of sarcomere length on the force—pCa 
% relation in fast‐and slow‐twitch skinned muscle fibres from the rat. 
% The Journal of Physiology. 1982 Dec 1;333(1):637-53.

fprintf('\n\n3. Find sources for lopt, fiso, and ltslk beyond Lemaire et al.');
fprintf(  '\n   for the rat soleus muscle.')

%%
% Plot configuration
%%
plotLayoutSettings = struct('numberOfHorizontalPlotColumns',  4,...
                            'numberOfVerticalPlotRows',       2,...
                            'flag_fixedPlotWidth',            1,...
                            'plotWidth',                      7,...
                            'plotHeight',                     7,...
                            'flag_usingOctave',               0);

numberOfHorizontalPlotColumns = plotLayoutSettings.numberOfHorizontalPlotColumns;
numberOfVerticalPlotRows      = plotLayoutSettings.numberOfVerticalPlotRows;
flag_fixedPlotWidth           = plotLayoutSettings.flag_fixedPlotWidth;
plotWidth                     = plotLayoutSettings.plotWidth;
plotHeight                    = plotLayoutSettings.plotHeight;
flag_usingOctave              = plotLayoutSettings.flag_usingOctave;

plotHorizMarginCm = 2;
plotVertMarginCm  = 2;

pageWidth   = (plotWidth+plotHorizMarginCm)*numberOfHorizontalPlotColumns...
                +plotHorizMarginCm;

pageHeight  = (plotHeight+plotVertMarginCm)*numberOfVerticalPlotRows...
                +plotVertMarginCm;

plotConfigGeneric;

%%%
%
% XE normalized parameters from the cat soleus with an elastic tendon
% from Table 1 
%
% Matthew Millard, David W Franklin, Walter Herzog (2024) 
% A three filament mechanistic model of musculotendon force and impedance 
% eLife 12:RP88344 https://doi.org/10.7554/eLife.88344.4    
%
%%%
normCrossBridgeStiffness    = 49.1;  %fiso/lopt
normCrossBridgeDamping      = 0.347; %fiso/(lopt/s)


%
% I haven't yet found any papers that report the molecular weight of
% titin from rat soleus titin.
%
titinMolecularWeightInkDDefault =[];

%%%
%
% From Figure 7 of Stephenson & Williams
%
% Using theoretical force-length relation:
%
%   lopt = 2.5um
%   fiso = 1
%
% Using a (manually identified) better fit to the data
%
%   lopt = 2.66 um
%   fiso = 1
%
% The lengths where fpeN = 0 and fpeN = 1 are
%
%   lceNFpeNZero  = 4.24um  / 2.5um = 1.696 lo
%   lceNFpeNOne   = 2.629um / 2.5um = 1.05 lo
%
% Using the manually identified lopt
%
%   lceNFpeNZero  = 4.24um  / 2.66um = 1.59 lo
%   lceNFpeNOne   = 2.629um / 2.66um = 0.988 lo
%
% Sticking with the theoretical force-length model, at least to start.
%
% Stephenson DG, Williams DA. Effects of sarcomere length on the force—pCa 
% relation in fast‐and slow‐twitch skinned muscle fibres from the rat. 
% The Journal of Physiology. 1982 Dec 1;333(1):637-53.
%
%%%
passiveForceKeyPoints = [1.05,   0;...
                         1.696,1.0];  

normFiberLengthAtOneNormPassiveForceRatSoleusFibril = passiveForceKeyPoints(2,1);

% Since these experiments use skinned fibers, the ECM is assumed to contribute
% nothing. 
ecmForceFractionRatSoleusFitted = 0.0;% 

%
% default value
%
normPevkToActinAttachmentPointRatSoleusFitted=0.5;

%
% The half relaxation time in Figure 1 of Tomalka et al. (a stretch of ~20%
% in 0.5 s) is 
%
% thalf = 0.0313 s
%
% The half relaxation time of the cat soleus in Herzog & Leonard 2002
% Figure 7C (a stretch of 21% lopt in 0.333 s) is
% 
% thalf = 0.111 s
%
% Scaling the value of normMaxActiveTitinToActinDamping used to simulate 
% Herzog & Leonard we have
% 
% normMaxActiveTitinToActinDamping = 71.9*(0.0313 / 0.111) 
%                                  = 20.3
%
%
% Tomalka A, Weidner S, Hahn D, Seiberl W, Siebert T. Power amplification 
% increases with contraction velocity during stretch-shortening cycles of 
% skinned muscle fibers. Frontiers in physiology. 2021 Mar 31;12:644981.
%
% Herzog W, Leonard TR. Force enhancement following stretching of skeletal 
% muscle: a new mechanism. Journal of Experimental Biology. 2002 
% May 1;205(9):1275-83.
%
normMaxActiveTitinToActinDamping = 20.3; %fo/(lo/s)


smallNumericallyNonZeroNumber = sqrt(sqrt(eps));


ratSoleusFibril = createRatSoleusModel(...
                              normCrossBridgeStiffness,...
                              normCrossBridgeDamping,...
                              normPevkToActinAttachmentPointRatSoleusFitted,...
                              normMaxActiveTitinToActinDamping,...
                              normFiberLengthAtOneNormPassiveForceRatSoleusFibril,...
                              ecmForceFractionRatSoleusFitted,...
                              titinMolecularWeightInkDDefault,...
                              useWlcTitinModel,...
                              useCalibratedCurves,...
                              useTwoSidedTitinCurves,...
                              smallNumericallyNonZeroNumber,...
                              flag_enableNumericallyNonZeroGradients,...
                              scaleOptimalFiberLengthRatSoleus,...
                              scaleMaximumIsometricTensionRatSoleus, ...
                              useElasticTendon,...
                              makeFibrilModel,...
                              projectFolders,...
                              flag_useOctave);

if(useWlcTitinModel==1)
    filePathRatSoleus = fullfile(projectFolders.output_structs_FittedModels,...
                                'ratSoleusFibrilWLC.mat');
    save(filePathRatSoleus,'ratSoleusFibril');
else
    filePathRatSoleus = fullfile(projectFolders.output_structs_FittedModels,...
                                'ratSoleusFibrilLinearTitin.mat');
    save(filePathRatSoleus,'ratSoleusFibril');
end
