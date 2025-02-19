%%
% SPDX-FileCopyrightText: 2023 Matthew Millard <millard.matthew@gmail.com>
%
% SPDX-License-Identifier: MIT
%
%%
clc;
close all;
clear all;

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
fitOptimalLengthToData=1;
fitPassiveForceLengthRelationToData=0;

indexReferenceDataSet=3;
% 1. Tomalka et al. 2017
% 2. Zuurbier et al. 1995
% 3. Stephenson & Williams 1982

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
% Plot configuration and data structs
%%

[plotDataConfig,...
 plotIndexes,... 
 plotSettings] = getRatSoleusModelPlottingStructures();


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
plotLayoutSettings = struct('numberOfHorizontalPlotColumns',  2,...
                            'numberOfVerticalPlotRows',       1,...
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


normFiberLengthAtZeroPassiveForce     = 0.8;%1.05;
normFiberLengthAtOneNormPassiveForce  = 1.9;

% Since these experiments use skinned fibers, the ECM is assumed to contribute
% nothing. 
ecmForceFractionRatSoleusFitted = 0.0;% 

%
% default value
%
normPevkToActinAttachmentPointRatSoleus=0;


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



%
% Load the reference data
%
[ratSoleusData, ratSoleusMetaData] = loadRatSoleusData(projectFolders);


%
% Select the reference data set
%

activeForceLengthData = [];
passiveForceLengthData = [];

if(fitOptimalLengthToData==1)
    referenceActiveForceLengthDataTable = ...
        ratSoleusData(indexReferenceDataSet).activeForceLengthData;
    for i=1:1:length(referenceActiveForceLengthDataTable)
        activeForceLengthData = ...
            [activeForceLengthData;...
             referenceActiveForceLengthDataTable(i).x, ...
             referenceActiveForceLengthDataTable(i).y];
    end
end

if(fitPassiveForceLengthRelationToData==1)
    referencePassiveForceLengthDataTable = ...
        ratSoleusData(indexReferenceDataSet).passiveForceLengthData;
    for i=1:1:length(referencePassiveForceLengthDataTable)
        passiveForceLengthData = ...
            [passiveForceLengthData;...
             referencePassiveForceLengthDataTable(i).x, ...
             referencePassiveForceLengthDataTable(i).y];            
    end
end   

%
%
% Rat soleus model with the titin-actin bond at the IgP-PEVK border (N2A)
%
%
ratSoleusFibrilActiveTitin = createRatSoleusModel(...
                      normCrossBridgeStiffness,...
                      normCrossBridgeDamping,...
                      normPevkToActinAttachmentPointRatSoleus,...
                      normMaxActiveTitinToActinDamping,...
                      normFiberLengthAtZeroPassiveForce,...
                      normFiberLengthAtOneNormPassiveForce,...
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
                      activeForceLengthData,...
                      passiveForceLengthData,...
                      projectFolders,...
                      flag_useOctave);

if(useWlcTitinModel==1)
    filePathRatSoleus = fullfile(projectFolders.output_structs_FittedModels,...
                                'ratSoleusFibrilActiveTitinWLC.mat');
    save(filePathRatSoleus,'ratSoleusFibrilActiveTitin');
else
    filePathRatSoleus = fullfile(projectFolders.output_structs_FittedModels,...
                                'ratSoleusFibrilActiveTitinLinearTitin.mat');
    save(filePathRatSoleus,'ratSoleusFibrilActiveTitin');
end

%
% Sample the experimental datda
%

%
% sample Tomalka et al. 2017
%

idx =  ratSoleusMetaData.index_TRSS2017;

plotDataConfig(plotIndexes.TRSS2017_fl).x =...
    ratSoleusData(idx).activeForceLengthData.x;
plotDataConfig(plotIndexes.TRSS2017_fl).y =...
    ratSoleusData(idx).activeForceLengthData.y;

plotDataConfig(plotIndexes.TRSS2017_fpe).x =...
    ratSoleusData(idx).passiveForceLengthData.x;
plotDataConfig(plotIndexes.TRSS2017_fpe).y =...
    ratSoleusData(idx).passiveForceLengthData.y;

%
% sample Zuurbier et al.
%

idx =  ratSoleusMetaData.index_ZHGL1995;

yNorm = 1/100;
xData=[];
yData=[];
for i=1:1:length(ratSoleusData(idx).activeForceLengthData)
    ratSoleusData(idx).activeForceLengthData(i).y = ...
        ratSoleusData(idx).activeForceLengthData(i).y.*yNorm;
    xData = [xData;ratSoleusData(idx).activeForceLengthData(i).x];
    yData = [yData;ratSoleusData(idx).activeForceLengthData(i).y];    
end

ratSoleusData(idx).passiveForceLengthData = [];
plotDataConfig(plotIndexes.ZHGL1995_fl).x = xData;
plotDataConfig(plotIndexes.ZHGL1995_fl).y = yData;


%
% sample Stephenson & Williams 
%

idx = ratSoleusMetaData.index_SW1982;

xData_fl=[];
yData_fl=[];
xData_fpe=[];
yData_fpe=[];
for i=1:1:length(ratSoleusData(idx).activeForceLengthData)
    xData_fl = [xData_fl;ratSoleusData(idx).activeForceLengthData(i).x];
    yData_fl = [yData_fl;ratSoleusData(idx).activeForceLengthData(i).y];       
end
for i=1:1:length(ratSoleusData(idx).passiveForceLengthData)
    xData_fpe= [xData_fpe;ratSoleusData(idx).passiveForceLengthData(i).x];
    yData_fpe= [yData_fpe;ratSoleusData(idx).passiveForceLengthData(i).y]; 
end

plotDataConfig(plotIndexes.SW1982_fl).x = xData_fl;
plotDataConfig(plotIndexes.SW1982_fl).y = yData_fl;
plotDataConfig(plotIndexes.SW1982_fpe).x = xData_fpe;
plotDataConfig(plotIndexes.SW1982_fpe).y = yData_fpe;



%
% Sample the model curves
%

% 
% fl
%
activeForceLengthCurveData   = ...
    calcBezierYFcnXCurveSampleVector( ...
        ratSoleusFibrilActiveTitin.curves.activeForceLengthCurve,...
        100,[]);
      
lsOpt = ratSoleusFibrilActiveTitin.sarcomere.optimalSarcomereLength;

plotDataConfig(plotIndexes.model_fl).x = activeForceLengthCurveData.x.*lsOpt;
plotDataConfig(plotIndexes.model_fl).y = activeForceLengthCurveData.y;

%
% fpe and titin detail
%
passiveForceLengthCurveData   = ...
    calcBezierYFcnXCurveSampleVector( ...
        ratSoleusFibrilActiveTitin.curves.fiberForceLengthCurve, ...
        100,[]);

plotDataConfig(plotIndexes.model_fpe).x = passiveForceLengthCurveData.x.*lsOpt;
plotDataConfig(plotIndexes.model_fpe).y = passiveForceLengthCurveData.y;

titinCurveSample = ...
  sampleTitinCurves20250217(...
    ratSoleusFibrilActiveTitin.curves,...
    ratSoleusFibrilActiveTitin.sarcomere,...
    100);




lambdaECM = ratSoleusFibrilActiveTitin.sarcomere.extraCellularMatrixPassiveForceFraction;


plotDataConfig(plotIndexes.model_titinPassive).x = ...
    titinCurveSample.curveSampleTitin.x.*(2*lsOpt);
plotDataConfig(plotIndexes.model_titinPassive).y = ...
    titinCurveSample.curveSampleTitin.y.*(1-lambdaECM) ...
   +titinCurveSample.curveSampleECMHalf.y.*(lambdaECM);

plotDataConfig(plotIndexes.model_titinActive).x = ...
    titinCurveSample.curveSampleTitinActive.x.*(2*lsOpt);
plotDataConfig(plotIndexes.model_titinActive).y = ...
    titinCurveSample.curveSampleTitinActive.y.*(1-lambdaECM) ...
   +titinCurveSample.curveSampleECMHalf.y.*(lambdaECM);



%
% fv
%

forceVelocityCurveData   = ...
    calcBezierYFcnXCurveSampleVector( ...
        ratSoleusFibrilActiveTitin.curves.fiberForceVelocityCurve, ...
        100,[]);

scaleVmax = ...
    ratSoleusFibrilActiveTitin.musculotendon.maximumNormalizedFiberVelocity;

plotDataConfig(plotIndexes.model_fv).x = forceVelocityCurveData.x .* scaleVmax;
plotDataConfig(plotIndexes.model_fv).y = forceVelocityCurveData.y;



%
% Plot the overview curves
%


normMaxActiveSarcomereLength = ...
         ratSoleusFibrilActiveTitin.sarcomere.normZLineLength ...
        +2*ratSoleusFibrilActiveTitin.sarcomere.normActinLength...
        +2*ratSoleusFibrilActiveTitin.sarcomere.normMyosinHalfLength;



plotSettings(1).xticks = [];
plotSettings(1).yticks = [];

plotSettings(1).xticks = [...
    ratSoleusFibrilActiveTitin.sarcomere.normSarcomereLengthZeroForce,...
    1,...
    normMaxActiveSarcomereLength];

plotSettings(1).xticks =...
    plotSettings(1).xticks...
    .*ratSoleusFibrilActiveTitin.sarcomere.optimalSarcomereLength;

plotSettings(1).xticks =...
    [plotSettings(1).xticks, ...
     max(plotDataConfig(plotIndexes.SW1982_fpe).x)];

plotSettings(1).yticks = [0,1];

plotSettings(2).xticks = [];
plotSettings(2).yticks = [];

plotSettings(2).xticks = [...
    -ratSoleusFibrilActiveTitin.musculotendon.maximumNormalizedFiberVelocity,...
    -0.5*ratSoleusFibrilActiveTitin.musculotendon.maximumNormalizedFiberVelocity,...
    0,...
    ratSoleusFibrilActiveTitin.musculotendon.maximumNormalizedFiberVelocity];

plotSettings(2).yticks = [...
    0.00,...
    ratSoleusFibrilActiveTitin.musculotendon.forceVelocityMultiplierAtHalfMaximumFiberVelocity,...
    1.00,...
    ratSoleusFibrilActiveTitin.musculotendon.forceVelocityMultiplierAtLowEccentricFiberVelocity,...
    ratSoleusFibrilActiveTitin.musculotendon.forceVelocityMultiplierAtMaximumEccentricFiberVelocity];

figModelCurves = figure;
    

for i=1:1:length(plotDataConfig)
    if(plotDataConfig(i).enablePlot==1)
        row=plotDataConfig(i).row;
        col=plotDataConfig(i).col;
        subplot('Position', reshape(subPlotPanel(row,col,:),1,4));
        plot(plotDataConfig(i).x,...
             plotDataConfig(i).y,...
             plotDataConfig(i).Mark,...
             'Color',plotDataConfig(i).LineColor,...
             'MarkerFaceColor',plotDataConfig(i).MarkerFaceColor,...
             'MarkerEdgeColor',plotDataConfig(i).MarkerEdgeColor,...
             'MarkerSize',plotDataConfig(i).MarkerSize,...
             'DisplayName',plotDataConfig(i).DisplayName,...
             'HandleVisibility',plotDataConfig(i).HandleVisibility);
        hold on;
    end
end

for i=1:1:length(plotSettings)
    figure(figModelCurves);
    row = plotSettings(i).row;
    col = plotSettings(i).col;
    subplot('Position', reshape(subPlotPanel(row,col,:),1,4));
    xlim(plotSettings(i).xlim);
    ylim(plotSettings(i).ylim);
    if(isempty(plotSettings(i).xticks)==0)
        xticks(round(plotSettings(i).xticks,2));
    end
    if(isempty(plotSettings(i).yticks)==0)
        yticks(round(plotSettings(i).yticks,2));
    end

    legend('Location',plotSettings(i).legendLocation);
    legend('boxoff');
    box off;
    xlabel(plotSettings(i).xlabel);
    ylabel(plotSettings(i).ylabel);            
    title(plotSettings(i).title);            
end


figure(figModelCurves);
configPlotExporter;
filePath = fullfile(projectFolders.output_plots_MuscleCurves,...
                    'fig_Pub_MusleCurves_RatSoleus.pdf');
print('-dpdf', filePath); 