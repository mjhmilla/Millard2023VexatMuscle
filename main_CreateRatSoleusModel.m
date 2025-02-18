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
indexReferenceDataSet=1;
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

[dataToPlot,...
 dataIndexes,... 
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
passiveForceKeyPoints = [1.05,   0;...
                         1.6,1.0];  

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

%
% Load the reference data
%
for i=1:1:3
    referenceData(i).activeForceLengthData=[];
    referenceData(i).passiveForceLengthData=[];
end
%
% Tomalka et al. 2017
%

fileActiveForceLength_TRSS2017 = ...
    fullfile(projectFolders.experiments_TRSS2017,...
             'data','fig_TomalkaRodeSchumacherSiebert2017_Fig2_fl.csv');

filePassiveForceLength_TRSS2017 = ...
    fullfile(projectFolders.experiments_TRSS2017,...
             'data','fig_TomalkaRodeSchumacherSiebert2017_Fig2_fpe.csv');

referenceData(1).activeForceLengthData =...
        loadDigitizedData( ...
        fileActiveForceLength_TRSS2017,...
        'Length (um)','Norm. Force (f/fo)',...
        {'Exp.'},...
         'Rat fibril (EDL) $$f^{PE}$$');

referenceData(1).passiveForceLengthData =...
        loadDigitizedData( ...
        filePassiveForceLength_TRSS2017,...
        'Length (um)','Norm. Force (f/fo)',...
        {'Exp.'},...
         'Rat fibril (EDL) $$f^{L}$$'); 

dataToPlot(dataIndexes.TRSS2017_fl).x =...
    referenceData(1).activeForceLengthData.x;
dataToPlot(dataIndexes.TRSS2017_fl).y =...
    referenceData(1).activeForceLengthData.y;

dataToPlot(dataIndexes.TRSS2017_fpe).x =...
    referenceData(1).passiveForceLengthData.x;
dataToPlot(dataIndexes.TRSS2017_fpe).y =...
    referenceData(1).passiveForceLengthData.y;

%
% Zuurbier et al.
%
fileActiveForceLength_ZHGL1995 = ...
    fullfile(projectFolders.experiments_ZHGL1995,...
    'data','ZuurbierHeslingaGrootLaarse1995.csv'     );
        
referenceData(2).activeForceLengthData =...
            loadDigitizedData( ...
            fileActiveForceLength_ZHGL1995,...
            'Length (um)','Norm. Force (f/fo)',...
            {'ZHGL1995 Exp: $$\mu$$',...
             'ZHGL1995 Exp: $$\mu+\sigma$$',...
             'ZHGL1995 Exp: $$\mu-\sigma$$',...
             'Model (GM)','Model (EDL)'},...
             'Rat fibril (GM + EDL) $$f^{L}$$');


yNorm = 1/100;
xData=[];
yData=[];
for i=1:1:length(referenceData(2).activeForceLengthData)
    referenceData(2).activeForceLengthData(i).y = ...
        referenceData(2).activeForceLengthData(i).y.*yNorm;
    xData = [xData;referenceData(2).activeForceLengthData(i).x];
    yData = [yData;referenceData(2).activeForceLengthData(i).y];    
end

referenceData(2).passiveForceLengthData = [];

dataToPlot(dataIndexes.ZHGL1995_fl).x = xData;
dataToPlot(dataIndexes.ZHGL1995_fl).y = yData;


%
% Stephenson & Williams 
%
fileActiveForceLength_SW1982 = ...
    fullfile(projectFolders.experiments_SW1982,...
    'data','StephensonWilliams1982_Fig7_fl.csv'     );

filePassiveForceLength_SW1982 = ...
    fullfile(projectFolders.experiments_SW1982,...
    'data','StephensonWilliams1982_Fig7_fpe.csv'     );

referenceData(3).activeForceLengthData =...
    loadDigitizedData(fileActiveForceLength_SW1982,...
        'Length (um)','Norm. Force (f/fo)',...
        {'EDL (low temp)','EDL (norm. temp)',...
         'SOL (low temp)','SOL (norm. temp)'},...
         'Rat fibril (SOL + EDL) $$f^{L}$$');

referenceData(3).passiveForceLengthData = ...
    loadDigitizedData(filePassiveForceLength_SW1982,...
        'Length (um)','Norm. Force (f/fo)',...
        {'EDL and SOL'},...
        'Rat fibril (SOL + EDL) $$f^{PE}$$');

xData_fl=[];
yData_fl=[];
xData_fpe=[];
yData_fpe=[];
for i=1:1:length(referenceData(3).activeForceLengthData)
    xData_fl = [xData_fl;referenceData(3).activeForceLengthData(i).x];
    yData_fl = [yData_fl;referenceData(3).activeForceLengthData(i).y];       
end
for i=1:1:length(referenceData(3).passiveForceLengthData)
    xData_fpe= [xData_fpe;referenceData(3).passiveForceLengthData(i).x];
    yData_fpe= [yData_fpe;referenceData(3).passiveForceLengthData(i).y]; 
end

dataToPlot(dataIndexes.SW1982_fl).x = xData_fl;
dataToPlot(dataIndexes.SW1982_fl).y = yData_fl;
dataToPlot(dataIndexes.SW1982_fpe).x = xData_fpe;
dataToPlot(dataIndexes.SW1982_fpe).y = yData_fpe;


%
% Select the reference data set
%

referenceActiveForceLengthDataTable = ...
    referenceData(indexReferenceDataSet).activeForceLengthData;

referencePassiveForceLengthDataTable = ...
    referenceData(indexReferenceDataSet).passiveForceLengthData;
    
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
                      referenceActiveForceLengthDataTable,...
                      referencePassiveForceLengthDataTable,...
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

%
% Sample the curve data
%

% 
% fl
%
activeForceLengthCurveData   = ...
    calcBezierYFcnXCurveSampleVector( ...
        ratSoleusFibril.curves.activeForceLengthCurve,...
        100,[]);
      
lsOpt = ratSoleusFibril.sarcomere.optimalSarcomereLength;

dataToPlot(dataIndexes.model_fl).x = activeForceLengthCurveData.x.*lsOpt;
dataToPlot(dataIndexes.model_fl).y = activeForceLengthCurveData.y;

%
% fpe
%
passiveForceLengthCurveData   = ...
    calcBezierYFcnXCurveSampleVector( ...
        ratSoleusFibril.curves.fiberForceLengthCurve, ...
        100,[]);

dataToPlot(dataIndexes.model_fpe).x = passiveForceLengthCurveData.x.*lsOpt;
dataToPlot(dataIndexes.model_fpe).y = passiveForceLengthCurveData.y;

titinCurveSample = ...
  sampleTitinCurves20250217(...
    ratSoleusFibril.curves,...
    ratSoleusFibril.sarcomere,...
    100);

%dataToPlot(dataIndexes.model_fpe).x = passiveForceLengthCurveData.x;
%dataToPlot(dataIndexes.model_fpe).y = passiveForceLengthCurveData.y;

lambdaECM = ratSoleusFibril.sarcomere.extraCellularMatrixPassiveForceFraction;

% dataToPlot(dataIndexes.model_fpe).x = ...
%     titinCurveSample.curveSampleTitin.x.*(2*lsOpt);
% 
% dataToPlot(dataIndexes.model_fpe).y = ...
%     titinCurveSample.curveSampleTitin.y.*(1-lambdaECM) ...
%    +titinCurveSample.curveSampleECMHalf.y.*(lambdaECM);

%
% fv
%

forceVelocityCurveData   = ...
    calcBezierYFcnXCurveSampleVector( ...
        ratSoleusFibril.curves.fiberForceVelocityCurve, ...
        100,[]);

dataToPlot(dataIndexes.model_fv).x = forceVelocityCurveData.x;
dataToPlot(dataIndexes.model_fv).y = forceVelocityCurveData.y;

%
% Titin detail
%


dataToPlot(dataIndexes.model_titinPassive).x = ...
    titinCurveSample.curveSampleTitin.x.*(2*lsOpt);
dataToPlot(dataIndexes.model_titinPassive).y = ...
    titinCurveSample.curveSampleTitin.y.*(1-lambdaECM) ...
   +titinCurveSample.curveSampleECMHalf.y.*(lambdaECM);

dataToPlot(dataIndexes.model_titinActive).x = ...
    titinCurveSample.curveSampleTitinActive.x.*(2*lsOpt);
dataToPlot(dataIndexes.model_titinActive).y = ...
    titinCurveSample.curveSampleTitinActive.y.*(1-lambdaECM) ...
   +titinCurveSample.curveSampleECMHalf.y.*(lambdaECM);



%
% Plot the overview curves
%


figModelCurves = figure;
    

for i=1:1:length(dataToPlot)
    row=dataToPlot(i).row;
    col=dataToPlot(i).col;
    subplot('Position', reshape(subPlotPanel(row,col,:),1,4));
    plot(dataToPlot(i).x,...
         dataToPlot(i).y,...
         dataToPlot(i).Mark,...
         'Color',dataToPlot(i).LineColor,...
         'MarkerFaceColor',dataToPlot(i).MarkerFaceColor,...
         'MarkerEdgeColor',dataToPlot(i).MarkerEdgeColor,...
         'MarkerSize',dataToPlot(i).MarkerSize,...
         'DisplayName',dataToPlot(i).DisplayName,...
         'HandleVisibility',dataToPlot(i).HandleVisibility);
    hold on;
end

for i=1:1:length(plotSettings)
    figure(figModelCurves);
    row = plotSettings(i).row;
    col = plotSettings(i).col;
    subplot('Position', reshape(subPlotPanel(row,col,:),1,4));
    xlim(plotSettings(i).xlim);
    ylim(plotSettings(i).ylim);   
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