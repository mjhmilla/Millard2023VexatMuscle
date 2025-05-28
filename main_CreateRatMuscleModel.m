%%
% SPDX-FileCopyrightText: 2025 Matthew Millard <millard.matthew@gmail.com>
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
validMuscles = {'SOL','EDL'};
muscleName = validMuscles{2};

%Least-squares fit of the optimal CE length and the passive-force-length
%relation
indexOfDataSetToFitOptCELength          =0;
indexOfDataSetToPassiveForceLengthCurve =0;

%Sets the initial stiffness and length at which the curve develops
%one norm force
indexOfDataSetForPassiveCurveParameters = 3;
% 1. Tomalka et al. 2017
% 2. Zuurbier et al. 1995
% 3. Stephenson & Williams 1982

expDataSetFittingData(3)=struct('optimalSarcomereLength',0,...
                               'minLengthWhereFpeIsLinear',0);

expDataSetFittingData(1).optimalSarcomereLength=2.525;
expDataSetFittingData(2).optimalSarcomereLength=2.525;
expDataSetFittingData(3).optimalSarcomereLength=2.525;

expDataSetFittingData(1).minLengthWhereFpeIsLinear=nan;
expDataSetFittingData(2).minLengthWhereFpeIsLinear=nan;
expDataSetFittingData(3).minLengthWhereFpeIsLinear=0.3;



%%
% Load the reference data
%%

[ratMuscleData, ratMuscleMetaData] = ...
        loadRatSkeletalMuscleData(projectFolders);
%
% Select the reference data set
%

activeForceLengthData = [];
passiveForceLengthData = [];

if(indexOfDataSetToFitOptCELength > 0)
    referenceActiveForceLengthDataTable = ...
        ratMuscleData(indexOfDataSetToFitOptCELength).activeForceLengthData;
    for i=1:1:length(referenceActiveForceLengthDataTable)
        activeForceLengthData = ...
            [activeForceLengthData;...
             referenceActiveForceLengthDataTable(i).x, ...
             referenceActiveForceLengthDataTable(i).y];
    end
end

if(indexOfDataSetToPassiveForceLengthCurve>0)
    referencePassiveForceLengthDataTable = ...
        ratMuscleData(indexOfDataSetToPassiveForceLengthCurve).passiveForceLengthData;
    for i=1:1:length(referencePassiveForceLengthDataTable)
        passiveForceLengthData = ...
            [passiveForceLengthData;...
             referencePassiveForceLengthDataTable(i).x, ...
             referencePassiveForceLengthDataTable(i).y];            
    end
end  

%%
% Manually set sarcomere properties 
%%
makeSkinnedFibrilModel = 1;

setSarcomereProperties.normCrossBridgeStiffness         = 75;%49.1;  %fiso/lopt
setSarcomereProperties.normCrossBridgeDamping           = 0.347*(75/49.1); %fiso/(lopt/s)
% Manually set to produce the force transient in TRSS2017 Figure 2A
% at the beginning of lengthening


setSarcomereProperties.useWLCTitinModel                 = 0;
% 0: Using a linear titin model
% 1: Using the WLC titin model
% As described in Millard, Franklin, Herzog 2024

setSarcomereProperties.titinMolecularWeightInkD         = []; 
% I haven't been able to find a report of the molecular weight of titin
% in rat skeletal muscle.

setSarcomereProperties.ecmForceFraction                 = (1-makeSkinnedFibrilModel)*0.56;% 
% No ECM for a skinned fibril
% Otherwise use 0.56, the average amount of ECM reported in Prado et al.
% across 5 rabbit skeletal muscles

setSarcomereProperties.normPevkToActinAttachmentPoint   = 0.5; 
%0 : Prox-Ig/PEVK boundary
%1 : PEVK/Dist-Ig boundary

setSarcomereProperties.normMaxActiveTitinToActinDamping = 20.3; %fo/(lo/s)
% setSarcomereProperties.normMaxActiveTitinToActinDamping
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
%
% Tomalka A, Weidner S, Hahn D, Seiberl W, Siebert T. Power amplification 
% increases with contraction velocity during stretch-shortening cycles of 
% skinned muscle fibers. Frontiers in physiology. 2021 Mar 31;12:644981.
%
% Herzog W, Leonard TR. Force enhancement following stretching of skeletal 
% muscle: a new mechanism. Journal of Experimental Biology. 2002 
% May 1;205(9):1275-83.

%%
% Manually set musculotendon properties
%%

setMusculotendonProperties.normFiberLengthAtZeroPassiveForce        = 0.6;
setMusculotendonProperties.normFiberLengthAtOneNormPassiveForce     = 1.9;
setMusculotendonProperties.normFiberStiffnessAtOneNormPassiveForce  = nan;
setMusculotendonProperties.scaleOptimalFiberLength                  = 1;
setMusculotendonProperties.scaleMaximumIsometricTension             = 1;
setMusculotendonProperties.makeSkinnedFibrilModel = makeSkinnedFibrilModel;
setMusculotendonProperties.useElasticTendon                         = ...
        1 && ~setMusculotendonProperties.makeSkinnedFibrilModel;


%%
% Manually set curves
%%


setCurveProperties.useCalibratedCurves                              = 1;
setCurveProperties.useTwoSidedTitinCurves                           = 0;
setCurveProperties.smallNumericallyNonZeroNumber      = sqrt(sqrt(eps));
setCurveProperties.flag_enableNumericallyNonZeroGradients           = 1;
setCurveProperties.useForceVelocityCurveWithSlopeDiscontinuity      = 1;
setCurveProperties.fitTitinToTRSS2017Data         =  ratMuscleData(1);
setCurveProperties.activeForceLengthData          =  [];
setCurveProperties.passiveForceLengthData         =  [];


specimenTemperature     = 12; %As in 12 degrees centrigrade

flag_useOctave                          = 0;



%%
% Plot configuration and data structs
%%

[plotDataConfig,...
 plotIndexes,... 
 plotSettings] = getRatMusculotendonModelPlottingStructures(muscleName);


%%
% Rat soleus fibril Model
%%

fprintf('\n\nCreating: rat EDL model \n');
fprintf('  used to simulate Tomalka, Rode, Schumacher, Siebert 2017.\n\n');

fprintf('\n\nTo do:');
fprintf('\n\n1. Look at Prado again: there are a lot of references related to');
fprintf(  '\n   rat muscle.');
fprintf(  '\n');


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

 
%%
% Fit the passive force-length relation
%   Even for fibril models this is useful as it serves as a fitting
%   reference for the passive-force-length relation for the entire 
%   titin segment.
%%

if(indexOfDataSetForPassiveCurveParameters>0)
    
    referencePassiveForceLengthDataTable = ...
        ratMuscleData(indexOfDataSetForPassiveCurveParameters).passiveForceLengthData;

    fittingFpeNMinForce = ...
            expDataSetFittingData(...
            indexOfDataSetForPassiveCurveParameters).minLengthWhereFpeIsLinear;

    fittingFpeNOptSarcomereLengthInUm = ...
            expDataSetFittingData(...
            indexOfDataSetForPassiveCurveParameters).optimalSarcomereLength; 

    fittingDataFpeN=[];
    for i=1:1:length(referencePassiveForceLengthDataTable)
        for j=1:1:length(referencePassiveForceLengthDataTable(i).x)
            if(referencePassiveForceLengthDataTable(i).y(j) > fittingFpeNMinForce)
                fittingDataFpeN = ...
                    [fittingDataFpeN;...
                     referencePassiveForceLengthDataTable(i).x(j), ...
                     referencePassiveForceLengthDataTable(i).y(j)];            
            end
        end
    end
    
    %Fit a line to the data
    A = [fittingDataFpeN(:,1),ones(size(fittingDataFpeN(:,1)))];
    b = fittingDataFpeN(:,2);
    x = (A'*A)\(A'*b);
    c = x(1,1);
    y0 = x(2,1);
    %y = c*x + x0
    
    normFiberStiffnessAtOneNormPassiveForce = ...
        c/(1/fittingFpeNOptSarcomereLengthInUm);
    normFiberLengthAtOneNormPassiveForce = ...
        ((1-y0)/c)/fittingFpeNOptSarcomereLengthInUm;

    here=1;

end


ratMuscleModelParameters = createRatSkeletalMuscleModel(...
                              setSarcomereProperties,...
                              setMusculotendonProperties,...
                              setCurveProperties,...
                              specimenTemperature,...                   
                              muscleName,...
                              projectFolders,...
                              flag_useOctave);

wlcStr = '';
if(setSarcomereProperties.useWLCTitinModel==1)
    wlcStr='WLC';
end
fibrilStr='';
if(setMusculotendonProperties.makeSkinnedFibrilModel==1)
    fibrilStr='Fibril';
end

fileName = ['rat',muscleName,fibrilStr,'ActiveTitin',wlcStr,'.mat'];
filePathRatMuscle = fullfile(projectFolders.output_structs_FittedModels,...
                             fileName);
save(filePathRatMuscle,'ratMuscleModelParameters');


%
% Sample the experimental datda
%

idx =  ratMuscleMetaData.index_TRSS2017;

plotDataConfig(plotIndexes.TRSS2017_fl).x =...
    ratMuscleData(idx).activeForceLengthData.x;
plotDataConfig(plotIndexes.TRSS2017_fl).y =...
    ratMuscleData(idx).activeForceLengthData.y;

plotDataConfig(plotIndexes.TRSS2017_fpe).x =...
    ratMuscleData(idx).passiveForceLengthData.x;
plotDataConfig(plotIndexes.TRSS2017_fpe).y =...
    ratMuscleData(idx).passiveForceLengthData.y;

%
% sample Zuurbier et al.
%

idx =  ratMuscleMetaData.index_ZHGL1995;

yNorm = 1/100;
xData=[];
yData=[];
for i=1:1:length(ratMuscleData(idx).activeForceLengthData)
    ratMuscleData(idx).activeForceLengthData(i).y = ...
        ratMuscleData(idx).activeForceLengthData(i).y.*yNorm;
    xData = [xData;ratMuscleData(idx).activeForceLengthData(i).x];
    yData = [yData;ratMuscleData(idx).activeForceLengthData(i).y];    
end

ratMuscleData(idx).passiveForceLengthData = [];
plotDataConfig(plotIndexes.ZHGL1995_fl).x = xData;
plotDataConfig(plotIndexes.ZHGL1995_fl).y = yData;


%
% sample Stephenson & Williams 
%

idx = ratMuscleMetaData.index_SW1982;

xData_fl=[];
yData_fl=[];
xData_fpe=[];
yData_fpe=[];
for i=1:1:length(ratMuscleData(idx).activeForceLengthData)
    xData_fl = [xData_fl;ratMuscleData(idx).activeForceLengthData(i).x];
    yData_fl = [yData_fl;ratMuscleData(idx).activeForceLengthData(i).y];       
end
for i=1:1:length(ratMuscleData(idx).passiveForceLengthData)
    xData_fpe= [xData_fpe;ratMuscleData(idx).passiveForceLengthData(i).x];
    yData_fpe= [yData_fpe;ratMuscleData(idx).passiveForceLengthData(i).y]; 
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
        ratMuscleModelParameters.curves.activeForceLengthCurve,...
        100,[]);
      
lsOpt = ratMuscleModelParameters.sarcomere.optimalSarcomereLength;

plotDataConfig(plotIndexes.model_fl).x = activeForceLengthCurveData.x.*lsOpt;
plotDataConfig(plotIndexes.model_fl).y = activeForceLengthCurveData.y;

%
% fpe and titin detail
%
passiveForceLengthCurveData   = ...
    calcBezierYFcnXCurveSampleVector( ...
        ratMuscleModelParameters.curves.fiberForceLengthCurve, ...
        100,[]);

plotDataConfig(plotIndexes.model_fpe).x = passiveForceLengthCurveData.x.*lsOpt;
plotDataConfig(plotIndexes.model_fpe).y = passiveForceLengthCurveData.y;

titinCurveSample = ...
  sampleTitinCurves20250217(...
    ratMuscleModelParameters.curves,...
    ratMuscleModelParameters.sarcomere,...
    100);




lambdaECM = ratMuscleModelParameters.sarcomere.extraCellularMatrixPassiveForceFraction;


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
        ratMuscleModelParameters.curves.fiberForceVelocityCurve, ...
        100,[]);

scaleVmax = ...
    ratMuscleModelParameters.musculotendon.maximumNormalizedFiberVelocity;

plotDataConfig(plotIndexes.model_fv).x = forceVelocityCurveData.x .* scaleVmax;
plotDataConfig(plotIndexes.model_fv).y = forceVelocityCurveData.y;



%
% Plot the overview curves
%

%
% Force-length ticks
%
maxActiveSarcomereLengthInUm = ...
         2*ratMuscleModelParameters.sarcomere.zLineLengthInUm ...
        +2*ratMuscleModelParameters.sarcomere.actinLengthInUm...
        +ratMuscleModelParameters.sarcomere.myosinLengthInUm;

optimalSarcomereLengthInUm = ...
         2*ratMuscleModelParameters.sarcomere.zLineLengthInUm ...
        +2*ratMuscleModelParameters.sarcomere.actinLengthInUm...
        +ratMuscleModelParameters.sarcomere.myosinBareLengthInUm;

shortSarcomereLengthInUm = ...
         2*ratMuscleModelParameters.sarcomere.zLineLengthInUm ...
        +ratMuscleModelParameters.sarcomere.myosinLengthInUm;

zeroForceSarcomereLengthInUm = ...
    ratMuscleModelParameters.sarcomere.zeroForceSarcomereLengthInUm;


plotSettings(1).xticks = [];
plotSettings(1).yticks = [];

plotSettings(1).xticks = [...
    zeroForceSarcomereLengthInUm,...
    shortSarcomereLengthInUm,...
    optimalSarcomereLengthInUm,...
    maxActiveSarcomereLengthInUm];



plotSettings(1).xticks =...
    [plotSettings(1).xticks, ...
     max(plotDataConfig(plotIndexes.SW1982_fpe).x)];

%Evaluate the max. isometric force at the length when the myosin tip
%touches the z-line. These expressions require a diagram to understand 
%but I've justed checked them and these expressions are correct.

%The length at which the two actin filaments overlap with
%the active part of  myosin
shallowPlateauInterference = ...
    2*(      ratMuscleModelParameters.sarcomere.actinLengthInUm ...
       - 0.5*ratMuscleModelParameters.sarcomere.myosinLengthInUm ...
       - 0.5*ratMuscleModelParameters.sarcomere.myosinBareLengthInUm);

%The length at which the actin filaments can interact with the 
%active part of myosin without over lap
shallowPlateauOverlap      = ...
    2*(ratMuscleModelParameters.sarcomere.myosinLengthInUm ...
     - ratMuscleModelParameters.sarcomere.actinLengthInUm);

%The maximum possible length at which myosin and actin can
% actively interact 
maxOverlap                 = ...
    ratMuscleModelParameters.sarcomere.myosinLengthInUm ...
  - ratMuscleModelParameters.sarcomere.myosinBareLengthInUm;

%With half of the cross-bridges pulling in one direction and the other half 
%pulling in the opposite direction the interference section contributes no force

interferenceTension     = 0.0; 

maxNormForceAtShortLength =...
    ( shallowPlateauInterference*interferenceTension ...
      + shallowPlateauOverlap )/maxOverlap;

plotSettings(1).yticks = [0,maxNormForceAtShortLength,1];

if(indexOfDataSetForPassiveCurveParameters>0)
    fittingFpeNMinForce = ...
        expDataSetFittingData(...
        indexOfDataSetForPassiveCurveParameters).minLengthWhereFpeIsLinear;
    
    plotSettings(1).yticks = ...
        [0,fittingFpeNMinForce,maxNormForceAtShortLength,1];
    plotSettings(1).yticks = sort(plotSettings(1).yticks);
end




%
% Force-velocity ticks
%

plotSettings(2).xticks = [];
plotSettings(2).yticks = [];

plotSettings(2).xticks = [...
    -ratMuscleModelParameters.musculotendon.maximumNormalizedFiberVelocity,...
    -0.5*ratMuscleModelParameters.musculotendon.maximumNormalizedFiberVelocity,...
    0,...
    ratMuscleModelParameters.musculotendon.maximumNormalizedFiberVelocity];

plotSettings(2).yticks = [...
    0.00,...
    ratMuscleModelParameters.musculotendon.forceVelocityMultiplierAtHalfMaximumFiberVelocity,...
    1.00,...
    ratMuscleModelParameters.musculotendon.forceVelocityMultiplierAtLowEccentricFiberVelocity,...
    ratMuscleModelParameters.musculotendon.forceVelocityMultiplierAtMaximumEccentricFiberVelocity];

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

    titleStr = plotSettings(i).title{:};

    title(titleStr);            
end


figure(figModelCurves);

filePath = fullfile(projectFolders.output_plots_MuscleCurves,...
                    ['fig_Pub_MuscleCurves_Rat',muscleName,'.pdf']);
print('-dpdf', filePath); 