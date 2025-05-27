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
mapToEDLModel = 1;
%Will map the rat soleus model to an EDL model. This option exists
%specifically for a rat edl fiber model: the literature on rat soleus
%fiber properties is much more complete than edl fibers. And so, to make
%an EDL model I use all of the normalized properties of the soleus model
%and change the fiber velocity to be consistent with an EDL (and also
%change the muscle name);

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


makeFibrilModel         = 1;
useElasticTendon        = 1 && ~makeFibrilModel;
useWlcTitinModel        = 0;
useCalibratedCurves     = 0;
useTwoSidedTitinCurves  = 0;

specimenTemperature     = 12; %As in 12 degrees centrigrade

flag_enableNumericallyNonZeroGradients  = 1;
scaleOptimalFiberLengthRatSoleus        = 1;
scaleMaximumIsometricTensionRatSoleus   = 1;
flag_useOctave                          = 0;



%%
% Plot configuration and data structs
%%

[plotDataConfig,...
 plotIndexes,... 
 plotSettings] = getRatSoleusModelPlottingStructures(mapToEDLModel);


%%
% Rat soleus fibril Model
%%

fprintf('\n\nCreating: default rat soleus fibril model\n');
fprintf('  used to simulate Tomalka, Weider, Hahn, Seiberl, Siebert 2020.\n\n');

fprintf('\n\nTo do:');
fprintf('\n1. Make a routine to solve for the curviness in fpeN, fTiPN, fTiDN');
fprintf('\n   that minimizes the squared errors w.r.t. the experimental data.');
fprintf('\n   the fitting method should probably be set explicitly using a flag.');
fprintf('\n\n2. Update the active-force-length relation to fit the data better?');
fprintf('\n   Perhaps use Guenter and Rockenfellers model');
fprintf('\n\n3. Look at Prado again: there are a lot of references related to');
fprintf(  '\n   rat muscle.');
fprintf(  '\n');
% Stephenson DG, Williams DA. Effects of sarcomere length on the force—pCa 
% relation in fast‐and slow‐twitch skinned muscle fibres from the rat. 
% The Journal of Physiology. 1982 Dec 1;333(1):637-53.

fprintf('\n\n4. Find sources for lopt, fiso, and ltslk beyond Lemaire et al.');
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

%
% Load the reference data
%
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
normCrossBridgeStiffness    = 75;%49.1;  %fiso/lopt
normCrossBridgeDamping      = 0.347*(75/49.1); %fiso/(lopt/s)


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


normFiberLengthAtZeroPassiveForce       = 0.6;%1.05;
normFiberLengthAtOneNormPassiveForce    = 1.9;
normFiberStiffnessAtOneNormPassiveForce = nan;



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

% Since these experiments use skinned fibers, the ECM is assumed to contribute
% nothing. 
ecmForceFractionRatSoleusFitted = 0.0;% 

%
% default value
%
normPevkToActinAttachmentPointRatSoleus= 0.5;


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
% Rat soleus model with the titin-actin bond at the IgP-PEVK border (N2A)
%
useSharpForceVelocityCurve          = 1;
%This makes a force-velocity curve that has a slope discontinuity 
%at vce=0 so that there is a sharp difference between the concentric
%and eccentric curves. This is useful to match the transient response of
%muscle during sharp length changes.
%
useModifiedPassiveForceLengthCurve  = 1;

ratMuscleModelParameters = createRatSoleusModel(...
                      normCrossBridgeStiffness,...
                      normCrossBridgeDamping,...
                      normPevkToActinAttachmentPointRatSoleus,...
                      normMaxActiveTitinToActinDamping,...
                      normFiberLengthAtZeroPassiveForce,...                      
                      normFiberLengthAtOneNormPassiveForce,...
                      normFiberStiffnessAtOneNormPassiveForce,...                      
                      ecmForceFractionRatSoleusFitted,...
                      titinMolecularWeightInkDDefault,...
                      useWlcTitinModel,...
                      useCalibratedCurves,...
                      useTwoSidedTitinCurves,...
                      smallNumericallyNonZeroNumber,...
                      flag_enableNumericallyNonZeroGradients,...
                      scaleOptimalFiberLengthRatSoleus,...
                      scaleMaximumIsometricTensionRatSoleus, ...
                      specimenTemperature,...
                      useSharpForceVelocityCurve,...
                      useModifiedPassiveForceLengthCurve,...                      
                      useElasticTendon,...
                      makeFibrilModel,...
                      activeForceLengthData,...
                      passiveForceLengthData,...
                      mapToEDLModel,...
                      projectFolders,...
                      flag_useOctave);

wlcStr = '';
if(useWlcTitinModel==1)
    wlcStr='WLC';
end
fibrilStr='';
if(makeFibrilModel==1)
    fibrilStr='Fibril';
end
muscleStr='Soleus';
if(mapToEDLModel==1)
    muscleStr='EDL';
end

fileName = ['rat',muscleStr,fibrilStr,'ActiveTitin',wlcStr,'.mat'];
filePathRatMuscle = fullfile(projectFolders.output_structs_FittedModels,...
                             fileName);
if(mapToEDLModel==1)
    save(filePathRatMuscle,'ratMuscleModelParameters');
else
    save(filePathRatMuscle,'ratMuscleModelParameters');
end


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

    if(mapToEDLModel==1)
        idx = strfind(titleStr,'Soleus');
        titleStr = [titleStr(1,1:(idx-1)),'EDL',titleStr(1,(idx+6):end)];
    end

    title(titleStr);            
end


figure(figModelCurves);

if(mapToEDLModel==1)
    filePath = fullfile(projectFolders.output_plots_MuscleCurves,...
                        'fig_Pub_MuscleCurves_RatSoleusMappedToEDL.pdf');
end
print('-dpdf', filePath); 