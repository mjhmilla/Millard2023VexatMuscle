%%
% SPDX-FileCopyrightText: 2025 Matthew Millard <millard.matthew@gmail.com>
%
% SPDX-License-Identifier: MIT
%
%%
flag_OuterLoopMode=0;

if(flag_OuterLoopMode==0)
    %clc;
    close all;
    clear all;

    simConfigInput.runFitting              = 1; 
    simConfigInput.generatePlots           = 1;
    simConfigInput.fitToIndividualTrials   = 0; 
    simConfigInput.manuallySetTimeConstant = 0;

end


%%
% Parameters
%%
simConfig.runFitting              = simConfigInput.runFitting; 
simConfig.generatePlots           = simConfigInput.generatePlots;
simConfig.fitToIndividualTrials   = simConfigInput.fitToIndividualTrials; 
% 0: A single parameter that best fits all trials will be solved for 
% 1: The parameters that best fit each trial will be solved for
% Note: this only applies to a subset of parameters
simConfig.runSandbox              = 0;
simConfig.runSimUpdateResults     = 0;

%Simulations for testing

simConfig.trials                  = [1,2,3];
simConfig.defaultTrialId          = 0; %Set to 1,2,3 to fit to just this trial
simConfig.useDefaultModel         = 1;
simConfig.flag_debugFitting       = 1;
simConfig.numericalDiffLowPassFilterFrequencyHz=30; %



modelConfig.fibrilOption     = 'Fibril'; %Fibril or ''
modelConfig.wlcOption        = ''; %WLC or ''6
modelConfig.muscleName       = 'EDL';
modelConfig.experimentName   = 'TRSS2017';

fittingConfig.fittingEvaluationMethod = 'approximationByStateCalculationMinimal';
% 'approximationByStateCalculationMinimal'
% 'approximationByStateCalculation'
% 'approximationByInitialization'
% 'simulateFullModel'

fittingConfig.fitFl             =1;
fittingConfig.fitFv             =1;
fittingConfig.fv.scaleEnhancement=1;
if(simConfigInput.manuallySetTimeConstant==1)
    fittingConfig.fv.scaleEnhancement=0.75;
end

fittingConfig.fitTimeConstant   =0; 
fittingConfig.manuallySetTimeConstant = simConfigInput.manuallySetTimeConstant;
assert(~(fittingConfig.fitTimeConstant ...
    && fittingConfig.manuallySetTimeConstant), ...
    'Error: it does not make sense to fit something that is set manually');
%Lengthening time constant in Eqn. 16 of Millard, Franklin, Herzog

fittingConfig.fitKx             =0;
fittingConfig.fitQToF           =1;
fittingConfig.fitQToK           =0;
fittingConfig.fitf1HNPreload    =0;
fittingConfig.fitl1HNOffset     =0;

assert(~(fittingConfig.fitf1HNPreload && fittingConfig.fitl1HNOffset),...
    'Error: using both fitf1HNPreload and fitl1HNOffset does not make sense');
assert(~(fittingConfig.fitQToF && fittingConfig.fitQToK),...
  'Error: fitting Q to force and also Q to stiffness does not make sense');

fittingConfig.numberOfBisections = 10; %Gives us an parameter error of < 0.00098 or < 0.098%
fittingConfig.idxFvKey = 3;
fittingConfig.idxFlKey = 1;
fittingConfig.titin.trials = simConfig.trials; 
% [1]    : will fit just trial 1
% [1,2,3]: will fit to trials 1,2,3
fittingConfig.titin.individuallyFit = simConfig.fitToIndividualTrials; 
% 0: a set of titin parameters that best fits all trials will be solved
% 1: a set of titin parameters the best fits each individual trial will be solved  
%
fittingConfig.titin.applyToAllTrials = 1;
    %If a fittingConfig.titin.trials is a single trial, then setting
    %this flag to 1 will apply the fitted parameters to all other trials

pubPlotOptions.plotSmoothedStiffnessData    = 1;
pubPlotOptions.plotRawStiffnessData         = 0;
pubPlotOptions.stiffnessLowerForceBound     = 0.05;

pubPlotOptions.fceNMax   = 3.1;
pubPlotOptions.kceNMax   = 6.5;
pubPlotOptions.lceNMin   = 0.6;
pubPlotOptions.lceNMax   = 1.45;

%%
% Analysis outputs
%%
nTrials = length(simConfig.trials);
titinAnalysis(nTrials)=struct('lce',[],'fceN',[],'fxeN',[],'f2N',[],'k2N',[]);
for idx=simConfig.trials
    titinAnalysis(idx)= struct('lce',[],'fceN',[],'fxeN',[],'f2N',[],'k2N',[]);
end

%%
% Set up the files
rootDir         = getRootProjectDirectory();
projectFolders  = getProjectFolders(rootDir);

addpath( genpath(projectFolders.parameters)     );
addpath( genpath(projectFolders.curves)         );
addpath( genpath(projectFolders.experiments)    );
addpath( genpath(projectFolders.simulation)     );
addpath( genpath(projectFolders.models)         );
addpath( genpath(projectFolders.postprocessing) );




%%
% Plot parameters
%%
%
length135mm = 5.25;
length160mm = 6.5;
plotLayoutSettings = struct('numberOfHorizontalPlotColumns',  2,...
                            'numberOfVerticalPlotRows',       6,...
                            'flag_fixedPlotWidth',            1,...
                            'plotWidth',                      length160mm,...
                            'plotHeight',                     length160mm,...
                            'flag_usingOctave',               0);

numberOfHorizontalPlotColumns = plotLayoutSettings.numberOfHorizontalPlotColumns;
numberOfVerticalPlotRows      = plotLayoutSettings.numberOfVerticalPlotRows;
flag_fixedPlotWidth           = plotLayoutSettings.flag_fixedPlotWidth;
plotWidth                     = plotLayoutSettings.plotWidth;
plotHeight                    = plotLayoutSettings.plotHeight;
flag_usingOctave              = plotLayoutSettings.flag_usingOctave;

plotHorizMarginCm = 2.0;
plotVertMarginCm  = 2.0;

pageWidth   = (plotWidth+plotHorizMarginCm)*numberOfHorizontalPlotColumns...
                +plotHorizMarginCm + (plotWidth+plotHorizMarginCm)*2;

pageHeight  = (plotHeight+plotVertMarginCm)*numberOfVerticalPlotRows...
                +plotVertMarginCm;

plotConfigGeneric;

subPlotPanel(:,:,1) = subPlotPanel(:,:,1)+0.25;

plotConfig.subPlotPanel = subPlotPanel;
plotConfig.pageWidth    = pageWidth;
plotConfig.pageHeight   = pageHeight;


%%
% Line colors
%% 
cs = getPaulTolColourSchemes('vibrant');

n=0.75;

lineColors.orig(1,:)=[0,0,0];
lineColors.orig(2,:)=[44,133,229]./255;
lineColors.orig(3,:)=[63,181,175]./255;
lineColors.timeSeriesL = cs.blue;
lineColors.timeSeriesF = [0,0,0];%cs.red;
lineColors.annotationBox = [1,1,1].*0.75;

for i=1:1:3
    n0 = 0.25*(i-1);
    n1 = n0;%n0+0.25;
    lineColors.exp(i,:)         = [0,0,0].*(1-n1)+[1,1,1].*n1;  
    lineColors.simXE(i,:)       = cs.grey;
    lineColors.calcTitinF(i,:)  = lineColors.exp(i,:);%cs.blue.*(1-n0)+[1,1,1].*n0;
    lineColors.calcTitinK(i,:)  = lineColors.exp(i,:);%lineColors.calcTitinF(i,:);
    lineColors.simF(i,:)        = cs.blue.*(1-n0)+[1,1,1].*n0;%cs.red.*(1-n0)+[1,1,1].*n0;
    lineColors.simTitinF(i,:)   = cs.red.*(1-n0)+[1,1,1].*n0;%cs.magenta.*(1-n0)+[1,1,1].*n0;
    lineColors.simTitinK(i,:)   = lineColors.simTitinF(i,:);
end


plotConfig.lineColors = lineColors;

%%
% Load the models parameters
%%

nTrials = length(simConfig.trials);

ratFibrilModelsDefault(nTrials) = struct('musculotendon',[],...
                            'sarcomere',[],...
                            'falData',[],...
                            'fpeData',[],...
                            'curves',[],...
                            'fitting',[]);


for idxTrial = simConfig.trials

    ratFibrilModelsDefault(idxTrial) = struct('musculotendon',[],...
                            'sarcomere',[],...
                            'falData',[],...
                            'fpeData',[],...
                            'curves',[],...
                            'fitting',[]);
    
    if(simConfig.useDefaultModel==1)
        fileName = [    'rat',modelConfig.experimentName,...
                        modelConfig.muscleName,...
                        modelConfig.fibrilOption,...
                        'ActiveTitin',...
                        modelConfig.wlcOption,...
                        '_',num2str(simConfig.defaultTrialId),'.mat'];        
    else
        fileName = [    'rat',modelConfig.experimentName,...
                        modelConfig.muscleName,...
                        modelConfig.fibrilOption,...
                        'ActiveTitin',...
                        modelConfig.wlcOption,...
                        '_',num2str(idxTrial),'.mat'];
    end
    filePathRatMuscle = ...
        fullfile(projectFolders.output_structs_FittedModels,fileName);

    tmpModel   = load(filePathRatMuscle);
    
    modelFields = fields(ratFibrilModelsDefault(idxTrial));

    for idxField = 1:1:length(modelFields)
        ratFibrilModelsDefault(idxTrial).(modelFields{idxField}) = ...
            tmpModel.ratMuscleModelParameters.(modelFields{idxField});
    end

    ratFibrilModelsDefault(idxTrial).sarcomere.l1HNOffset=0;
    %
    % Update the default settings to be consistent with a skinned fiber
    %
    ratFibrilModelsDefault(idxTrial).sarcomere.scaleECM            =0;
    ratFibrilModelsDefault(idxTrial).sarcomere.scaleTitinProximal  =1;
    ratFibrilModelsDefault(idxTrial).sarcomere.scaleTitinDistal    =1;
    
    ratFibrilModelsDefault(idxTrial).sarcomere.normLengthTitinActinBondMinimum = 0.;
    normLengthTitinActinBondMinimum = ...
        ratFibrilModelsDefault(idxTrial).sarcomere.normLengthTitinActinBondMinimum;
    %fprintf('%1.3e norm-length-titin-actin-bond-minimum\n',...
    %        normLengthTitinActinBondMinimum);
    
    %
    % Use the fast-shortening, slow-lengthening model
    %    
    responseTimeScaling=1;
    if(fittingConfig.fitTimeConstant==1 ...
            || fittingConfig.manuallySetTimeConstant==1)
        responseTimeScaling=30;
    end
    
    ratFibrilModelsDefault(idxTrial).sarcomere.slidingTimeConstantBlendingParameter = 0.01;
    
    ratFibrilModelsDefault(idxTrial).sarcomere.slidingTimeConstantLengthening= ...
        ratFibrilModelsDefault(idxTrial).sarcomere.slidingTimeConstant...
        *responseTimeScaling;
    
    ratFibrilModelsDefault(idxTrial).sarcomere.slidingTimeConstantShortening= ...
        ratFibrilModelsDefault(idxTrial).sarcomere.slidingTimeConstant;
    
    ratFibrilModelsDefault(idxTrial).sarcomere.useVariableSlidingTimeConstant = 1;

end

%%
% Reference data
%%
expTRSS2017 = loadTRSS2017Data(projectFolders);

ratFibrilModelsFitted = ratFibrilModelsDefault;
fidFitting = [];

% Add the fitting options to the name
fittingOptions = fields(fittingConfig);
fittingNames = [];
fittingAbbr = [];
for i=1:1:length(fittingOptions)
    if(contains(fittingOptions{i},'fit'))
        name = fittingOptions{i};
        abbr = name(4:end);
        if(length(abbr)>4)
            abbr=abbr(1,1:4);
        end
        
        if(isempty(fittingNames))
            fittingNames = [{name}];
            fittingAbbr = [{abbr}];
        else
            fittingNames = [fittingNames,{name}];
            fittingAbbr = [fittingAbbr,{abbr}];
        end
    end
    if(contains(fittingOptions{i},'manuallySet'))
            name = fittingOptions{i};
            abbr = name(12:end);
            if(length(abbr)>4)
                abbr=abbr(1,1:4);
            end
            
            if(isempty(fittingNames))
                fittingNames = [{name}];
                fittingAbbr = [{['m-',abbr]}];
            else
                fittingNames = [fittingNames,{name}];
                fittingAbbr = [fittingAbbr,{['m-',abbr]}];
            end
        end    
end



fittingTrialsStr = '_';
for i=1:1:length(fittingConfig.titin.trials)
    fittingTrialsStr = ...
        [fittingTrialsStr,num2str(fittingConfig.titin.trials(1,i))];
end

for i=1:1:length(fittingNames)
    if(fittingConfig.(fittingNames{i})==1)
        fittingTrialsStr = [fittingTrialsStr,...
            '_',fittingAbbr{i}];
    end
end

if(fittingConfig.titin.individuallyFit==1)
    fittingTrialsStr = [fittingTrialsStr,'_i'];
end


fittingConfig.trialStr = fittingTrialsStr;

%%
% Sandbox
%%
if(simConfig.runSandbox==1)

    figSandbox=figure;

    tmp = load(fullfile(projectFolders.output_structs_FittedModels,...
            ['ratTRSS2017EDLFibrilActiveTitinFitted',...
             fittingConfig.trialStr,'.mat']));

    ratFibrilModelsFitted=tmp.ratFibrilModelsFitted;

    success=runParameterStudyRatFibrilTRSS2017_01(...
                             ratFibrilModelsFitted,...
                             expTRSS2017,...
                             simConfig,...
                             fittingConfig,...
                             plotConfig);
end
%%
% Run simulations and update the results
%%
if(simConfig.runSimUpdateResults==1)
    tmp = load(fullfile(projectFolders.output_structs_FittedModels,...
            ['ratTRSS2017EDLFibrilActiveTitinFitted',...
             fittingConfig.trialStr,'.mat']));

    ratFibrilModelsFitted=tmp.ratFibrilModelsFitted;
    optParams.name  = 'simulate';
    optParams.value = nan;
    optParams.expX  = [];
    optParams.expY  = [];
    optParams.expDYDX = nan;
    
    simConfigTmp=simConfig;
    simConfigTmp.trials =simConfig.trials;
    simConfigTmp.flag_debugFitting=0;
    fittingFraction=1;
    npts=500;
    
    figDebugFitting=figure;
    
    [optError,optErrorValues,figDebugFitting,...
        ratFibrilModelsFittedUpd,benchRecordFitted] =...
        calcErrorTRSS2017RampFraction(...
            optParams,...
            fittingFraction,...
            npts,...
            ratFibrilModelsFitted,...
            expTRSS2017Digitized,...
            simConfigTmp,...
            figDebugFitting,...
            plotConfig.subPlotPanel,...
            plotConfig.lineColors.simTitinK,...
            'simulateFullModel',...
            projectFolders);


    save(fullfile(projectFolders.output_structs_FittedModels,...
        ['ratTRSS2017EDLFibrilActiveTitinFitted',...
        fittingConfig.trialStr,'.mat']),...
        'ratFibrilModelsFitted');

    save(fullfile(projectFolders.output_structs_TRSS2017,...
         ['benchRecordVexat_TRSS2017_fitted',...
          fittingConfig.trialStr,'.mat']),...
         'benchRecordFitted');

%     save(fullfile(projectFolders.output_structs_TRSS2017,...
%          ['fittingInfo_ratTRSS2017EDLFibrilActiveTitinFitted',...
%           fittingConfig.trialStr,'.mat']),...
%          'fitInfo');
    
end

%%
% Fitting
%%

if(simConfig.runFitting==1)


    [ratFibrilModelsFitted, ...
     benchRecordFitted,...
     fitInfo] = ...
        fitRatFibrilTRSS2017(ratFibrilModelsDefault,...
                             expTRSS2017,...
                             simConfig,...
                             fittingConfig,...
                             plotConfig,...
                             projectFolders);

    save(fullfile(projectFolders.output_structs_FittedModels,...
        ['ratTRSS2017EDLFibrilActiveTitinFitted',...
        fittingConfig.trialStr,'.mat']),...
        'ratFibrilModelsFitted');

    save(fullfile(projectFolders.output_structs_TRSS2017,...
         ['benchRecordVexat_TRSS2017_fitted',...
          fittingConfig.trialStr,'.mat']),...
         'benchRecordFitted');

    save(fullfile(projectFolders.output_structs_TRSS2017,...
         ['fittingInfo_ratTRSS2017EDLFibrilActiveTitinFitted',...
          fittingConfig.trialStr,'.mat']),...
         'fitInfo');
    
end

%%
% Plotting
%%

if(simConfig.generatePlots==1)
    load(fullfile(projectFolders.output_structs_TRSS2017,...
            ['fittingInfo_ratTRSS2017EDLFibrilActiveTitinFitted',...
             fittingConfig.trialStr,'.mat']));    
    load(fullfile(projectFolders.output_structs_TRSS2017,...
            ['benchRecordVexat_TRSS2017_fitted',...
             fittingConfig.trialStr,'.mat']));
    tmp = load(fullfile(projectFolders.output_structs_FittedModels,...
            ['ratTRSS2017EDLFibrilActiveTitinFitted',...
             fittingConfig.trialStr,'.mat']));

    ratFibrilModelsFitted=tmp.ratFibrilModelsFitted;


    if(simConfig.fitToIndividualTrials==1)
      filePathTable = fullfile(projectFolders.output_tables_TRSS2017,...
                      ['fittingInfo_ratTRSS2017EDLFibrilActiveTitinFitted',...
                        fittingConfig.trialStr,'.tex']);    

    else
      filePathTable = fullfile(projectFolders.output_tables_TRSS2017,...
                      ['fittingInfo_ratTRSS2017EDLFibrilActiveTitinFitted',...
                        fittingConfig.trialStr,'.tex']);    
    end

    success = writeTRSS2017ErrorTable(filePathTable,fitInfo,...
                                   simConfig,projectFolders);

    figPub=figure;

    figPub = plotRatFibrilTRSS2017( ...
                figPub,...
                ratFibrilModelsDefault,...
                ratFibrilModelsFitted,...
                fitInfo,...
                benchRecordFitted,...
                expTRSS2017,...
                simConfig,...
                fittingConfig,...
                plotConfig,...
                pubPlotOptions);


    pause(0.1);
    figure(figPub);   
    pause(0.1);

    configPlotExporter;

    filePath = fullfile(projectFolders.output_plots_TRSS2017,...
                        ['fig_Sim_TRSS2017',...
                        fittingConfig.trialStr,'.pdf']);
    
    print('-dpdf', filePath); 
    pause(0.1);

end
