%%
% SPDX-FileCopyrightText: 2025 Matthew Millard <millard.matthew@gmail.com>
%
% SPDX-License-Identifier: MIT
%
%%
clc;
close all;
clear all;

%%
% Parameters
%%
mapToEDLModel           = 1;
makeFibrilModel         = 1;
useElasticTendon        = 1 && ~makeFibrilModel;
useWlcTitinModel        = 0;

runSimulations          = 1;
specimenTemperature     = 12; %As in 12 degrees centrigrade


validMuscles = {'SOL','EDL'};
muscleName = validMuscles{1};

validExperiments = {'TRSS2017','TWHSS2021','WTRS2024'};
experimentName = validExperiments{2};

simulationFileName   = sprintf('benchRecordVexat_%s.mat',experimentName);


%%
% Set up the files
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
% Load the model parameters
%%

wlcStr = '';
if(useWlcTitinModel==1)
    wlcStr='WLC';
end
fibrilStr='';
if(makeFibrilModel==1)
    fibrilStr='Fibril';
end

fileName = ['rat',experimentName,muscleName,...
             fibrilStr,'ActiveTitin',wlcStr,'.mat'];
%ratTWHSS2021SOLFibrilActiveTitin
filePathRatMuscle = fullfile(projectFolders.output_structs_FittedModels,...
                             fileName);
load(filePathRatMuscle);

assert(ratMuscleModelParameters.musculotendon.temperature ...
        == specimenTemperature,...
       'Error: model temperature does not match desired temperature');

if(makeFibrilModel==1)
    assert(ratMuscleModelParameters.musculotendon.tendonSlackLength==0,...
           ['Error: fibril model is desired but the model has',...
           ' a non-zero tendon slack length']);
end

%%
% Plot parameters
%%
plotLayoutSettings = struct('numberOfHorizontalPlotColumns',  2,...
                            'numberOfVerticalPlotRows',       3,...
                            'flag_fixedPlotWidth',            1,...
                            'plotWidth',                      6,...
                            'plotHeight',                     6,...
                            'flag_usingOctave',               0);

numberOfHorizontalPlotColumns = plotLayoutSettings.numberOfHorizontalPlotColumns;
numberOfVerticalPlotRows      = plotLayoutSettings.numberOfVerticalPlotRows;
flag_fixedPlotWidth           = plotLayoutSettings.flag_fixedPlotWidth;
plotWidth                     = plotLayoutSettings.plotWidth;
plotHeight                    = plotLayoutSettings.plotHeight;
flag_usingOctave              = plotLayoutSettings.flag_usingOctave;

plotHorizMarginCm = 3;
plotVertMarginCm  = 3;

pageWidth   = (plotWidth+plotHorizMarginCm)*numberOfHorizontalPlotColumns...
                +plotHorizMarginCm;

pageHeight  = (plotHeight+plotVertMarginCm)*numberOfVerticalPlotRows...
                +plotVertMarginCm;

plotConfigGeneric;


%%
% Reference data
%%
[ratMuscleData, ratMuscleMetaData] = ...
        loadRatSkeletalMuscleData(projectFolders);

index_SW1982 = ratMuscleMetaData.index_SW1982;
index_TWHSS2021 = ratMuscleMetaData.index_TWHSS2021;

lceOptMdl   = ratMuscleModelParameters.musculotendon.optimalFiberLength;
lceOptData  = lceOptMdl;

%%
% Simulation parameters
%%
simSoln         = [];
folderTWHSS2021  = projectFolders.output_structs_TWHSS2021;

musculotendon   = ratMuscleModelParameters.musculotendon;
sarcomere       = ratMuscleModelParameters.sarcomere;
curves          = ratMuscleModelParameters.curves;

%
%This is a skinned fiber
%
sarcomere.scaleECM=0;
sarcomere.scaleTitinProximal=1;
sarcomere.scaleTitinDistal=1;

%
%Adjust the time-constants of lengthening and shortening
%

sarcomere.forceVelocityCalibrationFactor = 1;

responseTimeScaling=20;

sarcomere.slidingTimeConstantBlendingParameter = 0.5;

sarcomere.slidingTimeConstantLengthening= ...
    sarcomere.slidingTimeConstant*responseTimeScaling;

sarcomere.slidingTimeConstantShortening= ...
    sarcomere.slidingTimeConstant;

sarcomere.useVariableSlidingTimeConstant = 1;

sarcomere.normCrossBridgeCyclingDamping = ...
    sarcomere.normCrossBridgeCyclingDamping;%*0.1;

disp('Using variable sliding time constant:');
fprintf('%1.3f\ttau-shortening\n', ...
    sarcomere.slidingTimeConstantShortening);
fprintf('%1.3f\ttau-lengthening\n'  ,...
    sarcomere.slidingTimeConstantLengthening);

fprintf('%1.3f\tMax. Titin-Actin Damping\n'  ,...
    sarcomere.normMaxActiveTitinToActinDamping);

sarcomere.normMaxActiveTitinToActinDamping = 5;

disp('Note: updateMillard2023VexatCache.m update near line 1198');
disp('  to modulate the time constant so that it is smaller when');
disp('  shortening and larger when lengthening');

curves.useCalibratedCurves=1;


%%
% Run simulations
%%


if(runSimulations==1)

    [success] = runTomalkaWeidnerHahnSieberlSiebert2021SimulationsVEXAT( ...
                              ratMuscleData(index_TWHSS2021),...
                              musculotendon,...
                              sarcomere,...
                              curves,...
                              fullfile(folderTWHSS2021,simulationFileName));

    simSoln = load(fullfile(folderTWHSS2021,simulationFileName));
else
    simSoln = load(fullfile(folderTWHSS2021,simulationFileName));
end
%%
% Plot results
%%

plotColors = getPaulTolColourSchemes('vibrant');

seriesColors =zeros(3,3);
seriesColors(1,:) = [236,179, 41]./255;
seriesColors(2,:) = [220, 89, 36]./255;
seriesColors(3,:) = [ 18,123,193]./255;

expSeriesNames = {  'Exp: SSC 85\% $$v_o$$','Exp: SSC 60\% $$v_o$$',...
                    'Exp: SSC 30\% $$v_o$$','Exp: Pure SHO 85\% $$v_o$$',...
                    'Exp: $$P_1$$','Exp: $$P_2$$'};

mdlSeriesNames = {  'Sim: SSC 85\% $$v_o$$','Sim: SSC 60\% $$v_o$$',...
                    'Sim: SSC 30\% $$v_o$$','Sim: Pure SHO 85\% $$v_o$$',...
                    'Sim: $$P_1$$','Sim: $$P_2$$'};


greyLevel = 0.67;

lineColors =zeros(9,3);
lineColors(1,:) = seriesColors(1,:).*(1-greyLevel) + ([1,1,1].*0.5).*greyLevel;
lineColors(2,:) = seriesColors(2,:).*(1-greyLevel) + ([1,1,1].*0.5).*greyLevel;
lineColors(3,:) = seriesColors(3,:).*(1-greyLevel) + ([1,1,1].*0.5).*greyLevel;
lineColors(4,:) = [1,1,1].*greyLevel;
lineColors(5,:) = [0,0,0];
lineColors(6,:) = [0,0,0];

lineColors(7,:) = seriesColors(1,:);
lineColors(8,:) = seriesColors(2,:);
lineColors(9,:) = seriesColors(3,:);
lineColors(10,:)= [0,0,0];

lineColors(11,:) = [1,1,1].*.7;
lineColors(12,:) = [1,1,1].*.9;

idxExp = [1,2,3,4,5,6];
idxMdl = [7,8,9,10];

idxGreys = [11,12];

expLineWidth=1.5;
mdlLineWidth=0.5;

falSample = calcBezierYFcnXCurveSampleVector(...
            ratMuscleModelParameters.curves.activeForceLengthCurve,...
            100,...
            [1,4]./lceOptData);

fvSample = calcBezierYFcnXCurveSampleVector(...
            ratMuscleModelParameters.curves.fiberForceVelocityCurve,...
            100,...
            [-1,1]);

titinCurveSample = ...
  sampleTitinCurves20250217(...
    curves,...
    sarcomere,...
    100);

falFill.x = [falSample.x;...
             fliplr(falSample.x(2:end-1)')'];
falFill.x = falFill.x.*lceOptMdl;
falFill.y = [zeros(size(falSample.y));...
             fliplr(falSample.y(2:end-1)')'];

fvSample.x = fvSample.x*sarcomere.forceVelocityCalibrationFactor;
fvFill.x = [fvSample.x;...
             fliplr(fvSample.x(1:end)')'];
fvFill.x = fvFill.x.*lceOptMdl;
fvFill.y = [zeros(size(fvSample.y));...
             fliplr(fvSample.y(1:end)')'];

ftpFill.x = [titinCurveSample.curveSampleTitin.x;...
             fliplr(titinCurveSample.curveSampleTitin.x(2:end-1)')'];
ftpFill.x = ftpFill.x.*(2*lceOptMdl);
ftpFill.y = [zeros(size(titinCurveSample.curveSampleTitin.y));...
             fliplr(titinCurveSample.curveSampleTitin.y(2:end-1)')'];


index_TWHSS2021 = ratMuscleMetaData.index_TWHSS2021;

figH = figure;

subplot(subplot('Position', reshape(subPlotPanel(1,1,:),1,4)));
    fill(falFill.x,falFill.y,lineColors(idxGreys(2),:),...
        'EdgeColor','none', 'HandleVisibility','off');
    hold on;
    fill(ftpFill.x,ftpFill.y,lineColors(idxGreys(1),:),...
        'EdgeColor','none', 'HandleVisibility','off');
    hold on;    



    for i=1:1:4
        plot(ratMuscleData(index_TWHSS2021).activeStretchShorteningData(i).x,...
             ratMuscleData(index_TWHSS2021).activeStretchShorteningData(i).y,...
             '-','Color',lineColors(idxExp(i),:),...
             'LineWidth',expLineWidth,...
             'DisplayName',expSeriesNames{i});
        hold on;
    end

    for i=1:1:4
        plot(simSoln.benchRecord.normFiberLength(:,i).*lceOptMdl,...
             simSoln.benchRecord.fiberForce(:,i),...
             '-','Color',lineColors(idxMdl(i),:),...
             'DisplayName',mdlSeriesNames{i},...
             'LineWidth',mdlLineWidth);
        hold on;
    end
    box off;
    
    legend('Location','NorthWest');
    legend boxoff;


    xlim([1.9,2.5]);
    ylim([0,1.6]);
    xticks([2,2.4]);
    yticks([0,0.5,1,1.5]);

    xlabel('Length ($$\mu$$m)');
    ylabel('Norm. Force ($$f/f_o^M$$)');
    title('A. Simulation of Tomalka et al. 2021 Fig. 4');

subplot(subplot('Position', reshape(subPlotPanel(1,2,:),1,4)));

    fill(falFill.x,falFill.y,lineColors(idxGreys(2),:),...
        'EdgeColor','none', 'HandleVisibility','off');
    hold on;
    fill(ftpFill.x,ftpFill.y,lineColors(idxGreys(1),:),...
        'EdgeColor','none', 'HandleVisibility','off');
    hold on;    

    for i=1:1:3
        plot(ratMuscleData(index_TWHSS2021).activeStretchShorteningBLEData(i).x,...
             ratMuscleData(index_TWHSS2021).activeStretchShorteningBLEData(i).y,...
             '-','Color',lineColors(idxExp(i),:),...
             'LineWidth',expLineWidth,...
             'DisplayName',[expSeriesNames{i},' BLE']);
        hold on;
    end

    for i=1:1:4
        plot(simSoln.benchRecord.normFiberLength(:,i).*lceOptMdl,...
             simSoln.benchRecord.normDistalTitinForce(:,i),...
             '-','Color',lineColors(idxMdl(i),:),...
             'DisplayName',[mdlSeriesNames{i},' D.Titin'],...
             'LineWidth',mdlLineWidth);
        hold on;
    end
    box off;

    legend('Location','NorthWest');
    legend boxoff;


    xlim([1.9,2.5]);
    ylim([0,1.6]);
    xticks([2,2.4]);
    yticks([0,0.5,1,1.5]);

    xlabel('Length ($$\mu$$m)');
    ylabel('Norm. Force ($$f/f_o^M$$)');
    title('B. Simulation of Tomalka et al. 2021 Fig. 5');


subplot(subplot('Position', reshape(subPlotPanel(2,1,:),1,4)));
    lengthTicks =[];

    expTimeOffset=100-0.2;

    fceMdl0 = simSoln.benchRecord.fiberForce(1,1);
    expFceScale = ...
        fceMdl0 ...
        /ratMuscleData(index_TWHSS2021).activeStretchShorteningForceData(1).y(1);

    yyaxis left;  
    plot(ratMuscleData(index_TWHSS2021).activeStretchShorteningForceData(1).x-expTimeOffset,...
         ratMuscleData(index_TWHSS2021).activeStretchShorteningForceData(1).y.*expFceScale,...
         '-','Color',lineColors(idxGreys(1),:),...
         'LineWidth',expLineWidth,...
         'DisplayName','TWHSS2021: Fig. 2C');
    hold on;    

    yyaxis right;
    plot(ratMuscleData(index_TWHSS2021).activeStretchShorteningLengthData(1).x-expTimeOffset,...
         ratMuscleData(index_TWHSS2021).activeStretchShorteningLengthData(1).y,...
         '-','Color',lineColors(idxGreys(2),:),...
         'LineWidth',expLineWidth,...
         'DisplayName','TWHSS2021: Fig. 2D');
    hold on;    


    for idx=1:1:4
        yyaxis left;    
        plot(simSoln.benchRecord.time(:,idx),...
             simSoln.benchRecord.fiberForce(:,idx),...
             '-','Color',lineColors(idxMdl(idx),:),...
             'LineWidth',mdlLineWidth,...
             'DisplayName',mdlSeriesNames{i});
        hold on;

    
        yyaxis right;


        plot(simSoln.benchRecord.time(:,idx),...
             simSoln.benchRecord.pathLength(:,idx),...
             '--','Color',lineColors(idxMdl(idx),:),...
             'LineWidth',mdlLineWidth,...         
             'DisplayName',sprintf('Sim %i',idx));
        hold on;

        lengthTicks = [lengthTicks,...
            simSoln.benchRecord.pathLength(1,idx),...
            simSoln.benchRecord.pathLength(end,idx)];
    end

    yyaxis left;
    xlabel('Time (s)');
    ylabel('Norm. Force ($$f/f_o^M$$)');
    ylim([-1,1.5]);

    yyaxis right;

    lengthTicks= sort(lengthTicks);
    lengthTicks = unique(lengthTicks);
    lengthTicks = round(lengthTicks,2);

    dl = (max(lengthTicks)-min(lengthTicks))*0.1;
    ylim([min(lengthTicks)-dl,max(lengthTicks)*2]);
    yticks([lengthTicks]);
    ylabel('Norm. Path length ($$\ell / \ell_o^M$$)');
    
    box off;
    title('C. Tomalka et al. 2021 (time domain)');

subplot(subplot('Position', reshape(subPlotPanel(2,2,:),1,4)));
    for idx=1:1:4
        plot(simSoln.benchRecord.time(:,idx),...
             simSoln.benchRecord.normProximalTitinLength(:,idx),...
             '-','Color',lineColors(idxMdl(idx),:),...
             'LineWidth',mdlLineWidth,...         
             'DisplayName',mdlSeriesNames{idx});
        hold on;
    end
    box off;
    xlabel('Time (s)');
    ylabel('Norm. Length ($$\ell^1/\ell_o^M$$)');
    title('D. Titin-Actin bond length');


figure(figH);    
configPlotExporter;
filePath = fullfile(projectFolders.output_plots_TWHSS2021,...
                    'fig_Sim_TWHSS2021.pdf');
print('-dpdf', filePath); 