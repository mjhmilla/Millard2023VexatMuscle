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

runSimulations          = 0;
specimenTemperature     = 12; %As in 12 degrees centrigrade

simulationFileName      = 'benchRecordVexat_TRSS2017.mat';

validMuscles = {'SOL','EDL'};
muscleName = validMuscles{2};

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

fileName = ['rat',muscleName,fibrilStr,'ActiveTitin',wlcStr,'.mat'];
filePathRatMuscle = fullfile(projectFolders.output_structs_FittedModels,...
                             fileName);
load(filePathRatMuscle);

assert(ratMuscleModelParameters.musculotendon.temperature==specimenTemperature,...
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


lceOptMdl   = ratMuscleModelParameters.musculotendon.optimalFiberLength;
lceOptData  = min(ratMuscleData(1).activeLengtheningData(3).x);

%%
% Simulation parameters
%%
simSoln         = [];
folderTRSS2017  = projectFolders.output_structs_TRSS2017;

musculotendon   = ratMuscleModelParameters.musculotendon;
sarcomere       = ratMuscleModelParameters.sarcomere;
curves          = ratMuscleModelParameters.curves;

lceNTA_00 = sarcomere.normLengthTitinActinBondMinimum;
lceNTA_01 = sarcomere.normLengthTitinActinBondMinimum;

betaA_TiA_00 = sarcomere.normMaxActiveTitinToActinDamping;
betaA_TiA_01 = sarcomere.normMaxActiveTitinToActinDamping*10;

disp('Changing normLengthTitinActinBondMinimum:');
fprintf('%1.3f\tfrom\n', lceNTA_00);
fprintf('%1.3f\tto\n'  , lceNTA_01);

sarcomere.normLengthTitinActinBondMinimum = lceNTA_01;

disp('Changing normActivePevkDamping:');
fprintf('%1.3f\tfrom\n', betaA_TiA_00);
fprintf('%1.3f\tto\n'  , betaA_TiA_01);

sarcomere.normMaxActiveTitinToActinDamping = betaA_TiA_01;


sarcomere.scaleECM=0;
sarcomere.scaleTitinProximal=1;
sarcomere.scaleTitinDistal=1;



responseTimeScaling=20;

sarcomere.slidingTimeConstantBlendingParameter = 0.01;

sarcomere.slidingTimeConstantLengthening= ...
    sarcomere.slidingTimeConstant*responseTimeScaling;

sarcomere.slidingTimeConstantShortening= ...
    sarcomere.slidingTimeConstant;

sarcomere.useVariableSlidingTimeConstant = 1;



disp('Note: updateMillard2023VexatCache.m update near line 1198');
disp('  to modulate the time constant so that it is smaller when');
disp('  shortening and larger when lengthening');

curves.useCalibratedCurves=1;


%%
% Run simulations
%%

if(runSimulations==1)

    [success] = runTomalkaRodeSchumacherSiebert2017SimulationsVEXAT( ...
                              ratMuscleData(1),...
                              musculotendon,...
                              sarcomere,...
                              curves,...
                              fullfile(folderTRSS2017,simulationFileName));
    simSoln = load(fullfile(folderTRSS2017,simulationFileName));
else
    simSoln = load(fullfile(folderTRSS2017,simulationFileName));
end
%%
% Plot results
%%

plotColors = getPaulTolColourSchemes('vibrant');


lineColors =zeros(3,3);
lineColors(1,:) = [0,0,0];
lineColors(2,:) = [44,133,229]./255;
lineColors(3,:) = [63,181,175]./255;

lineColors(4,:) = [0,0,0];
lineColors(5,:) = [1,0,0];
lineColors(6,:) = [0,1,0];
lineColors(7,:) = [44,133,229]./255;

lineColors(8,:) = [1,1,1].*0.4;
lineColors(9,:) = [1,1,1].*0.4;

for i=1:1:7
    lineColors(i,:) = lineColors(i,:).*0.25 + [1,1,1].*0.75;
end

lineColors(10,:) = [0,0,0];
lineColors(11,:) = plotColors.orange;
lineColors(12,:) = [1,1,1].*0.9;
lineColors(13,:) = [1,1,1].*0.85;

lineColors(14,:) = [0,0,0];
lineColors(15,:) = [44,133,229]./255;
lineColors(16,:) = [63,181,175]./255;

lineColors(17,:) = [0,0,0];
lineColors(18,:) = [1,0,0];
lineColors(19,:) = [0,1,0];
lineColors(20,:) = [44,133,229]./255;

idxExpFig2A = [1:3];
idxExpFig3A = [4:7];
idxExpMark  = [8,9];
idxMdl      = [10,11,12,13];

idxMdlFig2A = [14:16];
idxMdlFig3A = [17:20];

idxGreys =[12,13];

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



figH = figure;
subplot(subplot('Position', reshape(subPlotPanel(1,1,:),1,4)));
    lengthTicks =[];
    for idx=1:1:3
        yyaxis left;    
        plot(simSoln.benchRecord.time(:,idx),...
             simSoln.benchRecord.fiberForce(:,idx),...
             '-','Color',lineColors(idxMdlFig2A(idx),:),...
             'LineWidth',mdlLineWidth,...
             'DisplayName',sprintf('Sim %i',idx));
        hold on;

    
        yyaxis right;
        n = (idx-1)/2;
        lpColor = [0.5,0.3,0].*(1-n) + [1,0.6,0].*(n);

        plot(simSoln.benchRecord.time(:,idx),...
             simSoln.benchRecord.pathLength(:,idx)./lceOptData,...
             '--','Color',lpColor,...
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
    ylim([-0.5,2]);

    yyaxis right;
    lengthTicks = lengthTicks ./ lceOptData;
    lengthTicks= sort(lengthTicks);
    lengthTicks = unique(lengthTicks);
    lengthTicks = round(lengthTicks,2);

    dl = (max(lengthTicks)-min(lengthTicks))*0.1;
    ylim([min(lengthTicks)-dl,max(lengthTicks)*2]);
    yticks([lengthTicks]);
    ylabel('Norm. Path length ($$\ell / \ell_o^M$$)');
    
    box off;
    title('Tomalka et al. 2017 Fig 2A (time domain)');

subplot(subplot('Position', reshape(subPlotPanel(1,2,:),1,4)));
    for idx=1:1:3
        plot(simSoln.benchRecord.time(:,idx),...
             simSoln.benchRecord.normProximalTitinLength(:,idx),...
             '-','Color',lineColors(idxMdlFig2A(idx),:),...
             'LineWidth',mdlLineWidth,...         
             'DisplayName',sprintf('Sim %i',idx));
        hold on;
    end
    box off;
    xlabel('Time (s)');
    ylabel('Norm. Length ($$\ell^1/\ell_o^M$$)');
    title('Titin-Actin bond length');

subplot(subplot('Position', reshape(subPlotPanel(2,1,:),1,4)));
    fill(falFill.x,falFill.y,lineColors(idxMdl(3),:),...
        'EdgeColor','none', 'HandleVisibility','off');
    hold on;
    fill(ftpFill.x,ftpFill.y,lineColors(idxMdl(4),:),...
        'EdgeColor','none', 'HandleVisibility','off');
    hold on;    

    plot(ratMuscleData(1).activeForceLengthData.x,...
         ratMuscleData(1).activeForceLengthData.y,...
         '.','Color',lineColors(idxExpMark(1),:),...
         'MarkerFaceColor',lineColors(idxExpMark(1),:),...
         'MarkerSize',3,'HandleVisibility','off');
    hold on;
    plot(ratMuscleData(1).passiveForceLengthData.x,...
         ratMuscleData(1).passiveForceLengthData.y,...
         'x','Color',lineColors(idxExpMark(2),:),...
         'MarkerFaceColor',lineColors(idxExpMark(2),:),...
         'MarkerSize',2,'HandleVisibility','off');
    hold on;
    for i=1:1:length(idxExpFig2A)
        plot(ratMuscleData(1).activeLengtheningData(i).x,...
             ratMuscleData(1).activeLengtheningData(i).y,...
             '-','Color',lineColors(idxExpFig2A(i),:),...
             'LineWidth',expLineWidth,...
             'DisplayName',sprintf('Exp %i',i));
        hold on;
    end
    for idx=1:1:3
        plot(simSoln.benchRecord.normFiberLength(:,idx).*lceOptMdl,...
             simSoln.benchRecord.fiberForce(:,idx),...
             '-','Color',lineColors(idxMdlFig2A(idx),:),...
             'DisplayName',sprintf('Sim %i',idx),...
             'LineWidth',mdlLineWidth);
        hold on;
    end
    box off;
    
    legend('Location','NorthWest');
    legend boxoff;


    xlim([1,4]);
    ylim([0,3]);
    xticks([1,1.5,2,2.5,3,3.5,4]);
    yticks([0,0.5,1,1.5,2,2.5]);

    xlabel('Norm. Length ($$\ell/\ell_o^M$$)');
    ylabel('Norm. Force ($$f/f_o^M$$)');
    title('Simulation of Tomalka et al. 2017 Fig. 2A');

subplot(subplot('Position', reshape(subPlotPanel(2,2,:),1,4)));
    fill(falFill.x,falFill.y,lineColors(idxMdl(3),:),...
        'EdgeColor','none','HandleVisibility','off');
    hold on;
    fill(ftpFill.x,ftpFill.y,lineColors(idxMdl(4),:),...
        'EdgeColor','none','HandleVisibility','off');
    hold on;    


    plot(ratMuscleData(1).activeForceLengthData.x,...
         ratMuscleData(1).activeForceLengthData.y,...
         '.','Color',lineColors(idxExpMark(1),:),...
         'MarkerFaceColor',lineColors(idxExpMark(1),:),...
         'MarkerSize',3,'HandleVisibility','off');
    hold on;
    plot(ratMuscleData(1).passiveForceLengthData.x,...
         ratMuscleData(1).passiveForceLengthData.y,...
         'x','Color',lineColors(idxExpMark(2),:),...
         'MarkerFaceColor',lineColors(idxExpMark(2),:),...
         'MarkerSize',2,'HandleVisibility','off');
    hold on;

    for i=1:1:length(idxExpFig3A)
        plot(ratMuscleData(1).activeLengtheningBDMData(i).x,...
             ratMuscleData(1).activeLengtheningBDMData(i).y,...
             '-','Color',lineColors(idxExpFig3A(i),:),...
             'LineWidth',expLineWidth,...
             'DisplayName',sprintf('BDM Exp %i',i));
        hold on;
    end

    plot(simSoln.benchRecord.normFiberLength(:,2).*lceOptMdl,...
         simSoln.benchRecord.normDistalTitinForce(:,2),...
         '-','Color',lineColors(idxMdlFig2A(2),:),...
         'DisplayName','Sim 2',...
         'LineWidth',mdlLineWidth);


    legend('Location','NorthWest');
    legend boxoff;
    
    xlim([1,4]);
    ylim([0,3]);
    xticks([1,1.5,2,2.5,3,3.5,4]);
    yticks([0,0.5,1,1.5,2,2.5]);
    hold on;
    box off;
    xlabel('Norm. Length ($$\ell/\ell_o^M$$)');
    ylabel('Norm. Force ($$f/f_o^M$$)');
    title('Simulation of Tomalka et al. 2017 Fig. 3A');

subplot(subplot('Position', reshape(subPlotPanel(3,1,:),1,4)));
    fill(fvFill.x,fvFill.y,lineColors(idxGreys(2),:),'EdgeColor','none');
    hold on;

    vmaxNorm = musculotendon.maximumNormalizedFiberVelocity;
    fvNSim = [];
    for idx=4:1:6
        vceN = simSoln.benchRecord.pathVelocity(end,idx) / (lceOptMdl);
        fceN = simSoln.benchRecord.normFiberForce(end,idx);
        falN = simSoln.benchRecord.fiberActiveForceLengthMultiplier(end,idx);
        lceN = simSoln.benchRecord.normFiberLength(end,idx);
        fpeN = simSoln.benchRecord.normPassiveFiberForce(end,idx);

        fvN = (fceN-fpeN)/falN;
        fvNSim = [fvNSim;...
                  vceN,fvN];
    end
    plot(fvNSim(:,1),fvNSim(:,2),'x','Color',[0,0,0],...
         'MarkerFaceColor',[0,0,0],'MarkerSize',4,'DisplayName','Sim Fv');
    hold on;
    xticks([-vmaxNorm,0,vmaxNorm]);
    yticks([0,1]);
    axis tight;
    box off;

    xlabel('Norm. Velocity ($$v/\ell_o^M$$)');
    ylabel('Norm. Force ($$f/f_o^M$$');
    title('Check: Simulated Force-Velocity Relation');

subplot(subplot('Position', reshape(subPlotPanel(3,2,:),1,4)));
    colorA =[1,0,0];
    lpMin = inf;
    lpMax = -inf;
    for idx=4:1:6
        n = (idx-4)/2;
        yyaxis left;
        fvColor = [0.25,0,0.25].*(1-n) + [1,0,1].*(n);
        lpColor = [0.5,0.3,0].*(1-n) + [1,0.6,0].*(n);

        plot(simSoln.benchRecord.time(:,idx),...
             simSoln.benchRecord.normFiberForce(:,idx),...
             '-','Color',fvColor);        
        hold on;

        plot(simSoln.benchRecord.time(end,idx),...
             simSoln.benchRecord.normFiberForce(end,idx),...
             'x','Color',[0,0,0],...
             'MarkerFaceColor',[0,0,0],...
             'MarkerSize',4,'DisplayName','Sim Fv');        
        hold on;

        yyaxis right;
        plot(simSoln.benchRecord.time(:,idx),...
             simSoln.benchRecord.pathLength(:,idx)./lceOptMdl,...
             '--','Color',lpColor);
        hold on;

        lpMin = min(lpMin, min(simSoln.benchRecord.pathLength(:,idx)./lceOptMdl));
        lpMax = max(lpMax, max(simSoln.benchRecord.pathLength(:,idx)./lceOptMdl));
    end    
    box off;
    yyaxis left
    ylim([-0.5,1]);
    xlabel('Time (s)');
    ylabel('Norm. Force ($$f/f_o^M$$)');

    yyaxis right;
    ylabel('Norm. Path Length ($$\ell_p/\ell_o^M$$)');
    yticks(round( [min(simSoln.benchRecord.pathLength(:,4)),...
                   max(simSoln.benchRecord.pathLength(:,4))]./lceOptMdl,2));

    dpDelta = (lpMax-lpMin);
    ylim([lpMin, lpMax+2*dpDelta]);
    title('Check: Simulated Force-Velocity Relation');

figure(figH);    
configPlotExporter;
filePath = fullfile(projectFolders.output_plots_TRSS2017,...
                    'fig_Sim_TRSS2017_Fig2A_3A.pdf');
print('-dpdf', filePath); 