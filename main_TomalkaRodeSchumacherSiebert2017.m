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

simulationFileName      = 'benchRecordVexat_TRSS2017.mat';


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
muscleStr='Soleus';
if(mapToEDLModel==1)
    muscleStr='EDL';
end

fileName = ['rat',muscleStr,fibrilStr,'ActiveTitin',wlcStr,'.mat'];
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


%%
% Reference data
%%
[ratMuscleData, ratMuscleMetaData] = ...
        loadRatSkeletalMuscleData(projectFolders);



%%
% Simulation: length change function
%%
lceOptMdl   = ratMuscleModelParameters.musculotendon.optimalFiberLength;
lceOptData  = min(ratMuscleData(1).activeLengtheningData(3).x);

vmax = ...
  ratMuscleModelParameters.musculotendon.maximumNormalizedFiberVelocity;

lmin = min(ratMuscleData(1).activeLengtheningData(2).x);
lmax = max(ratMuscleData(1).activeLengtheningData(2).x);

%Starting
e0 = 0;
l0 = lmin;
t0 = 0;

%Excitation 1
l1 = lmin;
e1 = 1;
t1 = t0+1;

%Ramp start
l2 = lmin;
e2 = 1;
t2 = t1+3;

%Ramp End
rampSlope = 0.11*vmax*lceOptData;
l3 = lmax;
t3 = t2 + (l3-l2)/(rampSlope);
e3 = 1;

timeSpan = [t0,t3];

%Setup functions
activation=1;
excitationFcn = @(argT)calcStepFunction(argT,...
                  t1,...
                  inf,...
                  activation);

pathLengthFcn = @(argT)calcRampStateSharp(...
                       argT,t2,t3,l2,rampSlope);

activationFcn = ...
    @(argU,argA)calcFirstOrderActivationDerivative(argU,argA, ...
                ratMuscleModelParameters.sarcomere.activationTimeConstant,...
                ratMuscleModelParameters.sarcomere.deactivationTimeConstant,0);

% Fitting points
lA = (lmin+(1/3)*(lmax-lmin));
lB =  lmax;
npts = 10;
lkp  = [lA:((lB-lA)/(npts-1)):lB]';

fkp  = interp1(ratMuscleData(1).activeLengtheningData(2).x,...
               ratMuscleData(1).activeLengtheningData(2).y,...
               lkp);

fittingKeyPoints = [lkp, fkp];


%%
% Simulation parameters
%%
simSoln         = [];
folderTRSS2017  = projectFolders.output_structs_TRSS2017;

if(runSimulations==1)
    flag_fitTitinModel  = 0;
    
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

    responseTimeScaling=15; %20

    sarcomere.slidingTimeConstant = ...
        sarcomere.slidingTimeConstant*responseTimeScaling;

    sarcomere.normCrossBridgeCyclingDamping = ...
        sarcomere.normCrossBridgeCyclingDamping...
        *responseTimeScaling*0.5;

    curves.useCalibratedCurves=1;
    
    [success] = runTomalkaRodeSchumacherSiebert2017SimulationsVEXAT( ...
                              timeSpan,...
                              excitationFcn,...
                              activationFcn,...
                              pathLengthFcn,...
                              musculotendon,...
                              sarcomere,...
                              curves,...
                              flag_fitTitinModel,...
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

idxExpFig2A = [1:3];
idxExpFig3A = [4:7];
idxExpMark  = [8,9];
idxMdl      = [10,11,12,13];

expLineWidth=1.5;
mdlLineWidth=0.5;

falSample = calcBezierYFcnXCurveSampleVector(...
            ratMuscleModelParameters.curves.activeForceLengthCurve,...
            100,...
            [1,4]./lceOptData);

fpeSample = calcBezierYFcnXCurveSampleVector(...
            ratMuscleModelParameters.curves.fiberForceLengthCurve,...
            100,...
            [1,4]./lceOptData);

falFill.x = [falSample.x;...
             fliplr(falSample.x(2:end-1)')'];
falFill.x = falFill.x.*lceOptMdl;
falFill.y = [zeros(size(falSample.y));...
             fliplr(falSample.y(2:end-1)')'];

fpeFill.x = [fpeSample.x;...
             fliplr(fpeSample.x(2:end-1)')'];
fpeFill.x = fpeFill.x.*lceOptMdl;
fpeFill.y = [zeros(size(fpeSample.y));...
             fliplr(fpeSample.y(2:end-1)')'];

figH = figure;
subplot(subplot('Position', reshape(subPlotPanel(1,1,:),1,4)));
    yyaxis left;
    plot(simSoln.benchRecord.time(:,1),...
         simSoln.benchRecord.fiberForce(:,1),...
         '-','Color',lineColors(idxMdl(1),:),...
         'LineWidth',mdlLineWidth,...
         'DisplayName','VEXAT');
    hold on;
    xlabel('Time (s)');
    ylabel('Norm. Force ($$f/f_o^M$$)');

    yyaxis right;
    plot(simSoln.benchRecord.time(:,1),...
         simSoln.benchRecord.pathLength(:,1),...
         '-','Color',lineColors(idxMdl(2),:),...
         'LineWidth',mdlLineWidth,...         
         'DisplayName','VEXAT');
    hold on;
    dl = ( max(simSoln.benchRecord.pathLength(:,1)) ...
          -min(simSoln.benchRecord.pathLength(:,1)))*0.1;

    ylim([0,...
          max(simSoln.benchRecord.pathLength(:,1))*3]);
    yticks([0,...
            round(simSoln.benchRecord.pathLength(1,1)  ,2),...
            round(simSoln.benchRecord.pathLength(end,1),2)]);
    ylabel('Path length ($$\mu m$$)');
    
    box off;
    title('Rat EDL fiber force');

subplot(subplot('Position', reshape(subPlotPanel(1,2,:),1,4)));
    plot(simSoln.benchRecord.time(:,1),...
         simSoln.benchRecord.normProximalTitinLength(:,1),...
         '-','Color',lineColors(idxMdl(1),:),...
         'LineWidth',mdlLineWidth,...         
         'DisplayName','VEXAT')
    hold on;
    box off;
    xlabel('Time (s)');
    ylabel('Norm. Length ($$\ell^1/\ell_o^M$$)');
    title('Prox. Titin Length');
subplot(subplot('Position', reshape(subPlotPanel(2,1,:),1,4)));
    fill(falFill.x,falFill.y,lineColors(idxMdl(3),:),'EdgeColor','none');
    hold on;
    fill(fpeFill.x,fpeFill.y,lineColors(idxMdl(4),:),'EdgeColor','none');
    hold on;    

    plot(ratMuscleData(1).activeForceLengthData.x,...
         ratMuscleData(1).activeForceLengthData.y,...
         '.','Color',lineColors(idxExpMark(1),:),...
         'MarkerFaceColor',lineColors(idxExpMark(1),:),...
         'MarkerSize',3);
    hold on;
    plot(ratMuscleData(1).passiveForceLengthData.x,...
         ratMuscleData(1).passiveForceLengthData.y,...
         'x','Color',lineColors(idxExpMark(2),:),...
         'MarkerFaceColor',lineColors(idxExpMark(2),:),...
         'MarkerSize',2);
    hold on;
    for i=1:1:length(idxExpFig2A)
        plot(ratMuscleData(1).activeLengtheningData(i).x,...
             ratMuscleData(1).activeLengtheningData(i).y,...
             '-','Color',lineColors(idxExpFig2A(i),:),...
             'LineWidth',expLineWidth);
        hold on;
    end
    plot(simSoln.benchRecord.normFiberLength(:,1).*lceOptMdl,...
         simSoln.benchRecord.fiberForce(:,1),...
         '-','Color',lineColors(idxMdl(1),:),...
         'DisplayName','VEXAT',...
         'LineWidth',mdlLineWidth);
    hold on;
    box off;
    
    xlim([1,4]);
    ylim([0,3]);
    xticks([1,1.5,2,2.5,3,3.5,4]);
    yticks([0,0.5,1,1.5,2,2.5]);

    xlabel('Norm. Length ($$\ell/\ell_o^M$$)');
    ylabel('Norm. Force ($$f/f_o^M$$)');
    title('Rat EDL fiber force');

subplot(subplot('Position', reshape(subPlotPanel(2,2,:),1,4)));
    fill(falFill.x,falFill.y,lineColors(idxMdl(3),:),'EdgeColor','none');
    hold on;
    fill(fpeFill.x,fpeFill.y,lineColors(idxMdl(4),:),'EdgeColor','none');
    hold on;    


    plot(ratMuscleData(1).activeForceLengthData.x,...
         ratMuscleData(1).activeForceLengthData.y,...
         '.','Color',lineColors(idxExpMark(1),:),...
         'MarkerFaceColor',lineColors(idxExpMark(1),:),...
         'MarkerSize',3);
    hold on;
    plot(ratMuscleData(1).passiveForceLengthData.x,...
         ratMuscleData(1).passiveForceLengthData.y,...
         'x','Color',lineColors(idxExpMark(2),:),...
         'MarkerFaceColor',lineColors(idxExpMark(2),:),...
         'MarkerSize',2);
    hold on;

    for i=1:1:length(idxExpFig3A)
        plot(ratMuscleData(1).activeLengtheningBDMData(i).x,...
             ratMuscleData(1).activeLengtheningBDMData(i).y,...
             '-','Color',lineColors(idxExpFig3A(i),:),...
             'LineWidth',expLineWidth);
        hold on;
    end

    plot(simSoln.benchRecord.normFiberLength(:,1).*lceOptMdl,...
         simSoln.benchRecord.normDistalTitinForce(:,1),...
         '-','Color',lineColors(idxMdl(1),:),...
         'DisplayName','VEXAT',...
         'LineWidth',mdlLineWidth);

    xlim([1,4]);
    ylim([0,3]);
    xticks([1,1.5,2,2.5,3,3.5,4]);
    yticks([0,0.5,1,1.5,2,2.5]);
    hold on;
    box off;
    xlabel('Norm. Length ($$\ell/\ell_o^M$$)');
    ylabel('Norm. Force ($$f/f_o^M$$)');
    title('Rat EDL active lengthening');

figure(figH);    
configPlotExporter;
filePath = fullfile(projectFolders.output_plots_TRSS2017,...
                    'fig_Sim_TRSS2017_Fig2A_3A.pdf');
print('-dpdf', filePath); 