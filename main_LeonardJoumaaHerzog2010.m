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

%Simulation of
%Leonard TR, Joumaa V, Herzog W. An activatable molecular spring reduces 
%muscle tearing during extreme stretching. Journal of biomechanics. 
%2010 Nov 16;43(15):3063-6.

flag_outerLoopMode = 1;

rootDir         = getRootProjectDirectory();
projectFolders  = getProjectFolders(rootDir);

if(flag_outerLoopMode == 0)
  clc;
  close all;
  clear all;

  rootDir         = getRootProjectDirectory();
  projectFolders  = getProjectFolders(rootDir);

  flag_useTunedRabbitPsoasModel     = 0;

  flag_simulateHillModel            = 0; 
  flag_simulateVexatModel           = 0;

  flag_useCalibratedVexatCurves     = 1;
  flag_useTwoSidedTitinCurves       = 0;

  flag_simulateActiveStretch  = 1;
  flag_simulatePassiveStretch = 1;
  
  flag_plotData         = 1;  
  flag_savePlotsToFile  = 1;
  
  nominalNormalizedFiberLength = 1.0;

%  flag_fitTitin                   = 0;  
  normActiveTitinToActinDamping   = [];% take the default. Prev. 6; 
  normPassiveTitinToActinDamping  = []; %[]: take the default. Prev. 1.5; 
     
  figLeonardJoumaaHerzog2010Fig2Comparision = figure;
  flag_useOctave                    = 0;
  %WLCVexatResults = 'experiments/LeonardJoumaaHerzog2010/benchRecordVexat_RigidTendon_K44p31D0p50Tau_LJH2010__TiAD1000p00_TiPD1p50_NomLen1p00_90Hz_TiAdj.mat';

  flag_useFig3KirchBoskovRymer1994              = 0; 
  flag_useElasticTendon                         = 0; 
  flag_useFiberDamping                          = 1;
  fiberDampingCoefficient = 0.01;
   

  
end 

assert(flag_useElasticTendon ==0,...
       'This is a skinned fibril experiment: there is no tendon!');


addpath( genpath(projectFolders.parameters)     );
addpath( genpath(projectFolders.curves)         );
addpath( genpath(projectFolders.experiments)    );
addpath( genpath(projectFolders.simulation)     );
addpath( genpath(projectFolders.models)         );
addpath( genpath(projectFolders.postprocessing) );

plotLayoutSettings = struct('numberOfHorizontalPlotColumns',  1,...
                            'numberOfVerticalPlotRows',       3,...
                            'flag_fixedPlotWidth',            1,...
                            'plotWidth',                      7,...
                            'plotHeight',                     4,...
                            'flag_usingOctave',               0);
% Single column plots of:
% length vs time
% force vs. time
% stiffness vs. time
% damping vs time.
numberOfHorizontalPlotColumns = plotLayoutSettings.numberOfHorizontalPlotColumns;
numberOfVerticalPlotRows      = plotLayoutSettings.numberOfVerticalPlotRows;
flag_fixedPlotWidth           = plotLayoutSettings.flag_fixedPlotWidth;
plotWidth                     = plotLayoutSettings.plotWidth;
plotHeight                    = plotLayoutSettings.plotHeight;
flag_usingOctave              = plotLayoutSettings.flag_usingOctave;
plotConfig;




%%
% Configure the muscle models:  load the configurations generated from 
%                               main_KirschBoskovRymer1994.m
%%

musculotendonProperties = [];
sarcomereProperties     = [];
normMuscleCurves        = [];
fitting                 = [];



outputFileEndingHill = 'Default';


%Basic parameters for the Hill model
if(flag_useTunedRabbitPsoasModel==0)
    tmp=load(fullfile(projectFolders.output_structs_FittedModels,...
                      'rabbitPsoasFibrilWLC.mat'));

    musculotendonProperties   = tmp.rabbitPsoasFibrilWLC.musculotendon;
    sarcomereProperties       = tmp.rabbitPsoasFibrilWLC.sarcomere;
    normMuscleCurves          = tmp.rabbitPsoasFibrilWLC.curves;
    fitting                   = tmp.rabbitPsoasFibrilWLC.fitting;    
    outputFileEndingVexat = 'WLC-Titin';

else
    tmp=load(fullfile(projectFolders.output_structs_FittedModels,...
                     'tunedRabbitPsoasFibril.mat'));
    musculotendonProperties   = tmp.tunedRabbitPsoasFibril.musculotendon;
    sarcomereProperties       = tmp.tunedRabbitPsoasFibril.sarcomere;
    normMuscleCurves          = tmp.tunedRabbitPsoasFibril.curves;
    fitting                   = tmp.tunedRabbitPsoasFibril.fitting;    
    outputFileEndingVexat = 'Linear-Titin';
    
end


assert(flag_useTwoSidedTitinCurves==normMuscleCurves.useTwoSidedTitinCurves,...
       'Error: curves struct does not contain the desired sided curves');


normMuscleCurves.useCalibratedCurves    = flag_useCalibratedVexatCurves;
normMuscleCurves.useTwoSidedTitinCurves = flag_useTwoSidedTitinCurves;

normTendonDampingConstant = ...
    musculotendonProperties.normTendonDampingConstant;
normTendonDampingLinear = ...
    musculotendonProperties.normTendonDampingLinear;

figNameGainPhase = 'Fig12';
if(flag_useFig3KirchBoskovRymer1994==1)
  figNameGainPhase = 'Fig3';  
end


%%
% Meta configuration properties: Do not touch.  
%%

dataFolder = [projectFolders.experiments_LJH2010,filesep];
structsFolder= [projectFolders.output_structs_LJH2010,filesep];
plotFolder = [projectFolders.output_plots_LJH2010,filesep];


%%
% Create the text tags that capture the configuration of the model
% and simulation: used the automatically generate unique informative file
% names based on the properties of the muscle
%
%%

%%
% Retreive the experiment information
%%

fileFig2Inset = fullfile(projectFolders.experiments_LJH2010,...
                        'data','fig_LeonardJoumaaHerzog2010_Fig2_Inset.csv');
fileFig2Main  = fullfile(projectFolders.experiments_LJH2010,...
                         'data','fig_LeonardJoumaaHerzog2010_Fig2_Main.csv');

%fileFig2Inset = ['experiments/LeonardJoumaaHerzog2010/data/',...
%                 'fig_LeonardJoumaaHerzog2010_Fig2_Inset.csv'];
%fileFig2Main  = ['experiments/LeonardJoumaaHerzog2010/data/',...
%                 'fig_LeonardJoumaaHerzog2010_Fig2_Main.csv'];

dataLJH2010Fig2Inset = ...
    loadDigitizedData( fileFig2Inset,...
        'Sarcomere Length (um)','Stress (kPa)',...
        {'Active LB','Active UB','Passive LB','Passive UB'},'');
                      
dataLJH2010Fig2Main = ...
    loadDigitizedData( fileFig2Main,...
        'Sarcomere Length (um)','Stress (kPa)',...
        {'Active Len. Std','Active Force Std',...
         'Passive Len. Std','Passive Force Std'},'');

optSarcomereLength = 2.2;
optStress = 200;

for i=1:1:length(dataLJH2010Fig2Inset)
  dataLJH2010Fig2Inset(i).x = dataLJH2010Fig2Inset(i).x./optSarcomereLength;
  dataLJH2010Fig2Inset(i).y = dataLJH2010Fig2Inset(i).y./optStress;
end
for i=1:1:length(dataLJH2010Fig2Main)
  dataLJH2010Fig2Main(i).x = dataLJH2010Fig2Main(i).x./optSarcomereLength;
  dataLJH2010Fig2Main(i).y = dataLJH2010Fig2Main(i).y./optStress;
end


lengthStart = 1-nominalNormalizedFiberLength;
lengthEnd   = mean(dataLJH2010Fig2Main(1).x)-nominalNormalizedFiberLength;

%In the paper it is reported that the stretch is performed at 0.1 um/sec.
%Scaling this rate up to a larger fiber
timeStretch = (lengthEnd-lengthStart) / (0.1/optSarcomereLength); 

lengthRampKeyPoints = [          0, lengthStart;...
                     (timeStretch), lengthEnd];
               
stimulationKeyTimes = [          0, 1;...
                     (timeStretch), 1];

timeSpan = [0,(timeStretch+2)];

passiveForceKeyPoints = [1.0,    0;...
                         2.86,1.31];  

activeForceKeyPoints = [1.0, 1.0;...
                        3.38, 5.14];
%%
% Simulate
%%

if(flag_simulateVexatModel==1)
  [success] = runLeonardJoumaaHerzog2010SimulationsVexat(...
                  nominalNormalizedFiberLength,...
                  activeForceKeyPoints,...                  
                  passiveForceKeyPoints,...
                  timeSpan,...
                  lengthRampKeyPoints,...
                  stimulationKeyTimes,...
                  flag_useElasticTendon,...
                  musculotendonProperties,...
                  sarcomereProperties,...
                  normMuscleCurves,...
                  outputFileEndingVexat, ...
                  structsFolder,...
                  flag_simulateActiveStretch,...
                  flag_simulatePassiveStretch,...
                  flag_useOctave);
end

if(flag_simulateHillModel==1)

  [success] = runLeonardJoumaaHerzog2010SimulationsDampedEquilibrium(...
                nominalNormalizedFiberLength,...
                passiveForceKeyPoints,...
                timeSpan,...
                lengthRampKeyPoints,...
                stimulationKeyTimes,...
                flag_useElasticTendon,...
                flag_useFiberDamping,...
                fiberDampingCoefficient,...
                musculotendonProperties,...
                sarcomereProperties,...
                normMuscleCurves,...
                outputFileEndingHill, ...
                structsFolder,...
                flag_simulateActiveStretch,...
                flag_simulatePassiveStretch,...
                flag_useOctave);
end

               







%%
% Plot the results
%%
if(flag_plotData == 1)



  figure(figLeonardJoumaaHerzog2010Fig2Comparision);

  subPlotSquareA = subPlotSquare;
  subPlotSquareB = subPlotSquareA;
  subPlotSquareB(1,2) = subPlotSquareB(1,2) - subPlotSquareB(1,4)*1.3;  
  subPlotSquareC = subPlotSquareB;
  subPlotSquareC(1,2) = subPlotSquareC(1,2) - subPlotSquareC(1,4)*1.3;  

  
  dataVexatLinearTitin = load([structsFolder,'benchRecordVexat_Linear-Titin.mat']);  
  dataVexatWLCTitin = load([structsFolder,'benchRecordVexat_WLC-Titin']);
  dataDampedEq = load([structsFolder,'benchRecordHill_Default.mat']);

  expActiveColor  = [1,1,1].*0.75;
  expPassiveColor = [1,1,1].*0.75;

  falSample = calcBezierYFcnXCurveSampleVector(...
                normMuscleCurves.activeForceLengthCurve, 200, [0,2]);
  fpeSample = calcBezierYFcnXCurveSampleVector(...
                normMuscleCurves.fiberForceLengthCurve, 200, [1,2.86]);

  %%
  % Force-length relation of the 3 segments of titin
  %%
  figH = subplot('Position',subPlotSquareA);


  flag_addWLCCurve=1;
  flag_addAnnotation=1;

  [figH,ligpV_00,figpV_00,figpWLCV_00] ...
        = addWLCPlot(figH, ...
            dataVexatWLCTitin.normMuscleCurves.forceLengthIgPTitinCurve,...
            dataVexatWLCTitin.normMuscleCurves.forceLengthIgPTitinCurve.xpts(6,1),...
            dataVexatWLCTitin.sarcomereProperties.IGPContourLengthNorm,...
            mean(dataLJH2010Fig2Main(2).y),...
            [1,0,0],'IgP',flag_addWLCCurve,flag_addAnnotation); 

  [figH,lpevkV_00,fpevkV_00,fpevkWLCV_00] ...
        = addWLCPlot(figH, ...
            dataVexatWLCTitin.normMuscleCurves.forceLengthPevkTitinCurve,...
            dataVexatWLCTitin.normMuscleCurves.forceLengthPevkTitinCurve.xpts(6,1),...
            dataVexatWLCTitin.sarcomereProperties.PEVKContourLengthNorm,...
            mean(dataLJH2010Fig2Main(2).y),...
            [85/255,0,0],'PEVK',flag_addWLCCurve,flag_addAnnotation); 

  [figH,ligdV_00,figdV_00,figdWLCV_00] ...
        = addWLCPlot(figH, ...
            dataVexatWLCTitin.normMuscleCurves.forceLengthIgDTitinCurve,...
            dataVexatWLCTitin.normMuscleCurves.forceLengthIgDTitinCurve.xpts(6,1),...
            dataVexatWLCTitin.sarcomereProperties.IGDFreeContourLengthNorm,...
            mean(dataLJH2010Fig2Main(2).y),...
            [1,0,0.5],'IgD',flag_addWLCCurve,flag_addAnnotation); 

  %Linear model
  flag_addWLCCurve =0;
  flag_addAnnotation=0;
  [figH,ligpV_01,figpV_01,figpWLCV_01] ...
        = addWLCPlot(figH, ...
            dataVexatLinearTitin.normMuscleCurves.forceLengthIgPTitinCurve,...
            dataVexatLinearTitin.normMuscleCurves.forceLengthIgPTitinCurve.xpts(6,1),...
            dataVexatLinearTitin.sarcomereProperties.IGPContourLengthNorm,...
            mean(dataLJH2010Fig2Main(2).y),...
            [1,0,0],'IgP-Lin',flag_addWLCCurve,flag_addAnnotation); 

  [figH,lpevkV_01,fpevkV_01,fpevkWLCV_01] ...
        = addWLCPlot(figH, ...
            dataVexatLinearTitin.normMuscleCurves.forceLengthPevkTitinCurve,...
            dataVexatLinearTitin.normMuscleCurves.forceLengthPevkTitinCurve.xpts(6,1),...
            dataVexatLinearTitin.sarcomereProperties.PEVKContourLengthNorm,...
            mean(dataLJH2010Fig2Main(2).y),...
            [85/255,0,0],'PEVK-Lin',flag_addWLCCurve,flag_addAnnotation); 

  [figH,ligdV_01,figdV_01,figdWLCV_01] ...
        = addWLCPlot(figH, ...
            dataVexatLinearTitin.normMuscleCurves.forceLengthIgDTitinCurve,...
            dataVexatLinearTitin.normMuscleCurves.forceLengthIgDTitinCurve.xpts(6,1),...
            dataVexatLinearTitin.sarcomereProperties.IGDFreeContourLengthNorm,...
            mean(dataLJH2010Fig2Main(2).y),...
            [1,0,0.5],'IgD-Lin',flag_addWLCCurve,flag_addAnnotation); 


  contourLengths = ...
      [dataVexatWLCTitin.sarcomereProperties.IGPContourLengthNorm;...
       dataVexatWLCTitin.sarcomereProperties.PEVKContourLengthNorm;...
       dataVexatWLCTitin.sarcomereProperties.IGDFreeContourLengthNorm];
  contourLengths=sort(contourLengths);

  %Just double check that these two titin models really are the same
  assert(abs(dataVexatWLCTitin.sarcomereProperties.IGPContourLengthNorm ...
            -dataVexatLinearTitin.sarcomereProperties.IGPContourLengthNorm) ...
            < 1e-6);
  assert(abs(dataVexatWLCTitin.sarcomereProperties.PEVKContourLengthNorm ...
            -dataVexatLinearTitin.sarcomereProperties.PEVKContourLengthNorm) ...
            < 1e-6);
  assert(abs(dataVexatWLCTitin.sarcomereProperties.IGDFreeContourLengthNorm ...
            -dataVexatLinearTitin.sarcomereProperties.IGDFreeContourLengthNorm) ...
            < 1e-6);



  xticks(round(contourLengths,3));
  
  titinForceFraction = 1-dataVexatWLCTitin.sarcomereProperties.extraCellularMatrixPassiveForceFraction;

  yticks(round([0,titinForceFraction,mean(dataLJH2010Fig2Main(2).y)],2));


  subplot('Position',subPlotSquareB);
  %  fill(fpeSample.x,...
  %        fpeSample.y.*scaleHillFpe,...
  %        [1,1,1].*0.85,'EdgeColor','none');
  %  hold on;

  fill( falSample.x,...
        falSample.y,...
        [1,1,1].*0.75,'EdgeColor','none');
  hold on;


  
  
  idx=1;
  lineExpActive = ...
    plot( dataLJH2010Fig2Main(idx).x,...
          dataLJH2010Fig2Main(idx).y,...
          'Color',expActiveColor,'LineWidth',2);
  hold on;
  idx=2;
  plot( dataLJH2010Fig2Main(idx).x,...
        dataLJH2010Fig2Main(idx).y,...
        'Color',expActiveColor,'LineWidth',2);
  
  text( mean(dataLJH2010Fig2Main(idx).x),...
        max(dataLJH2010Fig2Main(idx).y),...
        'Active Lengthening','FontSize',8,...
        'VerticalAlignment','bottom',...
        'HorizontalAlignment','center');
  hold on;
      
  idx=3;
  lineExpPassive = ...
    plot( dataLJH2010Fig2Main(idx).x,...
          dataLJH2010Fig2Main(idx).y,...
          '-','Color',expActiveColor,'LineWidth',2);
  hold on;
  idx=4;
  plot( dataLJH2010Fig2Main(idx).x,...
        dataLJH2010Fig2Main(idx).y,...
        '-','Color',expActiveColor,'LineWidth',2);
      
  text( mean(dataLJH2010Fig2Main(idx).x),...
        min(dataLJH2010Fig2Main(idx).y),...
        'Passive Lengthening','FontSize',8,...
        'VerticalAlignment','top',...
        'HorizontalAlignment','left');
  hold on;
      
  lineColorVexatWLCAct = [1,0,1];
  lineColorVexatWLCPas = [1,0,1].*0.5 + [1,0.75,1].*0.5;
      
  lineColorVexatAct = [0,0,1];
  lineColorVexatPas = [0,0,1].*0.5 + [0.75,0.75,1].*0.5;
  lineColorDampedEqAct = [1,0,0];
  lineColorDampedEqPas = [1,0,0].*0.5 + [0.75,0.75,1].*0.5;  
  
  lineVexatPas =plot(dataVexatLinearTitin.benchRecord.normFiberLength(:,2),...
                      dataVexatLinearTitin.benchRecord.normFiberForce(:,2),...
                      '-','Color',lineColorVexatPas,'LineWidth',1.0);
  hold on;
  
  if(isempty(dataVexatWLCTitin)==0)
      lineVexatWLCPas =plot(dataVexatWLCTitin.benchRecord.normFiberLength(:,2),...
                               dataVexatWLCTitin.benchRecord.normFiberForce(:,2),...
                          '-','Color',lineColorVexatWLCPas,'LineWidth',0.5);
      hold on;        
  end
  lineDampedEqPas =plot(dataDampedEq.benchRecord.normFiberLength(:,2),...
                        dataDampedEq.benchRecord.normFiberForce(:,2),...
                        '-','Color',lineColorDampedEqPas,'LineWidth',0.5);
  hold on;

  
  lineVexatAct =plot(dataVexatLinearTitin.benchRecord.normFiberLength(:,1),...
                      dataVexatLinearTitin.benchRecord.normFiberForce(:,1),...
                      'Color',lineColorVexatAct,'LineWidth',1);
  hold on;
  if(isempty(dataVexatWLCTitin)==0)
      lineVexatWLCAct =plot(dataVexatWLCTitin.benchRecord.normFiberLength(:,1),...
                          dataVexatWLCTitin.benchRecord.normFiberForce(:,1),...
                          'Color',lineColorVexatWLCAct,'LineWidth',1);
      hold on;
  end

  n = length(dataVexatLinearTitin.benchRecord.time);

  if(isempty(dataVexatWLCTitin)==0)
      assert(n==length(dataVexatWLCTitin.benchRecord.time));
      %Annotate the length at which the titin tip passes over the end
      %of actin
      lActN = sarcomereProperties.normActinLength;
      lTitinActinN = reshape(dataVexatLinearTitin.benchRecord.extra(:,1,1),n,1) ...
          +  sarcomereProperties.ZLineToT12NormLengthAtOptimalFiberLength;
      lTitinActinWLCN = reshape(dataVexatWLCTitin.benchRecord.extra(:,1,1),n,1) ...
          +  sarcomereProperties.ZLineToT12NormLengthAtOptimalFiberLength;

      idxSO = find( lActN-lTitinActinN <= 0, 1 );
      idxSO=idxSO-10; 
      %The slip off is smoothed - I'm going back to get an index before 
      %the slip off starts to reduce titin's force

      idxSOWLC = find( lActN-lTitinActinWLCN <= 0, 1 );
      idxSOWLC=idxSOWLC-10;

      lceNSO = dataVexatLinearTitin.benchRecord.normFiberLength(idxSO,1);
      fceNSO = dataVexatLinearTitin.benchRecord.normFiberForce(idxSO,1);

      lceNSOWLC = dataVexatWLCTitin.benchRecord.normFiberLength(idxSOWLC,1);
      fceNSOWLC = dataVexatWLCTitin.benchRecord.normFiberForce(idxSOWLC,1);

      plot(lceNSO,fceNSO,'*','MarkerSize',5,'Color',[0,0,0]);
      hold on;
      plot(lceNSOWLC,fceNSOWLC,'*','MarkerSize',5,'Color',[0,0,0]);
      hold on;

      text(lceNSOWLC*0.8,fceNSOWLC,'Titin attachment',...
          'HorizontalAlignment','right','FontSize',8);
      hold on;
      text(lceNSOWLC*0.8,fceNSOWLC*0.9,'slips off of actin',...
          'HorizontalAlignment','right','FontSize',8);
      hold on;

      xArrow = [0.8,1].*lceNSOWLC;
      yArrow = [1,1].*fceNSOWLC;
      plot(xArrow,yArrow,'Color',[0,0,0]);
      hold on;

      xArrow(1,2) = lceNSOWLC;
      yArrow(1,2) = fceNSOWLC;
      %plot(xArrow,yArrow,'Color',[0,0,0]);
      %hold on;

      %Annotate the length at which titin fails
      normLengthCEAtTitinContour = ...
         2*(sarcomereProperties.normContourLengthTitinProximal ...
         +sarcomereProperties.normContourLengthTitinDistal ...
         +sarcomereProperties.normLengthTitinFixed);

      plot([1;1].*normLengthCEAtTitinContour,...
           [0;1].*max(dataVexatWLCTitin.benchRecord.normFiberForce(:,1)),...
           '--','Color',lineColorVexatWLCAct.*0.5+[1,1,1].*0.5,'LineWidth',1);
      hold on;
      th= text( normLengthCEAtTitinContour-0.1,...
            round(mean(dataLJH2010Fig2Main(2).y),2)*0.7,...
            'Contour Length Reached',...
            'FontSize',8,...
            'HorizontalAlignment','left',...
            'Color',lineColorVexatWLCAct.*0.5+[1,1,1].*0.5);
      set(th,'Rotation',90);
      hold on;
  end

%   plot(dataDampedEq.benchRecord.normFiberLength(:,1),...
%                         dataDampedEq.benchRecord.normFiberForce(:,1),...
%                         'Color',[1,1,1],'LineWidth',1);
%   hold on;
  
  
  lineDampedEqAct =plot(dataDampedEq.benchRecord.normFiberLength(:,1),...
                        dataDampedEq.benchRecord.normFiberForce(:,1),...
                        'Color',lineColorDampedEqAct,'LineWidth',1);
  hold on;

  

  if(isempty(dataVexatWLCTitin)==0)
  
      legend([lineExpActive,lineVexatAct,lineVexatWLCAct,lineDampedEqAct], ...
              'Exp.','Model: Linear Titin','Model: WLC Titin','Hill-type ',...
              'Location','NorthWest');
  else
      legend([lineExpActive,lineVexatAct,lineDampedEqAct], ...
              'Exp.','Model: Linear Titin','Hill-type ',...
              'Location','NorthWest');

  end
  legend boxoff;
  
  title({'Simulation of Leonard, Joumaa, \& Herzog 2010','(skinned fibril)'});
  
  
  
  box off;
  
  if(isempty(dataVexatWLCTitin)==0)
      xticks(round(...
              [1,...
              1.65,...
              lceNSOWLC,...
              normLengthCEAtTitinContour,...
              mean(dataLJH2010Fig2Main(3).x),...         
              mean(dataLJH2010Fig2Main(1).x)], ...
              2));
  else
      xticks(round(...
              [1,...
              1.65,...
              mean(dataLJH2010Fig2Main(3).x),...         
              mean(dataLJH2010Fig2Main(1).x)], ...
              2));
  end

    if(isempty(dataVexatWLCTitin)==0)
    
      
      yticks(round(...
              [0, ...
              1,...
              mean(dataLJH2010Fig2Main(4).y),...
              fceNSOWLC,...
              mean(dataLJH2010Fig2Main(2).y)],...
              2));      
    else
      yticks(round(...
              [0, ...
              1,...
              mean(dataLJH2010Fig2Main(4).y),...
              mean(dataLJH2010Fig2Main(2).y)],...
              2));
    end        
        
  xlim([0, max(dataLJH2010Fig2Main(1).x*1.1)]);
  ylim([0, max(dataLJH2010Fig2Main(2).y*1.1)]);
   
  xlabel('Norm. Length ($\ell^{M}/\ell^{M}_{\circ}$)');    
  ylabel('Norm. Force ($f^{M}/ f^{M}_{\circ}$)');
  

  
  subplot('Position',subPlotSquareC);
  
  
  lceOpt = musculotendonProperties.optimalFiberLength;
  lenIGp = dataVexatLinearTitin.benchRecord.state(:,1,3)./lceOpt;
  lce    = dataVexatLinearTitin.benchRecord.normFiberLength(:,1);
  lenIGd = 0.5.*lce - lenIGp;
  
  plot( dataVexatLinearTitin.benchRecord.normFiberLength(:,1),...
        lenIGd,...
        'Color',[0,0,1],'LineWidth',1);
  xlabel('Norm. Length ($\ell^{M}/\ell^{M}_{\circ}$)');    
  ylabel('Norm. Length ($\ell^{IGd}/\ell^{M}_{\circ}$)');    
  
  box off;
  
%maximumNormalizedFiberVelocity  = 3.5*2; 
%maximumPennationAngle           = 85*(pi/180); 
%tendonStrainAtOneNormForce      = 0.049;

  if(flag_savePlotsToFile==1)
    set(figLeonardJoumaaHerzog2010Fig2Comparision,'Units','centimeters',...
    'PaperUnits','centimeters',...
    'PaperSize',[pageWidth pageHeight],...
    'PaperPositionMode','manual',...
    'PaperPosition',[0 0 pageWidth pageHeight]);     
    %set(findall(figList(i).h,'-property','FontSize'),'FontSize',10);     
    set(figLeonardJoumaaHerzog2010Fig2Comparision,'renderer','painters');     
    set(gcf,'InvertHardCopy','off')

    print('-dpdf', [plotFolder,'fig_Pub_LeonardJoumaaHerzog2010.pdf']);   
  end

  %% Titin's segments
  % Have a look at titin's segments and see how close they come to critical
  % points:
  %   - Titin-Actin attachment slipping off of actin?
  %   - Prox. segment vs. contour length
  %   - Dist. segment vs. contour length

  if(exist('figTitinSegments','builtin')==0)
    figTitinSegments=figure;
  end
  figure(figTitinSegments);



  lActN = sarcomereProperties.normActinLength;
  lTitinActinN = reshape(dataVexatLinearTitin.benchRecord.extra(:,1,1),n,1) ...
      +  sarcomereProperties.ZLineToT12NormLengthAtOptimalFiberLength;
  lTitinActinNWLC = reshape(dataVexatWLCTitin.benchRecord.extra(:,1,1),n,1) ...
      +  sarcomereProperties.ZLineToT12NormLengthAtOptimalFiberLength;


  subplot(2,3,1);
      plot(dataVexatLinearTitin.benchRecord.time,lTitinActinN,...
          'Color',[0,0,1],'LineWidth',1);
      hold on;
      plot(dataVexatWLCTitin.benchRecord.time,lTitinActinNWLC,...
          'Color',[1,0,1],'LineWidth',1);
      hold on;
    
      timeMin = min(dataVexatLinearTitin.benchRecord.time);
      timeMax = max(dataVexatLinearTitin.benchRecord.time);
      plot([timeMin;timeMax],[lActN;lActN],'--','Color',[1,0,0]);
      hold on;
    
      xlabel('Time (s)');
      ylabel('Norm. Length');
      title('Titin-Actin Attachment vs. Actin Tip');
    
      box off;

  subplot(2,3,2);
    
      plot(dataVexatLinearTitin.benchRecord.time,...
          reshape(dataVexatLinearTitin.benchRecord.extra(:,1,1),n,1),...
          'Color',[0,0,1],'LineWidth',1,'DisplayName',...
          [dataVexatLinearTitin.benchRecord.extraLabels{1}]);
      hold on;
      plot(dataVexatWLCTitin.benchRecord.time,...
           reshape(dataVexatWLCTitin.benchRecord.extra(:,1,1),n,1),...
          'Color',[1,0,1],'LineWidth',1,'DisplayName',...
          [dataVexatWLCTitin.benchRecord.extraLabels{1}]);
      hold on;
      plot([timeMin;timeMax],...
           [1;1].*sarcomereProperties.normContourLengthTitinProximal,...
           '--','Color',[1,0,0]);
      hold on;
      legend;
      legend boxoff;
      box off;
    
      xlabel('Time (s)');
      ylabel('Norm. Length');
      title('Prox. Titin Seg. vs Contour Length');
  

  subplot(2,3,3);
      plot(dataVexatLinearTitin.benchRecord.time,...
           reshape(dataVexatLinearTitin.benchRecord.extra(:,1,5),n,1),...
          'Color',[0,0,1],'LineWidth',1,'DisplayName',...
          [dataVexatLinearTitin.benchRecord.extraLabels{5}]);
      hold on;
      plot(dataVexatWLCTitin.benchRecord.time,...
           reshape(dataVexatWLCTitin.benchRecord.extra(:,1,5),n,1),...
          'Color',[1,0,1],'LineWidth',1,'DisplayName',...
          [dataVexatWLCTitin.benchRecord.extraLabels{5}]);
      hold on;
      plot([timeMin;timeMax],...
           [1;1].*sarcomereProperties.normContourLengthTitinDistal,...
           '--','Color',[1,0,0]);
      hold on;
      box off;

      legend;
      legend boxoff;      
    
      xlabel('Time (s)');
      ylabel('Norm. Length');
      title('Distal. Titin Seg. vs Contour Length');


  subplot(2,3,4);
      plot(reshape(dataVexatLinearTitin.benchRecord.extra(:,1,1),n,1),...
           reshape(dataVexatLinearTitin.benchRecord.extra(:,1,2),n,1),...
          'Color',[0,0,1],'LineWidth',1,'DisplayName',...
          [dataVexatLinearTitin.benchRecord.extraLabels{1},...
          '-',dataVexatLinearTitin.benchRecord.extraLabels{2}]);
      hold on;
      plot(reshape(dataVexatWLCTitin.benchRecord.extra(:,1,1),n,1),...
           reshape(dataVexatWLCTitin.benchRecord.extra(:,1,2),n,1),...
          'Color',[1,0,1],'LineWidth',1,'DisplayName',...
          [dataVexatWLCTitin.benchRecord.extraLabels{1},...
          '-',dataVexatWLCTitin.benchRecord.extraLabels{2}]);
      hold on;  
      box off;
      legend;
      legend boxoff;      
      xlabel('Norm. Length');
      ylabel('Norm. Force')
      title('Proximal Titin Seg');
  
  subplot(2,3,5);
      plot(reshape(dataVexatLinearTitin.benchRecord.extra(:,1,5),n,1),...
           reshape(dataVexatLinearTitin.benchRecord.extra(:,1,6),n,1),...
          'Color',[0,0,1],'LineWidth',1,'DisplayName',...
          [dataVexatLinearTitin.benchRecord.extraLabels{5},...
          '-',dataVexatLinearTitin.benchRecord.extraLabels{6}]);
      hold on;
      plot(reshape(dataVexatWLCTitin.benchRecord.extra(:,1,5),n,1),...
           reshape(dataVexatWLCTitin.benchRecord.extra(:,1,6),n,1),...
          'Color',[1,0,1],'LineWidth',1,'DisplayName',...
          [dataVexatWLCTitin.benchRecord.extraLabels{5},...
          '-',dataVexatWLCTitin.benchRecord.extraLabels{6}]);
      hold on;
      legend;
      legend boxoff;      
      box off;
      xlabel('Norm. Length');
      ylabel('Norm. Force')
      title('Distal Titin Seg');

  here=1;

  %dataVexatLinearTitin.benchRecord.normFiberLength(:,1),...
  %dataVexatLinearTitin.benchRecord.normFiberForce(:,1)
  %dataVexatWLCTitin.benchRecord.normFiberLength(:,1),...
  %dataVexatWLCTitin.benchRecord.normFiberForce(:,1)

  %subplot(1,3,1);

  %subplot(1,2,1);
  %plot()

end




