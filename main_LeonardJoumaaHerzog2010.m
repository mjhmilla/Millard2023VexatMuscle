%Replication of
%Leonard TR, Joumaa V, Herzog W. An activatable molecular spring reduces 
%muscle tearing during extreme stretching. Journal of biomechanics. 
%2010 Nov 16;43(15):3063-6.

flag_outerLoopMode = 1;

if(flag_outerLoopMode == 0)
  clc;
  close all;
  clear all;

  flag_simulateHillModel            = 0; 
  flag_simulateOpus31Model          = 0;

  flag_useCalibratedOpus31Curves    = 1;
  flag_useTwoSidedTitinCurves       = 1;

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
  %tunedOpus31Results = 'experiments/LeonardJoumaaHerzog2010/benchRecordOpus31_RigidTendon_K44p31D0p50Tau_LJH2010__TiAD1000p00_TiPD1p50_NomLen1p00_90Hz_TiAdj.mat';

  flag_useFig3KirchBoskovRymer1994              = 0; 
  flag_useElasticTendon                         = 0; 
  flag_useFiberDamping                          = 1;
  fiberDampingCoefficient = 0.01;
   

  
end 

assert(flag_useElasticTendon ==0,...
       'This is a skinned fibril experiment: there is no tendon!');


parametersDirectoryTree       = genpath('parameters');
curvesDirectoryTree           = genpath('curves');
experimentsDirectoryTree      = genpath('experiments');
simulationDirectoryTree       = genpath('simulation');
modelDirectoryTree            = genpath('models');
postprocessingDirectoryTree   = genpath('postprocessing');

addpath(parametersDirectoryTree       );
addpath(curvesDirectoryTree           );
addpath(experimentsDirectoryTree      );
addpath(simulationDirectoryTree       );
addpath(modelDirectoryTree            );
addpath(postprocessingDirectoryTree   );

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

%Basic parameters for the Hill model
load('output/structs/defaultFelineSoleus.mat')
musculotendonProperties   = defaultFelineSoleus.musculotendon;
sarcomereProperties       = defaultFelineSoleus.sarcomere;
normMuscleCurves          = defaultFelineSoleus.curves;

normMuscleCurves.useCalibratedCurves    = flag_useCalibratedOpus31Curves;
normMuscleCurves.useTwoSidedTitinCurves = flag_useTwoSidedTitinCurves;

normTendonDampingConstant = ...
    musculotendonProperties.normTendonDampingConstant;
normTendonDampingLinear = ...
    musculotendonProperties.normTendonDampingLinear;

figNameGainPhase = 'Fig12';
if(flag_useFig3KirchBoskovRymer1994==1)
  figNameGainPhase = 'Fig3';  
end

%Basic parameters + extra parameters needed for Opus 31: 
% cross-bridge stiffness and damping 
% tendon damping
% titin's multiple segments

load(['output/structs/felineSoleusRigidTendonKBR1994',figNameGainPhase,'.mat']);
musculotendonPropertiesOpus31_RT = felineSoleusRigidTendonKBR1994.musculotendon;
sarcomerePropertiesOpus31_RT     = felineSoleusRigidTendonKBR1994.sarcomere;
normMuscleCurves_RT              = felineSoleusRigidTendonKBR1994.curves;


sarcomerePropertiesOpus31_ET      = [];
musculotendonPropertiesOpus31_ET  = [];


if(flag_useElasticTendon==1)
  assert(flag_useElasticTendon ==0,...
         'This is a skinned fibril experiment: there is no tendon!');  
  %sarcomerePropertiesOpus31     = sarcomerePropertiesOpus31_ET;
  %musculotendonPropertiesOpus31 = musculotendonPropertiesOpus31_ET;
else
  sarcomerePropertiesOpus31     = sarcomerePropertiesOpus31_RT;  
  musculotendonPropertiesOpus31 = musculotendonPropertiesOpus31_RT;
  normMuscleCurvesOpus31 = normMuscleCurves_RT;
end

normMuscleCurvesOpus31.useCalibratedCurves    = flag_useCalibratedOpus31Curves;
normMuscleCurvesOpus31.useTwoSidedTitinCurves = flag_useTwoSidedTitinCurves;

if(isempty(normPassiveTitinToActinDamping)==0)
  sarcomereProperties.normPassiveTitinToActinDamping          = normPassiveTitinToActinDamping;
  sarcomerePropertiesOpus31.normPassiveTitinToActinDamping    = normPassiveTitinToActinDamping;
  sarcomerePropertiesOpus31_RT.normPassiveTitinToActinDamping = normPassiveTitinToActinDamping;
else
  normPassiveTitinToActinDamping = sarcomereProperties.normPassiveTitinToActinDamping;
end

if(isempty(normActiveTitinToActinDamping)==0)
  sarcomereProperties.normMaxActiveTitinToActinDamping          = normActiveTitinToActinDamping;
  sarcomerePropertiesOpus31.normMaxActiveTitinToActinDamping    = normActiveTitinToActinDamping;
  sarcomerePropertiesOpus31_RT.normMaxActiveTitinToActinDamping = normActiveTitinToActinDamping;
else
  normActiveTitinToActinDamping = sarcomereProperties.normMaxActiveTitinToActinDamping;
end

%%
% Meta configuration properties: Do not touch.  
%%

dataFolder = 'experiments/LeonardJoumaaHerzog2010/';
plotFolder = 'output/plots/LeonardJoumaaHerzog2010/';


%%
% Create the text tags that capture the configuration of the model
% and simulation: used the automatically generate unique informative file
% names based on the properties of the muscle
%
%%

generateModelSpecificKeyWords;    

tendonTag = '_ElasticTendon';
if(flag_useElasticTendon==0)
  tendonTag = '_RigidTendon';
end



outputFileEndingOpus31    = '';
outputFileEndingOpus31_RT = '';


outputFileEndingOpus31_RT = sprintf('_K%sD%sTau%s_LJH2010_%s_%s_%s_%s',...
  kScaleStr_RT, dScaleStr_RT, tScaleStr,...
  ['_TiAD',titinActiveDampingStr],...
  ['TiPD',titinPassiveDampingStr],...    
  ['NomLen',nominalNormalizedFiberLengthStr],...
  strFittingBandwidth);

outputFileEndingOpus31 = outputFileEndingOpus31_RT;

outputFileEndingHill = sprintf('_D%i_LJH2010_%s',flag_useFiberDamping);

outputFileEndingHill_RT = ...
  sprintf('_D%i_LJH2010',0);
  
%%
% Retreive the experiment information
%%
fileFig2Inset = ['experiments/LeonardJoumaaHerzog2010/data/',...
                 'fig_LeonardJoumaaHerzog2010_Fig2_Inset.csv'];
fileFig2Main  = ['experiments/LeonardJoumaaHerzog2010/data/',...
                 'fig_LeonardJoumaaHerzog2010_Fig2_Main.csv'];

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
                       timeStretch, lengthEnd];
               
stimulationKeyTimes = [          0, 1;...
                       timeStretch, 1];

timeSpan = [0,timeStretch];

passiveForceKeyPoints = [1.0,    0;...
                         2.86,1.31];  

activeForceKeyPoints = [1.0, 1.0;...
                        3.38, 5.14];
%%
% Simulate
%%

if(flag_simulateOpus31Model==1)
  [success] = runLeonardJoumaaHerzog2010SimulationsOpus31(...
                  nominalNormalizedFiberLength,...
                  activeForceKeyPoints,...                  
                  passiveForceKeyPoints,...
                  timeSpan,...
                  lengthRampKeyPoints,...
                  stimulationKeyTimes,...
                  flag_useElasticTendon,...
                  musculotendonPropertiesOpus31,...
                  sarcomerePropertiesOpus31,...
                  normMuscleCurvesOpus31,...
                  outputFileEndingOpus31, ...
                  dataFolder,...
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
                dataFolder,...
                flag_simulateActiveStretch,...
                flag_simulatePassiveStretch,...
                flag_useOctave);
end

               







%%
% Plot the results
%%
if(flag_plotData == 1)
  figure(figLeonardJoumaaHerzog2010Fig2Comparision);

  nameModification = '';
  if(flag_useElasticTendon == 1)
    nameModification = 'ElasticTendon';
  else
    nameModification = 'RigidTendon';          
  end
  
  
  fileNameOpus31 = [dataFolder,'benchRecordOpus31_',...
                    nameModification,outputFileEndingOpus31,'_TiDefault.mat'];
  dataOpus31 = load(fileNameOpus31);
  
  dataFolderContents = dir(dataFolder);

  strA = sprintf('_TiAD%1.2f',tunedNormActiveTitinToActinDamping);
  strA(strfind(strA,'.'))='p';
  tunedFileDateNum = 0;
  tunedFileName = '';

  for idxFile=1:1:length(dataFolderContents)
    if(contains(dataFolderContents(idxFile).name,strA))
        if (dataFolderContents(idxFile).datenum > tunedFileDateNum)
            tunedFileName = dataFolderContents(idxFile).name;
            tunedFileDateNum = dataFolderContents(idxFile).datenum;
        end
    end
  end
  
  dataOpus31Tuned=[];
  if(isempty(tunedFileName)==0)
    dataOpus31Tuned = load([dataFolder,tunedFileName]);
  end

  fileNameDampedEq = [dataFolder,'benchRecordHill_',...
                      nameModification,outputFileEndingHill,'.mat'];
  dataDampedEq = load(fileNameDampedEq);
  
  
  

  expActiveColor  = [1,1,1].*0.75;
  expPassiveColor = [1,1,1].*0.75;

  falSample = calcBezierYFcnXCurveSampleVector(...
                normMuscleCurves.activeForceLengthCurve, 200, [0,2]);
  fpeSample = calcBezierYFcnXCurveSampleVector(...
                normMuscleCurves.fiberForceLengthCurve, 200, [1,2.86]);


  subplot('Position',subPlotSquare);
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
      
  lineColorOpus31TunedAct = [1,0,1];
  lineColorOpus31TunedPas = [1,0,1].*0.5 + [1,0.75,1].*0.5;
      
  lineColorOpus31Act = [0,0,1];
  lineColorOpus31Pas = [0,0,1].*0.5 + [0.75,0.75,1].*0.5;
  lineColorDampedEqAct = [1,0,0];
  lineColorDampedEqPas = [1,0,0].*0.5 + [0.75,0.75,1].*0.5;  
  
  lineOpus31Pas =plot(dataOpus31.benchRecord.normFiberLength(:,2),...
                      dataOpus31.benchRecord.normFiberForce(:,2),...
                      '-','Color',lineColorOpus31Pas,'LineWidth',1.0);
  hold on;
  
  if(isempty(tunedFileName)==0)
      lineOpus31TunedPas =plot(dataOpus31Tuned.benchRecord.normFiberLength(:,2),...
                               dataOpus31Tuned.benchRecord.normFiberForce(:,2),...
                          '-','Color',lineColorOpus31TunedPas,'LineWidth',0.5);
      hold on;  
  end
  lineDampedEqPas =plot(dataDampedEq.benchRecord.normFiberLength(:,2),...
                        dataDampedEq.benchRecord.normFiberForce(:,2),...
                        '-','Color',lineColorDampedEqPas,'LineWidth',0.5);
  hold on;


%   plot(dataOpus31.benchRecord.normFiberLength(:,1),...
%                       dataOpus31.benchRecord.normFiberForce(:,1),...
%                       'Color',[1,1,1],'LineWidth',1);
%   hold on;  
  
  lineOpus31Act =plot(dataOpus31.benchRecord.normFiberLength(:,1),...
                      dataOpus31.benchRecord.normFiberForce(:,1),...
                      'Color',lineColorOpus31Act,'LineWidth',1);
  hold on;
  if(isempty(tunedFileName)==0)
      lineOpus31TunedAct =plot(dataOpus31Tuned.benchRecord.normFiberLength(:,1),...
                          dataOpus31Tuned.benchRecord.normFiberForce(:,1),...
                          'Color',lineColorOpus31TunedAct,'LineWidth',1);
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

  

  if(isempty(tunedFileName)==0)
  
      legend([lineExpActive,lineOpus31Act,lineOpus31TunedAct,lineDampedEqAct], ...
              'Exp.','Model','Model Adjusted','Hill-type ',...
              'Location','NorthWest');
  else
      legend([lineExpActive,lineOpus31Act,lineDampedEqAct], ...
              'Exp.','Model','Hill-type ',...
              'Location','NorthWest');

  end
  legend boxoff;
  
  title({'Simulation of Leonard, Joumaa, \& Herzog 2010','(skinned fibril)'});
  
  
  
  box off;
  
  
  xticks(round(...
          [1,...
          1.65,...
          mean(dataLJH2010Fig2Main(3).x),...         
          mean(dataLJH2010Fig2Main(1).x)], ...
          2));

  yticks(round(...
          [0, ...
          1,...
          mean(dataLJH2010Fig2Main(4).y),...
          mean(dataLJH2010Fig2Main(2).y)],...
          2));
        
        
  xlim([0, max(dataLJH2010Fig2Main(1).x*1.1)]);
  ylim([0, max(dataLJH2010Fig2Main(2).y*1.1)]);
   
  xlabel('Norm. Length ($\ell^{M}/\ell^{M}_{\circ}$)');    
  ylabel('Norm. Force ($f^{M}/ f^{M}_{\circ}$)');
  
  subPlotSquareBelow = subPlotSquare;
  subPlotSquareBelow(1,2) = subPlotSquareBelow(1,2) - subPlotSquare(1,4)*1.2;
  
  subplot('Position',subPlotSquareBelow);
  
  
  lceOpt = musculotendonProperties.optimalFiberLength;
  lenIGp = dataOpus31.benchRecord.state(:,1,3)./lceOpt;
  lce    = dataOpus31.benchRecord.normFiberLength(:,1);
  lenIGd = 0.5.*lce - lenIGp;
  
  plot( dataOpus31.benchRecord.normFiberLength(:,1),...
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
end




