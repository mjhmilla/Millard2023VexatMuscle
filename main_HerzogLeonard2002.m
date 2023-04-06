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

%
% A simulation of Herzog & Leonard's 2002 experiment using matlab 
% implementations of the proposed model and a Hill model (the Matlab
% version of the Millard2012 model that is a part of OpenSim)
%
%   Proposed model: rigid and elastic tendon versions
%   Hill model    : rigid and elastic tendon versions 
%
% Herzog W, Leonard TR. Force enhancement following stretching of skeletal 
% muscle: a new mechanism. Journal of Experimental Biology. 
% 2002 May 1;205(9):1275-83.
%
flag_outerLoopMode = 1;

rootDir         = getRootProjectDirectory();
projectFolders  = getProjectFolders(rootDir);

if(flag_outerLoopMode == 0)% && flag_buildCombinedPlot == 1)
  clc;
  close all;
  clear all;

  rootDir         = getRootProjectDirectory();
  projectFolders  = getProjectFolders(rootDir);

  figDescendingCombined = figure;
  flag_buildCombinedPlot = 1;

  
  flag_simulateHillModel                        = 0; 
  flag_simulateVexatModel                      = 0;
  flag_useCalibratedVexatCurves                = 1;
  flag_useTitinCurvesWithRigidIgDSegment        = 0;
  flag_useTwoSidedTitinCurves                   = 0;

  flag_useFig3KirchBoskovRymer1994              = 0; 
  flag_useElasticTendon                         = 0;  
  flag_useFiberDamping                          = 1;
  
  flag_simulateActiveStretch  = 0;
  flag_simulatePassiveStretch = 0;
  flag_simulateStatic         = 0;
  

  flag_plotData = 1;
  flag_savePlotsToFile=1;
  pubOutputFolder = '';       
  

  figureNumber       = 7;
  subFigureNumber    = 1;
  trialNumber        = 2;  


  fileSimParamsHL2002 = fullfile( projectFolders.experiments_HL2002,... 
                                  'simulationParametersHerzogLeonard2002.csv');

  simParamsHL2002 = csvread(fileSimParamsHL2002);
  nominalNormalizedFiberLength = simParamsHL2002(1,1);



  scaleMaximumIsometricForce     = [];
  scaleOptimalFiberLength        = [];%[]: take the default. Prev. 1.0;
  scaleSlidingTimeConstant       = [];%[]: take the default. Prev. 0.2;
  normActiveTitinToActinDamping  = [];%[]: take the default. Prev. 6; 
  normPassiveTitinToActinDamping = [];%[]: take the default. Prev. 1.5; 
  
   
  figHerzogLeonard2002Fig7Comparision = figure;
  
  flag_testVexatDerivativeFunction = 1;
  flag_useOctave                    = 0;
%else
%  flag_useElasticTendon                         = 1;
end


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

%subPlotPanel(1,1,4) = 0.5*subPlotPanel(1,1,4);                %height
%subPlotPanel(1,1,2) = subPlotPanel(1,1,2)-subPlotPanel(1,1,4);%y top left


flag_runSimulations = (flag_simulateActiveStretch ...
                  ||  flag_simulatePassiveStretch ...
                  ||  flag_simulateStatic);

%%
% Configure the muscle models:  load the configurations generated from 
%                               main_KirschBoskovRymer1994.m
%%

%Basic parameters for the Hill model
fileDefaultFelineSoleus = fullfile( projectFolders.output_structs_FittedModels,...
                                    'defaultFelineSoleus.mat');

load(fileDefaultFelineSoleus);
musculotendonProperties   = defaultFelineSoleus.musculotendon;
sarcomereProperties       = defaultFelineSoleus.sarcomere;
normMuscleCurves          = defaultFelineSoleus.curves;

normMuscleCurves.useCalibratedCurves = flag_useCalibratedVexatCurves;
normMuscleCurves.useTitinCurvesWithRigidIgDSegment=...
    flag_useTitinCurvesWithRigidIgDSegment;
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

fileFittedFelineSoleusRT = fullfile(projectFolders.output_structs_FittedModels,...
    ['fittedFelineSoleusHL2002KBR1994',figNameGainPhase,'_RT.mat']);

fileFittedFelineSoleusET = fullfile(projectFolders.output_structs_FittedModels,...
    ['fittedFelineSoleusHL2002KBR1994',figNameGainPhase,'_ET.mat']);


tmp=load(fileFittedFelineSoleusRT);
musculotendonPropertiesVexat_RT = tmp.fittedFelineSoleus.musculotendon;
sarcomerePropertiesVexat_RT     = tmp.fittedFelineSoleus.sarcomere;
normMuscleCurvesVexat_RT        = tmp.fittedFelineSoleus.curves;
fittingVexat_RT                 = tmp.fittedFelineSoleus.fitting;


tmp=load(fileFittedFelineSoleusET);
musculotendonPropertiesVexat_ET = tmp.fittedFelineSoleus.musculotendon;
sarcomerePropertiesVexat_ET     = tmp.fittedFelineSoleus.sarcomere;
normMuscleCurvesVexat_ET        = tmp.fittedFelineSoleus.curves;
fittingVexat_ET                 = tmp.fittedFelineSoleus.fitting;

sarcomerePropertiesVexat     = [];
musculotendonPropertiesVexat = [];

if(flag_useElasticTendon==1)
  sarcomerePropertiesVexat     = sarcomerePropertiesVexat_ET;
  musculotendonPropertiesVexat = musculotendonPropertiesVexat_ET;
  normMuscleCurvesVexat        = normMuscleCurvesVexat_ET;
  fitting = fittingVexat_ET;
else
  sarcomerePropertiesVexat     = sarcomerePropertiesVexat_RT;  
  musculotendonPropertiesVexat = musculotendonPropertiesVexat_RT;
  normMuscleCurvesVexat        = normMuscleCurvesVexat_RT;  
  fitting = fittingVexat_RT;
  
end

assert(flag_useTwoSidedTitinCurves==normMuscleCurvesVexat.useTwoSidedTitinCurves,...
       'Error: curves struct does not contain the desired sided curves');

normMuscleCurvesVexat.useCalibratedCurves = flag_useCalibratedVexatCurves;
normMuscleCurvesVexat.useTitinCurvesWithRigidIgDSegment=...
    flag_useTitinCurvesWithRigidIgDSegment;
normMuscleCurvesVexat.useTwoSidedTitinCurves = flag_useTwoSidedTitinCurves;

if(isempty(scaleOptimalFiberLength)==0)
  musculotendonProperties.optimalFiberLength = ...
    musculotendonProperties.optimalFiberLength*scaleOptimalFiberLength;
  musculotendonPropertiesVexat.optimalFiberLength = ...
    musculotendonPropertiesVexat.optimalFiberLength*scaleOptimalFiberLength;
  musculotendonPropertiesVexat_RT.optimalFiberLength = ...
    musculotendonPropertiesVexat_RT.optimalFiberLength*scaleOptimalFiberLength;
  musculotendonPropertiesVexat_ET.optimalFiberLength = ...
    musculotendonPropertiesVexat_ET.fiso*scaleOptimalFiberLength;
else
  scaleOptimalFiberLength=1;
end

if(isempty(scaleMaximumIsometricForce)==0)
  musculotendonProperties.fiso = ...
    musculotendonProperties.fiso*scaleMaximumIsometricForce;
  musculotendonPropertiesVexat.fiso = ...
    musculotendonPropertiesVexat.fiso*scaleMaximumIsometricForce;
  musculotendonPropertiesVexat_RT.fiso = ...
    musculotendonPropertiesVexat_RT.fiso*scaleMaximumIsometricForce;
  musculotendonPropertiesVexat_ET.fiso = ...
    musculotendonPropertiesVexat_ET.fiso*scaleMaximumIsometricForce;
else
  scaleMaximumIsometricForce = 1.0;
end

if(isempty(scaleSlidingTimeConstant)==0)

  sarcomereProperties.slidingTimeConstant = ...
    sarcomereProperties.slidingTimeConstant*scaleSlidingTimeConstant;
  sarcomerePropertiesVexat.slidingTimeConstant = ...
    sarcomerePropertiesVexat.slidingTimeConstant*scaleSlidingTimeConstant;
  sarcomerePropertiesVexat_RT.slidingTimeConstant = ...
    sarcomerePropertiesVexat_RT.slidingTimeConstant*scaleSlidingTimeConstant;
  sarcomerePropertiesVexat_ET.slidingTimeConstant = ...
    sarcomerePropertiesVexat_ET.slidingTimeConstant*scaleSlidingTimeConstant;
  
  
else
  scaleSlidingTimeConstant = 1.0;
end

if(isempty(normPassiveTitinToActinDamping)==0)
  sarcomereProperties.normPassiveTitinToActinDamping          = normPassiveTitinToActinDamping;
  sarcomerePropertiesVexat.normPassiveTitinToActinDamping    = normPassiveTitinToActinDamping;
  sarcomerePropertiesVexat_RT.normPassiveTitinToActinDamping = normPassiveTitinToActinDamping;
  sarcomerePropertiesVexat_ET.normPassiveTitinToActinDamping = normPassiveTitinToActinDamping;
else
  normPassiveTitinToActinDamping = sarcomereProperties.normPassiveTitinToActinDamping;
end

if(isempty(normActiveTitinToActinDamping)==0)
  sarcomereProperties.normMaxActiveTitinToActinDamping          = normActiveTitinToActinDamping;
  sarcomerePropertiesVexat.normMaxActiveTitinToActinDamping    = normActiveTitinToActinDamping;
  sarcomerePropertiesVexat_RT.normMaxActiveTitinToActinDamping = normActiveTitinToActinDamping;
  sarcomerePropertiesVexat_ET.normMaxActiveTitinToActinDamping = normActiveTitinToActinDamping;
else
  normActiveTitinToActinDamping = sarcomereProperties.normMaxActiveTitinToActinDamping;
end





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

subFigTitle = '';
switch subFigureNumber
  case 1
    subFigTitle = 'A';
  case 2
    subFigTitle = 'B';
  case 3
    subFigTitle = 'C';
  otherwise
    assert(0);
end

outputFileEndingVexat    = '';
outputFileEndingVexat_ET = '';
outputFileEndingVexat_RT = '';

outputFileEndingVexat_ET = sprintf('_K%sD%sTau%s_KTC%s_KTL%s_HL2002_%i%i%i_%s_%s_%s_%s',...
  kScaleStr_ET, dScaleStr_ET, tScaleStr, kTConstStr, kTLinearStr,...
  figureNumber,subFigureNumber,trialNumber,...
  ['_TiAD',titinActiveDampingStr],...
  ['TiPD',titinPassiveDampingStr],...  
  ['NomLen',nominalNormalizedFiberLengthStr],...
  strFittingBandwidth);

outputFileEndingVexat_RT = sprintf('_K%sD%sTau%s_HL2002_%i%i%i_%s_%s_%s_%s',...
  kScaleStr_RT, dScaleStr_RT, tScaleStr,...
  figureNumber,subFigureNumber,...
  trialNumber,['_TiAD',titinActiveDampingStr],...
  ['TiPD',titinPassiveDampingStr],...    
  ['NomLen',nominalNormalizedFiberLengthStr],...
  strFittingBandwidth);

if(flag_useElasticTendon==1)
  outputFileEndingVexat = outputFileEndingVexat_ET;
else
  outputFileEndingVexat = outputFileEndingVexat_RT;
end

%scaleHillFpeStr = sprintf('%1.2f',scaleHillFpe);
%z=strfind(scaleHillFpeStr,'.');
%scaleHillFpeStr = [scaleHillFpeStr(1,1:(z-1)),'p',scaleHillFpeStr(1,(z+1))];

outputFileEndingHill = sprintf('_D%i_HL2002_%i%i%i_%s',...
  flag_useFiberDamping,...%scaleHillFpeStr,...
  figureNumber,subFigureNumber,trialNumber,...
  ['NomLen',nominalNormalizedFiberLengthStr]);

outputFileEndingHill_RT = ...
  sprintf('_D%i__HL2002_%i%i%i_%s',...
    0,...%scaleHillFpeStr,...
    figureNumber,subFigureNumber,trialNumber,...
    ['NomLen',nominalNormalizedFiberLengthStr]);
  
outputFileEndingHill_ET = ...
  sprintf('_D%i_HL2002_%i%i%i_%s',...
    1,...%scaleHillFpeStr,...
    figureNumber,subFigureNumber,trialNumber,...
    ['NomLen',nominalNormalizedFiberLengthStr]);

%%
% Retreive the experiment information
%%

%maximumNormalizedFiberVelocity  = 3.5*2; 
%maximumPennationAngle           = 85*(pi/180); 
%tendonStrainAtOneNormForce      = 0.049;

expConfigHerzogLeonard2002 =...
 getHerzogLeonard2002Configuration( figureNumber,...
                                    subFigureNumber, ...
                                    trialNumber,...
                                    projectFolders);




if(flag_testVexatDerivativeFunction==1)
    
  a     = 1;
  dadt  = 0;
  
  lceOpt    = musculotendonPropertiesVexat.optimalFiberLength; 
  alphaOpt  = musculotendonPropertiesVexat.pennationAngle;
  ltSlk     = musculotendonPropertiesVexat.tendonSlackLength;
  etIso     = musculotendonPropertiesVexat.tendonStrainAtOneNormForce;
  ltIso     = ltSlk*(1+etIso);
  
  lceN    = 1;
  lceHN   = 0.5*lceN;
  lce     = lceN*lceOpt;
  
  l1 = 0.5*normMuscleCurves.forceLengthIgPTitinCurve.xEnd(1)*lceOpt;

  lx = 0;
  dlx = 0;
  la = lce*0.5-sarcomerePropertiesVexat.normMyosinHalfLength*lceOpt;
  dla = 0;


  lp = lceOpt*cos(alphaOpt) + ltIso;
  dlp= 0.1;
  


  activationState = [dadt;a];  
  pathState       = [dlp;lp];
  muscleState     = [];
  if(flag_useElasticTendon==1)
    muscleState     = [lce;dla;la;l1];    
  else
    muscleState     = [dla;la;l1];        
  end
  
  modelConfig = struct('useElasticTendon',flag_useElasticTendon,...
                       'initializeState',0,...
                       'iterMax', 100,...
                       'tol', 1e-8);
  
  mtInfo = calcMillard2023VexatMuscleInfo( activationState,...
                                            pathState,...
                                            muscleState,...
                                            musculotendonPropertiesVexat,...
                                            sarcomerePropertiesVexat,...
                                            normMuscleCurvesVexat,...
                                            modelConfig);
                                                
end

if(flag_simulateVexatModel==1 && flag_runSimulations == 1)
  %expConfigHerzogLeonard2002.timeSpan  
  %tspanAct = [0,expConfigHerzogLeonard2002.lengthRampKeyPoints(1,1)];  
  [success] = runHerzogLeonard2002SimulationsVexat(...
                            nominalNormalizedFiberLength,...
                            expConfigHerzogLeonard2002.nominalForce,...                            
                            expConfigHerzogLeonard2002.timeSpan,...
                            expConfigHerzogLeonard2002.lengthRampKeyPoints,...
                            expConfigHerzogLeonard2002.stimulationKeyTimes,...
                            flag_useElasticTendon,...
                            musculotendonPropertiesVexat,...
                            sarcomerePropertiesVexat,...
                            normMuscleCurvesVexat,...
                            outputFileEndingVexat, ...
                            projectFolders.output_structs_HL2002,...
                            flag_simulateActiveStretch,...
                            flag_simulatePassiveStretch,...
                            flag_simulateStatic,...
                            flag_useOctave);
end

if(flag_simulateHillModel==1 && flag_runSimulations==1)

  [success] = runHerzogLeonard2002SimulationsDampedEquilibrium(...
                            nominalNormalizedFiberLength,...
                            expConfigHerzogLeonard2002.nominalForce,...                            
                            expConfigHerzogLeonard2002.timeSpan,...
                            expConfigHerzogLeonard2002.lengthRampKeyPoints,...
                            expConfigHerzogLeonard2002.stimulationKeyTimes,...
                            flag_useElasticTendon,...
                            flag_useFiberDamping,...
                            musculotendonProperties,...
                            sarcomereProperties,...
                            normMuscleCurves,...
                            outputFileEndingHill, ...
                            projectFolders.output_structs_HL2002,...
                            flag_simulateActiveStretch,...
                            flag_simulatePassiveStretch,...
                            flag_simulateStatic,...
                            flag_useOctave);
end

if(flag_plotData == 1)
  nameModification = '';
  if(flag_useElasticTendon == 1)
    nameModification = 'ElasticTendon';
  else
    nameModification = 'RigidTendon';          
  end

  fileNameVexat = fullfile( projectFolders.output_structs_HL2002,...
    ['benchRecordVexat_',nameModification,outputFileEndingVexat,'.mat']);
  dataVexat = load(fileNameVexat);

  fileNameDampedEq = fullfile(projectFolders.output_structs_HL2002,...
    ['benchRecordHill_',nameModification,outputFileEndingHill,'.mat']);
  dataDampedEq = load(fileNameDampedEq);
  

  figHerzogLeonard2002Fig7Comparision = ...
    plotHerzogLeonardFig7Comparision(...
                  expConfigHerzogLeonard2002,dataVexat, dataDampedEq,  ...
                  figureNumber,subFigureNumber,trialNumber,...
                  figHerzogLeonard2002Fig7Comparision,subPlotPairPanel); 
  
  figDesLimbName= ['fig_Pub_HerzogLeonard2002_DescendingLimbStability_',...
                    num2str(figureNumber), subFigTitle,num2str(trialNumber),...
                    'RigidElastic','_',tScaleStr,'_TiAD',titinActiveDampingStr,...
                    '_TiPD',titinPassiveDampingStr,...  
                    '_NomLen',nominalNormalizedFiberLengthStr,...
                    '_',strFittingBandwidth,'.pdf'];

  if(flag_buildCombinedPlot > 0 ...
          && isempty(dataVexat.benchRecord)==0 ...
          && isempty(dataDampedEq.benchRecord)==0)
    
    
    figDescendingCombined = ...
      plotHerzogLeonardDescendingLimb(figDescendingCombined,...
                      expConfigHerzogLeonard2002,dataVexat, dataDampedEq,  ...
                      nominalNormalizedFiberLength,flag_useElasticTendon,...
                      figureNumber,subFigureNumber,trialNumber,...
                      subPlotHerzogLeonard2000Stability,...
                      projectFolders.output_plots_HL2002,...
                      figDesLimbName,pageWidth,pageHeight);
  end
  %axis square;
  
  if(flag_savePlotsToFile==1)
    
    figure(figHerzogLeonard2002Fig7Comparision);
    
    set(figHerzogLeonard2002Fig7Comparision,'Units','centimeters',...
    'PaperUnits','centimeters',...
    'PaperSize',[pageWidth pageHeight],...
    'PaperPositionMode','manual',...
    'PaperPosition',[0 0 pageWidth pageHeight]);     
    %set(findall(figList(i).h,'-property','FontSize'),'FontSize',10);     
    set(figHerzogLeonard2002Fig7Comparision,'renderer','painters');     
    set(gcf,'InvertHardCopy','off')


    filePlotName = fullfile(projectFolders.output_plots_HL2002,...
                    ['fig_Pub_HerzogLeonard2002_',...
                    num2str(figureNumber), subFigTitle,...
                    tendonTag,'_',tScaleStr,'_TiAD',titinActiveDampingStr,...
                    '_TiPD',titinPassiveDampingStr,...  
                    '_NomLen',nominalNormalizedFiberLengthStr,...
                    '_',strFittingBandwidth,'.pdf']);

    
    print('-dpdf', filePlotName);   


  end
end       