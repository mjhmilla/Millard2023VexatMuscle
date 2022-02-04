
flag_outerLoopMode = 1;


% Run once set to 1
% Run a second time set to 2 

if(flag_outerLoopMode == 0)
  clc;  
  close all;      
  clear all;

  flag_outerLoopMode      = 0;
  flag_buildCombinedPlot  = 0;
  
  figCombined = figure;
  
  flag_simulateHillModel            = 0; 
  flag_simulateOpus31Model          = 0;

  flag_useCrossbridgeStiffnessCalibratedCurves = 0;

  flag_plotData                     = 1;  
  flag_savePlotsToFile              = 1;

  flag_useElasticTendon             = 0; 
  flag_useFiberDamping              = 1;
  fiberDampingCoefficient           = 0.01;
  
  flag_useFig3KirchBoskovRymer1994  = 0; 
  flag_useOctave                    = 0;
  
  flag_activeForceLengthSimulations  = 0;
  flag_passiveForceLengthSimulations = 0;  
  flag_forceVelocitySimulations      = 0;
  
  normFiberLengthAtForceVelocitySample  = 1.;
  
  maxShorteningVelocity = 4.5;
  forceVelocityNormFiberHalfLength = 0.05;
  
  numberOfLengthSteps   = 20;
  numberOfVelocitySteps = 20;
  
  musculotendonPropertyAdjOpus31(4) = struct('name',[],'value',[]);
  
  musculotendonPropertyAdjOpus31(1).name = 'normTendonDampingConstant';
  musculotendonPropertyAdjOpus31(1).value = NaN;%0.0;
  
  musculotendonPropertyAdjOpus31(2).name = 'normTendonDampingLinear';
  musculotendonPropertyAdjOpus31(2).value = NaN;%0.1;
  
  musculotendonPropertyAdjOpus31(3).name = 'maximumNormalizedFiberVelocity';
  musculotendonPropertyAdjOpus31(3).value = NaN;%4.5;%4.5; 
  
  musculotendonPropertyAdjOpus31(4).name = 'normNumericalDamping';
  musculotendonPropertyAdjOpus31(4).value = NaN;%1e-4; 
  
            
  sarcomerePropertyAdjOpus31(5) = struct('name',[],'value',[]);
  
  
  sarcomerePropertyAdjOpus31(1).name = 'normCrossBridgeCyclingDamping';
  sarcomerePropertyAdjOpus31(1).value = NaN;%60;
  
  sarcomerePropertyAdjOpus31(2).name = 'normMaxActiveTitinToActinDamping';
  sarcomerePropertyAdjOpus31(2).value = NaN;
  
  sarcomerePropertyAdjOpus31(3).name = 'normPassiveTitinToActinDamping';  
  sarcomerePropertyAdjOpus31(3).value = NaN;
  
  sarcomerePropertyAdjOpus31(4).name = 'slidingTimeConstant';  
  sarcomerePropertyAdjOpus31(4).value = NaN;%0.01;  
 
  
  
  
  
  figClassicFl = figure;
  figClassicFv = figure;
  figClassicFv_Opus31AccelerationFactors = figure;
  figClassicFl_Opus31AccelerationFactors = figure;

  

end

dataFolder            = 'experiments/StandardTests/';
plotFolder            = 'output/plots/StandardTests/';

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
plotConfig;

%%
% Load the latest cat soleus configuration
%%
%Basic parameters for the Hill model
load('output/structs/defaultFelineSoleus.mat')
musculotendonProperties   = defaultFelineSoleus.musculotendon;
sarcomereProperties       = defaultFelineSoleus.sarcomere;
normMuscleCurves          = defaultFelineSoleus.curves;

normTendonDampingConstant = ...
    musculotendonProperties.normTendonDampingConstant;
normTendonDampingLinear = ...
    musculotendonProperties.normTendonDampingLinear;

figNameGainPhase = 'Fig12';
if(flag_useFig3KirchBoskovRymer1994==1)
  figNameGainPhase = 'Fig3';  
end

load(['output/structs/felineSoleusRigidTendonKBR1994',figNameGainPhase,'.mat']);
musculotendonPropertiesOpus31_RT = felineSoleusRigidTendonKBR1994.musculotendon;
sarcomerePropertiesOpus31_RT     = felineSoleusRigidTendonKBR1994.sarcomere;

load(['output/structs/felineSoleusElasticTendonKBR1994',figNameGainPhase,'.mat']);
musculotendonPropertiesOpus31_ET = felineSoleusElasticTendonKBR1994.musculotendon;
sarcomerePropertiesOpus31_ET     = felineSoleusElasticTendonKBR1994.sarcomere;

sarcomerePropertiesOpus31       = [];
musculotendonPropertiesOpus31   = [];

if(exist('musculotendonPropertyAdjOpus31','var')==1)
  for i=1:1:length(musculotendonPropertyAdjOpus31)
    if(isnan(musculotendonPropertyAdjOpus31(i).value)==0)
      musculotendonPropertiesOpus31_RT.( ...
        musculotendonPropertyAdjOpus31(i).name) = ...
        musculotendonPropertyAdjOpus31(i).value;
      musculotendonPropertiesOpus31_ET.( ...
        musculotendonPropertyAdjOpus31(i).name) = ...
        musculotendonPropertyAdjOpus31(i).value;    
    end
  end
end

if(exist('sarcomerePropertyAdjOpus31','var')==1)
  for i=1:1:length(sarcomerePropertyAdjOpus31)
    if(isnan(sarcomerePropertyAdjOpus31(i).value)==0)
      sarcomerePropertiesOpus31_RT.( ...
        sarcomerePropertyAdjOpus31(i).name) = ...
        sarcomerePropertyAdjOpus31(i).value;
      sarcomerePropertiesOpus31_ET.( ...
        sarcomerePropertyAdjOpus31(i).name) = ...
        sarcomerePropertyAdjOpus31(i).value;    
    end
  end
end

if(flag_useElasticTendon==1)
  sarcomerePropertiesOpus31       = sarcomerePropertiesOpus31_ET;
  musculotendonPropertiesOpus31   = musculotendonPropertiesOpus31_ET;
  
else
  sarcomerePropertiesOpus31       = sarcomerePropertiesOpus31_RT;
  musculotendonPropertiesOpus31   = musculotendonPropertiesOpus31_RT;
  
end

%%
% Generate model & parameter specific names
%%

generateModelSpecificKeyWords;    

tendonTag = '_ElasticTendon';
if(flag_useElasticTendon==0)
  tendonTag = '_RigidTendon';
end


outputFileEndingOpus31    = '';
outputFileEndingOpus31_ET = '';
outputFileEndingOpus31_RT = '';

outputFileEndingOpus31_ET = ...
  sprintf('_K%sD%sTau%s_KTC%s_KTL%s_StandardTests_%s_%s_%s_%s',...
  kScaleStr_ET, dScaleStr_ET, tScaleStr, kTConstStr, kTLinearStr,...
  ['TiAD',titinActiveDampingStr],...
  ['TiPD',titinPassiveDampingStr],...  
  strFittingBandwidth,...
  ['useCalibCurves_',num2str(flag_useCrossbridgeStiffnessCalibratedCurves)]);

outputFileEndingOpus31_RT = ...
 sprintf('_K%sD%sTau%s_StandardTests_%s_%s_%s_%s',...
  kScaleStr_RT, dScaleStr_RT, tScaleStr, ...
  ['TiAD',titinActiveDampingStr],...
  ['TiPD',titinPassiveDampingStr],...    
  strFittingBandwidth,...
  ['useCalibCurves_',num2str(flag_useCrossbridgeStiffnessCalibratedCurves)]);

if(flag_useElasticTendon==1)
  outputFileEndingOpus31 = outputFileEndingOpus31_ET;
else
  outputFileEndingOpus31 = outputFileEndingOpus31_RT;
end

outputFileEndingHill = sprintf('_D%i_StandardTests',flag_useFiberDamping);

outputFileEndingHill_RT = sprintf('_D%i_StandardTests',0);
outputFileEndingHill_ET = sprintf('_D%i_StandardTests',1);

%%
% Simulations
%%

if(flag_simulateHillModel == 1)

  [success] = runStandardForceLengthVelocitySimulationsDampedEquilibrium(...
                normFiberLengthAtForceVelocitySample,...    
                flag_useElasticTendon,...
                flag_useFiberDamping,...
                fiberDampingCoefficient,...
                musculotendonProperties,...
                sarcomereProperties,...
                normMuscleCurves,...
                outputFileEndingHill, ...
                dataFolder,...
                flag_activeForceLengthSimulations,...
                flag_passiveForceLengthSimulations,...
                flag_forceVelocitySimulations,...
                numberOfLengthSteps,...
                numberOfVelocitySteps,...   
                maxShorteningVelocity,...
                forceVelocityNormFiberHalfLength,...
                flag_useOctave);
end
if(flag_simulateOpus31Model == 1)

  normMuscleCurves.useCalibratedCurves = flag_useCrossbridgeStiffnessCalibratedCurves;

  [success] = runStandardForceLengthVelocitySimulationsOpus31(...
                normFiberLengthAtForceVelocitySample,...
                flag_useElasticTendon,...
                musculotendonPropertiesOpus31,...
                sarcomerePropertiesOpus31,...
                normMuscleCurves,...
                outputFileEndingOpus31, ...
                dataFolder,...
                flag_activeForceLengthSimulations,...
                flag_passiveForceLengthSimulations,...
                flag_forceVelocitySimulations,...
                numberOfLengthSteps,...
                numberOfVelocitySteps,...
                maxShorteningVelocity,...
                forceVelocityNormFiberHalfLength,...
                flag_useOctave);
end


%%
% Force-length experiments
%%
if(flag_plotData == 1)

  nameModification = '';
  if(flag_useElasticTendon == 1)
   nameModification = 'ElasticTendon';
  else
   nameModification = 'RigidTendon';          
  end
  
  dataActiveForceLengthHill = ...
    load([dataFolder,'activeForceLengthHill_',...
          nameModification,outputFileEndingHill,'.mat']);

  dataPassiveForceLengthHill = ...
    load([dataFolder,'passiveForceLengthHill_',...
          nameModification,outputFileEndingHill,'.mat']);

  dataForceVelocityHill = ...
    load([dataFolder,'forceVelocityHill_',...
          nameModification,outputFileEndingHill,'.mat']);

  dataActiveForceLength31 = ...
    load([dataFolder,'activeForceLengthOpus31_',...
          nameModification,outputFileEndingOpus31,'.mat']);

  dataPassiveForceLength31 = ...
    load([dataFolder,'passiveForceLengthOpus31_',...
          nameModification,outputFileEndingOpus31,'.mat']);

  dataForceVelocity31 = ...
    load([dataFolder,'forceVelocityOpus31_',...
          nameModification,outputFileEndingOpus31,'.mat']);        
     
        
   

  lineColorOpus31 = [0,0,1];
  lineColorHill   = [1,0,0];  
  nameOpus31 = 'Model: ET';
  nameHill   = 'Hill: ET';
  lineWidthHill = 0.5;
  lineWidthOpus31 = 0.5;
  
  
  markerSize = 2;
  if(flag_buildCombinedPlot==1)
    lineColorOpus31 = lineColorOpus31.*0.25+[1,1,1].*0.75;
    lineColorHill   = lineColorHill.*0.25+[1,1,1].*0.75;  
    markerSize = 4;    
    nameOpus31 = 'Model: RT';
    nameHill   = 'Hill: RT';
    lineWidthHill = 2;
    lineWidthOpus31 = 2;
    
  end


  
  subPlotPanel(1,1,:)=subPlotHerzogLeonard2000Stability(2,1,:);
  subPlotPanel(2,1,:)=subPlotHerzogLeonard2000Stability(3,1,:);
  subPlotPanel(3,1,:)=subPlotPanel(1,1,:);
  subPlotPanel(3,1,1)=subPlotPanel(3,1,1)+subPlotPanel(3,1,3)*1.2;
  %%
  % Force-Length Plots: Plot the reference curves
  %%
  figure(figClassicFl);

  activeForceLengthCurveSample = ...
    calcBezierYFcnXCurveSampleVector(...
    normMuscleCurves.activeForceLengthCurve, 100,[]);

  passiveForceLengthCurveSample = ...
    calcBezierYFcnXCurveSampleVector(...
    normMuscleCurves.fiberForceLengthCurve, 100,[]);

  figure(figClassicFl);

  subplot('Position',reshape(subPlotPanel(1,1,:),1,4));    
  
  if(flag_buildCombinedPlot==1)  
    plot(activeForceLengthCurveSample.x,...
         activeForceLengthCurveSample.y,...
         '-','Color',[0.75,0.75,0.75].*1,'LineWidth',3);
    hold on;  
    plot(passiveForceLengthCurveSample.x,...
         passiveForceLengthCurveSample.y,...
         '-','Color',[0.75,0.75,0.75].*1,'LineWidth',3);
    hold on;
  end

  for k=1:1:size(dataActiveForceLengthHill.benchRecord.time,2)  
    subplot('Position',reshape(subPlotPanel(1,1,:),1,4));
            
    fpeATHill = interpLocallyMonotonicCurve(...
      dataPassiveForceLengthHill.benchRecord.normFiberLength,...
      dataPassiveForceLengthHill.benchRecord.normFiberForce,...
      dataActiveForceLengthHill.benchRecord.normFiberLength(end,k),...
      5);
        
    plot(dataActiveForceLengthHill.benchRecord.normFiberLength(end,k), ...
         dataActiveForceLengthHill.benchRecord.normFiberForce(end,k)...
         -fpeATHill,...
         'x',...
         'Color',lineColorHill,'LineWidth',0.5,...
         'MarkerSize',markerSize,...
         'MarkerFaceColor',lineColorHill,...
         'DisplayName',nameHill);
    hold on;
    box off;
    
    if( abs(dataActiveForceLengthHill.benchRecord.normFiberLength(end,k)-1)<0.025)
      subplot('Position',reshape(subPlotPanel(2,1,:),1,4));
        plot(dataActiveForceLengthHill.benchRecord.time(:,k), ...
             dataActiveForceLengthHill.benchRecord.normFiberForce(:,k),...
              '-','Color',lineColorHill,...
              'LineWidth',lineWidthHill,...
              'DisplayName',nameHill);
        hold on;        
        plot(dataActiveForceLengthHill.benchRecord.time(end,k), ...
             dataActiveForceLengthHill.benchRecord.normFiberForce(end,k),...
             'x','Color',lineColorHill,'LineWidth',0.5,...
             'MarkerSize',markerSize,...
             'MarkerFaceColor',lineColorHill,...
             'DisplayName',nameHill);
        box off;
        hold on;
        
        subplot('Position',reshape(subPlotPanel(1,1,:),1,4));
        plot(dataPassiveForceLengthHill.benchRecord.normFiberLength(:,1), ...
             dataPassiveForceLengthHill.benchRecord.normFiberForce(:,1),...
             '-','Color',lineColorHill,...
             'LineWidth',lineWidthHill,...
             'DisplayName',nameHill);
        hold on;
        box off;
    end     
  end
    
  for k=1:1:size(dataActiveForceLength31.benchRecord.time,2)  
    
    fpeATOpus31 = interpLocallyMonotonicCurve(...
      dataPassiveForceLength31.benchRecord.normFiberLength,...
      dataPassiveForceLength31.benchRecord.normFiberForce,...
      dataActiveForceLength31.benchRecord.normFiberLength(end,k),...
      5);
            
    subplot('Position',reshape(subPlotPanel(1,1,:),1,4));
    plot( dataActiveForceLength31.benchRecord.normFiberLength(end,k), ...
          dataActiveForceLength31.benchRecord.normFiberForce(end,k)...
          -fpeATOpus31,...
          'o',...
          'Color',lineColorOpus31,'LineWidth',0.5,...
          'MarkerSize',markerSize,...
          'MarkerFaceColor',lineColorOpus31,...
          'DisplayName',nameOpus31);
    hold on;

    if( abs(dataActiveForceLength31.benchRecord.normFiberLength(end,k)-1)<0.025)
      subplot('Position',reshape(subPlotPanel(2,1,:),1,4));
      
        plot( dataActiveForceLength31.benchRecord.time(:,k), ...
              dataActiveForceLength31.benchRecord.normFiberForce(:,k),...
              '-','Color',lineColorOpus31,...
              'LineWidth',lineWidthOpus31,...
              'DisplayName',nameOpus31);
        hold on;    
        
        plot( dataActiveForceLength31.benchRecord.time(end,k), ...
              dataActiveForceLength31.benchRecord.normFiberForce(end,k),...
              'o','Color',lineColorOpus31,'LineWidth',0.5,...
              'MarkerSize',markerSize,...
              'MarkerFaceColor',lineColorOpus31,...
              'DisplayName',nameOpus31);
        hold on;  
        box off;
        
        subplot('Position',reshape(subPlotPanel(1,1,:),1,4));
        plot(dataPassiveForceLength31.benchRecord.normFiberLength(:,1), ...
             dataPassiveForceLength31.benchRecord.normFiberForce(:,1),...
             '-','Color',lineColorOpus31,...
             'LineWidth',lineWidthOpus31,...
             'DisplayName',nameOpus31);
        hold on;
        box off;        

    end      
  end

  subplot('Position',reshape(subPlotPanel(1,1,:),1,4));
  ylim([0,1.1]);
  if(flag_buildCombinedPlot==2)  
    xlabel('Norm. Length ($\ell^{M} / \ell^M_\circ$)');
    ylabel('Norm. Force ($f^{M} / f^M_\circ$)');
    
    %legend('Location','NorthWest');
    %legend boxoff;
    colorHillRT  = [1,0,0].*0.25+[1,1,1].*0.75;
    colorModelRT = [0,0,1].*0.25+[1,1,1].*0.75;
    colorHillET  = [1,0,0];%.*0.75+[1,1,1].*0.25;
    colorModelET = [0,0,1];%.*0.75+[1,1,1].*0.25;
    
    
    dy = 1*0.05;
    dx = 2*0.05;
    
    x0 = dx*2;    
    
    plot(  x0,1.,'o','Color',colorModelRT,'MarkerFaceColor',colorModelRT,...
          'MarkerSize',4,'LineWidth',0.5);
    hold on;
    text(  x0+dx,1,'Model: RT');
    hold on;

    plot( x0,1-dy,'x','Color',colorHillRT,'MarkerFaceColor',colorHillRT,...
          'MarkerSize',4,'LineWidth',0.5);
    hold on;    
    text( x0+dx,1-dy,'Hill: RT');
    hold on;
    
    plot( x0,1-2*dy,'o','Color',colorModelET,'MarkerFaceColor',colorModelET,...
          'MarkerSize',2,'LineWidth',0.5);
    hold on;        
    text( x0+dx,1-2*dy,'Model: ET');
    hold on;

    plot( x0,1-3*dy,'o','Color',colorHillET,'MarkerFaceColor',colorHillET,...
          'MarkerSize',2,'LineWidth',0.5);
    hold on;            
    text( x0+dx,1-3*dy,'Hill: ET');
    hold on;
    
  end
  subplot('Position',reshape(subPlotPanel(2,1,:),1,4));
  ylim([0,1.1]);
  if(flag_buildCombinedPlot==2)  
    xlabel('Time (s)');
    ylabel('Norm. Force ($f^{M} / f^M_\circ$)');
  end
  %%
  % Force Velocity Plots
  %%

  figure(figClassicFv);

  normPassiveForceAtFvSampleLengthHill = ...
    interpLocallyMonotonicCurve(...
      dataPassiveForceLengthHill.benchRecord.normFiberLength,...
      dataPassiveForceLengthHill.benchRecord.normFiberForce,...
      normFiberLengthAtForceVelocitySample,...
      5);
  
          
  normPassiveForceAtFvSampleLengthOpus31 = ...
    interpLocallyMonotonicCurve(...
      dataPassiveForceLength31.benchRecord.normFiberLength,...
      dataPassiveForceLength31.benchRecord.normFiberForce,...
      normFiberLengthAtForceVelocitySample,...
      5);   
  
  subPlotAppended = [];

  forceVelocityCurveSample = ...
    calcBezierYFcnXCurveSampleVector(...
    normMuscleCurves.fiberForceVelocityCurve, 100,[]);

  figure(figClassicFv);
  subplot('Position',reshape(subPlotPanel(1,1,:),1,4));    
  
  if(flag_buildCombinedPlot~=0)  
    plot(forceVelocityCurveSample.x,...
         forceVelocityCurveSample.y,...
         '-','Color',[0.75,0.75,0.75].*1,'LineWidth',3);
    hold on;  
  end
  


  halfTrials = size(dataForceVelocityHill.benchRecord.time,2)/2;
  idxA = round(halfTrials*0.5);
  idxB = round(halfTrials*0.5)+halfTrials;

  
  for k=1:1:size(dataForceVelocityHill.benchRecord.time,2)  
    


      subplot('Position',reshape(subPlotPanel(1,1,:),1,4));

      plot(dataForceVelocityHill.benchRecord.eventNormFiberVelocity(1,k), ...
           dataForceVelocityHill.benchRecord.eventTendonForce(1,k)./musculotendonProperties.fiso,...
           'x',...
           'Color',lineColorHill,'LineWidth',0.5,...
           'MarkerSize',markerSize,...
           'MarkerFaceColor',lineColorHill,...
           'DisplayName',nameHill);
      hold on;
      box off;

  %     subPlotAppended = reshape(subPlotPanel(2,1,:),1,4);
  %     subPlotAppended(1,2) = subPlotAppended(1,2) ... 
  %                             - 0.5*subPlotAppended(1,4)...
  %                             - plotVertMargin;
  %     subPlotAppended(1,4) = 0.5*subPlotAppended(1,2);
  %     subplot('Position',subPlotAppended);
  %       plot(dataForceVelocity31.benchRecord.time(:,k),...
  %            dataForceVelocity31.benchRecord.pathLength(:,k).*1000,...
  %            '-','Color',[1,1,1].*0.75, 'LineWidth',3);
  %       hold on;

      %if( k== 2)% round(0.5*(size(dataForceVelocity31.benchRecord.time,2)-1)) )
  %         subplot('Position',subPlotAppended);
  %         plot(dataForceVelocity31.benchRecord.time(:,k),...
  %              dataForceVelocity31.benchRecord.pathLength(:,k),...
  %              '-','Color',[1,1,1].*0., 'LineWidth',3);
  %         hold on;
    if(k==idxA || k==idxB)
          subplot('Position',reshape(subPlotPanel(2,1,:),1,4));

          plot( dataForceVelocityHill.benchRecord.time(:,k), ...
                dataForceVelocityHill.benchRecord.tendonForce(:,k)./musculotendonProperties.fiso, ...
               '-','Color',lineColorHill,...
               'LineWidth',lineWidthHill,...
               'DisplayName',nameHill);
          hold on;        
          plot( dataForceVelocityHill.benchRecord.eventTime(1,k), ...
                dataForceVelocityHill.benchRecord.eventTendonForce(1,k)./musculotendonProperties.fiso, ...
                'x','Color',lineColorHill,'LineWidth',0.5,...
                'MarkerSize',markerSize,...
                'MarkerFaceColor',lineColorHill,...
                'LineWidth',lineWidthHill,...
                'DisplayName',nameHill);
          box off;
          hold on; 

         
      %end     
    end
  end
    
  assert(  size(dataForceVelocity31.benchRecord.time,2) ...
        == size(dataForceVelocityHill.benchRecord.time,2) );

  for k=1:1:size(dataForceVelocity31.benchRecord.time,2)  
    
    subplot('Position',reshape(subPlotPanel(1,1,:),1,4));
    plot(dataForceVelocity31.benchRecord.eventNormFiberVelocity(1,k), ...
         dataForceVelocity31.benchRecord.eventTendonForce(1,k)./musculotendonPropertiesOpus31.fiso,...
         'o',...
         'Color',lineColorOpus31,'LineWidth',0.5,...
         'MarkerSize',markerSize,...
         'MarkerFaceColor',lineColorOpus31,...
         'DisplayName',nameOpus31);
    hold on;

    %if( k == 2) %round(0.5*(size(dataForceVelocity31.benchRecord.time,2)-1)) )
    if(k==idxA || k==idxB)
      subplot('Position',reshape(subPlotPanel(2,1,:),1,4));
        plot( dataForceVelocity31.benchRecord.time(:,k), ...
              dataForceVelocity31.benchRecord.tendonForce(:,k)./musculotendonPropertiesOpus31.fiso,...
             '-','Color',lineColorOpus31,...
             'LineWidth',lineWidthOpus31,...
             'DisplayName',nameOpus31);
        hold on;        
        plot( dataForceVelocity31.benchRecord.eventTime(1,k), ...
              dataForceVelocity31.benchRecord.eventTendonForce(1,k)./musculotendonPropertiesOpus31.fiso,...
             'o','Color',lineColorOpus31,'LineWidth',0.5,...
             'MarkerSize',markerSize,...
             'MarkerFaceColor',lineColorOpus31,...
             'LineWidth',lineWidthOpus31,...
             'DisplayName',nameOpus31);
        hold on;  
        box off;
        

        
        subplot('Position',reshape(subPlotPanel(3,1,:),1,4));
          plot( dataForceVelocity31.benchRecord.time(:,k), ...
              dataForceVelocity31.benchRecord.normFiberLength(:,k),...
             '-','Color',lineColorOpus31,...
             'LineWidth',lineWidthOpus31,...
             'DisplayName',nameOpus31);
          hold on;
          plot( dataForceVelocity31.benchRecord.eventTime(:,k), ...              
              dataForceVelocity31.benchRecord.eventNormFiberLength(1,k),...
             'o','Color',lineColorOpus31,'LineWidth',0.5,...
             'MarkerSize',markerSize,...
             'MarkerFaceColor',lineColorOpus31,...
             'LineWidth',lineWidthOpus31,...
             'DisplayName',nameOpus31);
          hold on;
        xlabel('Time (s)');
        ylabel('Norm. Fiber Length');

        box off;
        subplot('Position',reshape(subPlotPanel(2,1,:),1,4));
    end
   
    %end      
  end

  subplot('Position',reshape(subPlotPanel(1,1,:),1,4));
  ylim([-0.05,1.6]);
  if(flag_buildCombinedPlot==2)    
    xlabel('Norm. Velocity ($v^{M} / v^M_\circ$)');
    ylabel('Norm. Force ($f^{M} / f^M_\circ$)');
    box off;
    
    dy = 1.5*0.05;
    dx = 2.4*0.05;
    
    x0 = -1.1;
    
    plot(x0,1.45,'o','Color',colorModelRT,'MarkerSize',4,...
         'MarkerFaceColor',colorModelRT,'LineWidth',0.5);
    hold on;
    text(x0+dx,1.45,'Model: RT');
    hold on;

    plot(x0,1.45-dy,'x','Color',colorHillRT,'MarkerSize',4,...
         'MarkerFaceColor',colorHillRT,'LineWidth',0.5);
    hold on;    
    text(x0+dx,1.45-dy,'Hill: RT');
    hold on;
    
    plot(x0,1.45-2*dy,'o','Color',colorModelET,'MarkerSize',2,...
         'MarkerFaceColor',colorModelET,'LineWidth',0.5);
    hold on;        
    text(x0+dx,1.45-2*dy,'Model: ET');
    hold on;

    plot(x0,1.45-3*dy,'o','Color',colorHillET,'MarkerSize',2,...
          'MarkerFaceColor',colorHillET,'LineWidth',0.5);
    hold on;            
    text(x0+dx,1.45-3*dy,'Hill: ET');
    hold on;

  end

  subplot('Position',reshape(subPlotPanel(2,1,:),1,4));
  ylim([0.0,1.75]);
  if(flag_buildCombinedPlot==2)  
    xlabel('Time (s)');
    ylabel('Norm. Force ($f^{M} / f^M_\circ$)');
    box off;
    
  end
  %subplot('Position',subPlotAppended);
  %xlabel('Time (s)');
  %ylabel('Length (mm)');
  %box off;

  %%
  % Write the plots to file
  %% 

  figure(figClassicFl);
  pause(0.1);

  set(figClassicFl,'Units','centimeters',...
  'PaperUnits','centimeters',...
  'PaperSize',[pageWidth pageHeight],...
  'PaperPositionMode','manual',...
  'PaperPosition',[0 0 pageWidth pageHeight]);     
  set(figClassicFl,'renderer','painters');     
  set(gcf,'InvertHardCopy','off')


  if(flag_buildCombinedPlot > 0)
    tendonTag = 'ElasticRigid';
  end
  
  print('-dpdf', [plotFolder,'fig_Pub_ForceLength_',...
                  tendonTag,'_',tScaleStr,'_TiAD',titinActiveDampingStr,...
                  '_TiPD',titinPassiveDampingStr,...  
                  '_',strFittingBandwidth,'.pdf']);   
  pause(0.1);


  figure(figClassicFv);
  pause(0.1);

  set(figClassicFv,'Units','centimeters',...
  'PaperUnits','centimeters',...
  'PaperSize',[pageWidth pageHeight],...
  'PaperPositionMode','manual',...
  'PaperPosition',[0 0 pageWidth pageHeight]);     
  set(figClassicFv,'renderer','painters');     
  set(gcf,'InvertHardCopy','off')


  print('-dpdf', [plotFolder,'fig_Pub_ForceVelocity_',...
                  tendonTag,'_',tScaleStr,'_TiAD',titinActiveDampingStr,...
                  '_TiPD',titinPassiveDampingStr,...  
                  '_',strFittingBandwidth,'.pdf']);   

  pause(0.1);
  
end
