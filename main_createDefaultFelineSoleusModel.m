flag_outerLoopMode = 0;

if(flag_outerLoopMode == 0)
  clc;
  close all;  
  clear all;

  fitCrossBridgeStiffnessDampingToKirch199490Hz=1;
  flag_makeAndSavePubPlots                         = 1;
  flag_fitToFig3KirchBoskovRymer1994               = 0;
  fitCrossBridgeStiffnessDampingToKirch199490Hz    = 1;
  flag_fitActiveTitinProperties                    = 0;

  rigidTendonReferenceModel             = [];
  elasticTendonReferenceModel           = [];
  normPevkToActinAttachmentPoint        = 0.5;
  normFiberLengthAtOneNormPassiveForce  = 1.367732948060934e+00;

end

pubOutputFolder                 = 'output/plots/MuscleCurves/';
postprocessingDirectoryTree     = genpath('postprocessing');
addpath(postprocessingDirectoryTree   );

flag_useOctave = 0;
flag_enableNumericallyNonZeroGradients    = 1;
flag_plotEveryDefaultCurve = 0;

% 1. Active force length curve vs. data
% Solution: There were some initial descrepencies between the experimental force
%length data and a theoretical curve. These errors almost completely go
%away if it is assumed that the experimental recordings are of total 
%path length, rather than fiber length. In this case, when the elastiticy
%of the tendon is taken into account the theoretical active-force-length 
%curve and the transformed data nicely align.

%Failed attempt:
%This creates a cat soleus with an optimal fiber length of 58 mm: this
%is simply way too big to be realistic (given the data I'm working from)
flag_solveForOptimalFiberLengthOfBestFit  = 0; 

%Failed attempt:
shiftLengthActiveForceLengthCurveDescendingCurve = 0.;%...
%  (1/3)*( (1.154-1.087) + (1.23-1.162) + (1.077-1.039) );

smallNumericallyNonZeroNumber           = sqrt(sqrt(eps));

%%
% Add the directories needed to run this script
%%
parametersDirectoryTreeMTParams     = genpath('parameters');
parametersDirectoryTreeExperiments  = genpath('experiments');
parametersDirectoryTreeModels       = genpath('models');
parametersDirectoryTreeCurves       = genpath('curves');
parametersDirectoryTreeSimulation   = genpath('simulation');

addpath(parametersDirectoryTreeMTParams);
addpath(parametersDirectoryTreeExperiments);
addpath(parametersDirectoryTreeModels);
addpath(parametersDirectoryTreeCurves);
addpath(parametersDirectoryTreeSimulation);

scaleOptimalFiberLength      = 1.0; 
scaleMaximumIsometricTension = 1;

if(exist('fitCrossBridgeStiffnessDampingToKirch199490Hz','var')==0)
  fitCrossBridgeStiffnessDampingToKirch199490Hz = 1;
end

[felineSoleusMusculotendonProperties, ...
 felineSoleusSarcomereProperties,...
 felineSoleusActiveForceLengthData,...
 felineSoleusPassiveForceLengthData] = ...
 createFelineSoleus(    scaleOptimalFiberLength,...
                        scaleMaximumIsometricTension,...
                        normFiberLengthAtOneNormPassiveForce,...
                        normPevkToActinAttachmentPoint,...
                        fitCrossBridgeStiffnessDampingToKirch199490Hz,...
                        flag_useOctave);

createMusculoTendonFcn = ...
  @(argScaleFiberLength,argScaleFiso)createFelineSoleus(...
                                        argScaleFiberLength,...
                                        argScaleFiso,...
                                        normFiberLengthAtOneNormPassiveForce,...
                                        normPevkToActinAttachmentPoint,...
                                        fitCrossBridgeStiffnessDampingToKirch199490Hz,...
                                        flag_useOctave); 
                                        
[felineSoleusNormMuscleCurvesDefault,...
 felineSoleusMusculotendonPropertiesDefault,...
 felineSoleusSarcomerePropertiesDefault,... 
 activeForceLengthCurveAnnotationPoints,...
 felineSoleusActiveForceLengthDataDefault,...
 felineSoleusPassiveForceLengthDataDefault,...
 felineSoleusPassiveForceLengthCurveSettings]= ...
    createFittedMuscleCurves( ...
      felineSoleusMusculotendonProperties,...
      felineSoleusSarcomereProperties,...
      felineSoleusActiveForceLengthData,...
      felineSoleusPassiveForceLengthData,...
      shiftLengthActiveForceLengthCurveDescendingCurve,...
      flag_enableNumericallyNonZeroGradients,...
      smallNumericallyNonZeroNumber,...
      flag_solveForOptimalFiberLengthOfBestFit,...
      createMusculoTendonFcn,...
      flag_useOctave);


%%
%Check to make sure that
% The normFiberLengthAtOneNormPassiveForce used to create keypoints for the
% titin model actually matches the norm fiber length at which the passive
% curve develops 1 normalized force.
%%
fpe2 = felineSoleusNormMuscleCurvesDefault.fiberForceLengthCurve.yEnd(1,2);
lpe2 = felineSoleusNormMuscleCurvesDefault.fiberForceLengthCurve.xEnd(1,2);
assert(abs(fpe2-1)<1e-3);
assert(abs(lpe2-felineSoleusSarcomerePropertiesDefault.normFiberLengthAtOneNormPassiveForce)<1e-3);



defaultFelineSoleus = struct('musculotendon',...
                            felineSoleusMusculotendonPropertiesDefault,...
                            'sarcomere',...
                            felineSoleusSarcomerePropertiesDefault,...
                            'falData',...
                            felineSoleusActiveForceLengthDataDefault,...
                            'fpeData',...
                            felineSoleusPassiveForceLengthDataDefault,...
                            'curves',...
                            felineSoleusNormMuscleCurvesDefault);
                      
save('output/structs/defaultFelineSoleus.mat',...
     'defaultFelineSoleus');  

%%
%% Update the active titin properties of the model to fit 
%%   a single trial of Herzog & Leonard 2002
%%

%Fit the elastic tendon model.
if(flag_fitActiveTitinProperties== 0 )
    if(isempty(elasticTendonReferenceModel)==0)
        disp('Using reference model');        
        tmp=load(elasticTendonReferenceModel);
        modelName = fields(tmp);
        felineSoleusSarcomerePropertiesUpd_ET= tmp.(modelName{1}).sarcomere;
        felineSoleusNormMuscleCurvesUpd_ET= tmp.(modelName{1}).curves;
    else
        %assert(0,['Error: flag_fitActiveTitinProperties is disabled ',...
        %    'but no elasticTendonReferenceModel is provided']);
        disp('Using default model');
        felineSoleusSarcomerePropertiesUpd_ET = felineSoleusSarcomerePropertiesDefault;
        felineSoleusNormMuscleCurvesUpd_ET    = felineSoleusNormMuscleCurvesDefault;
    end

else
    %The parameters updated are the 
    %  :normPevkToActinAttachmentPoint (default of 0.5)
    %  :normMaxActiveTitinToActinDamping (default of 20)
    %
    %The hand-tuned default values are quite good, but fitting is required to
    %minimize the error. Since the process of both force development and 
    %relaxation are nonlinear, there is not an elegant and fast way to find 
    %these parameters without simulating the model directly. 
    
    figureNumber       = 7;
    subFigureNumber    = 2;
    trialNumber        = 3;  
    
    expConfigHerzogLeonard2002 =...
     getHerzogLeonard2002Configuration( figureNumber,...
                                        subFigureNumber, ...
                                        trialNumber);
    
    dataFolder = 'experiments/HerzogLeonard2002/fitting/';
    
    if(felineSoleusSarcomerePropertiesDefault.titinModelType==1)
        felineSoleusNormMuscleCurvesDefault.useTwoSidedTitinCurves=0;
    end
    felineSoleusNormMuscleCurvesDefault.useCalibratedCurves=1;
    flag_useElasticTendon       = 1;
    
    
    [felineSoleusSarcomerePropertiesUpd_ET,...
        felineSoleusNormMuscleCurvesUpd_ET] = ...
        updateActiveTitinParameters(felineSoleusMusculotendonPropertiesDefault, ...
                                 felineSoleusSarcomerePropertiesDefault,...
                                 felineSoleusNormMuscleCurvesDefault,...
                                 felineSoleusPassiveForceLengthCurveSettings,...
                                 expConfigHerzogLeonard2002,...
                                 flag_useElasticTendon,...
                                 dataFolder,...
                                 flag_fitActiveTitinProperties,...
                                 flag_useOctave);
    
        disp('felineSoleusSarcomerePropertiesUpd_ET');
        fprintf('\t%e\t%s\n',felineSoleusSarcomerePropertiesUpd_ET.normPevkToActinAttachmentPoint,... 
                            'normPevkToActinAttachmentPoint');
        fprintf('\t%e\t%s\n',felineSoleusSarcomerePropertiesUpd_ET.normMaxActiveTitinToActinDamping,... 
                            'normMaxActiveTitinToActinDamping');
end

%Fit the rigid tendon model
if(flag_fitActiveTitinProperties== 0 )
    if(isempty(rigidTendonReferenceModel)==0)
        disp('Using default model');        
        tmp=load(rigidTendonReferenceModel);
        modelName = fields(tmp);
        felineSoleusSarcomerePropertiesUpd_RT= tmp.(modelName{1}).sarcomere;
        felineSoleusNormMuscleCurvesUpd_RT= tmp.(modelName{1}).curves;
    else
        disp('Using reference model');
        felineSoleusSarcomerePropertiesUpd_RT = felineSoleusSarcomerePropertiesDefault;
        felineSoleusNormMuscleCurvesUpd_RT    = felineSoleusNormMuscleCurvesDefault;

    end

else

    if(felineSoleusSarcomerePropertiesDefault.titinModelType==0)
        felineSoleusNormMuscleCurvesDefault.useTwoSidedTitinCurves=1;
    end    
    felineSoleusNormMuscleCurvesDefault.useCalibratedCurves=1;
    flag_useElasticTendon       = 0;
    
    [felineSoleusSarcomerePropertiesUpd_RT,...
        felineSoleusNormMuscleCurvesUpd_RT] = ...
        updateActiveTitinParameters(felineSoleusMusculotendonPropertiesDefault, ...
                                 felineSoleusSarcomerePropertiesDefault,...
                                 felineSoleusNormMuscleCurvesDefault,...
                                 felineSoleusPassiveForceLengthCurveSettings,...
                                 expConfigHerzogLeonard2002,...
                                 flag_useElasticTendon,...
                                 dataFolder,...
                                 flag_fitActiveTitinProperties,...
                                 flag_useOctave);
    
    disp('felineSoleusSarcomerePropertiesUpd_RT');
    fprintf('\t%e\t%s\n',felineSoleusSarcomerePropertiesUpd_RT.normPevkToActinAttachmentPoint,... 
                        'normPevkToActinAttachmentPoint');
    fprintf('\t%e\t%s\n',felineSoleusSarcomerePropertiesUpd_RT.normMaxActiveTitinToActinDamping,... 
                        'normMaxActiveTitinToActinDamping');
end
%%
%% Update the cross-bridge properties of the Opus 31 model to fit the 
%% Frequency response of Kirch, Boskov, & Rymer 1994.
%%

%Since the force response of the model to small perturbations is well 
%approximated as linear we can directly evaluate the visco elastic
%properties of the lumped-cross bridge that best fit the data of Kirsch,
%Boskov, and Rymer. 

%%
% Fitting Data: Kirsch, Boskov, & Rymer 1994
%%
fittingFilesGain      = 'experiments/KirschBoskovRymer1994/data/fig_KirschBoskovRymer1994_Fig3_gain.csv';
fittingFilesPhase     = 'experiments/KirschBoskovRymer1994/data/fig_KirschBoskovRymer1994_Fig3_phase.csv';
fittingFilesCoherence = 'experiments/KirschBoskovRymer1994/data/fig_KirschBoskovRymer1994_Fig3_coherence.csv';
fittingFilesK = {'experiments/KirschBoskovRymer1994/data/fig_KirschBoskovRymer1994_Fig9A.csv',...
                 'experiments/KirschBoskovRymer1994/data/fig_KirschBoskovRymer1994_Fig9B.csv',...
                 'experiments/KirschBoskovRymer1994/data/fig_KirschBoskovRymer1994_Fig12_K.csv'}; 
fittingFilesD = {'experiments/KirschBoskovRymer1994/data/fig_KirschBoskovRymer1994_Fig10.csv',...
                 'experiments/KirschBoskovRymer1994/data/fig_KirschBoskovRymer1994_Fig12_D.csv'}; 



dataKBR1994Fig3Force = 5; %N, as mentioned in the caption

dataKBR1994Fig3Gain = loadDigitizedData(fittingFilesGain,...
                        'Frequency (Hz)','Stiffness (N/mm)',...
                        {'1.6mm, 90Hz','1.6mm, 15Hz'},'');
dataKBR1994Fig3Phase = loadDigitizedData(fittingFilesPhase,...
                         'Frequency (Hz)','Phase (deg)',...
                         {'1.6mm, 90Hz','1.6mm, 15Hz'},'');
dataKBR1994Fig3Coherence = loadDigitizedData(fittingFilesCoherence,...
                          'Frequency (Hz)','Coherence$$^2$$',...
                          {'1.6mm, 90Hz','1.6mm, 15Hz'},'');

dataKBR1994Fig9A = loadDigitizedData(fittingFilesK{1},...
                      'Force (N)','K (N/mm)',...
                      {'0.4mm 15Hz','0.8mm 15Hz','1.6mm 15Hz'},'');
                    
dataKBR1994Fig9B = loadDigitizedData(fittingFilesK{2},...
                      'Force (N)','K (N/mm)',...
                      {'0.4mm 15Hz','0.4mm 35Hz','0.4mm 90Hz',...
                       '1.6mm 15Hz','1.6mm 35Hz','1.6mm 90Hz'},'');

dataKBR1994Fig10 = loadDigitizedData(fittingFilesD{1},...
                      'Force (N)','K (N/mm/s)',...
                      {'15Hz','35Hz','90Hz'},'');
                     
dataKBR1994Fig12K = loadDigitizedData(fittingFilesK{3},...
                          'Force (N)','K (N/mm)',...
                          {'Soleus','MG'},'0.8mm 35Hz');
                        
dataKBR1994Fig12D = loadDigitizedData(fittingFilesD{2},...
                          'Force (N)','B (N/mm/s)',...
                          {'Soleus','MG'},'0.8mm 35Hz');

%%

% Load the muscle and sarcomere properties of the cat soleus
% 
% Kirsch, Boskov & Rymer mention that the soleus rest length is
% 13 mm to 6 mm short of the 'maximum physiological length'. I have
% no idea how the 'maximum physiological length' is defined. I'm going
% to start the simulations at the optimal fiber length for now.
%%


nominalNormFiberLengthAtSlack  = 1.0;
nominalForceKDFit              = 5; %Stiffness & damping fit done at 5N tension
scaleSlidingTimeConstant       = 1;
scaleCrossBridgeCyclingDamping = 1;

%tmp = load('output/structs/defaultFelineSoleus.mat');
%musculotendonProperties   = tmp.defaultFelineSoleus.musculotendon;
%sarcomereProperties       = tmp.defaultFelineSoleus.sarcomere;
%felineSoleusNormMuscleCurvesDefault = tmp.defaultFelineSoleus.curves;

musculotendonProperties = felineSoleusMusculotendonPropertiesDefault;

sarcomereProperties_ET     = felineSoleusSarcomerePropertiesUpd_ET;
normMuscleCurves_ET        = felineSoleusNormMuscleCurvesUpd_ET;
sarcomereProperties_RT     = felineSoleusSarcomerePropertiesUpd_RT;
normMuscleCurves_RT        = felineSoleusNormMuscleCurvesUpd_RT;


normTendonDampingConstant = ...
    musculotendonProperties.normTendonDampingConstant;
normTendonDampingLinear = ...
    musculotendonProperties.normTendonDampingLinear;

flag_updateNormFiberLengthByTendonStretch = 1;

                                  
[musculotendonPropertiesOpus31_RT,sarcomerePropertiesOpus31_RT] = ...
  updateOpus31CrossBridgeParameters(nominalForceKDFit,...
                                    nominalNormFiberLengthAtSlack,...
                                    flag_fitToFig3KirchBoskovRymer1994,...
                                    dataKBR1994Fig3Gain,...
                                    dataKBR1994Fig3Phase,...
                                    dataKBR1994Fig12K,...
                                    dataKBR1994Fig12D,...
                                    normTendonDampingConstant,...
                                    normTendonDampingLinear,...
                                    scaleSlidingTimeConstant,...
                                    scaleCrossBridgeCyclingDamping,...
                                    0,...
                                    musculotendonProperties,...
                                    sarcomereProperties_RT,...
                                    normMuscleCurves_RT,...
                                    flag_useOctave);                                  
                                  
[musculotendonPropertiesOpus31_ET,sarcomerePropertiesOpus31_ET] = ...
  updateOpus31CrossBridgeParameters(nominalForceKDFit,...
                                    nominalNormFiberLengthAtSlack,...
                                    flag_fitToFig3KirchBoskovRymer1994,...
                                    dataKBR1994Fig3Gain,...
                                    dataKBR1994Fig3Phase,...
                                    dataKBR1994Fig12K,...
                                    dataKBR1994Fig12D,...
                                    normTendonDampingConstant,...
                                    normTendonDampingLinear,...
                                    scaleSlidingTimeConstant,...
                                    scaleCrossBridgeCyclingDamping,...
                                    1,...
                                    musculotendonProperties,...
                                    sarcomereProperties_ET,...
                                    normMuscleCurves_ET,...
                                    flag_useOctave);   


felineSoleusRigidTendonKBR1994 = defaultFelineSoleus;
felineSoleusRigidTendonKBR1994.curves        = normMuscleCurves_RT;
felineSoleusRigidTendonKBR1994.musculotendon = musculotendonPropertiesOpus31_RT;
felineSoleusRigidTendonKBR1994.sarcomere     = sarcomerePropertiesOpus31_RT;

figNameGainPhase = 'Fig12';
if(flag_fitToFig3KirchBoskovRymer1994==1)
  figNameGainPhase = 'Fig3';  
end

save(['output/structs/felineSoleusRigidTendonKBR1994',figNameGainPhase,'.mat'],...
      'felineSoleusRigidTendonKBR1994');

felineSoleusElasticTendonKBR1994              = defaultFelineSoleus;
felineSoleusElasticTendonKBR1994.curves       = normMuscleCurves_ET;
felineSoleusElasticTendonKBR1994.musculotendon= musculotendonPropertiesOpus31_ET;
felineSoleusElasticTendonKBR1994.sarcomere    = sarcomerePropertiesOpus31_ET;

save(['output/structs/felineSoleusElasticTendonKBR1994',figNameGainPhase,'.mat'],...
      'felineSoleusElasticTendonKBR1994');


%%
% Note the average offset between the active-force-length curve and
% the transformed data
%%

xExp = felineSoleusActiveForceLengthDataDefault(2:end,1);
yExp = felineSoleusActiveForceLengthDataDefault(2:end,2);
xCurve = zeros(size(xExp));

for i=1:1:length(xExp)
xCurve(i,1) = calcBezierFcnXGivenY(yExp(i,1), ...
  felineSoleusNormMuscleCurvesUpd_ET.activeForceLengthCurve,... 
  xExp(i,1));
end                                    
dx = mean(xCurve-xExp);
%felineSoleusActiveForceLengthDataDefault(:,1)=...
%  felineSoleusActiveForceLengthDataDefault(:,1)+dx;

disp('Ajusted optimal fiber length');
fprintf('%1.6f lce/lopt \n',felineSoleusActiveForceLengthDataDefault(1,1));
disp('Average error on the descending limb');
fprintf('%1.6f lce/lopt \n',dx);
fprintf('%1.6f mm \n',dx*(musculotendonPropertiesOpus31_ET.optimalFiberLength*1000));

lceNStart = felineSoleusActiveForceLengthDataDefault(1,1);
save('output/structs/normalizedFiberLengthStartHerzogLeonard2002.mat',...
     'lceNStart');


    figH = plotStructOfBezierSplines( felineSoleusNormMuscleCurvesUpd_ET,...
                                      {'Inverse','use'});                          

%%
% Plot the derivative of the tendon force length curve on top of the
% stiffness curve
%% 
figure(figH.tendonStiffnessCurve);

    curveSample = calcBezierYFcnXCurveSampleVector(...
                    felineSoleusNormMuscleCurvesUpd_ET.('tendonForceLengthCurve'), 200,[]);

    xmin = min(curveSample.x);
    xmax = max(curveSample.x);
    ymin = min(curveSample.y);
    ymax = max(curveSample.y);

subplot(2,2,1);
  plot(curveSample.x, curveSample.dydx,...
    '--','Color',[1,1,1].*0.5,'LineWidth',2);
  hold on;
%%
%Plot experimental data over top of the curves where it is available.
%%
figure(figH.activeForceLengthCurve);
  subplot(2,2,1);  
  plot(  felineSoleusActiveForceLengthDataDefault(:,1),...
       felineSoleusActiveForceLengthDataDefault(:,2),'xb');
  hold on;          

figure(figH.fiberForceLengthCurve);
  subplot(2,2,1);
  plot(   felineSoleusPassiveForceLengthDataDefault(:,1),...
          felineSoleusPassiveForceLengthDataDefault(:,2),'xb');
  hold on;          
  
%%
% Plot the passive force length curve of the 2 segment titin model and
% compare it to the passive force length curve
%%

lceN0 = calcBezierFcnXGivenY(0, ...
          felineSoleusNormMuscleCurvesUpd_ET.fiberForceLengthCurve);
lceN1 = calcBezierFcnXGivenY(1, ...
          felineSoleusNormMuscleCurvesUpd_ET.fiberForceLengthCurve);
        
npts = 100;
lceNSeries = [lceN0:((lceN1-lceN0)/(npts-1)):lceN1]';

fNSeries = zeros(npts,4); % fpe, fecm, fIgp, fPevkIgd,
lNSeries = zeros(npts,4); % lpe, lecm, lIgp, lPevkIgd,

normLengthIgdFixed = ...
  sarcomerePropertiesOpus31_ET.IGDFixedNormLengthAtOptimalFiberLength;

normLengthT12ToZ = ...
  sarcomerePropertiesOpus31_ET.ZLineToT12NormLengthAtOptimalFiberLength;


for i=1:1:npts
  lceN  = lceNSeries(i,1);
  fpeN  = calcBezierYFcnXDerivative(lceN,...
            felineSoleusNormMuscleCurvesUpd_ET.fiberForceLengthCurve,0);
  fecmN = calcBezierYFcnXDerivative(lceN*0.5,...
            felineSoleusNormMuscleCurvesUpd_ET.forceLengthECMHalfCurve,0);
          
  lIgpPevkN = lceN*0.5 - normLengthIgdFixed - normLengthT12ToZ; 
          
  [lPN, lDN, fTiN] = calcSeriesSpringStretch(lIgpPevkN,...
            felineSoleusNormMuscleCurvesUpd_ET.forceLengthProximalTitinCurve,...
            felineSoleusNormMuscleCurvesUpd_ET.forceLengthProximalTitinInverseCurve, ...
            felineSoleusNormMuscleCurvesUpd_ET.forceLengthDistalTitinCurve,...
            felineSoleusNormMuscleCurvesUpd_ET.forceLengthDistalTitinInverseCurve);          

  fNSeries(i,:) = [fpeN,fecmN,fTiN, fTiN];
  lNSeries(i,:) = [lceN, lceN,lPN*2, lDN*2];
          
end


fig_forceLength = figure;
  plot(lNSeries(:,1),fNSeries(:,1),'k');
  hold on;
  plot(lNSeries(:,1), fNSeries(:,2)+fNSeries(:,3),'b');
  hold on;
  plot(lNSeries(:,1), fNSeries(:,2),'r');
  hold on;
  plot(felineSoleusPassiveForceLengthDataDefault(:,1),...
       felineSoleusPassiveForceLengthDataDefault(:,2),'xb');
  legend('fpe','fecm+fti','fecm','data');
  xlabel('Norm. Length')
  ylabel('Norm. Force');
  


                   

%%
% Generate publication quality plots
%%
if(flag_makeAndSavePubPlots==1)
  plotMuscleCurves( felineSoleusNormMuscleCurvesUpd_ET,...
                      activeForceLengthCurveAnnotationPoints,...
                      musculotendonPropertiesOpus31_ET,...
                      sarcomerePropertiesOpus31_ET,...
                      felineSoleusActiveForceLengthDataDefault,...
                      felineSoleusPassiveForceLengthDataDefault,...
                      pubOutputFolder);
end
   
%%
% Remove the directories ...
%%
rmpath(parametersDirectoryTreeMTParams);
rmpath(parametersDirectoryTreeExperiments);
rmpath(parametersDirectoryTreeModels);
rmpath(parametersDirectoryTreeCurves);
