flag_outerLoopMode = 0;

if(flag_outerLoopMode == 0)
  clc;
  close all;  
  clear all;


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



[humanSoleusMusculotendonProperties, ...
 humanSoleusSarcomereProperties] = ...
    createHumanSoleus(  scaleOptimalFiberLength,...
                        scaleMaximumIsometricTension,...
                        normFiberLengthAtOneNormPassiveForce,...
                        normPevkToActinAttachmentPoint,...
                        flag_useOctave);

createMusculoTendonFcn = ...
  @(argScaleFiberLength,argScaleFiso)createHumanSoleus(...
                                        argScaleFiberLength,...
                                        argScaleFiso,...
                                        normFiberLengthAtOneNormPassiveForce,...
                                        normPevkToActinAttachmentPoint,...
                                        flag_useOctave); 
                                        
humanSoleusActiveForceLengthData  = [];
humanSoleusPassiveForceLengthData = [];

%We have no data to fit to, and so these options are not used
flag_solveForOptimalFiberLengthOfBestFit         = 0; 
shiftLengthActiveForceLengthCurveDescendingCurve = 0.;

[humanSoleusNormMuscleCurvesDefault,...
 humanSoleusMusculotendonPropertiesDefault,...
 humanSoleusSarcomerePropertiesDefault,... 
 activeForceLengthCurveAnnotationPoints,...
 humanSoleusActiveForceLengthDataDefault,...
 humanSoleusPassiveForceLengthDataDefault,...
 humanSoleusPassiveForceLengthCurveSettings]= ...
    createFittedMuscleCurves( ...
      humanSoleusMusculotendonProperties,...
      humanSoleusSarcomereProperties,...
      humanSoleusActiveForceLengthData,...
      humanSoleusPassiveForceLengthData,...
      shiftLengthActiveForceLengthCurveDescendingCurve,...
      flag_enableNumericallyNonZeroGradients,...
      smallNumericallyNonZeroNumber,...
      flag_solveForOptimalFiberLengthOfBestFit,...
      createMusculoTendonFcn,...
      flag_useOctave);



defaultHumanSoleus = struct('musculotendon',...
                            humanSoleusMusculotendonPropertiesDefault,...
                            'sarcomere',...
                            humanSoleusSarcomerePropertiesDefault,...
                            'falData',...
                            humanSoleusActiveForceLengthDataDefault,...
                            'fpeData',...
                            humanSoleusPassiveForceLengthDataDefault,...
                            'curves',...
                            humanSoleusNormMuscleCurvesDefault);
                      
save('output/structs/defaultHumanSoleus.mat',...
     'defaultHumanSoleus');  



if(isempty(elasticTendonReferenceModel)==0)
    disp('Using reference model');        
    tmp=load(elasticTendonReferenceModel);
    modelName = fields(tmp);
    humanSoleusSarcomerePropertiesUpd_ET = tmp.(modelName{1}).sarcomere;
    humanSoleusNormMuscleCurvesUpd_ET    = tmp.(modelName{1}).curves;
else
    disp('Using default model');
    humanSoleusSarcomerePropertiesUpd_ET = humanSoleusSarcomerePropertiesDefault;
    humanSoleusNormMuscleCurvesUpd_ET    = humanSoleusNormMuscleCurvesDefault;
end


if(isempty(rigidTendonReferenceModel)==0)
    disp('Using default model');        
    tmp=load(rigidTendonReferenceModel);
    modelName = fields(tmp);
    humanSoleusSarcomerePropertiesUpd_RT = tmp.(modelName{1}).sarcomere;
    humanSoleusNormMuscleCurvesUpd_RT    = tmp.(modelName{1}).curves;
else
    disp('Using reference model');
    humanSoleusSarcomerePropertiesUpd_RT = humanSoleusSarcomerePropertiesDefault;
    humanSoleusNormMuscleCurvesUpd_RT    = humanSoleusNormMuscleCurvesDefault;

end



figH = plotStructOfBezierSplines( humanSoleusNormMuscleCurvesUpd_ET,...
                                  {'Inverse','use'});                          


   
%%
% Remove the directories ...
%%
rmpath(parametersDirectoryTreeMTParams);
rmpath(parametersDirectoryTreeExperiments);
rmpath(parametersDirectoryTreeModels);
rmpath(parametersDirectoryTreeCurves);
