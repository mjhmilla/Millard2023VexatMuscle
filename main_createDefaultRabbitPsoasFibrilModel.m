flag_outerLoopMode = 0;

if(flag_outerLoopMode == 0)
  clc;
  close all;  
  clear all;


  rigidTendonReferenceModel             = [];
  elasticTendonReferenceModel           = [];
  normPevkToActinAttachmentPoint        = 0.5;

%Taken from the data for a cat soleus. I have no passive-force-length
%measurements from a rabbit psoas that I can find in the literature, so
%for now I'm using this.  
  normFiberLengthAtOneNormPassiveForce  = 1.367732948060934e+00;
  fitCrossBridgeStiffnessDampingToKirch199490Hz = 1;

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



[rabbitPsoasFibrilMusculotendonProperties, ...
 rabbitPsoasFibrilSarcomereProperties] = ...
    createRabbitPsoasFibril(  scaleOptimalFiberLength,...
                        scaleMaximumIsometricTension,...
                        normFiberLengthAtOneNormPassiveForce,...
                        normPevkToActinAttachmentPoint,...
                        fitCrossBridgeStiffnessDampingToKirch199490Hz,...
                        flag_useOctave);

createMusculoTendonFcn = ...
  @(argScaleFiberLength,argScaleFiso)createRabbitPsoasFibril(...
                                        argScaleFiberLength,...
                                        argScaleFiso,...
                                        normFiberLengthAtOneNormPassiveForce,...
                                        normPevkToActinAttachmentPoint,...
                                        flag_useOctave); 
                                        
rabbitPsoasFibrilActiveForceLengthData  = [];
rabbitPsoasFibrilPassiveForceLengthData = [];

%We have no data to fit to, and so these options are not used
flag_solveForOptimalFiberLengthOfBestFit         = 0; 
shiftLengthActiveForceLengthCurveDescendingCurve = 0.;

[rabbitPsoasFibrilNormMuscleCurvesDefault,...
 rabbitPsoasFibrilMusculotendonPropertiesDefault,...
 rabbitPsoasFibrilSarcomerePropertiesDefault,... 
 activeForceLengthCurveAnnotationPoints,...
 rabbitPsoasFibrilActiveForceLengthDataDefault,...
 rabbitPsoasFibrilPassiveForceLengthDataDefault,...
 rabbitPsoasFibrilPassiveForceLengthCurveSettings]= ...
    createFittedMuscleCurves( ...
      rabbitPsoasFibrilMusculotendonProperties,...
      rabbitPsoasFibrilSarcomereProperties,...
      rabbitPsoasFibrilActiveForceLengthData,...
      rabbitPsoasFibrilPassiveForceLengthData,...
      shiftLengthActiveForceLengthCurveDescendingCurve,...
      flag_enableNumericallyNonZeroGradients,...
      smallNumericallyNonZeroNumber,...
      flag_solveForOptimalFiberLengthOfBestFit,...
      createMusculoTendonFcn,...
      flag_useOctave);



defaultRabbitPsoasFibril = struct('musculotendon',...
                            rabbitPsoasFibrilMusculotendonPropertiesDefault,...
                            'sarcomere',...
                            rabbitPsoasFibrilSarcomerePropertiesDefault,...
                            'falData',...
                            rabbitPsoasFibrilActiveForceLengthDataDefault,...
                            'fpeData',...
                            rabbitPsoasFibrilPassiveForceLengthDataDefault,...
                            'curves',...
                            rabbitPsoasFibrilNormMuscleCurvesDefault);
                      
save('output/structs/defaultRabbitPsoasFibril.mat',...
     'defaultRabbitPsoasFibril');  



if(isempty(elasticTendonReferenceModel)==0)
    disp('Using reference model');        
    tmp=load(elasticTendonReferenceModel);
    modelName = fields(tmp);
    rabbitPsoasFibrilSarcomerePropertiesUpd_ET = tmp.(modelName{1}).sarcomere;
    rabbitPsoasFibrilNormMuscleCurvesUpd_ET    = tmp.(modelName{1}).curves;
else
    disp('Using default model');
    rabbitPsoasFibrilSarcomerePropertiesUpd_ET = rabbitPsoasFibrilSarcomerePropertiesDefault;
    rabbitPsoasFibrilNormMuscleCurvesUpd_ET    = rabbitPsoasFibrilNormMuscleCurvesDefault;
end


if(isempty(rigidTendonReferenceModel)==0)
    disp('Using default model');        
    tmp=load(rigidTendonReferenceModel);
    modelName = fields(tmp);
    rabbitPsoasFibrilSarcomerePropertiesUpd_RT = tmp.(modelName{1}).sarcomere;
    rabbitPsoasFibrilNormMuscleCurvesUpd_RT    = tmp.(modelName{1}).curves;
else
    disp('Using reference model');
    rabbitPsoasFibrilSarcomerePropertiesUpd_RT = rabbitPsoasFibrilSarcomerePropertiesDefault;
    rabbitPsoasFibrilNormMuscleCurvesUpd_RT    = rabbitPsoasFibrilNormMuscleCurvesDefault;

end



figH = plotStructOfBezierSplines( rabbitPsoasFibrilNormMuscleCurvesUpd_ET,...
                                  {'Inverse','use'});                          


   
%%
% Remove the directories ...
%%
rmpath(parametersDirectoryTreeMTParams);
rmpath(parametersDirectoryTreeExperiments);
rmpath(parametersDirectoryTreeModels);
rmpath(parametersDirectoryTreeCurves);
