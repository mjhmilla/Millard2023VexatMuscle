clc;
close all;
clear all;

%%
% Global model parameters
%%
flag_useOctave                        = 0; 
%It's been a long time since this script has been tested on Octave: it 
%probably no longer works.

flag_makeAndSavePubPlots              = 1;
plotOutputFolder                      = '../output/plots/MuscleCurves/';

normPevkToActinAttachmentPointDefault = 0.5;
%0: Use previously computed values
%1: Solve for the titin-actin attachment point and active damping
%   coefficent (this is time consuming ~ 30 min)

normFiberLengthAtOneNormPassiveForceDefault = 1.367732948060934e+00;
% This the normalized fiber length at which the passive-force-length curve 
% used in this work develops 1 maximum isometric force passively when fit to 
% passive-force length data from Fig. 7 of Herzog & Leonard. This is used as 
% a default value for the other models used in this work for which the 
% passive force length properties are unknown.
%
% Herzog W, Leonard TR. Force enhancement following stretching of skeletal 
% muscle: a new mechanism. Journal of Experimental Biology. 
% 2002 May 1;205(9):1275-83.

smallNumericallyNonZeroNumber           = sqrt(sqrt(eps));
flag_enableNumericallyNonZeroGradients  = 1;
%0: Parts of curves that go to zero in reality will go to zero in the curves.
%   For example, when the tendon is slack the force-length curve in this segment
%   will be a perfectly horizonal time.
%1: Parts of curves that go to zero in reality will have a slope of 
%   smallNumericallyNonZeroNumber. This is useful if the muscle model is 
%   being used in a gradient based optimization routine.


useCalibratedCurves     = 1;
useTwoSidedTitinCurves  = 0;

%%
% Rabbit Psoas Model parameters
%%
flag_plotAllRabbitPsoasFibrilCurves     = 0;
scaleOptimalFiberLengthRabbitPsoas      = 1;
scaleMaximumIsometricTensionRabbitPsoas = 1;

%%
% Human Soleus Model parameters
%%
scaleOptimalFiberLengthHumanSoleus      = 1.0; 
scaleMaximumIsometricTensionHumanSoleus = 1;
flag_plotAllHumanSoleusCurves           = 0;


%%
% Cat Soleus Model Parameters
%%
flag_fitFromScratchCatSoleusModel                       = 1;

flag_fitFelineSoleusActiveTitinProperties               = 1;
%Takes 10-20 minutes, but must be done once.
%Numerically identifies the point of attachement between the PEVK segment
%and actin that produces simulated forces that most closely matches 
%Herzog & Leonard 2002.

flag_fitToFig3KirchBoskovRymer1994                      = 0;
flag_fitCrossBridgeStiffnessDampingTo90HzDataKirch1994  = 1;

scaleOptimalFiberLengthCatSoleus        = 1.0; 
scaleMaximumIsometricTensionCatSoleus   = 1;

rigidTendonReferenceModelCatSoleus      = [];
elasticTendonReferenceModelCatSoleus    = [];

%if(flag_fitFromScratchCatSoleusModel==0)
%    rigidTendonReferenceModelCatSoleus =   ...
%        'output/structs/felineSoleusRigidTendonKBR1994Fig12.mat';
%    elasticTendonReferenceModelCatSoleus= ...
%        'output/structs/felineSoleusElasticTendonKBR1994Fig12.mat';
%end

%%
% Paths
%%

parametersDirectoryTreeMTParams     = genpath('parameters');
parametersDirectoryTreeExperiments  = genpath('experiments');
parametersDirectoryTreeModels       = genpath('models');
parametersDirectoryTreeCurves       = genpath('curves');
parametersDirectoryTreeSimulation   = genpath('simulation');
postprocessingDirectoryTree         = genpath('postprocessing');

addpath(parametersDirectoryTreeMTParams     );
addpath(parametersDirectoryTreeExperiments  );
addpath(parametersDirectoryTreeModels       );
addpath(parametersDirectoryTreeCurves       );
addpath(parametersDirectoryTreeSimulation   );
addpath(postprocessingDirectoryTree         );


%%
% Plot configuration
%%
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end



%%
% Rabbit psoas fibril Model
%%
disp('Creating: default rabbit-psoas fibril model');
disp('  used to simulate Leonard, Joumaa, Herzog 2010.');

defaultRabbitPsoasFibril = createRabbitPsoasFibrilModel(...
                              normPevkToActinAttachmentPointDefault,...
                              normFiberLengthAtOneNormPassiveForceDefault,...
                              useCalibratedCurves,...
                              useTwoSidedTitinCurves,...
                              smallNumericallyNonZeroNumber,...
                              flag_enableNumericallyNonZeroGradients,...
                              scaleOptimalFiberLengthRabbitPsoas,...
                              scaleMaximumIsometricTensionRabbitPsoas, ...
                              flag_useOctave);

save('output/structs/defaultRabbitPsoasFibril.mat',...
     'defaultRabbitPsoasFibril');  


if(flag_plotAllRabbitPsoasFibrilCurves==1)
    figRabbitPsoasFibrilCurves = ...
    plotStructOfBezierSplines( defaultRabbitPsoasFibril.curves,...
                                      {'Inverse','use'});       
end

%%
% Human soleus and achilles tendon model
%%
disp('Creating: default human soleus model');
disp('  used to simulate the Ig and PEVK kinematics from Trombitas et al.');
defaultHumanSoleus = createHumanSoleusModel(...
                        normPevkToActinAttachmentPointDefault,...
                        normFiberLengthAtOneNormPassiveForceDefault,... 
                        useCalibratedCurves,...
                        useTwoSidedTitinCurves,...
                        smallNumericallyNonZeroNumber,...
                        flag_enableNumericallyNonZeroGradients,...
                        scaleOptimalFiberLengthHumanSoleus,...
                        scaleMaximumIsometricTensionHumanSoleus,...
                        flag_useOctave);

save('output/structs/defaultHumanSoleus.mat',...
     'defaultHumanSoleus');  

if(flag_plotAllHumanSoleusCurves==1)
    figHumanSoleusCurves = ...
    plotStructOfBezierSplines( defaultHumanSoleus.curves,...
                                      {'Inverse','use'});       
end



%%
% Cat soleus and tendon model
%%


disp('Creating: default feline soleus model')


[ defaultFelineSoleus,...
  activeForceLengthCurveAnnotationPoints,...
  felineSoleusActiveForceLengthDataDefault,...
  felineSoleusPassiveForceLengthDataDefault,...
  felineSoleusPassiveForceLengthCurveSettings ] = ...
        createFelineSoleusModel(...
                normPevkToActinAttachmentPointDefault,...
                normFiberLengthAtOneNormPassiveForceDefault,... 
                useCalibratedCurves,...
                useTwoSidedTitinCurves,...
                smallNumericallyNonZeroNumber,...
                flag_enableNumericallyNonZeroGradients,...
                scaleOptimalFiberLengthCatSoleus,... 
                scaleMaximumIsometricTensionCatSoleus,...
                flag_useOctave);

save('output/structs/defaultFelineSoleus.mat',...
     'defaultFelineSoleus');  


disp('You are here');

fittedFelineSoleusHL2002_ET = [];
fittedFelineSoleusHL2002_RT = [];

if(flag_fitFelineSoleusActiveTitinProperties==1)

    disp([' Running: fitFelineSoleusPevkActinBondLocation (10-20 min)']);
    disp(['   Numerically solving for the titin-actin bond location']);
    disp(['   that best fits Herzog and Leonard 2002']);

    %%
    % Numerically solve for the point of attachment within the PEVK
    % segment which produces active lengthening profiles that 
    % closely fit Herzog and Leonard 2002. This requires numerical
    % simulations that take around 10-20 minutes to complete.
    %%

    flag_useElasticTendon = 1;

    fittedFelineSoleusHL2002_ET = ...
        fitFelineSoleusPevkActinBondLocation( ...
            defaultFelineSoleus,...
            flag_useElasticTendon,...
            useCalibratedCurves,...
            useTwoSidedTitinCurves,...            
            felineSoleusPassiveForceLengthCurveSettings);

    save(['output/structs/fittedFelineSoleusHL2002_ET',...
            fittedFelineSoleusHL2002_ET]);

    flag_useElasticTendon = 0;

    fittedFelineSoleusHL2002_RT = ...
        fitFelineSoleusPevkActinBondLocation( ...
            defaultFelineSoleus,...
            flag_useElasticTendon,...
            useCalibratedCurves,...
            useTwoSidedTitinCurves,...
            felineSoleusPassiveForceLengthCurveSettings);

    save(['output/structs/fittedFelineSoleusHL2002_RT',...
            fittedFelineSoleusHL2002_RT]);



else 

    disp([' Loading: elastic and rigid tendon models previously fit using ']);
    disp(['   fitFelineSoleusPevkActinBondLocation']);


    tmp=load('output/structs/fittedFelineSoleusHL2002_ET');
    fittedFelineSoleusHL2002_ET=tmp.fittedFelineSoleusHL2002_ET;

    tmp=load('output/structs/fittedFelineSoleusHL2002_RT');
    fittedFelineSoleusHL2002_RT=tmp.fittedFelineSoleusHL2002_RT;

end


%%
% Numerically solve for the cross-bridge viscoelastic stiffness and
% damping so that the frequency response of the model matches Kirsch
% et al.'s analysis as closely as possible. This process does not 
% require simulation but just instead solving of a QP and completes in
% seconds.
%%
fittedFelineSoleusHL2002KBR1994Fig3_ET  = [];
fittedFelineSoleusHL2002KBR1994Fig3_RT  = [];
fittedFelineSoleusHL2002KBR1994Fig12_ET = [];
fittedFelineSoleusHL2002KBR1994Fig12_RT = [];


gainPhaseTypeNames  = {'Fig3','Fig12'};
gainPhaseTypeValues = [3,12];
tendonTypeNames     = {'_RT','_ET'};
tendonTypeValues    = [0,1];

for indexTendon = 1:1:length(tendonTypeNames)

    for indexGainPhase = 1:1:length(gainPhaseTypeNames)
        flag_useElasticTendon           = tendonTypeValues(1,indexTendon);
        figNameGainPhase                = gainPhaseTypeNames{indexGainPhase};

        fittedFelineSoleusHL2002KBR1994 = ...
            fitFelineSoleusCrossbridgeViscoelasticity( ...
                                        fittedFelineSoleusHL2002_ET,...
                                        flag_useElasticTendon,...
                                        figNameGainPhase,...
                                        felineSoleusActiveForceLengthDataDefault,...
                                        felineSoleusPassiveForceLengthCurveSettings);


        save(['output/structs/fittedFelineSoleusHL2002KBR1994',...
            figNameGainPhase,tendonTypeNames{indexTendon},'.mat'],...
              'fittedFelineSoleusHL2002KBR1994')

        typeNumber = (tendonTypeValues(1,indexTendon))*100 ...
                     + gainPhaseTypeValues(1,indexGainPhase);

        switch typeNumber
            case 3
                fittedFelineSoleusHL2002KBR1994Fig3_RT      = fittedFelineSoleusHL2002KBR1994;
            case 12
                fittedFelineSoleusHL2002KBR1994Fig12_RT     = fittedFelineSoleusHL2002KBR1994;
            case 103
                fittedFelineSoleusHL2002KBR1994Fig3_ET      = fittedFelineSoleusHL2002KBR1994;
            case 112
                fittedFelineSoleusHL2002KBR1994Fig12_ET     = fittedFelineSoleusHL2002KBR1994;

            otherwise
                assert(0,'Error: invalid model configuration');
        end
    end

end


if(flag_makeAndSavePubPlots==1)
  [success] = plotMuscleCurves( fittedFelineSoleusHL2002KBR1994Fig12_ET,...
                                defaultHumanSoleus,...
                                activeForceLengthCurveAnnotationPoints,...
                                felineSoleusActiveForceLengthDataDefault,...
                                felineSoleusPassiveForceLengthDataDefault,...
                                plotOutputFolder);
end 

%close all;

%flag_fitToFig3KirchBoskovRymer1994               = 1;
%fitCrossBridgeStiffnessDampingToKirch199490Hz    = 0;
%flag_fitFelineSoleusActiveTitinProperties        = 0; %Takes previously generated versions 
%flag_makeAndSavePubPlots                         = 0; 
%%

%%
% Same as before but now:
%   1. ... same as before
%   2. ... best match the frequency response of Figure 12 of Kirsch, 
%      Boskov, & Rymer as closely as possible
%
% Generate rigid-tendon and elastic-tendon cat soleus model with
%   1. The titin-actin attachment point fitted to match Herzog & Leonard
%      as closely as possible.
%   2. The lumped cross-bridge stiffness and damping tuned to best match
%      the frequency response of Figure 12 of Kirsch, Boskov, & Rymer 
%      as closely as possible.
%%
%main_createDefaultFelineSoleusModel;