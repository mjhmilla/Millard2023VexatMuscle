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
%The default value for the point of attachment between the PEVK segment
%and actin. This point of attachment is expressed as a normalized length
%where 0 corresponds to the start of the PEVK segment (at the prox Ig/PEVK
% boundary), 0.5 would be the middle of the PEVK segment, and 1.0 would 
% be the distal end of the PEVK segment (at the PEVK/distal Ig border).

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

ecmForceFractionDefault = 0.56;
% This is the average contribution of the ECM to the passive force length
% curve. Prado et al. reports the average contribution of titin to the 
% passive stiffness of 5 muscles in a rabbit. For a default value of the
% ECM's contribution we use 1-mean(titinContribution) = 0.56:
%
%Page 472 column 2 last paragraph of Prado et al. reports:
%
%     These results show that titin’s relative contribution to
%     total passive stiffness is much higher in some muscles,
%     like psoas (57%) and diaphragm (56%), than in oth-
%     ers, like soleus (24%), EDL (42%), and gastrocne-
%     mius (41%)
%
% Prado LG, Makarenko I, Andresen C, Krüger M, Opitz CA, Linke WA. Isoform 
% diversity of giant proteins in relation to passive and active contractile 
% properties of rabbit skeletal muscles. The Journal of general physiology. 
% 2005 Nov;126(5):461-80.

ecmForceFractionRabbitPsoas     = 1-(0.728-0.158); %Fig 8A, Prado et al.
ecmForceFractionHumanSoleus     = ecmForceFractionDefault;
ecmForceFractionFelineSoleus    = ecmForceFractionDefault;

smallNumericallyNonZeroNumber           = sqrt(sqrt(eps));
flag_enableNumericallyNonZeroGradients  = 1;
%0: Parts of curves that go to zero in reality will go to zero in the curves.
%   For example, when the tendon is slack the force-length curve in this segment
%   will be a perfectly horizonal time.
%1: Parts of curves that go to zero in reality will have a slope of 
%   smallNumericallyNonZeroNumber. This is useful if the muscle model is 
%   being used in a gradient based optimization routine.


useCalibratedCurves     = 1;
% Here calibrated curves refers to an active-force-length curve and a 
% force velocity curve that have been adjusted so that opus 31 (the
% proposed model) can reproduce the desired active-force-length and force
% velocity curves. This is necessary because the deformation of the 
% viscoelastic cross-bridge element is not taken into consideration 
% in the formulation of these curves, but does affect the output.

useTwoSidedTitinCurves  = 0;
% useTwoSidedTitinCurves = 0? 
% The force length curve of the prox. IG and PEVK segments has a slack 
% length: below the slack length the force goes to zero just like a tendon.
% 
% useTwoSidedTitinCurves = 1?
% The elastic titin curves have a shape, broadly speaking, like a tan 
% function. In the case, that the titin-actin bond is formed at a long 
% sarcomere length and then the sarcomere is rapidly shortened, the distal
% section of titin may become more proximal to the Z-line than the 
% titin-actin bond (a figure is really needed to properly explain this):
% in these circumstances the PEVK segment would have a negative length
% and would generate an elastic forces that impedes further shortening.
% If this came to pass, titin would reduce the tension a muscle can generate
% during shortening. I found this not to be the case during the ramp 
% shortening simulations used in this work and so by default I set 
% useTwoSidedTitinCurves=0. It migth be the case that a shortening ramp 
% beginning at a longer length would see some additional force reduction
% due to this effect.

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
flag_fitFelineSoleusActiveTitinProperties               = 0;
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

fprintf('\n\nCreating: default rabbit psoas fibril model\n');
fprintf('  used to simulate Leonard, Joumaa, Herzog 2010.\n\n');


defaultRabbitPsoasFibril = createRabbitPsoasFibrilModel(...
                              normPevkToActinAttachmentPointDefault,...
                              normFiberLengthAtOneNormPassiveForceDefault,...
                              ecmForceFractionRabbitPsoas,...
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
fprintf('\n\nCreating: default human soleus model\n');
fprintf('\tused to simulate the Ig and PEVK kinematics from Trombitas et al.\n\n');

defaultHumanSoleus = createHumanSoleusModel(...
                        normPevkToActinAttachmentPointDefault,...
                        normFiberLengthAtOneNormPassiveForceDefault,... 
                        ecmForceFractionHumanSoleus,...
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

fprintf('\n\nCreating: default feline soleus model\n\n');


[ defaultFelineSoleus,...
  activeForceLengthCurveAnnotationPoints,...
  felineSoleusActiveForceLengthDataDefault,...
  felineSoleusPassiveForceLengthDataDefault,...
  felineSoleusPassiveForceLengthCurveSettings ] = ...
        createFelineSoleusModel(...
                normPevkToActinAttachmentPointDefault,...
                normFiberLengthAtOneNormPassiveForceDefault,... 
                ecmForceFractionFelineSoleus,...
                useCalibratedCurves,...
                useTwoSidedTitinCurves,...
                smallNumericallyNonZeroNumber,...
                flag_enableNumericallyNonZeroGradients,...
                scaleOptimalFiberLengthCatSoleus,... 
                scaleMaximumIsometricTensionCatSoleus,...
                flag_useOctave);

save('output/structs/defaultFelineSoleus.mat',...
     'defaultFelineSoleus');  



fittedFelineSoleusHL2002_ET = [];
fittedFelineSoleusHL2002_RT = [];

if(flag_fitFelineSoleusActiveTitinProperties==1)

    fprintf('\n\nRunning: fitFelineSoleusPevkActinBondLocation (10-20 min)\n');
    fprintf('\tNumerically solving for the titin-actin bond location\n');
    fprintf('\tthat best fits Herzog and Leonard 2002\n\n');

    %%
    % Numerically solve for the point of attachment within the PEVK
    % segment which produces active lengthening profiles that 
    % closely fit Herzog and Leonard 2002. This requires numerical
    % simulations that take around 10-20 minutes to complete.
    %%

    flag_useElasticTendon = 0;

    fittedFelineSoleusHL2002_RT = ...
        fitFelineSoleusPevkActinBondLocation( ...
            defaultFelineSoleus,...
            flag_useElasticTendon,...
            felineSoleusPassiveForceLengthCurveSettings,...
            flag_useOctave);

    save('output/structs/fittedFelineSoleusHL2002_RT',...
            'fittedFelineSoleusHL2002_RT');

    flag_useElasticTendon = 1;

    fittedFelineSoleusHL2002_ET = ...
        fitFelineSoleusPevkActinBondLocation( ...
            defaultFelineSoleus,...
            flag_useElasticTendon,...          
            felineSoleusPassiveForceLengthCurveSettings,...
            flag_useOctave);

    save('output/structs/fittedFelineSoleusHL2002_ET',...
            'fittedFelineSoleusHL2002_ET');





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