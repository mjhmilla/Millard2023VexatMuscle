%This flag allows us to avoid the memory clearing functions so that
%this can be timed using tic and tock from within main_OuterLoop
flag_OuterOuterLoopMode =1;
if(flag_OuterOuterLoopMode ==0)
    clc;
    close all;
    clear all;
end

rootDir         = getRootProjectDirectory();
projectFolders  = getProjectFolders(rootDir);

disp('----------------------------------------');
disp(' running main_CreateModels_OuterLoop');
disp('----------------------------------------');
disp('   :run-time: ~26  minutes with titin fitting ');
disp('             2-3 minutes without titin fitting');
disp('            *Intel i7-3630QM @ 2.40 GHz, Ubuntu 22');
disp '             8 GB ram, SSD harddrive');
disp('----------------------------------------');


%%
% Global model parameters
%%
flag_loadPreviouslyOptimizedParameters = 0;
% 0: A lengthy optimization is done to find the best point within the
%    PEVK segment to bond to actin. This value is saved to file for later
%    use.
% 1: Previously calculated parameters are loaded. 

%It has been a long time since this script has been tested on Octave: it 
%probably no longer works.
flag_useOctave            = 0; 

flag_makeAndSavePubPlots  = 1;
plotOutputFolder          = [projectFolders.output_plots_MuscleCurves,filesep];

normMaxActiveTitinToActinDamping = 65;

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

%
%0: Parts of curves that go to zero in reality will go to zero in the curves.
%   For example, when the tendon is slack the force-length curve in this segment
%   will be a perfectly horizonal time.
%1: Parts of curves that go to zero in reality will have a slope of 
%   smallNumericallyNonZeroNumber. This is useful if the muscle model is 
%   being used in a gradient based optimization routine.

wlcTitinModel = 1;
% This titin model type will extend the force length curve of each titin
% segment to fit the WLC model up to a high force without going to a 
% singularity.

linearTitinModel=0;
% This titin model type will linear extrapolate the force length curve
% from the length at which the muscle develops a passive force equivalent
% to one maximum sometric force.

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
% Feline Soleus (additional flags)
%%

flag_plotAllDefaultFelineSoleusCurves = 0;

%%
% Rabbit Psoas Model parameters
%%
flag_plotAllRabbitPsoasFibrilCurves     = 0;
scaleOptimalFiberLengthRabbitPsoas      = 1;
scaleMaximumIsometricTensionRabbitPsoas = 1;

%%
% Human Soleus Model parameters
%%
flag_plotAllHumanSoleusCurves           = 0;
scaleOptimalFiberLengthHumanSoleus      = 1.0; 
scaleMaximumIsometricTensionHumanSoleus = 1;


%%
% Cat Soleus Model Parameters
%%

flag_fitFelineSoleusActiveTitinProperties               = nan; 
flag_loadFittedFelineSoleusActiveTitinProperties        = nan;
flag_fitFelineCrossbridgeProperties                     = nan;
flag_loadFittedFelineCrossbridgeProperties              = nan;


if(flag_loadPreviouslyOptimizedParameters==1)
    flag_fitFelineSoleusActiveTitinProperties               = 0; 
    flag_loadFittedFelineSoleusActiveTitinProperties        = 1;
    flag_fitFelineCrossbridgeProperties                     = 0;
    flag_loadFittedFelineCrossbridgeProperties              = 1;

else
    flag_fitFelineSoleusActiveTitinProperties               = 1; 
    flag_loadFittedFelineSoleusActiveTitinProperties        = 0;
    flag_fitFelineCrossbridgeProperties                     = 1;
    flag_loadFittedFelineCrossbridgeProperties              = 0;

end

scaleOptimalFiberLengthCatSoleus        = 1.0; 
scaleMaximumIsometricTensionCatSoleus   = 1;

rigidTendonReferenceModelCatSoleus      = [];
elasticTendonReferenceModelCatSoleus    = [];


%%
% Paths
%%
addpath( genpath(projectFolders.parameters)     );
addpath( genpath(projectFolders.curves)         );
addpath( genpath(projectFolders.experiments)    );
addpath( genpath(projectFolders.simulation)     );
addpath( genpath(projectFolders.models)         );
addpath( genpath(projectFolders.postprocessing) );



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
                normMaxActiveTitinToActinDamping,...
                normFiberLengthAtOneNormPassiveForceDefault,... 
                ecmForceFractionFelineSoleus,...
                linearTitinModel,...
                useCalibratedCurves,...
                useTwoSidedTitinCurves,...
                smallNumericallyNonZeroNumber,...
                flag_enableNumericallyNonZeroGradients,...
                scaleOptimalFiberLengthCatSoleus,... 
                scaleMaximumIsometricTensionCatSoleus,...
                projectFolders,...
                flag_useOctave);


filePathDefault = fullfile(   projectFolders.output_structs_FittedModels,...
                                    'defaultFelineSoleus.mat');

save(filePathDefault,'defaultFelineSoleus');  

if(flag_plotAllDefaultFelineSoleusCurves==1)
    figHumanSoleusCurves = ...
    plotStructOfBezierSplines( defaultFelineSoleus.curves,...
                                      {'Inverse','use'});       
end

fittedFelineSoleus = [];
fittingTag = '';
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
            projectFolders,...
            flag_useOctave);

    fittedFelineSoleus = fittedFelineSoleusHL2002_RT;
    fittingTag = 'HL2002';
    
    filePathHL2002RT = fullfile(   projectFolders.output_structs_FittedModels,...
                                    ['fittedFelineSoleus',fittingTag,'_RT']);

    save(filePathHL2002RT,'fittedFelineSoleus');
    %save(['output/structs/fittedFelineSoleus',fittingTag,'_RT'],...
    %        'fittedFelineSoleus');

    flag_useElasticTendon = 1;
    fittedFelineSoleusHL2002_ET = ...
        fitFelineSoleusPevkActinBondLocation( ...
            defaultFelineSoleus,...
            flag_useElasticTendon,...          
            felineSoleusPassiveForceLengthCurveSettings,...
            projectFolders,...
            flag_useOctave);

    fittedFelineSoleus = fittedFelineSoleusHL2002_ET;
    fittingTag = 'HL2002';

    filePathHL2002ET = fullfile(   projectFolders.output_structs_FittedModels,...
                                    ['fittedFelineSoleus',fittingTag,'_ET']);
    save(filePathHL2002ET,'fittedFelineSoleus');

    %save(['output/structs/fittedFelineSoleus',fittingTag,'_ET'],...
    %        'fittedFelineSoleus');

end

if(flag_loadFittedFelineSoleusActiveTitinProperties==1)
    disp([' Loading: elastic and rigid tendon models previously fit using ']);
    disp(['   fitFelineSoleusPevkActinBondLocation']);

    fittingTag = 'HL2002';
    filePathHL2002RT = fullfile(   projectFolders.output_structs_FittedModels,...
                                    ['fittedFelineSoleus',fittingTag,'_RT']);
    tmp=load(filePathHL2002RT);    
    %tmp=load(['output/structs/fittedFelineSoleus',fittingTag,'_RT']);

    fittedFelineSoleusHL2002_RT=tmp.fittedFelineSoleus;
    

    filePathHL2002ET = fullfile(   projectFolders.output_structs_FittedModels,...
                                    ['fittedFelineSoleus',fittingTag,'_ET']);


    tmp=load(filePathHL2002ET);
    %tmp=load(['output/structs/fittedFelineSoleus',fittingTag,'_ET']);

    fittedFelineSoleusHL2002_ET=tmp.fittedFelineSoleus;
    
else
    disp([' Using default feline solues ']);        
end


%%
% Numerically solve for the cross-bridge viscoelastic stiffness and
% damping so that the frequency response of the model matches Kirsch
% et al.'s analysis as closely as possible. This process does not 
% require simulation but just instead solving of a QP and completes in
% seconds.
%%
fittedFelineSoleusKBR1994Fig3_ET  = [];
fittedFelineSoleusKBR1994Fig3_RT  = [];
fittedFelineSoleusKBR1994Fig12_ET = [];
fittedFelineSoleusKBR1994Fig12_RT = [];

flag_alreadyFitted = 0;
fittingTag = [fittingTag,'KBR1994'];


perturbationFrequencyName = '90Hz';
% Kirsch et al. report the frequency response of muscle to 15 Hz and 90 Hz
% stochastic perturbations in Figure 3. In Figure 9 and 10 the stiffness 
% and damping coefficients are reported for the 15Hz, 35Hz, and 90Hz 
% perturbations. In Figure 12 the stiffness and damping coefficients that 
% best fit the 35 Hz

gainPhaseTypeNames    = {'Fig3','Fig12'};
gainPhaseTypeValues   = [3,12];
tendonTypeNames       = {'_RT','_ET'};
tendonTypeValues      = [0,1];

for indexGainPhase = 1:1:length(gainPhaseTypeNames)


    for indexTendon = 1:1:length(tendonTypeNames)

        switch tendonTypeValues(1,indexTendon)
            case 0
                if(flag_fitFelineSoleusActiveTitinProperties || ...
                        flag_loadFittedFelineSoleusActiveTitinProperties)
                    fittedFelineSoleusStarting=fittedFelineSoleusHL2002_RT;
                else
                    fittedFelineSoleusStarting=defaultFelineSoleus;
                end
            case 1
                if(flag_fitFelineSoleusActiveTitinProperties || ...
                        flag_loadFittedFelineSoleusActiveTitinProperties)
                    fittedFelineSoleusStarting=fittedFelineSoleusHL2002_ET;
                else
                    fittedFelineSoleusStarting=defaultFelineSoleus;
                end
    
            otherwise assert(0,'Error: Incorrect tendon type');
        end

        flag_useElasticTendon           = tendonTypeValues(1,indexTendon);
        figNameGainPhase                = gainPhaseTypeNames{indexGainPhase};

        fittedFelineSoleus = ...
            fitFelineSoleusCrossbridgeViscoelasticity( ...
                fittedFelineSoleusStarting,...
                flag_useElasticTendon,...
                figNameGainPhase,...
                perturbationFrequencyName,...
                felineSoleusActiveForceLengthDataDefault,...
                felineSoleusPassiveForceLengthCurveSettings,...
                flag_useOctave,...
                projectFolders);

        
        filePathKBR1994 = fullfile(projectFolders.output_structs_FittedModels,...
            ['fittedFelineSoleus',fittingTag,...
            figNameGainPhase,tendonTypeNames{indexTendon},'.mat']);
        save(filePathKBR1994,'fittedFelineSoleus');

        %save([['output/structs/fittedFelineSoleus',fittingTag],...
        %    figNameGainPhase,tendonTypeNames{indexTendon},'.mat'],...
        %      'fittedFelineSoleus')

        typeNumber = (tendonTypeValues(1,indexTendon))*100 ...
                     + gainPhaseTypeValues(1,indexGainPhase);

        switch typeNumber
            case 3
                fittedFelineSoleusKBR1994Fig3_RT      = fittedFelineSoleus;
            case 12
                fittedFelineSoleusKBR1994Fig12_RT     = fittedFelineSoleus;
            case 103
                fittedFelineSoleusKBR1994Fig3_ET      = fittedFelineSoleus;
            case 112
                fittedFelineSoleusKBR1994Fig12_ET     = fittedFelineSoleus;

            otherwise
                assert(0,'Error: invalid model configuration');
        end
    end

end


%%
% Rabbit psoas fibril Model
%%

fprintf('\n\nCreating: default rabbit psoas fibril model\n');
fprintf('  used to simulate Leonard, Joumaa, Herzog 2010.\n\n');

normCrossBridgeStiffness    = fittedFelineSoleusKBR1994Fig12_RT.sarcomere.normCrossBridgeStiffness;
normCrossBridgeDamping      = fittedFelineSoleusKBR1994Fig12_RT.sarcomere.normCrossBridgeDamping;

titinMolecularWeightInkDDefault =[];



% An adjusted value to bring the simulation data closer to the 
% experiments of Leonard, Joumaa, and Herzog 2010

ecmForceFractionRabbitPsoasFitted = 7.482620034367735e-01;
%titinMolecularWeightInkDLong = -1; 
%This is a long titin molecule at 3700 kD, equivalent to the titin
%that's found in human soleus.

normPevkToActinAttachmentPointRabbitPsoasFitted=0.675;


rabbitPsoasFibrilWLC = createRabbitPsoasFibrilModel(...
                              normCrossBridgeStiffness,...
                              normCrossBridgeDamping,...
                              normPevkToActinAttachmentPointRabbitPsoasFitted,...
                              normMaxActiveTitinToActinDamping,...
                              normFiberLengthAtOneNormPassiveForceDefault,...
                              ecmForceFractionRabbitPsoasFitted,...
                              titinMolecularWeightInkDDefault,...
                              wlcTitinModel,...
                              useCalibratedCurves,...
                              useTwoSidedTitinCurves,...
                              smallNumericallyNonZeroNumber,...
                              flag_enableNumericallyNonZeroGradients,...
                              scaleOptimalFiberLengthRabbitPsoas,...
                              scaleMaximumIsometricTensionRabbitPsoas, ...
                              projectFolders,...
                              flag_useOctave);

filePathRabbitPsoas = fullfile(projectFolders.output_structs_FittedModels,...
                                'rabbitPsoasFibrilWLC.mat');

save(filePathRabbitPsoas,'rabbitPsoasFibrilWLC');

%save('output/structs/rabbitPsoasFibrilWLC.mat',...
%     'rabbitPsoasFibrilWLC');  

if(flag_plotAllRabbitPsoasFibrilCurves==1)
    figRabbitPsoasFibrilCurves = ...
    plotStructOfBezierSplines( rabbitPsoasFibrilWLC.curves,...
                                      {'Inverse','use'});       
end




tunedRabbitPsoasFibril = createRabbitPsoasFibrilModel(...
                              normCrossBridgeStiffness,...
                              normCrossBridgeDamping,...
                              normPevkToActinAttachmentPointRabbitPsoasFitted,...
                              normMaxActiveTitinToActinDamping,...
                              normFiberLengthAtOneNormPassiveForceDefault,...
                              ecmForceFractionRabbitPsoasFitted,...
                              titinMolecularWeightInkDDefault,...
                              linearTitinModel,...
                              useCalibratedCurves,...
                              useTwoSidedTitinCurves,...
                              smallNumericallyNonZeroNumber,...
                              flag_enableNumericallyNonZeroGradients,...
                              scaleOptimalFiberLengthRabbitPsoas,...
                              scaleMaximumIsometricTensionRabbitPsoas, ...
                              projectFolders,...
                              flag_useOctave);

tunedRabbitPsoasFibril.sarcomere.normMaxActiveTitinToActinDamping = ...
    tunedRabbitPsoasFibril.sarcomere.normMaxActiveTitinToActinDamping*10;

filePathTunedRabbitPsoas = fullfile(projectFolders.output_structs_FittedModels,...
                                'tunedRabbitPsoasFibril.mat');
save(filePathTunedRabbitPsoas,'tunedRabbitPsoasFibril');

%save('output/structs/tunedRabbitPsoasFibril.mat',...
%     'tunedRabbitPsoasFibril');  

flag_plotAllTunedRabbitPsoasFibrilCurves=1;
if(flag_plotAllTunedRabbitPsoasFibrilCurves==1)
    figRabbitPsoasFibrilCurves = ...
    plotStructOfBezierSplines( tunedRabbitPsoasFibril.curves,...
                                      {'Inverse','use'});       
end

%%
% Human soleus and achilles tendon model
%%
fprintf('\n\nCreating: default human soleus model\n');
fprintf('\tused to simulate the Ig and PEVK kinematics from Trombitas et al.\n\n');



defaultHumanSoleus = createHumanSoleusModel(...
                        normPevkToActinAttachmentPointDefault,...
                        normMaxActiveTitinToActinDamping,...                        
                        normFiberLengthAtOneNormPassiveForceDefault,... 
                        ecmForceFractionHumanSoleus,...
                        linearTitinModel,...
                        useCalibratedCurves,...
                        useTwoSidedTitinCurves,...
                        smallNumericallyNonZeroNumber,...
                        flag_enableNumericallyNonZeroGradients,...
                        scaleOptimalFiberLengthHumanSoleus,...
                        scaleMaximumIsometricTensionHumanSoleus,...
                        projectFolders,...
                        flag_useOctave);

filePathHumanSoleus = fullfile(projectFolders.output_structs_FittedModels,...
                                'defaultHumanSoleus.mat');
save(filePathHumanSoleus,'defaultHumanSoleus');

%save('output/structs/defaultHumanSoleus.mat',...
%     'defaultHumanSoleus');  

if(flag_plotAllHumanSoleusCurves==1)
    figHumanSoleusCurves = ...
    plotStructOfBezierSplines( defaultHumanSoleus.curves,...
                                      {'Inverse','use'});       
end




if(flag_makeAndSavePubPlots==1)
  [success] = plotMuscleCurves( fittedFelineSoleusKBR1994Fig12_ET,...
                                defaultHumanSoleus,...
                                activeForceLengthCurveAnnotationPoints,...
                                felineSoleusActiveForceLengthDataDefault,...
                                felineSoleusPassiveForceLengthDataDefault,...
                                normFiberLengthAtOneNormPassiveForceDefault,...
                                plotOutputFolder,...
                                projectFolders);
end 


%main_createDefaultFelineSoleusModel;