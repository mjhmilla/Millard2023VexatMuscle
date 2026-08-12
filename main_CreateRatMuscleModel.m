%%
% SPDX-FileCopyrightText: 2025 Matthew Millard <millard.matthew@gmail.com>
%
% SPDX-License-Identifier: MIT
%
%%

flag_OuterLoopMode=0;

if(flag_OuterLoopMode==0)
    clc;
    close all;
    clear all;
    flag_OuterLoopMode=0;
end

rootDir         = getRootProjectDirectory();
projectFolders  = getProjectFolders(rootDir);

addpath( genpath(projectFolders.parameters)     );
addpath( genpath(projectFolders.curves)         );
addpath( genpath(projectFolders.experiments)    );
addpath( genpath(projectFolders.simulation)     );
addpath( genpath(projectFolders.models)         );
addpath( genpath(projectFolders.postprocessing) );

%%
% Load the reference data
%%

[ratMuscleData, ratMuscleMetaData] = ...
        loadRatSkeletalMuscleData(projectFolders);

%%
%Script configuration
%%
flag_makeAndSavePubPlots        =1;
flag_makeDetailedExpDataPlots   =1;

validExperiments = {'TRSS2017','TWHSS2021','WTRS2024'};
experimentName = validExperiments{1};

flag_saveIdenticalModels=1;
trialId = 0;


validMuscles = {'SOL','EDL'};

switch experimentName
    case 'TRSS2017'
        muscleName = validMuscles{2};
    case 'TWHSS2021'
        muscleName = validMuscles{1};        
    case 'WTRS2024'
        muscleName = validMuscles{1}; 
        disp('Note WTRS2024 used both SOL and EDL');
    otherwise
        assert(0,'Error: experimental series name not recognized');
end

flag_manuallySetTitinParameters = 1;
if(flag_OuterLoopMode==1)
  flag_manuallySetTitinParameters = 0;
end

manuallySet_forceLengthCurveSettings = ...
  struct('normLengthZero',0,...
         'normLengthToe',1.500101242582240,...
         'fToe',0.5,...
         'kZero',1.2207e-4,...
         'kToe',1.688112423563665,...
         'curviness',0.659629158856933,...
         'kLow',2.1791e-6);

%Least-squares fit of the optimal CE length and the passive-force-length
%relation
indexOfDataSetToFitOptCELength          =0;
indexOfDataSetToPassiveForceLengthCurve =0;

%Sets the initial stiffness and length at which the curve develops
%one norm force
%indexTRSS2017 = ratMuscleMetaData.index;


indexPassiveDataSetLinearStiffness = ratMuscleMetaData.index_SW1982;
indexPassiveDataSetToeRegion       = ratMuscleMetaData.index_TRSS2017;
% 1. Tomalka et al. 2017
% 2. Zuurbier et al. 1995
% 3. Stephenson & Williams 1982



expDataSetFittingData(3)=struct('optimalSarcomereLength',0,...
                               'minLengthWhereFpeIsLinear',0);

expDataSetFittingData(1).optimalSarcomereLength=2.525;
expDataSetFittingData(2).optimalSarcomereLength=2.525;
expDataSetFittingData(3).optimalSarcomereLength=2.525;

expDataSetFittingData(1).minLengthWhereFpeIsLinear=nan;
expDataSetFittingData(2).minLengthWhereFpeIsLinear=nan;
expDataSetFittingData(3).minLengthWhereFpeIsLinear=0.3;

expTRSS2017 = loadTRSS2017Data(projectFolders);

fittingDataSets.fpe = [];
fittingDataSets.fl  = [];
fittingDataSets.fv  = [];

expTRSS2017Sets ={'lN070','lN085','lN145'};

flN.x=[];
flN.y=[];

fvN.x=[];
fvN.y=[];

[b,a] = butter(2,30/500,'low');


for i=1:1:length(expTRSS2017Sets)

  vN = filtfilt(b,a,expTRSS2017.(expTRSS2017Sets{i}).lN);
  dt =1/expTRSS2017.sampleFrequencyHz;
  timeSeries = [0:1:(length(vN)-1)]' .* dt;

  vN = calcCentralDifferenceDataSeries(timeSeries,vN);

  idxA = expTRSS2017.keyIndices.(expTRSS2017Sets{i})(1);
  idxB = expTRSS2017.keyIndices.(expTRSS2017Sets{i})(3);
  if(i==1)
    flN.x = expTRSS2017.(expTRSS2017Sets{i}).lN(idxA); 
    flN.y = expTRSS2017.(expTRSS2017Sets{i}).fNavg(idxA);    
    fvN.x = vN(idxA);

    %This is an approximation but it works in this case 
    fvN.y = expTRSS2017.(expTRSS2017Sets{i}).fNavg(idxB) ...
           /expTRSS2017.(expTRSS2017Sets{i}).fNavg(idxA);
  else
    flN.x = [flN.x;expTRSS2017.(expTRSS2017Sets{i}).lN(idxA)]; 
    flN.y = [flN.y;expTRSS2017.(expTRSS2017Sets{i}).fNavg(idxA)];

    fvNVal= expTRSS2017.(expTRSS2017Sets{i}).fNavg(idxB) ...
           /expTRSS2017.(expTRSS2017Sets{i}).fNavg(idxA);
    fvN.x = [fvN.x;vN(idxA)];
    fvN.y = [fvN.y;fvNVal];    
  end

end

fittingDataSets.fl  = flN;
fittingDataSets.fv  = fvN;

%
% Select the reference data set
%

activeForceLengthData = [];
passiveForceLengthData = [];

if(indexOfDataSetToFitOptCELength > 0)
    referenceActiveForceLengthDataTable = ...
        ratMuscleData(indexOfDataSetToFitOptCELength).activeForceLengthData;
    for i=1:1:length(referenceActiveForceLengthDataTable)
        activeForceLengthData = ...
            [activeForceLengthData;...
             referenceActiveForceLengthDataTable(i).x, ...
             referenceActiveForceLengthDataTable(i).y];
    end
end

if(indexOfDataSetToPassiveForceLengthCurve>0)
    referencePassiveForceLengthDataTable = ...
        ratMuscleData(indexOfDataSetToPassiveForceLengthCurve).passiveForceLengthData;
    for i=1:1:length(referencePassiveForceLengthDataTable)
        passiveForceLengthData = ...
            [passiveForceLengthData;...
             referencePassiveForceLengthDataTable(i).x, ...
             referencePassiveForceLengthDataTable(i).y];            
    end
end  

%%
% Manually set sarcomere properties 
%%



makeSkinnedFibrilModel = 1;

setSarcomereProperties.normCrossBridgeStiffness         = 49.1;  %fiso/lopt
setSarcomereProperties.normCrossBridgeDamping           = 0.347; %fiso/(lopt/s)
% Manually set to produce the force transient in TRSS2017 Figure 2A
% at the beginning of lengthening


setSarcomereProperties.useWLCTitinModel                 = 0;
% 0: Using a linear titin model
% 1: Using the WLC titin model
% As described in Millard, Franklin, Herzog 2024

setSarcomereProperties.titinMolecularWeightInkD         = []; 
% I haven't been able to find a report of the molecular weight of titin
% in rat skeletal muscle.

setSarcomereProperties.ecmForceFraction                 = (1-makeSkinnedFibrilModel)*0.56;% 
% No ECM for a skinned fibril
% Otherwise use 0.56, the average amount of ECM reported in Prado et al.
% across 5 rabbit skeletal muscles

setSarcomereProperties.normPevkToActinAttachmentPoint   = 0.85; 
%0 : Prox-Ig/PEVK boundary
%1 : PEVK/Dist-Ig boundary

setSarcomereProperties.normLengthTitinActinBondMinimum = 1.0;
%This sets the shortest length at which a titin actin-bond is possible;
%

setSarcomereProperties.normMaxActiveTitinToActinDamping = 20.3; %fo/(lo/s)
% setSarcomereProperties.normMaxActiveTitinToActinDamping
%
% The half relaxation time in Figure 1 of Tomalka et al. (a stretch of ~20%
% in 0.5 s) is 
%
% thalf = 0.0313 s
%
% The half relaxation time of the cat soleus in Herzog & Leonard 2002
% Figure 7C (a stretch of 21% lopt in 0.333 s) is
% 
% thalf = 0.111 s
%
% Scaling the value of normMaxActiveTitinToActinDamping used to simulate 
% Herzog & Leonard we have
% 
% normMaxActiveTitinToActinDamping = 71.9*(0.0313 / 0.111) 
%
% Tomalka A, Weidner S, Hahn D, Seiberl W, Siebert T. Power amplification 
% increases with contraction velocity during stretch-shortening cycles of 
% skinned muscle fibers. Frontiers in physiology. 2021 Mar 31;12:644981.
%
% Herzog W, Leonard TR. Force enhancement following stretching of skeletal 
% muscle: a new mechanism. Journal of Experimental Biology. 2002 
% May 1;205(9):1275-83.

setSarcomereProperties.setTitinSlackLengthToZero = 0;
setSarcomereProperties.setT12ToZLineSegmentToZero = 0;
% In the present model the Z-line to T12 segment is modeled as a rigid
% rod. This means that as the sarcomere shortens, the titin element cannot
% attach to actin at a distance that is shorter than 0.1. However, the 
% simulations of Tomalka et al. 2017 (The force length relation is
% invisible ...) suggests that titin might be able to attach to actin
% at lengths that are shorter. By setting the Z-line to T12 length to be
% zero, this length becomes a part of the proximal Ig segment. This means
% that the proximal Ig segment will really shorten to the base of the
% actin filament when it is under a load of zero.

setSarcomereProperties.normFiberLengthAtZeroPassiveForce        = 0;
setSarcomereProperties.normFiberLengthAtOneNormPassiveForce     = 1.9;
setSarcomereProperties.normFiberStiffnessAtOneNormPassiveForce  = nan;
setSarcomereProperties.normFiberForceAtPassiveToe=1;

%%
% Manually set musculotendon properties
%%

setMusculotendonProperties.scaleOptimalFiberLength                  = 1;
setMusculotendonProperties.scaleMaximumIsometricTension             = 1;
setMusculotendonProperties.makeSkinnedFibrilModel = makeSkinnedFibrilModel;
setMusculotendonProperties.useElasticTendon                         = ...
        1 && ~setMusculotendonProperties.makeSkinnedFibrilModel;


%%
% Manually set curves
%%


setCurveProperties.useCalibratedCurves                              = 1;
setCurveProperties.useTwoSidedTitinCurves                           = 0;
setCurveProperties.smallNumericallyNonZeroNumber      = sqrt(sqrt(eps));
setCurveProperties.flag_enableNumericallyNonZeroGradients           = 1;
setCurveProperties.useForceVelocityCurveWithSlopeDiscontinuity      = 1;
setCurveProperties.fitTitinToTRSS2017Data         =  ratMuscleData(1);
setCurveProperties.activeForceLengthData          =  [];
setCurveProperties.passiveForceLengthData         =  [];
setCurveProperties.forceVelocityMultiplierAtLowEccentricFiberVelocity = 1.25;
setCurveProperties.forceVelocityMultiplierAtMaximumEccentricFiberVelocity = 1.35; 


setCurveProperties.useTitinModel2025           = 0;
% 0. Use the default PEVK segment as described in Millard, Frankling,
% Herzog 2024
% 1. Use a PEVK segment that has been modified to fit the simulations of
%    Tomalka, Weidner, Sieberl, and Siebert 2017
%



specimenTemperature     = 12; %As in 12 degrees centrigrade

flag_useOctave                          = 0;



%%
%
% Series-specific parameters
%
%%
shiftPassiveCurve = 0;

switch experimentName

    case 'TRSS2017'
        muscleName = validMuscles{2};
                
        setSarcomereProperties.normLengthTitinActinBondMinimum  = 1.70/2.525;
        setSarcomereProperties.normMaxActiveTitinToActinDamping = 200;

        setSarcomereProperties.fiberForceLengthCurviness = 0.9;

        setCurveProperties.useTitinModel2025                 = 1;

        setCurveProperties.forceVelocityMultiplierAtLowEccentricFiberVelocity = 1.35;
        setCurveProperties.forceVelocityMultiplierAtMaximumEccentricFiberVelocity = 1.45; 

        setSarcomereProperties.setT12ToZLineSegmentToZero=0;
        setSarcomereProperties.setTitinSlackLengthToZero =0;

        setCurveProperties.flC.adjustAscendingLimbNormalizedForce = 0.0;%0.2;
        setCurveProperties.flC.useRoundAscendingLimb = 0;
        %setCurveProperties.adjustActiveForceLengthCurviness = 0.5;

        setCurveProperties.shiftLengthActiveForceLengthCurveDescendingCurve=0;

        %setSarcomereProperties.normCrossBridgeStiffness         =75;%49.1;  %fiso/lopt
        %setSarcomereProperties.normCrossBridgeDamping           =0.347;%*(75/49.1); %fiso/(lopt/s)


        if(flag_saveIdenticalModels==0)
            switch trialId
                case 0
                    setSarcomereProperties.normPevkToActinAttachmentPoint   = 0.75;                 
                case 1
                    setSarcomereProperties.normPevkToActinAttachmentPoint   = 0.60; 
                case 2
                    setSarcomereProperties.normPevkToActinAttachmentPoint   = 0.6875; 
                case 3
                    setSarcomereProperties.normPevkToActinAttachmentPoint = 0.775; 
                otherwise assert(0,'Error: trialId not recognized');
            end
        end
    case 'TWHSS2021'
        muscleName = validMuscles{1};        
        setSarcomereProperties.normPevkToActinAttachmentPoint   = 0.95;         
        setSarcomereProperties.normLengthTitinActinBondMinimum  = 0.5;
        setSarcomereProperties.normMaxActiveTitinToActinDamping = 100;  
        setCurveProperties.forceVelocityMultiplierAtLowEccentricFiberVelocity ...
            = 1.025;
        setCurveProperties.forceVelocityMultiplierAtMaximumEccentricFiberVelocity ...
            = 1.05;         

        setSarcomereProperties.fiberForceLengthCurviness = 0.5;
        setCurveProperties.shiftLengthActiveForceLengthCurveDescendingCurve=0;

        setSarcomereProperties.normCrossBridgeStiffness         = 75;%49.1;  %fiso/lopt
        setSarcomereProperties.normCrossBridgeDamping           = 0.347*(75/49.1); %fiso/(lopt/s)

        shiftPassiveCurve = -0.1;

    case 'WTRS2024'
        muscleName = validMuscles{2}; 
        disp('Note WTRS2024 used both SOL and EDL');
        switch muscleName
            case 'SOL'
                setSarcomereProperties.normPevkToActinAttachmentPoint   = 0.85;         
                setSarcomereProperties.normLengthTitinActinBondMinimum  = 1.0;
                setSarcomereProperties.normMaxActiveTitinToActinDamping = 20.3;           
            case 'EDL'
                setSarcomereProperties.normPevkToActinAttachmentPoint   = 0.85;         
                setSarcomereProperties.normLengthTitinActinBondMinimum  = 1.0;
                setSarcomereProperties.normMaxActiveTitinToActinDamping = 20.3;                 
            otherwise assert(0,'Error: muscleName not recognized');
        end
    otherwise
        assert(0,'Error: experimental series name not recognized');
end


%%
% Plot configuration and data structs
%%

[plotDataConfig,...
 plotIndexes,... 
 plotSettings] = getRatMusculotendonModelPlottingStructures(...
                    experimentName,muscleName);


%%
% Rat soleus fibril Model
%%

fprintf('\n\nCreating: rat EDL model \n');
fprintf('  used to simulate %s.\n\n',experimentName);

fprintf('\n\nTo do:');
fprintf('\n\n1. Look at Prado again: there are a lot of references related to');
fprintf(  '\n   rat muscle.');
fprintf(  '\n');


%%
% Plot configuration
%%
plotLayoutSettings = struct('numberOfHorizontalPlotColumns',  3,...
                            'numberOfVerticalPlotRows',       1,...
                            'flag_fixedPlotWidth',            1,...
                            'plotWidth',                      3.75,...
                            'plotHeight',                     3.75,...
                            'flag_usingOctave',               0);

numberOfHorizontalPlotColumns = plotLayoutSettings.numberOfHorizontalPlotColumns;
numberOfVerticalPlotRows      = plotLayoutSettings.numberOfVerticalPlotRows;
flag_fixedPlotWidth           = plotLayoutSettings.flag_fixedPlotWidth;
plotWidth                     = plotLayoutSettings.plotWidth;
plotHeight                    = plotLayoutSettings.plotHeight;
flag_usingOctave              = plotLayoutSettings.flag_usingOctave;

plotHorizMarginCm = 1;
plotVertMarginCm  = 2;

pageWidth   = (plotWidth+plotHorizMarginCm)*numberOfHorizontalPlotColumns...
                +plotHorizMarginCm;

pageHeight  = (plotHeight+plotVertMarginCm)*numberOfVerticalPlotRows...
                +plotVertMarginCm;

plotConfigGeneric;

 



%%
% Human soleus and achilles tendon model
%
% This is used in calcTitinCurves2025 to make sure that the
% modified force-length relations that I'm making for 
% the PEVK segment don't break Trombitas et al.
%%

fprintf('\n\nCreating: default human soleus model\n');
fprintf('\tused to simulate the Ig and PEVK kinematics from Trombitas et al.\n\n');


useWLCTitinModel = 0;

scaleOptimalFiberLengthHumanSoleus = 1;
scaleMaximumIsometricTensionHumanSoleus = 1;

useTitinModel2025_humanReferenceModel=0;

defaultHumanSoleus = createHumanSoleusModel2025(...
                        setSarcomereProperties.normPevkToActinAttachmentPoint,...
                        setSarcomereProperties.normMaxActiveTitinToActinDamping,...                        
                        setSarcomereProperties.normFiberLengthAtOneNormPassiveForce,... 
                        setSarcomereProperties.ecmForceFraction,...
                        useWLCTitinModel,...
                        setCurveProperties.useCalibratedCurves,...
                        setCurveProperties.useTwoSidedTitinCurves,...
                        setCurveProperties.smallNumericallyNonZeroNumber,...
                        setCurveProperties.flag_enableNumericallyNonZeroGradients,...
                        scaleOptimalFiberLengthHumanSoleus,...
                        scaleMaximumIsometricTensionHumanSoleus,...
                        setSarcomereProperties.setT12ToZLineSegmentToZero,...
                        useTitinModel2025_humanReferenceModel,...
                        projectFolders,...
                        flag_useOctave);

filePathHumanSoleus = fullfile(projectFolders.output_structs_FittedModels,...
                                'defaultHumanSoleus.mat');
save(filePathHumanSoleus,'defaultHumanSoleus');

%%
%
% Rat skeletal muscle
%
%%

%
% default rat muscle model
%
[defaultRatMuscleModelParameters,...
 activeForceLengthCurveAnnotationPoints] ...
    = createRatSkeletalMuscleModel(...
          setSarcomereProperties,...
          setMusculotendonProperties,...
          setCurveProperties,...
          [],...
          specimenTemperature,...   
          experimentName,...
          muscleName,...
          projectFolders,...
          flag_useOctave);

%
% make the fitted rat muscle model
%

%%
% Fit the passive force-length relation
%   Even for fibril models this is useful as it serves as a fitting
%   reference for the passive-force-length relation for the entire 
%   titin segment.
%%

if(indexPassiveDataSetLinearStiffness > 0)
      
    referencePassiveForceLengthDataTable = ...
        ratMuscleData(indexPassiveDataSetLinearStiffness).passiveForceLengthData;

    fittingFpeNMinForce = ...
            expDataSetFittingData(...
            indexPassiveDataSetLinearStiffness).minLengthWhereFpeIsLinear;

    fittingFpeNOptSarcomereLengthInUm = ...
            expDataSetFittingData(...
            indexPassiveDataSetLinearStiffness).optimalSarcomereLength; 

    fittingDataFpeN=[];
    for i=1:1:length(referencePassiveForceLengthDataTable)
        for j=1:1:length(referencePassiveForceLengthDataTable(i).x)
            if(referencePassiveForceLengthDataTable(i).y(j) > fittingFpeNMinForce)
                fittingDataFpeN = ...
                    [fittingDataFpeN;...
                     referencePassiveForceLengthDataTable(i).x(j), ...
                     referencePassiveForceLengthDataTable(i).y(j)];            
            end
        end
    end
    
    %Fit a line to the data
    A = [fittingDataFpeN(:,1),ones(size(fittingDataFpeN(:,1)))];
    b = fittingDataFpeN(:,2);
    x = (A'*A)\(A'*b);
    c = x(1,1);
    y0 = x(2,1);
    %y = c*x + x0
    
    setSarcomereProperties.normFiberStiffnessAtOneNormPassiveForce = ...
        c/(1/fittingFpeNOptSarcomereLengthInUm);
    setSarcomereProperties.normFiberLengthAtOneNormPassiveForce = ...
        ((1-y0)/c)/fittingFpeNOptSarcomereLengthInUm;

    if(strcmp(experimentName,'TWHSS2021')==1)
        setSarcomereProperties.normFiberLengthAtOneNormPassiveForce = ...
            setSarcomereProperties.normFiberLengthAtOneNormPassiveForce ...
            +shiftPassiveCurve;
    
        setSarcomereProperties.normFiberStiffnessAtLowPassiveForce = ...
            setSarcomereProperties.normFiberStiffnessAtOneNormPassiveForce*0.4;
    end
        
    disp(setSarcomereProperties.normFiberStiffnessAtOneNormPassiveForce);
    disp(setSarcomereProperties.normFiberLengthAtOneNormPassiveForce);

    fittingDataSets.fpe(1).x = fittingDataFpeN(:,1) ...
      ./ defaultRatMuscleModelParameters.sarcomere.optimalSarcomereLength;
    fittingDataSets.fpe(1).y = fittingDataFpeN(:,2);
end

%%
% Solve for the shape of fpe to minimize the squared errors in the toe 
% region
%%
if(indexPassiveDataSetToeRegion > 0)
  assert(indexPassiveDataSetLinearStiffness > 0,...
         ['Error: this fitting function assumes that the linear',...
          ' portion of the curve has been fit to data.']);

  %%
  %
  %%

  forceLengthCurveSettings.fToe = 0.5;

  forceLengthCurveSettings.kToe = ...
    setSarcomereProperties.normFiberStiffnessAtOneNormPassiveForce;

  forceLengthCurveSettings.normLengthZero = ...
    setSarcomereProperties.normFiberLengthAtZeroPassiveForce;

  dl = (forceLengthCurveSettings.fToe-1)/forceLengthCurveSettings.kToe;
  forceLengthCurveSettings.normLengthToe  = ...
    setSarcomereProperties.normFiberLengthAtOneNormPassiveForce +dl;


  forceLengthCurveSettings.kZero = 0;
  if(setCurveProperties.flag_enableNumericallyNonZeroGradients==1)
    forceLengthCurveSettings.kZero = ...
      setCurveProperties.smallNumericallyNonZeroNumber;        
  end


  optParams = {'kLow','curviness'};
  
  optParamsLB = [1,1].*sqrt(eps);
  
  kLowMax     = 0.26; 
  % This upper bound has been found manually. It's not easy to establish
  % this upper bound because, ultimately, the value of kLow is scaled
  % and used to establish kLow for the curves that define titin's segments
  % and each of these segments has its own bounds. These bounds depend
  % on both the geometry of the target passive curve as well as the 
  % number of prox Ig domains, PEVK residues, and distal Ig domains. 
  % There is a way to evaluate this bound analytically, but its bound to
  % be an ugly difficult-to-derive expression that I don't have time for
  % at this moment.

  optParamsUB = [kLowMax,0.95];

  kLow      = 0.2;
  curviness = 0.7;
  optParamsDefault = [kLow;curviness];


  forceLengthCurveSettings.kLow       = nan;
  forceLengthCurveSettings.curviness  = nan;

  fittingFpeNOptSarcomereLengthInUm = ...
        expDataSetFittingData(...
        indexPassiveDataSetLinearStiffness).optimalSarcomereLength; 

  expDataFpeN.lceN = ...
    ratMuscleData(indexPassiveDataSetToeRegion).passiveForceLengthData.x ...
    ./ fittingFpeNOptSarcomereLengthInUm;

  expDataFpeN.fpeN = ...
    ratMuscleData(indexPassiveDataSetToeRegion).passiveForceLengthData.y;

  fittingDataSets.fpe(2).x = expDataFpeN.lceN;
  fittingDataSets.fpe(2).y = expDataFpeN.fpeN;


  flag_generateCurveSample=1;
  [errDefault,defaultCurveSample,testRatMuscleModelParameters] = ...
    calcErrorTitinForceLengthCurves(...
                        optParamsDefault,...
                        expDataFpeN,...
                        optParams,...
                        forceLengthCurveSettings,...
                        defaultRatMuscleModelParameters,...
                        projectFolders,...
                        flag_generateCurveSample);

  flag_generateCurveSample=0;
  errFcnFpeN = @(arg)calcErrorTitinForceLengthCurves(...
                        arg,...
                        expDataFpeN,...
                        optParams,...
                        forceLengthCurveSettings,...
                        defaultRatMuscleModelParameters,...
                        projectFolders,...
                        flag_generateCurveSample);

  flag_plotfpeNToeFit=1;

  if(flag_plotfpeNToeFit==1)

    figFpeNToeFit=figure;    
    plot(expDataFpeN.lceN,...
         expDataFpeN.fpeN,'xk');
    hold on;
    plot(defaultCurveSample.x,...
         defaultCurveSample.y,'-k');
    hold on;
    text(defaultCurveSample.x(end),...
         defaultCurveSample.y(end),...
         sprintf('%1.3e RMSE',sqrt(mean((errDefault).^2))),...
         'HorizontalAlignment','right',...
         'VerticalAlignment','top');
    hold on;      
    box off;
    xlabel('Norm. Length');
    ylabel('Norm. Force');
    title('Fitting of the toe region of the fpeN curve');

  end


  options = optimset('AlwaysHonorConstraints','bounds');
  [x,resnorm,residual,exitflag] = ...
    lsqnonlin(errFcnFpeN,optParamsDefault,optParamsLB,optParamsUB);
  
  assert(resnorm < norm(errDefault));
  %assert(exitflag==1);

  fprintf('\n\nOptimal fpeN (fTiN) and titin curve parameters');
  fprintf('\n\t%1.3e\tkLow',x(1));
  fprintf('\n\t%1.3e\tcurviness\n\n',x(2));

  setSarcomereProperties.normFiberStiffnessAtLowPassiveForce = x(1);
  setSarcomereProperties.fiberForceLengthCurviness           = x(2); 
  setSarcomereProperties.normFiberForceAtPassiveToe = ...
    forceLengthCurveSettings.fToe;

  flag_generateCurveSample=1;
  [errOpt,optCurveSample,updRatMuscleModelParameters] = ...
    calcErrorTitinForceLengthCurves(...
                        x,...
                        expDataFpeN,...
                        optParams,...
                        forceLengthCurveSettings,...
                        defaultRatMuscleModelParameters,...
                        projectFolders,...
                        flag_generateCurveSample);

  fprintf('\nFitted titin model parameters\n\n');
  fprintf('\t%1.4f\tlToe\tfTiN\n', ...
    updRatMuscleModelParameters.curves.fittingReferenceTitin.normLengthToe);
  fprintf('\t%1.4f\tfToe\tfTiN\n', ...
    updRatMuscleModelParameters.curves.fittingReferenceTitin.fToe);
  fprintf('\t%1.4f\tkToe\tfTiN\n\n', ...
    updRatMuscleModelParameters.curves.fittingReferenceTitin.kToe);

  fprintf('\t%1.4f\tlToe\tfIgPN\n', ...
    updRatMuscleModelParameters.curves.forceLengthIgPTitinCurve.xEnd(2));
  fprintf('\t%1.4f\tfToe\tfIgPN\n', ...
    updRatMuscleModelParameters.curves.forceLengthIgPTitinCurve.yEnd(2));
  fprintf('\t%1.4f\tkToe\tfIgPN\n\n', ...
    updRatMuscleModelParameters.curves.forceLengthIgPTitinCurve.dydxEnd(2));

  fprintf('\t%1.4f\tlToe\tfPevkN\n', ...
    updRatMuscleModelParameters.curves.forceLengthPevkTitinCurve.xEnd(2));
  fprintf('\t%1.4f\tfToe\tfPevkN\n', ...
    updRatMuscleModelParameters.curves.forceLengthPevkTitinCurve.yEnd(2));
  fprintf('\t%1.4f\tkToe\tfPevkN\n\n', ...
    updRatMuscleModelParameters.curves.forceLengthPevkTitinCurve.dydxEnd(2));

  fprintf('\t%1.4f\tlToe\tfIgDN\n', ...
    updRatMuscleModelParameters.curves.forceLengthIgDTitinCurve.xEnd(2));  
  fprintf('\t%1.4f\tfToe\tfIgDN\n', ...
    updRatMuscleModelParameters.curves.forceLengthIgDTitinCurve.yEnd(2));
  fprintf('\t%1.4f\tkToe\tfIgDN\n\n', ...
    updRatMuscleModelParameters.curves.forceLengthIgDTitinCurve.dydxEnd(2));
  
  pause(0.1);

  if(flag_plotfpeNToeFit==1)
    %sample the default curve
    plot(optCurveSample.x,...
         optCurveSample.y,'-b');
    hold on;
    text(optCurveSample.x(end),...
         optCurveSample.y(end),...
         sprintf('%1.3e RMSE',sqrt(mean((errOpt).^2))),...
         'HorizontalAlignment','left',...
         'VerticalAlignment','top',...
         'Color',[0,0,1]);
    hold on;     
    here=1;
  end

end

[ratMuscleModelParameters,...
 activeForceLengthCurveAnnotationPoints] ...
    = createRatSkeletalMuscleModel(...
          setSarcomereProperties,...
          setMusculotendonProperties,...
          setCurveProperties,...
          [],...
          specimenTemperature,...   
          experimentName,...
          muscleName,...
          projectFolders,...
          flag_useOctave);

ratMuscleModelParametersManual_fTiN=[];
if(flag_manuallySetTitinParameters==1)

 
  [ratMuscleModelParametersManual_fTiN,...
   activeForceLengthCurveAnnotationPointsManual_fTiN] ...
      = createRatSkeletalMuscleModel(...
            setSarcomereProperties,...
            setMusculotendonProperties,...
            setCurveProperties,...
            manuallySet_forceLengthCurveSettings,...
            specimenTemperature,...   
            experimentName,...
            muscleName,...
            projectFolders,...
            flag_useOctave);  
end
%
%
%

wlcStr = '';
if(setSarcomereProperties.useWLCTitinModel==1)
    wlcStr='WLC';
end
fibrilStr='';
if(setMusculotendonProperties.makeSkinnedFibrilModel==1)
    fibrilStr='Fibril';
end

if(flag_saveIdenticalModels==1)
    for i=1:1:3
        fileName = ['rat',experimentName,muscleName,...
                     fibrilStr,'ActiveTitin',wlcStr,'.mat'];
        
        if(strcmp(experimentName,'TRSS2017')==1)
            fileName = ['rat',experimentName,muscleName,...
                         fibrilStr,'ActiveTitin',wlcStr,'_',...
                         num2str(i-1),'.mat'];
        
        end
        
        
        
        
        filePathRatMuscle = fullfile(projectFolders.output_structs_FittedModels,...
                                     fileName);
        save(filePathRatMuscle,'ratMuscleModelParameters');
    end

else

    fileName = ['rat',experimentName,muscleName,...
                 fibrilStr,'ActiveTitin',wlcStr,'.mat'];
    
    if(strcmp(experimentName,'TRSS2017')==1)
        fileName = ['rat',experimentName,muscleName,...
                     fibrilStr,'ActiveTitin',wlcStr,'_',...
                     num2str(trialId),'.mat'];
    
    end
    
    
    
    
    filePathRatMuscle = fullfile(projectFolders.output_structs_FittedModels,...
                                 fileName);
    save(filePathRatMuscle,'ratMuscleModelParameters');
end






%%
% Plotting
%%
if(flag_makeAndSavePubPlots==1)
    previousPlotFullFilePathName=[];
    plotFullFilePathName = fullfile(...
       projectFolders.output_plots_MuscleCurves,...
       ['fig_Pub_RatMuscleCurves_',experimentName]);

    if(strcmp(experimentName,'TRSS2017')==1)
        plotFullFilePathName = fullfile(...
               projectFolders.output_plots_MuscleCurves,...
               ['fig_Pub_RatMuscleCurves_',experimentName,...
                '_',num2str(trialId)]);  
        if(trialId > 1)
            previousPlotFullFilePathName = fullfile(...
               projectFolders.output_plots_MuscleCurves,...
               ['fig_Pub_RatMuscleCurves_',experimentName,...
                '_',num2str(trialId-1)]); 
        end
    end


  updateTitinPlotsOnly=0;
  if(strcmp(experimentName,'TRSS2017') && trialId > 1)
      updateTitinPlotsOnly=1;
  end


 activeForceLengthDataSeriesPlot(2)=struct('x',[],'y',[],'label','');
 passiveForceLengthDataSeriesPlot(2)=struct('x',[],'y',[],'label','');

idxSeries = [ratMuscleMetaData.index_SW1982; ...
                       ratMuscleMetaData.index_TRSS2017];

labelSeries = {'SW1982 F (Exp)';'TRSS2017 F (Exp)'};
colorSeriesActive            = [0.5,0.5,0.5; ...
                                                    0,0,0];
colorSeriesActiveFace   = [1,1,1; ...
                                                    0,0,0];
colorSeriesPassive          = [0.5,0.5,0.5; ...
                                                    0,0,0];%[0.75,0.75,1; 0.25,0.25,1];
colorSeriesPassiveFace = [0.5,0.5,0.5; ...
                                                    1,1,1];%[0.75,0.75,1; 0.25,0.25,1];

activeMarkSeries    = {'.','d'};
passiveMarkSeries = {'+','o'};
seriesColors = [0.67,0.67,0.67;...
                0.33,0.33,0.33];

for i=1:1:2
  idx = idxSeries(i,1);    
  for j=1:1:length(ratMuscleData(idx).activeForceLengthData)
      activeForceLengthDataSeriesPlot(i).x = [...
          activeForceLengthDataSeriesPlot(i).x;...
          ratMuscleData(idx).activeForceLengthData(j).x];

      activeForceLengthDataSeriesPlot(i).y = [...
          activeForceLengthDataSeriesPlot(i).y;...
          ratMuscleData(idx).activeForceLengthData(j).y];

      activeForceLengthDataSeriesPlot(i).label = labelSeries{i};
      activeForceLengthDataSeriesPlot(i).color = seriesColors(i,:);%colorSeriesActive(i,:);
      activeForceLengthDataSeriesPlot(i).MarkerFaceColor = [1,1,1];%colorSeriesActiveFace(i,:);
      activeForceLengthDataSeriesPlot(i).mark = activeMarkSeries{i};
  end
  for j=1:1:length(ratMuscleData(idx).passiveForceLengthData)  
      passiveForceLengthDataSeriesPlot(i).x = [...
          passiveForceLengthDataSeriesPlot(i).x;...
          ratMuscleData(idx).passiveForceLengthData(j).x];

      passiveForceLengthDataSeriesPlot(i).y = [...
          passiveForceLengthDataSeriesPlot(i).y;...
          ratMuscleData(idx).passiveForceLengthData(j).y];

      passiveForceLengthDataSeriesPlot(i).label = labelSeries{i};      
      passiveForceLengthDataSeriesPlot(i).color = seriesColors(i,:);%colorSeriesPassive(i,:);
      passiveForceLengthDataSeriesPlot(i).MarkerFaceColor = [1,1,1];%colorSeriesPassiveFace(i,:);
      passiveForceLengthDataSeriesPlot(i).mark = passiveMarkSeries{i};

  end 
end



  flag_make2by2Plot=1;
  flag_addFpeCurves=0;
  [success] = plotRatMuscleCurvesTRSS2017(...
                ratMuscleModelParameters,...
                defaultHumanSoleus,...
                activeForceLengthCurveAnnotationPoints,...
                fittingDataSets,...
                activeForceLengthDataSeriesPlot,...
                passiveForceLengthDataSeriesPlot,...
                ratMuscleModelParameters.sarcomere.normFiberLengthAtOneNormPassiveForce,...
                trialId,...
                flag_make2by2Plot,...
                flag_addFpeCurves,...
                updateTitinPlotsOnly,...
                previousPlotFullFilePathName,...
                plotFullFilePathName,...
                projectFolders);


  if(~isempty(previousPlotFullFilePathName))
    previousPlotTitinFittingFullFilePathName = ...
        strrep(previousPlotFullFilePathName,...
        '_RatMuscleCurves_',...
        '_RatMuscleTitinFittingCurves_');
  else
    previousPlotTitinFittingFullFilePathName=[];
  end

  plotTitinFittingFullFilePathName = ...
    strrep(plotFullFilePathName,...
    '_RatMuscleCurves_',...
    '_RatMuscleTitinFittingCurves_');

  [success] = plotRatMuscleTitinFittingCurvesTRSS2017_Pub(...
              ratMuscleModelParameters,...
              ratMuscleModelParametersManual_fTiN,...
              defaultHumanSoleus,...
              activeForceLengthCurveAnnotationPoints,...
              fittingDataSets,...
              activeForceLengthDataSeriesPlot,...
              passiveForceLengthDataSeriesPlot,...
              ratMuscleModelParameters.sarcomere.normFiberLengthAtOneNormPassiveForce,...
              trialId,...
              previousPlotTitinFittingFullFilePathName,...
              plotTitinFittingFullFilePathName,...
              projectFolders);

  [success] = plotRatMuscleTitinFittingCurvesTRSS2017_Pres(...
              ratMuscleModelParameters,...              
              defaultHumanSoleus,...
              activeForceLengthCurveAnnotationPoints,...
              fittingDataSets,...
              activeForceLengthDataSeriesPlot,...
              passiveForceLengthDataSeriesPlot,...
              ratMuscleModelParameters.sarcomere.normFiberLengthAtOneNormPassiveForce,...
              trialId,...
              previousPlotTitinFittingFullFilePathName,...
              plotTitinFittingFullFilePathName,...
              projectFolders);
  
   here=1;
end 




if(flag_makeDetailedExpDataPlots==1)
    %
    % Sample the experimental datda
    %
    
    idx =  ratMuscleMetaData.index_TRSS2017;
    
    plotDataConfig(plotIndexes.TRSS2017_fl).x =...
        ratMuscleData(idx).activeForceLengthData.x;
    plotDataConfig(plotIndexes.TRSS2017_fl).y =...
        ratMuscleData(idx).activeForceLengthData.y;
    
    plotDataConfig(plotIndexes.TRSS2017_fpe).x =...
        ratMuscleData(idx).passiveForceLengthData.x;
    plotDataConfig(plotIndexes.TRSS2017_fpe).y =...
        ratMuscleData(idx).passiveForceLengthData.y;
    
    %
    % sample Zuurbier et al.
    %
    
    idx =  ratMuscleMetaData.index_ZHGL1995;
    
    yNorm = 1/100;
    xData=[];
    yData=[];
    for i=1:1:length(ratMuscleData(idx).activeForceLengthData)
        ratMuscleData(idx).activeForceLengthData(i).y = ...
            ratMuscleData(idx).activeForceLengthData(i).y.*yNorm;
        xData = [xData;ratMuscleData(idx).activeForceLengthData(i).x];
        yData = [yData;ratMuscleData(idx).activeForceLengthData(i).y];    
    end
    
    ratMuscleData(idx).passiveForceLengthData = [];
    plotDataConfig(plotIndexes.ZHGL1995_fl).x = xData;
    plotDataConfig(plotIndexes.ZHGL1995_fl).y = yData;
    
    
    %
    % sample Stephenson & Williams 
    %
    
    idx = ratMuscleMetaData.index_SW1982;
    
    xData_fl=[];
    yData_fl=[];
    xData_fpe=[];
    yData_fpe=[];
    for i=1:1:length(ratMuscleData(idx).activeForceLengthData)
        xData_fl = [xData_fl;ratMuscleData(idx).activeForceLengthData(i).x];
        yData_fl = [yData_fl;ratMuscleData(idx).activeForceLengthData(i).y];       
    end
    for i=1:1:length(ratMuscleData(idx).passiveForceLengthData)
        xData_fpe= [xData_fpe;ratMuscleData(idx).passiveForceLengthData(i).x];
        yData_fpe= [yData_fpe;ratMuscleData(idx).passiveForceLengthData(i).y]; 
    end
    
    plotDataConfig(plotIndexes.SW1982_fl).x = xData_fl;
    plotDataConfig(plotIndexes.SW1982_fl).y = yData_fl;
    plotDataConfig(plotIndexes.SW1982_fpe).x = xData_fpe;
    plotDataConfig(plotIndexes.SW1982_fpe).y = yData_fpe;
    
    
    
    %
    % Sample the model curves
    %
    
    % 
    % fl
    %
    activeForceLengthCurveData   = ...
        calcBezierYFcnXCurveSampleVector( ...
            ratMuscleModelParameters.curves.activeForceLengthCurve,...
            100,[]);
          
    lsOpt = ratMuscleModelParameters.sarcomere.optimalSarcomereLength;
    
    plotDataConfig(plotIndexes.model_fl).x = activeForceLengthCurveData.x.*lsOpt;
    plotDataConfig(plotIndexes.model_fl).y = activeForceLengthCurveData.y;
    
    %
    % fpe and titin detail
    %
    passiveForceLengthCurveData   = ...
        calcBezierYFcnXCurveSampleVector( ...
            ratMuscleModelParameters.curves.fiberForceLengthCurve, ...
            100,[]);
    
    plotDataConfig(plotIndexes.model_fpe).x = passiveForceLengthCurveData.x.*lsOpt;
    plotDataConfig(plotIndexes.model_fpe).y = passiveForceLengthCurveData.y;
    
    titinCurveSample = ...
      sampleTitinCurves20250217(...
        ratMuscleModelParameters.curves,...
        ratMuscleModelParameters.sarcomere,...
        100);
    
    
    
    
    lambdaECM = ratMuscleModelParameters.sarcomere.extraCellularMatrixPassiveForceFraction;
    
    
    plotDataConfig(plotIndexes.model_titinPassive).x = ...
        titinCurveSample.curveSampleTitin.x.*(2*lsOpt);
    plotDataConfig(plotIndexes.model_titinPassive).y = ...
        titinCurveSample.curveSampleTitin.y.*(1-lambdaECM) ...
       +titinCurveSample.curveSampleECMHalf.y.*(lambdaECM);
    
    plotDataConfig(plotIndexes.model_titinActive).x = ...
        titinCurveSample.curveSampleTitinActive.x.*(2*lsOpt);
    plotDataConfig(plotIndexes.model_titinActive).y = ...
        titinCurveSample.curveSampleTitinActive.y.*(1-lambdaECM) ...
       +titinCurveSample.curveSampleECMHalf.y.*(lambdaECM);
    
    
    
    %
    % fv
    %
    
    forceVelocityCurveData   = ...
        calcBezierYFcnXCurveSampleVector( ...
            ratMuscleModelParameters.curves.fiberForceVelocityCurve, ...
            100,[]);
    
    scaleVmax = ...
        ratMuscleModelParameters.musculotendon.maximumNormalizedFiberVelocity;
    
    plotDataConfig(plotIndexes.model_fv).x = forceVelocityCurveData.x .* scaleVmax;
    plotDataConfig(plotIndexes.model_fv).y = forceVelocityCurveData.y;
    
    
    
    %
    % Plot the overview curves
    %
    
    %
    % Force-length ticks
    %
    maxActiveSarcomereLengthInUm = ...
             2*ratMuscleModelParameters.sarcomere.zLineLengthInUm ...
            +2*ratMuscleModelParameters.sarcomere.actinLengthInUm...
            +ratMuscleModelParameters.sarcomere.myosinLengthInUm;
    
    optimalSarcomereLengthInUm = ...
             2*ratMuscleModelParameters.sarcomere.zLineLengthInUm ...
            +2*ratMuscleModelParameters.sarcomere.actinLengthInUm...
            +ratMuscleModelParameters.sarcomere.myosinBareLengthInUm;
    
    shortSarcomereLengthInUm = ...
             2*ratMuscleModelParameters.sarcomere.zLineLengthInUm ...
            +ratMuscleModelParameters.sarcomere.myosinLengthInUm;
    
    zeroForceSarcomereLengthInUm = ...
        ratMuscleModelParameters.sarcomere.zeroForceSarcomereLengthInUm;
    
    
    plotSettings(1).xticks = [];
    plotSettings(1).yticks = [];
    
    plotSettings(1).xticks = [...
        zeroForceSarcomereLengthInUm,...
        shortSarcomereLengthInUm,...
        optimalSarcomereLengthInUm,...
        maxActiveSarcomereLengthInUm];
    
    
    
    plotSettings(1).xticks =...
        [plotSettings(1).xticks, ...
         max(plotDataConfig(plotIndexes.SW1982_fpe).x)];
    
    %Evaluate the max. isometric force at the length when the myosin tip
    %touches the z-line. These expressions require a diagram to understand 
    %but I've justed checked them and these expressions are correct.
    
    %The length at which the two actin filaments overlap with
    %the active part of  myosin
    shallowPlateauInterference = ...
        2*(      ratMuscleModelParameters.sarcomere.actinLengthInUm ...
           - 0.5*ratMuscleModelParameters.sarcomere.myosinLengthInUm ...
           - 0.5*ratMuscleModelParameters.sarcomere.myosinBareLengthInUm);
    
    %The length at which the actin filaments can interact with the 
    %active part of myosin without over lap
    shallowPlateauOverlap      = ...
        2*(ratMuscleModelParameters.sarcomere.myosinLengthInUm ...
         - ratMuscleModelParameters.sarcomere.actinLengthInUm);
    
    %The maximum possible length at which myosin and actin can
    % actively interact 
    maxOverlap                 = ...
        ratMuscleModelParameters.sarcomere.myosinLengthInUm ...
      - ratMuscleModelParameters.sarcomere.myosinBareLengthInUm;
    
    %With half of the cross-bridges pulling in one direction and the other half 
    %pulling in the opposite direction the interference section contributes no force
    
    interferenceTension     = 0.0; 
    
    maxNormForceAtShortLength =...
        ( shallowPlateauInterference*interferenceTension ...
          + shallowPlateauOverlap )/maxOverlap;
    
    plotSettings(1).yticks = [0,maxNormForceAtShortLength,1];
    
    if(indexPassiveDataSetLinearStiffness>0)
        fittingFpeNMinForce = ...
            expDataSetFittingData(...
            indexPassiveDataSetLinearStiffness).minLengthWhereFpeIsLinear;
        
        plotSettings(1).yticks = ...
            [0,fittingFpeNMinForce,maxNormForceAtShortLength,1];
        plotSettings(1).yticks = sort(plotSettings(1).yticks);
    end
    
    
    
    
    %
    % Force-velocity ticks
    %
    
    plotSettings(2).xticks = [];
    plotSettings(2).yticks = [];
    
    plotSettings(2).xticks = [...
        -ratMuscleModelParameters.musculotendon.maximumNormalizedFiberVelocity,...
        -0.5*ratMuscleModelParameters.musculotendon.maximumNormalizedFiberVelocity,...
        0,...
        ratMuscleModelParameters.musculotendon.maximumNormalizedFiberVelocity];
    
    plotSettings(2).yticks = [...
        0.00,...
        ratMuscleModelParameters.musculotendon.forceVelocityMultiplierAtHalfMaximumFiberVelocity,...
        1.00,...
        ratMuscleModelParameters.musculotendon.forceVelocityMultiplierAtLowEccentricFiberVelocity,...
        ratMuscleModelParameters.musculotendon.forceVelocityMultiplierAtMaximumEccentricFiberVelocity];
    
    figModelCurves = figure;
        
    
    for i=1:1:length(plotDataConfig)
        if(plotDataConfig(i).enablePlot==1)
            row=plotDataConfig(i).row;
            col=plotDataConfig(i).col;
            subplot('Position', reshape(subPlotPanel(row,col,:),1,4));

            if(strcmp(plotDataConfig(i).Mark,'-') ...
                    || strcmp(plotDataConfig(i).Mark,'--') ...
                    || strcmp(plotDataConfig(i).Mark,'.-'))

                plot(plotDataConfig(i).x,...
                     plotDataConfig(i).y,...
                     '-',...
                     'LineWidth',2.0,...
                     'Color',[1,1,1],...
                     'HandleVisibility','off');
                hold on;
            end

            plot(plotDataConfig(i).x,...
                 plotDataConfig(i).y,...
                 plotDataConfig(i).Mark,...
                 'LineWidth',1.,...
                 'Color',plotDataConfig(i).LineColor,...
                 'MarkerFaceColor',plotDataConfig(i).MarkerFaceColor,...
                 'MarkerEdgeColor',plotDataConfig(i).MarkerEdgeColor,...
                 'MarkerSize',plotDataConfig(i).MarkerSize,...
                 'DisplayName',[plotDataConfig(i).DisplayName,' (',plotDataConfig(i).type,')'],...
                 'HandleVisibility',plotDataConfig(i).HandleVisibility);
            hold on;
        end
    end
    
    for i=1:1:length(plotSettings)
        figure(figModelCurves);
        row = plotSettings(i).row;
        col = plotSettings(i).col;
        subplot('Position', reshape(subPlotPanel(row,col,:),1,4));
        xlim(plotSettings(i).xlim);
        ylim(plotSettings(i).ylim);
        if(isempty(plotSettings(i).xticks)==0)
            xticks(round(plotSettings(i).xticks,2));
        end
        if(isempty(plotSettings(i).yticks)==0)
            yticks(round(plotSettings(i).yticks,2));
        end
    
        legend('Location',plotSettings(i).legendLocation);
        legend('boxoff');
        box off;
        xlabel(plotSettings(i).xlabel);
        ylabel(plotSettings(i).ylabel);    
    
    
        title(plotSettings(i).title);            
    end
    
    
    pause(0.1);
    figure(figModelCurves);
    pause(0.1);
    filePath = fullfile(projectFolders.output_plots_MuscleCurves,...
                        ['fig_Pub_MuscleCurves_Rat',experimentName,muscleName,'.pdf']);

    if(strcmp(experimentName,'TRSS2017')==1)
        filePath = fullfile(projectFolders.output_plots_MuscleCurves,...
                   ['fig_Pub_MuscleCurves_Rat',experimentName,...
                     muscleName,'_',num2str(trialId),'.pdf']);
    end
    print('-dpdf', filePath); 
    pause(0.1);
end