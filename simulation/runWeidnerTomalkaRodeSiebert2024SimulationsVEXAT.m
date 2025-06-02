function [success] = runWeidnerTomalkaRodeSiebert2024SimulationsVEXAT( ...
                          expDataTRSS2017,...
                          musculotendonProperties,...
                          sarcomereProperties,...
                          normMuscleCurves,...
                          outputFilePath)

disp('Running VEXAT model on Tomalka, Rode, Schumacher Siebert 2017');

%%
% Initialize the model's state appropriately
%%
assert(musculotendonProperties.tendonSlackLength==0,...
        ['Error: This is a simulation of a muscle',...
        ' fiber: tendonSlackLength should be 0']);
assert(musculotendonProperties.pennationAngle==0,...
        ['Error: This is a simulation of a muscle',...
        ' fiber: pennationAngle should be 0']);


numberOfSimulations     = 6;
benchRecord             = [];

vmax = musculotendonProperties.maximumNormalizedFiberVelocity;

lceOptData  = min(expDataTRSS2017.activeLengtheningData(3).x);
lceOptMdl   = musculotendonProperties.optimalFiberLength;

npts=400;

for idxSimulation = 1:1:numberOfSimulations

    %%          
    % Setup the simulation structs
    %%   



    timeStart       = 0;
    timeEnd         = 0;
    timeStimulation = 0;
    timeRampStart   = 0;
    timeRampEnd     = 0;
    rampLengthStart = 0;
    rampVelocity    = 0;

    timeWait = 0.5;
    
    if(idxSimulation <= 3)
            lstart  = expDataTRSS2017.activeLengtheningData(idxSimulation).x(1,1);
            lend    = expDataTRSS2017.activeLengtheningData(idxSimulation).x(end,1); 

            timeStart       = 0;
            timeStimulation = timeStart + timeWait;
            timeRampStart   = timeStimulation + timeWait;
            
            rampLengthStart = lstart;
            rampVelocity    = 0.11*vmax*lceOptData;
            timeRampEnd     = timeRampStart + (lend-lstart)/(rampVelocity);
            timeEnd         = timeRampEnd;
    end

    if(idxSimulation >= 4 && idxSimulation <= 6)
            lstart  = lceOptMdl*1.1;
            lend    = lceOptMdl; 

            timeStart       = 0;
            timeStimulation = timeStart + timeWait;
            timeRampStart   = timeStimulation + timeWait;
            
            vceNorm = 0;
            switch idxSimulation
                case 4
                    vceNorm = -0.1;
                case 5
                    vceNorm = -0.5;                    
                case 6
                    vceNorm = -0.8;                    
                otherwise assert(0,'Error: idxSimulation is not 4,5, or 6')
            end

            rampLengthStart = lstart;
            rampVelocity    = (vmax*lceOptMdl)*(vceNorm);
            timeRampEnd     = timeRampStart + (lend-lstart)/(rampVelocity);
            timeEnd         = timeRampEnd;
    end
    
    timeSpan = [timeStart,timeEnd];

    activation=1;
    excitationFcn = @(argT)calcStepFunction(argT,...
                      timeStimulation,...
                      inf,...
                      activation);
    
    pathLengthFcn = @(argT)calcRampStateSharp(argT,...
                            timeRampStart,timeRampEnd,...
                            rampLengthStart,rampVelocity);
    
    activationFcn = ...
        @(argU,argA)calcFirstOrderActivationDerivative(argU,argA, ...
            sarcomereProperties.activationTimeConstant,...
            sarcomereProperties.deactivationTimeConstant,0);
    

    
    %
    % Bench config
    %
    benchConfig.npts                  = npts;
    benchConfig.relTol                = 1e-6;
    benchConfig.absTol                = 1e-6;
    benchConfig.minActivation         = 0;
    benchConfig.color0                = [0,0,1].*0.5;
    benchConfig.color1                = [0,0,1];
    
    nStates     = 3;
    labelStates = {'$$\dot{\ell}_{a}$$', '$$\ell_{a}$$', '$$\ell_1$$'};       
    
    benchConfig.numberOfMuscleStates  = nStates;
    benchConfig.stateLabels           = labelStates;
    benchConfig.name                  = '';
    benchConfig.initialState          = [];
    benchConfig.initialActivation     = 0;
    benchConfig.pathFcn               = [];
    benchConfig.excitationFcn         = [];
    benchConfig.activationFcn         = activationFcn; 
    benchConfig.tspan                 = timeSpan;  
    
    benchConfig.useFiberDamping  = 1;
    benchConfig.useElasticTendon = 0;
    benchConfig.damping          = 0.1;
    benchConfig.iterMax          = 100;
    benchConfig.tol              = 1e-6;
    
    loopTolerance = min(benchConfig.relTol,benchConfig.absTol)/100;
    
    %
    % Initialization 
    %
    
    modelConfig = struct( ...
      'iterMax'                 , 100             , ...
      'tol'                     , loopTolerance   , ... 
      'tolInit'                 , sqrt(eps)       , ...
      'minActivation'           , 0.0             , ...
      'useElasticTendon'        , 0 , ...
      'initializeState'         , 0                     );  
    
    modelConfig.initializeState =1;
    activationState0    = [0;excitationFcn(0)];
    pathState0          = pathLengthFcn(0);
    muscleState0        = zeros(nStates,1);
    mtInfo = calcMillard2023VexatMuscleInfo(activationState0,...
                                            pathState0,...
                                            muscleState0,...
                                            musculotendonProperties,...
                                            sarcomereProperties,...
                                            normMuscleCurves,...
                                            modelConfig);
    muscleState0                = mtInfo.state.value;
    modelConfig.initializeState = 0;           
    
    benchConfig.numberOfMuscleStates = length(muscleState0);
    benchConfig.initialState         = muscleState0;
    
    benchConfig.minimumActivation    = 0;
    benchConfig.name                 = 'Vexat';
    benchConfig.eventFcn             = [];            
    
    calcMillard2023VexatMuscleInfoFcn = ...
         @(activationState1,pathState2,mclState3) ...
         calcMillard2023VexatMuscleInfo(    activationState1,...
                                            pathState2,...
                                            mclState3,...
                                            musculotendonProperties,...
                                            sarcomereProperties,...
                                            normMuscleCurves,...
                                            modelConfig);
    
    %
    % Run the simulation
    %
    idx                     = idxSimulation;
    flag_appendEnergetics   = 0;
    flag_useOctave          = 0;
    
    benchConfig.pathFcn               = pathLengthFcn;
    benchConfig.excitationFcn         = excitationFcn;             
    benchRecord = runPrescribedLengthActivationSimulation(...
                               calcMillard2023VexatMuscleInfoFcn,...
                               [],...
                               benchConfig,...
                               benchRecord,...
                               idx, ...
                               numberOfSimulations,...
                               flag_appendEnergetics,...
                               flag_useOctave);
    fprintf('%i / %i\n', idx, numberOfSimulations);                            


end

 save(outputFilePath,...
    'benchRecord',...
    'timeSpan',...
    'musculotendonProperties',...
    'sarcomereProperties',...
    'normMuscleCurves');

success=1;
