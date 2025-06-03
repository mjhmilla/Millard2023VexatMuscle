function [success] = runTomalkaWeidnerHahnSieberlSiebert2021SimulationsVEXAT( ...
                          expDataTWHSS2021,...
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


numberOfSimulations     = 4;
benchRecord             = [];

vmax = musculotendonProperties.maximumNormalizedFiberVelocity;

lceOptMdl   = musculotendonProperties.optimalFiberLength;
vceMaxMdl   = musculotendonProperties.maximumNormalizedFiberVelocity;
lceOptData  = lceOptMdl;

npts=400;

for idxSimulation = 1:1:numberOfSimulations

    %%          
    % Setup the simulation structs
    %%   


    lengthProfile.time      = [];
    lengthProfile.length    = [];


    excitationProfile.time = [];
    excitationProfile.value =[];

    timeWait         = 0.1;
    
    if(idxSimulation <= 3)
            lenA     = 2.0;
            lenB     = 2.4;
            lenC     = 2.0; 

            vceA = 0;
            vceB = 0;
            vceC = 0;
            switch idxSimulation
                case 1
                    vceA =  vceMaxMdl*lceOptMdl*0.85;
                    vceB = -vceMaxMdl*lceOptMdl*0.85;                    
                case 2
                    vceA =  vceMaxMdl*lceOptMdl*0.6;
                    vceB = -vceMaxMdl*lceOptMdl*0.6;                    
                case 3
                    vceA =  vceMaxMdl*lceOptMdl*0.3;
                    vceB = -vceMaxMdl*lceOptMdl*0.3;                    
                    
                otherwise assert(0,'Error: invalid idxSimulation');
            end

            timeStart       = 0;
            timeStimulation = timeStart + timeWait;
            timeA           = timeStimulation + timeWait;
            timeB           = timeA + (lenB-lenA)/vceA;
            timeC           = timeB + (lenC-lenB)/vceB;
            timeEnd         = timeC + timeWait;
            
            lengthProfile.time      = [timeStart,timeA,timeB,timeC,timeEnd]';
            lengthProfile.length    = [lenA,lenA,lenB,lenC,lenC]';
        
%             excitationProfile.time = [timeStart,(timeStimulation-1e-6),...
%                                       timeStimulation,timeEnd]';
%             excitationProfile.value =[0,0,1,1]';   
            excitationProfile.time = [timeStart,timeEnd]';
            excitationProfile.value =[1,1]';  

            timeLength     = [0:0.01:timeEnd]';
            valueLength   = zeros(length(timeLength),2);

            timeEx = timeLength;
            valueEx = zeros(length(timeLength),2);
                            
    end

    if(idxSimulation > 3)
            lenA     = 2.4;
            lenB     = 2.4;
            lenC     = 2.0; 

            vceA =  0;
            vceB = -vceMaxMdl*lceOptMdl*0.85;                    
            vceC =  0;

            timeStart       = 0;
            timeStimulation = timeStart + timeWait;
            timeA           = timeStimulation + timeWait;
            timeB           = timeA + abs(lenC-lenB)/abs(vceB);
            timeC           = timeB + (lenC-lenB)/(vceB);
           

            timeEnd         = timeC + timeWait;

            lengthProfile.time      = [timeStart,timeA,timeB,timeC,timeEnd]';
            lengthProfile.length    = [lenA,lenA,lenB,lenC,lenC]';
        
            %excitationProfile.time = [timeStart,(timeStimulation-1e-6),...
            %                          timeStimulation,timeEnd]';
            %excitationProfile.value =[0,0,1,1]';   

            excitationProfile.time = [timeStart,timeEnd]';
            excitationProfile.value =[1,1]';  

            timeLength     = [0:0.01:timeEnd]';
            valueLength   = zeros(length(timeLength),2);

            timeEx = timeLength;
            valueEx = zeros(length(timeLength),2);
                     
    end
    
    flag_debug=0;
    if(flag_debug == 1)
        for i=1:1:length(timeEx)
            tmp = calcLinearlyInterpolatedState(timeLength(i,1),...
                        lengthProfile.time,...
                        lengthProfile.length,1);

            valueLength(i,1)=tmp(2,1);
            valueLength(i,2)=tmp(1,1);

            tmp = calcLinearlyInterpolatedState(timeLength(i,1),...
                        excitationProfile.time,...
                        excitationProfile.value,0);

            valueEx(i,1) = tmp(1,1);
        end

        figTest=figure;
        subplot(1,2,1);
            yyaxis left;
            plot(timeLength,valueLength(:,1),'b');
            hold on;
            ylabel('Length');
            yyaxis right;
            plot(timeLength,valueLength(:,2),'r');
            hold on;
            xlabel('Time (s)');
            ylabel('Velocity');
            title('Length function');
        subplot(1,2,2);
            plot(timeEx,valueEx(:,1),'b');
            hold on;
            ylabel('Excitation');
            xlabel('Time (s)');
            ylabel('Excitation');
            title('Excitation function');
    end   

    timeSpan = [timeStart,timeEnd];

    pathLengthFcn = @(argT)calcLinearlyInterpolatedState(argT,...
                lengthProfile.time,lengthProfile.length,1);


    excitationFcn = @(argT)calcLinearlyInterpolatedState(argT,...
                excitationProfile.time,excitationProfile.value,0);

    activation=1;
    
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
    benchConfig.initialActivation     = excitationFcn(0);
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
