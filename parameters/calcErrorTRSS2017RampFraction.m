function [optError,optErrorValues, figDebugFitting,...
            ratFibrilModelsUpd,benchRecord] =...
            calcErrorTRSS2017RampFraction(...
                       optParams,...
                       optErrorValues,...
                       fittingFraction, ...
                       npts, ...
                       ratFibrilModels, ...
                       expTRSS2017Raw,...
                       simConfig,...
                       figDebugFitting,...
                       subPlotPanel,...
                       lineColorsSimTRSS2017,...
                       benchRecord,...                       
                       evaluationMethod,...
                       projectFolders)

    ratFibrilModelsUpd  = ratFibrilModels;
    errVec              = zeros(length(ratFibrilModelsUpd),npts);
    numberOfSimulations = simConfig.numberOfTrials;
    
    optError=0;

    for idxTrial = simConfig.trials
        
        %Update the model parameters        
        switch optParams.name

            case 'responseTimeScaling'
                ratFibrilModelsUpd(idxTrial).sarcomere.slidingTimeConstantLengthening= ...
                    ratFibrilModelsUpd(idxTrial).sarcomere.slidingTimeConstant...
                    *optParams.value;     
            case 'xeStiffnessDampingScaling'
                ratFibrilModelsUpd(idxTrial).sarcomere.normCrossBridgeStiffness = ...
                    ratFibrilModelsUpd(idxTrial).sarcomere.normCrossBridgeStiffness ...
                    *optParams.value;
                ratFibrilModelsUpd(idxTrial).sarcomere.normCrossBridgeDamping = ...
                    ratFibrilModelsUpd(idxTrial).sarcomere.normCrossBridgeDamping ...
                    *optParams.value;
            case 'QToF'
                ratFibrilModelsUpd(idxTrial).sarcomere.normPevkToActinAttachmentPoint = ...
                    optParams.value;

              [ratFibrilModelsUpd(idxTrial).curves.forceLengthProximalTitinCurve, ...
                  ratFibrilModelsUpd(idxTrial).curves.forceLengthProximalTitinInverseCurve,...
               ratFibrilModelsUpd(idxTrial).curves.forceLengthDistalTitinCurve, ...
                  ratFibrilModelsUpd(idxTrial).curves.forceLengthDistalTitinInverseCurve,...
               ratFibrilModelsUpd(idxTrial).curves.forceLengthIgPTitinCurve, ...
                  ratFibrilModelsUpd(idxTrial).curves.forceLengthIgPTitinInverseCurve,...
               ratFibrilModelsUpd(idxTrial).curves.forceLengthPevkTitinCurve, ...
                  ratFibrilModelsUpd(idxTrial).curves.forceLengthPevkTitinInverseCurve,...
               ratFibrilModelsUpd(idxTrial).curves.forceLengthIgDTitinCurve, ...
                  ratFibrilModelsUpd(idxTrial).curves.forceLengthIgDTitinInverseCurve,...
               ratFibrilModelsUpd(idxTrial).curves.fittingReferenceTitin] ...
                        = createTitinCurves2025( ...
                             ratFibrilModelsUpd(idxTrial).curves.fiberForceLengthCurve,...                                   
                             ratFibrilModelsUpd(idxTrial).curves.forceLengthCurveSettings,...
                             ratFibrilModelsUpd(idxTrial).curves.forceLengthECMHalfCurve,...
                             ratFibrilModelsUpd(idxTrial).sarcomere,...
                             ratFibrilModelsUpd(idxTrial).musculotendon.name,...
                             ratFibrilModelsUpd(idxTrial).curves.setTitinSlackLengthToZero,...
                             ratFibrilModelsUpd(idxTrial).curves.useWLCTitinModel,...
                             ratFibrilModelsUpd(idxTrial).curves.flag_createTwoSidedCurves,...
                             ratFibrilModelsUpd(idxTrial).curves.flag_computeCurveIntegrals,...
                             ratFibrilModelsUpd(idxTrial).curves.flag_useElasticIgD,...
                             ratFibrilModelsUpd(idxTrial).sarcomere.titinModelType,... 
                             projectFolders,...
                             ratFibrilModelsUpd(idxTrial).curves.flag_useOctave);

            case 'QToK'
                ratFibrilModelsUpd(idxTrial).sarcomere.normPevkToActinAttachmentPoint = ...
                    optParams.value;

              [ratFibrilModelsUpd(idxTrial).curves.forceLengthProximalTitinCurve, ...
                  ratFibrilModelsUpd(idxTrial).curves.forceLengthProximalTitinInverseCurve,...
               ratFibrilModelsUpd(idxTrial).curves.forceLengthDistalTitinCurve, ...
                  ratFibrilModelsUpd(idxTrial).curves.forceLengthDistalTitinInverseCurve,...
               ratFibrilModelsUpd(idxTrial).curves.forceLengthIgPTitinCurve, ...
                  ratFibrilModelsUpd(idxTrial).curves.forceLengthIgPTitinInverseCurve,...
               ratFibrilModelsUpd(idxTrial).curves.forceLengthPevkTitinCurve, ...
                  ratFibrilModelsUpd(idxTrial).curves.forceLengthPevkTitinInverseCurve,...
               ratFibrilModelsUpd(idxTrial).curves.forceLengthIgDTitinCurve, ...
                  ratFibrilModelsUpd(idxTrial).curves.forceLengthIgDTitinInverseCurve,...
               ratFibrilModelsUpd(idxTrial).curves.fittingReferenceTitin] ...
                        = createTitinCurves2025( ...
                             ratFibrilModelsUpd(idxTrial).curves.fiberForceLengthCurve,...                                   
                             ratFibrilModelsUpd(idxTrial).curves.forceLengthCurveSettings,...
                             ratFibrilModelsUpd(idxTrial).curves.forceLengthECMHalfCurve,...
                             ratFibrilModelsUpd(idxTrial).sarcomere,...
                             ratFibrilModelsUpd(idxTrial).musculotendon.name,...
                             ratFibrilModelsUpd(idxTrial).curves.setTitinSlackLengthToZero,...                             
                             ratFibrilModelsUpd(idxTrial).curves.useWLCTitinModel,...
                             ratFibrilModelsUpd(idxTrial).curves.flag_createTwoSidedCurves,...
                             ratFibrilModelsUpd(idxTrial).curves.flag_computeCurveIntegrals,...
                             ratFibrilModelsUpd(idxTrial).curves.flag_useElasticIgD,...
                             ratFibrilModelsUpd(idxTrial).sarcomere.titinModelType,... 
                             projectFolders,...
                             ratFibrilModelsUpd(idxTrial).curves.flag_useOctave);
              
            case 'f1HNPreload'
                ratFibrilModelsUpd(idxTrial).sarcomere.f1HNPreload= ...
                    optParams.value;

            case 'l1HNOffset'
                ratFibrilModelsUpd(idxTrial).sarcomere.l1HNOffset=...
                    optParams.value;

            case 'FpeQToF'
              ratFibrilModelsUpd(idxTrial).sarcomere.normPevkToActinAttachmentPoint = ...
                  optParams.value(1);
              ratFibrilModelsUpd(idxTrial).curves.forceLengthCurveSettings.kToe = ...
                  optParams.value(2);
              ratFibrilModelsUpd(idxTrial).curves.forceLengthCurveSettings.normLengthToe = ...
                  optParams.value(3);
              ratFibrilModelsUpd(idxTrial).curves.forceLengthCurveSettings.curviness = ...
                  optParams.value(4);
              if(optParams.value(4)>1 ||optParams.value(4)<0)
                here=1;
              end

              [ratFibrilModelsUpd(idxTrial).curves.forceLengthProximalTitinCurve, ...
                  ratFibrilModelsUpd(idxTrial).curves.forceLengthProximalTitinInverseCurve,...
               ratFibrilModelsUpd(idxTrial).curves.forceLengthDistalTitinCurve, ...
                  ratFibrilModelsUpd(idxTrial).curves.forceLengthDistalTitinInverseCurve,...
               ratFibrilModelsUpd(idxTrial).curves.forceLengthIgPTitinCurve, ...
                  ratFibrilModelsUpd(idxTrial).curves.forceLengthIgPTitinInverseCurve,...
               ratFibrilModelsUpd(idxTrial).curves.forceLengthPevkTitinCurve, ...
                  ratFibrilModelsUpd(idxTrial).curves.forceLengthPevkTitinInverseCurve,...
               ratFibrilModelsUpd(idxTrial).curves.forceLengthIgDTitinCurve, ...
                  ratFibrilModelsUpd(idxTrial).curves.forceLengthIgDTitinInverseCurve,...
               ratFibrilModelsUpd(idxTrial).curves.fittingReferenceTitin] ...
                        = createTitinCurves2025( ...
                             ratFibrilModelsUpd(idxTrial).curves.fiberForceLengthCurve,...                                   
                             ratFibrilModelsUpd(idxTrial).curves.forceLengthCurveSettings,...
                             ratFibrilModelsUpd(idxTrial).curves.forceLengthECMHalfCurve,...
                             ratFibrilModelsUpd(idxTrial).sarcomere,...
                             ratFibrilModelsUpd(idxTrial).musculotendon.name,...
                             ratFibrilModelsUpd(idxTrial).curves.setTitinSlackLengthToZero,...
                             ratFibrilModelsUpd(idxTrial).curves.useWLCTitinModel,...
                             ratFibrilModelsUpd(idxTrial).curves.flag_createTwoSidedCurves,...
                             ratFibrilModelsUpd(idxTrial).curves.flag_computeCurveIntegrals,...
                             ratFibrilModelsUpd(idxTrial).curves.flag_useElasticIgD,...
                             ratFibrilModelsUpd(idxTrial).sarcomere.titinModelType,... 
                             projectFolders,...
                             ratFibrilModelsUpd(idxTrial).curves.flag_useOctave);                
            case 'simulate'
                %nothing to do here.
            otherwise
                assert(0,'Error: invalid optParams.name');
        end

        lceOptMdl   = ratFibrilModelsUpd(idxTrial).musculotendon.optimalFiberLength;
        vmax        = ratFibrilModelsUpd(idxTrial).musculotendon.maximumNormalizedFiberVelocity;
        lceOptData  = lceOptMdl; 
        %min(expTRSS2017.activeLengtheningData(3).x);

        %rampLengthStart  = expTRSS2017.activeLengtheningData(idxTrial).x(1,1);
        %rampLengthEnd    = expTRSS2017.activeLengtheningData(idxTrial).x(end,1); 

        seriesName = expTRSS2017Raw.trials{idxTrial};
        rampLengthStart = expTRSS2017Raw.(seriesName).lN(1,1).*lceOptMdl;
        rampLengthEnd   = expTRSS2017Raw.(seriesName).lN(end,1).*lceOptMdl;

        lpfFreqHz = simConfig.numericalDiffLowPassFilterFrequencyHz;
        [b,a]     = butter(2,lpfFreqHz/(0.5*expTRSS2017Raw.sampleFrequencyHz),'low');
        vceNnum  = calcCentralDifferenceDataSeries(...
                      expTRSS2017Raw.(seriesName).time,...
                      filtfilt(b,a,expTRSS2017Raw.(seriesName).lN));


        timeStart       = 0;
        timeRampStart   = 0.1;
        rampVelocity    = mean(vceNnum)*lceOptData;%0.11*vmax*lceOptData;
        timeRampEnd     = timeRampStart + ...
                          fittingFraction*(...
                            rampLengthEnd-rampLengthStart)/rampVelocity;         
        timeEnd         = timeRampEnd;


        timeSpan = [timeStart,timeEnd];
    
        timeStimulation = timeStart;
        activation      = 1;
    
        excitationFcn = @(argT)calcStepFunction(argT,...
                          timeStimulation-1,...
                          inf,...
                          activation);
        
        pathLengthFcn = @(argT)calcRampStateSharp(argT,...
                                timeRampStart,timeRampEnd,...
                                rampLengthStart,rampVelocity);   

        activationFcn = ...
            @(argU,argA)calcFirstOrderActivationDerivative(argU,argA, ...
                ratFibrilModelsUpd(idxTrial).sarcomere.activationTimeConstant,...
                ratFibrilModelsUpd(idxTrial).sarcomere.deactivationTimeConstant,0);        

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
        
        %%
        % Evaluate the cost
        %%

        exactEvaluation               =0;
        approximationByInitialization =1;
        approximationByStateSetting   =2;
        approximationMinimal          =3;

        flag_compareEvaluationMethods=1;

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
        mtInfo = calcMillard2023VexatMuscleInfo(...
                    activationState0,...
                    pathState0,...
                    muscleState0,...
                    ratFibrilModelsUpd(idxTrial).musculotendon,...
                    ratFibrilModelsUpd(idxTrial).sarcomere,...
                    ratFibrilModelsUpd(idxTrial).curves,...
                    modelConfig);  

        %
        % Create the derivative evaluation function
        %
        muscleState0                = mtInfo.state.value;
        muscleState0(3,1)           = muscleState0(3,1) ...
            + ratFibrilModelsUpd(idxTrial).sarcomere.l1HNOffset;
        modelConfig.initializeState = 0;           
        
        benchConfig.numberOfMuscleStates = length(muscleState0);
        benchConfig.initialState         = muscleState0;
        
        benchConfig.minimumActivation    = 0;
        benchConfig.name                 = 'Vexat';
        benchConfig.eventFcn             = [];      
        benchConfig.pathFcn               = pathLengthFcn;
        benchConfig.excitationFcn         = excitationFcn;          
        
       
        %%
        %
        % Evalaute the model
        %
        %%
        flag_compareEvaluationMethods=0;
        benchRecordApprox= [];
        benchRecordExact = [];
        durationApprox=0;
        durationExact=0;


        if(strcmp(evaluationMethod,'approximationByStateCalculationMinimal'))

          tstart=tic;
          nPts = 100;
          timeVector = [0:(1/(nPts-1)):1]'.*(timeRampEnd-timeRampStart)...
                       + timeRampStart;  
          activationState0=[0,1];
          calcStateDerFcn = ...
            @(argt, argState)calcApproximationMillard2023VexatMuscleInfo(...
                              argt,...
                              argState,...
                              muscleState0,...           
                              activationState0,...                                            
                              pathLengthFcn, ... 
                              pathState0,...
                              ratFibrilModelsUpd(idxTrial).musculotendon,...
                              ratFibrilModelsUpd(idxTrial).sarcomere,...
                              ratFibrilModelsUpd(idxTrial).curves,...
                              modelConfig,...
                              'quasiStaticSkinnedFiberApproximation');
          [t,y]=ode15s(calcStateDerFcn,timeVector,muscleState0(3));

          if(isempty(benchRecord))
            benchRecord.normFiberForce=zeros(nPts,numberOfSimulations);
            benchRecord.normFiberLength=zeros(nPts,numberOfSimulations);
          end
          
          assert(sum(abs(t-timeVector))<1e-6);

          for idxTime=1:1:nPts

            [muscleStateDerivativeUpd, mtInfoUpd] ...
              = calcApproximationMillard2023VexatMuscleInfo(...
                    timeVector(idxTime),...
                    y(idxTime),...
                    muscleState0,...           
                    activationState0,...                                            
                    pathLengthFcn, ... 
                    pathState0,...
                    ratFibrilModelsUpd(idxTrial).musculotendon,...
                    ratFibrilModelsUpd(idxTrial).sarcomere,...
                    ratFibrilModelsUpd(idxTrial).curves,...
                    modelConfig,...
                    'quasiStaticSkinnedFiberApproximation');

            benchRecord.normFiberForce(idxTime,idxTrial)= ...
              mtInfoUpd.muscleDynamicsInfo.normFiberForce;
            benchRecord.normFiberLength(idxTime,idxTrial)= ...
              mtInfoUpd.muscleDynamicsInfo.normFiberLength;              
          end

          durationApprox=toc(tstart);
          if(flag_compareEvaluationMethods==1)
            benchRecordApprox=benchRecord;
            fprintf('\n%1.3fs\tApprox-by-minimal-state-calculation simulation time\n',...
                    durationApprox);
          end
        end

        if(strcmp(evaluationMethod,'approximationByStateCalculation'))
          tstart=tic;
          nPts = 100;
          timeVector = [0:(1/(nPts-1)):1]'.*(timeRampEnd-timeRampStart)...
                       + timeRampStart;           
          mtInfoUpd=mtInfo;        

          for idxTime =1:1:nPts
            activationStateUpd    = [0;excitationFcn(timeVector(idxTime))];
            pathStateUpd          = pathLengthFcn(timeVector(idxTime));
            muscleStateUpd        = zeros(nStates,1);
          

            muscleStateUpd = mtInfo.state.value;
            muscleStateUpd(1) = pathStateUpd(1)*0.5;
            muscleStateUpd(2) = muscleState0(2)...
                + (pathStateUpd(2)-pathState0(2))*0.5;
            modelConfig.initializeState =0;

            mtInfoUpd= calcMillard2023VexatMuscleInfo(    ...
                        activationStateUpd,...
                        pathStateUpd,...
                        muscleStateUpd,...
                        ratFibrilModelsUpd(idxTrial).musculotendon,...
                        ratFibrilModelsUpd(idxTrial).sarcomere,...
                        ratFibrilModelsUpd(idxTrial).curves,...
                        modelConfig);    

            benchRecord = updateMillard2023VexatMuscleSimulationRecord(...
                         mtInfoUpd,...
                         benchConfig,...
                         benchRecord,...
                         idxTime,...
                         idxTrial,...
                         numberOfSimulations); 

            falN = mtInfoUpd.muscleLengthInfo.fiberActiveForceLengthMultiplier;
            fvN  = mtInfoUpd.fiberVelocityInfo.fiberForceVelocityMultiplier;
            f2N = mtInfoUpd.muscleDynamicsInfo.normTitin2Force;

            benchRecord.normFiberForce(idxTime,idxTrial) = falN*fvN+f2N;
              
%             if(idxTime==1)
%             fprintf('idxTrial\tlp\tvp\tlaHN\tvceN\t\tflN\tfvN\tf2N\n')
%             end
%             fprintf('%i.\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t\t%1.3f\t%1.3f\t%1.3f\n',...
%               idxTime,...
%               mtInfoUpd.muscleLengthInfo.pathLength,...
%               mtInfoUpd.muscleLengthInfo.pathLengthDerivative,...
%               mtInfoUpd.muscleLengthInfo.normFiberSlidingLength,...
%               mtInfoUpd.muscleLengthInfo.normFiberSlidingVelocity,...
%               falN,fvN,f2N);

          end
          durationApprox=toc(tstart);
          if(flag_compareEvaluationMethods==1)
            benchRecordApprox=benchRecord;
            fprintf('\n%1.3fs\tApprox-by-state-calculation simulation time\n',...
                    durationApprox);
          end
        end

        if(strcmp(evaluationMethod,'approximationByInitialization'))  
          tstart=tic;
          nPts = 100;
          timeVector = [0:(1/(nPts-1)):1]'.*(timeRampEnd-timeRampStart)...
                       + timeRampStart;  
          mtInfoUpd=mtInfo;              
          
          for idxTime =1:1:nPts
            activationStateUpd    = [0;excitationFcn(timeVector(idxTime))];
            pathStateUpd          = pathLengthFcn(timeVector(idxTime));
            muscleStateUpd        = zeros(nStates,1);
          
            modelConfig.initializeState =1;


            mtInfoUpd = calcMillard2023VexatMuscleInfo(...
                        activationStateUpd,...
                        pathStateUpd,...
                        muscleStateUpd,...
                        ratFibrilModelsUpd(idxTrial).musculotendon,...
                        ratFibrilModelsUpd(idxTrial).sarcomere,...
                        ratFibrilModelsUpd(idxTrial).curves,...
                        modelConfig);  

            %Set the titin state back to the starting value
            muscleStateUpd = mtInfoUpd.state.value;
            muscleStateUpd(3,1) = muscleState0(3,1);

            modelConfig.initializeState =0;

            mtInfoUpd= calcMillard2023VexatMuscleInfo(    ...
                        activationStateUpd,...
                        pathStateUpd,...
                        muscleStateUpd,...
                        ratFibrilModelsUpd(idxTrial).musculotendon,...
                        ratFibrilModelsUpd(idxTrial).sarcomere,...
                        ratFibrilModelsUpd(idxTrial).curves,...
                        modelConfig);  

            

            benchRecord = updateMillard2023VexatMuscleSimulationRecord(...
                         mtInfoUpd,...
                         benchConfig,...
                         benchRecord,...
                         idxTime,...
                         idxTrial,...
                         numberOfSimulations);

          end
          durationApprox=toc(tstart);
          if(flag_compareEvaluationMethods==1)
            benchRecordApprox=benchRecord;
            fprintf('\n%1.3fs\tApprox-by-initialization simulation time\n',...
                    durationApprox);
          end          
        end 
        if(strcmp(evaluationMethod,'simulateFullModel') ...
            || flag_compareEvaluationMethods==1)
          tstart=tic;
          calcMillard2023VexatMuscleInfoFcn = ...
               @(activationState1,pathState2,mclState3) ...
               calcMillard2023VexatMuscleInfo(    ...
                  activationState1,...
                  pathState2,...
                  mclState3,...
                  ratFibrilModelsUpd(idxTrial).musculotendon,...
                  ratFibrilModelsUpd(idxTrial).sarcomere,...
                  ratFibrilModelsUpd(idxTrial).curves,...
                  modelConfig);           
            %
            % Run the simulation
            %

            flag_appendEnergetics   = 0;
            flag_useOctave          = 0;
                       
            tstart=tic;
            benchRecord = runPrescribedLengthActivationSimulation2025(...
                                       calcMillard2023VexatMuscleInfoFcn,...
                                       [],...
                                       benchConfig,...
                                       benchRecord,...
                                       idxTrial, ...
                                       numberOfSimulations,...
                                       flag_appendEnergetics,...
                                       flag_useOctave);
            durationExact=toc(tstart);
            if(flag_compareEvaluationMethods==1)
              benchRecordExact=benchRecord;
              fprintf('\n%1.3fs\tFull-model simulation time\n',durationExact);
            end
        end
        
               
        if(flag_compareEvaluationMethods==1 ...
            && ~isempty(benchRecordExact) ...
            && ~isempty(benchRecordApprox))

          figComp=figure;
            plot(benchRecordExact.normFiberLength,...
                 benchRecordExact.normFiberForce,...
                 '-k',...
                 'DisplayName','Exact');            
            hold on;
            plot(benchRecordApprox.normFiberLength,...
                 benchRecordApprox.normFiberForce,...
                 '-b',...
                 'DisplayName','Approximate');            
            hold on;
            xlabel('Norm. Length');
            ylabel('Norm. Force');
            here=1;
            
        end

        switch optParams.name
            case 'responseTimeScaling'
                seriesName = expTRSS2017Raw.trials{idxTrial};
                [expLceU,iq] = ...
                  unique(expTRSS2017Raw.(seriesName).lN*lceOptMdl);
                expfNU = expTRSS2017Raw.(seriesName).fNavg(iq);
                expfNStdU = expTRSS2017Raw.(seriesName).fNstd(iq);


        
                errV  = zeros(npts,1);
                nerrV = zeros(npts,1);                
                yV    = zeros(npts,1);
                yFitV = zeros(npts,1);
                yStdV = zeros(npts,1);
                for k=1:1:npts
                    lceN  = benchRecord.normFiberLength(k,idxTrial).*lceOptMdl;
                    fN    = benchRecord.normFiberForce(k,idxTrial);
                    
                    expfN = interp1(expLceU,...
                                    expfNU,...
                                    lceN); 
                    expStdfN   = interp1(expLceU,...
                                         expfNStdU,...
                                         lceN);                    
                    yV(k,1)    = expfN;    
                    yFitV(k,1) = fN;                           
                    yStdV(k,1) = expStdfN;
                    errV(k,1)  = (expfN-fN);
                    nerrV(k,1) = (expfN-fN)/expStdfN;

                end

                n = optParams.exp(1).x;
                for idxTrial=1:1:length(optParams.exp)
                  assert(length(optParams.exp(idxTrial).x)==n);
                end
            
                if(isempty(optErrorValues))
                  optErrorValues.x      = zeros(npts,numberOfSimulations);
                  optErrorValues.y      = zeros(npts,numberOfSimulations);
                  optErrorValues.yFit   = zeros(npts,numberOfSimulations);
                  optErrorValues.yErr   = zeros(npts,numberOfSimulations);    
                  optErrorValues.rmse   = zeros(1,numberOfSimulations);
                  optErrorValues.yStd   = zeros(npts,numberOfSimulations);
                  optErrorValues.yNErr  = zeros(npts,numberOfSimulations);
                  optErrorValues.nrmse  = zeros(1,numberOfSimulations);
                end

                optErrorValues.x(:,idxTrial) = ...
                    benchRecord.normFiberLength(:,idxTrial).*lceOptMdl;
                optErrorValues.y(:,idxTrial) = yV;
                optErrorValues.yFit(:,idxTrial) = yFitV;
                optErrorValues.yErr(:,idxTrial) = errV;
                optErrorValues.yStd(:,idxTrial) = yStdV;
                optErrorValues.yNErr(:,idxTrial)= nerrV;   
                optErrorValues.rmse(1,idxTrial) = ...
                  sqrt(mean(optErrorValues.yErr(:,idxTrial).^2));
                optErrorValues.nrmse(1,idxTrial) = ...
                  sqrt(mean(optErrorValues.yNErr(:,idxTrial).^2));
                optError = optError+optErrorValues.nrmse(1,idxTrial);                
            case 'xeStiffnessDampingScaling'
                %[expLceU,iq] = unique(expTRSS2017.activeLengtheningData(idxTrial).x);
                %expfNU = expTRSS2017.activeLengtheningData(idxTrial).y(iq);

                seriesName = expTRSS2017Raw.trials{idxTrial};
                [expLceU,iq] = ...
                  unique(expTRSS2017Raw.(seriesName).lN*lceOptMdl);
                expfNU = expTRSS2017Raw.(seriesName).fNavg(iq);
                expfNStdU = expTRSS2017Raw.(seriesName).fNstd(iq);
        
                errV = zeros(npts,1);
                nerrV= zeros(npts,1);
                yV = zeros(npts,1);
                yFitV = zeros(npts,1);
                yStdV = zeros(npts,1);
                
                for k=1:1:npts
                    lceN = benchRecord.normFiberLength(k,idxTrial).*lceOptMdl;
                    fN   = benchRecord.normFiberForce(k,idxTrial);
                    
                    expfN = interp1(expLceU,...
                                    expfNU,...
                                    lceN);

                    expStdfN   = interp1(expLceU,...
                                         expfNStdU,...
                                         lceN);
                    yV(k,1)=expfN;
                    yFitV(k,1)=fN;
                    yStdV(k,1)=expStdfN;
                    errV(k,1) = expfN-fN;
                    nerrV(k,1)= (expfN-fN)/expStdfN;
                end

                if(isempty(optErrorValues))
                  optErrorValues.x      = zeros(npts,numberOfSimulations);
                  optErrorValues.y      = zeros(npts,numberOfSimulations);
                  optErrorValues.yFit   = zeros(npts,numberOfSimulations);
                  optErrorValues.yErr   = zeros(npts,numberOfSimulations);    
                  optErrorValues.rmse   = zeros(1,numberOfSimulations);
                  optErrorValues.yStd   = zeros(npts,numberOfSimulations);
                  optErrorValues.yNErr  = zeros(npts,numberOfSimulations);
                  optErrorValues.nrmse  = zeros(1,numberOfSimulations);
                end                

                optErrorValues.x(:,idxTrial) = ...
                    benchRecord.normFiberLength(:,idxTrial).*lceOptMdl;
                optErrorValues.y(:,idxTrial)     = yV; 
                optErrorValues.yFit(:,idxTrial)  = yFitV; 
                optErrorValues.yErr(:,idxTrial)  = errV;
                optErrorValues.yStd(:,idxTrial)  = yStdV;
                optErrorValues.yNErr(:,idxTrial) = nerrV;
                optErrorValues.rmse(1,idxTrial)  = ...
                  sqrt(mean(optErrorValues.yErr(:,idxTrial).^2));
                optErrorValues.nrmse(1,idxTrial) = ...
                  sqrt(mean(optErrorValues.yNErr(:,idxTrial).^2));
                
                optError = optError+optErrorValues.nrmse(1,idxTrial);

            case 'QToF'
               lopt = ratFibrilModelsUpd(idxTrial).sarcomere.optimalSarcomereLength;

               x0N = optParams.exp(idxTrial).x(1,1)/lopt;
               x1N = optParams.exp(idxTrial).x(end)/lopt;
               idx0 = find(benchRecord.normFiberLength(:,idxTrial)<x0N,1,'last');
               idx1 = size(benchRecord.normFiberLength,1);

               %
               % Sample the simulation record at the same lengths as the
               % experimental data
               %
               mdlX = optParams.exp(idxTrial).x;               
               mdlY = zeros(size(optParams.exp(idxTrial).x));

                
               for i=1:1:length(mdlX)
                   mdlY(i,1)=interp1(benchRecord.normFiberLength(idx0:idx1,idxTrial).*lopt,...
                                     benchRecord.normFiberForce(idx0:idx1,idxTrial),...
                                     mdlX(i,1));                   
               end

               n=length(optParams.exp(idxTrial).x);
               for i=1:1:length(optParams.exp)
                 assert(length(optParams.exp(i).x)==n);
               end

               if(isempty(optErrorValues))
                 optErrorValues.x      = zeros(n,numberOfSimulations);
                 optErrorValues.y      = zeros(n,numberOfSimulations);
                 optErrorValues.yFit   = zeros(n,numberOfSimulations);
                 optErrorValues.yErr   = zeros(n,numberOfSimulations);    
                 optErrorValues.rmse   = zeros(1,numberOfSimulations);
                 optErrorValues.yStd   = zeros(n,numberOfSimulations);
                 optErrorValues.yNErr  = zeros(n,numberOfSimulations);
                 optErrorValues.nrmse  = zeros(1,numberOfSimulations);
               end


               optErrorValues.x(:,idxTrial)      = mdlX;
               optErrorValues.y(:,idxTrial)      = optParams.exp(idxTrial).y;
               optErrorValues.yStd(:,idxTrial)   = optParams.exp(idxTrial).yStd;
               optErrorValues.yFit(:,idxTrial)   = mdlY;
               optErrorValues.yErr(:,idxTrial)   = ...
                 (mdlY - optParams.exp(idxTrial).y); 
               optErrorValues.yNErr(:,idxTrial)  = ...
                 optErrorValues.yErr(:,idxTrial) ...
                 ./optParams.exp(idxTrial).yStd; 
                optErrorValues.rmse(1,idxTrial)  = ...
                  sqrt(mean(optErrorValues.yErr(:,idxTrial).^2));
                optErrorValues.nrmse(1,idxTrial) = ...
                  sqrt(mean(optErrorValues.yNErr(:,idxTrial).^2));               
               
               %meanStd  = mean(optParams.exp(idxTrial).yStd);
               %nrmse    = sqrt(mean((mdlY - optParams.exp(idxTrial).y).^2)) ...
               %           /meanStd;
               
               optError = optError+optErrorValues.nrmse(1,idxTrial);                       
            case 'QToK'

               assert(0,['Error: QToK is out of date. The slope of every',...
                         ' experimental trial would need to be calculated',...
                         ' so that the mean and standard deviation of the',...
                         ' slope can be evaluated. Only then can NRMSE be',...
                         ' evaluated']);

               lopt = ratFibrilModelsUpd(idxTrial).sarcomere.optimalSarcomereLength;

               x0N = optParams.exp(idxTrial).x(1,1)/lopt;
               x1N = optParams.exp(idxTrial).x(end)/lopt;
               idx0 = find(benchRecord.normFiberLength(:,idxTrial)<x0N,1,'last');
               idx1 = size(benchRecord.normFiberLength,1);

               %
               % Sample the simulation record at the same lengths as the
               % experimental data
               %
               mdlX = optParams.exp(idxTrial).x;               
               mdlY = zeros(size(optParams.exp(idxTrial).x));

               for i=1:1:length(mdlX)
                   mdlY(i,1)=interp1(benchRecord.normFiberLength(idx0:idx1,idxTrial).*lopt,...
                                     benchRecord.normFiberForce(idx0:idx1,idxTrial),...
                                     mdlX(i,1));                   
               end
                
               %
               % Fit a line to this data
               %
               A = [mdlX,ones(size(mdlX))];
               b= mdlY;
               p = (A'*A)\(A'*b);
               mdlSlope = p(1,1);

               %
               % Calculate the error as the squared difference between the 
               % two slopes
               %    
               n=length(optParams.exp(idxTrial).x);
               for i=1:1:length(optParams.exp)
                 assert(length(optParams.exp(i).x)==n);
               end
               
               if(isempty(optErrorValues))
                 optErrorValues.x      = zeros(n,numberOfSimulations);
                 optErrorValues.y      = zeros(n,numberOfSimulations);
                 optErrorValues.yFit   = zeros(n,numberOfSimulations);
                 optErrorValues.yErr   = zeros(n,numberOfSimulations);    
               end

               optErrorValues.x(:,idxTrial) = [min(mdlX);max(mdlX)];                              
               optErrorValues.y(:,idxTrial)      = optParams.exp(idxTrial).dydx;
               optErrorValues.yFit(:,idxTrial)   = mdlSlope;
               optErrorValues.yErr(:,idxTrial)   = (mdlSlope-optParams.exp(idxTrial).dydx);               

               optError = optError+(mdlSlope-optParams.exp(idxTrial).dydx).^2;
               
            case 'f1HNPreload'

              assert(0,'Error: f1HNPreload needs to have its logic updated');
               lopt = ratFibrilModelsUpd(idxTrial).sarcomere.optimalSarcomereLength;

               x0N = optParams.exp(idxTrial).x(1,1)/lopt;
               x1N = optParams.exp(idxTrial).x(end,1)/lopt;
               idx0 = find(benchRecord.normFiberLength(:,idxTrial)<x0N,1,'last');
               idx1 = size(benchRecord.normFiberLength,1);

               %
               % Sample the simulation record at the same lengths as the
               % experimental data
               %
               mdlX = optParams.exp(idxTrial).x;               
               mdlY = zeros(size(optParams.exp(idxTrial).x));

               for i=1:1:length(mdlX)
                   mdlY(i,1)=interp1(benchRecord.normFiberLength(idx0:idx1,idxTrial).*lopt,...
                                     benchRecord.normFiberForce(idx0:idx1,idxTrial),...
                                     mdlX(i,1));                   
               end

               n=length(optParams.exp(idxTrial).x);
               for i=1:1:length(optParams.exp)
                 assert(length(optParams.exp(i).x)==n);
               end
               
               if(isempty(optErrorValues))
                 optErrorValues.x      = zeros(n,numberOfSimulations);
                 optErrorValues.y      = zeros(n,numberOfSimulations);
                 optErrorValues.yFit   = zeros(n,numberOfSimulations);
                 optErrorValues.yErr   = zeros(n,numberOfSimulations);    
               end

               optErrorValues.x(:,idxTrial) = mdlX;
               optErrorValues.y(:,idxTrial)      = optParams.exp(idxTrial).y;               
               optErrorValues.yFit(:,idxTrial)   = mdlY;                              
               optErrorValues.yErr(:,idxTrial)   = optParams.exp(idxTrial).y-mdlY;               

               optError = optError+sqrt(mean((mdlY - optParams.exp(idxTrial).y).^2));

            case 'l1HNOffset'
              assert(0,'Error: l1HNOffset needs to have its logic updated');              
               lopt = ratFibrilModelsUpd(idxTrial).sarcomere.optimalSarcomereLength;

               x0N = optParams.exp(idxTrial).x(1,1)/lopt;
               x1N = optParams.exp(idxTrial).x(end,1)/lopt;
               idx0 = find(benchRecord.normFiberLength(:,idxTrial)<x0N,1,'last');
               idx1 = size(benchRecord.normFiberLength,1);

               %
               % Sample the simulation record at the same lengths as the
               % experimental data
               %
               mdlX = optParams.exp(idxTrial).x;               
               mdlY = zeros(size(optParams.exp(idxTrial).x));

               for i=1:1:length(mdlX)
                   mdlY(i,1)=interp1(benchRecord.normFiberLength(idx0:idx1,idxTrial).*lopt,...
                                     benchRecord.normFiberForce(idx0:idx1,idxTrial),...
                                     mdlX(i,1));                   
               end

               n=length(optParams.exp(idxTrial).x);
               for i=1:1:length(optParams.exp)
                 assert(length(optParams.exp(i).x)==n);
               end
               
               if(isempty(optErrorValues))
                 optErrorValues.x      = zeros(n,numberOfSimulations);
                 optErrorValues.y      = zeros(n,numberOfSimulations);
                 optErrorValues.yFit   = zeros(n,numberOfSimulations);
                 optErrorValues.yErr   = zeros(n,numberOfSimulations);    
               end

               optErrorValues.x(:,idxTrial) = mdlX;
               optErrorValues.y(:,idxTrial)      = optParams.exp(idxTrial).y;               
               optErrorValues.yFit(:,idxTrial)   = mdlY;                              
               optErrorValues.yErr(:,idxTrial)   = optParams.exp(idxTrial).y-mdlY;               

               optError = optError+sqrt(mean((mdlY - optParams.exp(idxTrial).y).^2));

            case 'FpeQToF'
               lopt = ratFibrilModelsUpd(idxTrial).sarcomere.optimalSarcomereLength;

               x0N = optParams.exp(idxTrial).x(1,1)/lopt;
               x1N = optParams.exp(idxTrial).x(end)/lopt;
               idx0 = find(benchRecord.normFiberLength(:,idxTrial)<x0N,1,'last');
               idx1 = size(benchRecord.normFiberLength,1);

               %
               % Sample the simulation record at the same lengths as the
               % experimental data
               %
               mdlX = optParams.exp(idxTrial).x;               
               mdlY = zeros(size(optParams.exp(idxTrial).x));

                
               for i=1:1:length(mdlX)
                   mdlY(i,1)=...
                     interp1(benchRecord.normFiberLength(idx0:idx1,idxTrial).*lopt,...
                             benchRecord.normFiberForce(idx0:idx1,idxTrial),...
                             mdlX(i,1));                   
               end

               n=length(optParams.exp(idxTrial).x);
               for i=1:1:length(optParams.exp)
                 assert(length(optParams.exp(i).x)==n);
               end
               
               if(isempty(optErrorValues))
                 optErrorValues.x      = zeros(n,numberOfSimulations);
                 optErrorValues.y      = zeros(n,numberOfSimulations);
                 optErrorValues.yFit   = zeros(n,numberOfSimulations);
                 optErrorValues.yErr   = zeros(n,numberOfSimulations);    
                 optErrorValues.rmse   = zeros(1,numberOfSimulations);
                 optErrorValues.yStd   = zeros(n,numberOfSimulations);
                 optErrorValues.yNErr  = zeros(n,numberOfSimulations);
                 optErrorValues.nrmse  = zeros(1,numberOfSimulations);
               end

               optErrorValues.x(:,idxTrial)      = mdlX;
               optErrorValues.y(:,idxTrial)      = optParams.exp(idxTrial).y;
               optErrorValues.yStd(:,idxTrial)   = optParams.exp(idxTrial).yStd;
               optErrorValues.yFit(:,idxTrial)   = mdlY;
               optErrorValues.yErr(:,idxTrial)   = ...
                 (mdlY - optParams.exp(idxTrial).y); 
               optErrorValues.yNErr(:,idxTrial)  = ...
                 optErrorValues.yErr(:,idxTrial) ...
                 ./optParams.exp(idxTrial).yStd; 
                optErrorValues.rmse(1,idxTrial)  = ...
                  sqrt(mean(optErrorValues.yErr(:,idxTrial).^2));
                optErrorValues.nrmse(1,idxTrial) = ...
                  sqrt(mean(optErrorValues.yNErr(:,idxTrial).^2));               
               
               
               optError = optError+optErrorValues.nrmse(1,idxTrial); 

            case 'simulate'
                optError=nan;
            otherwise
                assert(0,'Error: invalid optParams.name');
        end

        

        if(simConfig.flag_debugFitting==1)
    
            figure(figDebugFitting);
            subplot('Position',reshape(subPlotPanel(2,1,:),1,4));
    
    
            switch optParams.name
                case 'responseTimeScaling'
                    %txtName = expTRSS2017.activeLengtheningData(idxTrial).seriesName;
                    %i0=strfind(txtName,'Exp.');
                    %txtName(1,i0:4)='Sim.';
                    seriesName = expTRSS2017Raw.trials{idxTrial};
                    txtName = ['Sim. ', seriesName];

                    plot(benchRecord.normFiberLength(:,idxTrial).*lceOptMdl,...
                         benchRecord.normCrossBridgeForceAlongTendon(:,idxTrial),...
                         '-','Color',lineColorsSimTRSS2017(idxTrial,:),...
                         'DisplayName',...
                         txtName);
                    hold on;
                    
                case 'xeStiffnessDampingScaling'
                    %txtName = expTRSS2017.activeLengtheningData(idxTrial).seriesName;
                    %i0=strfind(txtName,'Exp.');
                    %txtName(1,i0:4)='Sim.';
                    seriesName = expTRSS2017Raw.trials{idxTrial};
                    txtName = ['Sim. ', seriesName];
            
                    plot(benchRecord.normFiberLength(:,idxTrial).*lceOptMdl,...
                         benchRecord.normCrossBridgeForceAlongTendon(:,idxTrial),...
                         '-','Color',lineColorsSimTRSS2017(idxTrial,:),...
                         'DisplayName',...
                         txtName);
                    hold on;
                    
                case 'QToF'
                    %txtName = expTRSS2017.activeLengtheningData(idxTrial).seriesName;
                    %i0=strfind(txtName,'Exp.');
                    %txtName(1,i0:4)='Sim.';
                    seriesName = expTRSS2017Raw.trials{idxTrial};
                    txtName = ['Sim. ', seriesName];            
                    
                    plot(benchRecord.normFiberLength(:,idxTrial).*lceOptMdl,...
                         benchRecord.normFiberForce(:,idxTrial),...
                         '-','Color',lineColorsSimTRSS2017(idxTrial,:),...
                         'DisplayName',...
                         txtName);
                    hold on;

                case 'QToK'
                    %txtName = expTRSS2017.activeLengtheningData(idxTrial).seriesName;
                    %i0=strfind(txtName,'Exp.');
                    %txtName(1,i0:4)='Sim.';
                    seriesName = expTRSS2017Raw.trials{idxTrial};
                    txtName = ['Sim. ', seriesName];                    
            
                    plot(benchRecord.normFiberLength(:,idxTrial).*lceOptMdl,...
                         benchRecord.normFiberForce(:,idxTrial),...
                         '-','Color',lineColorsSimTRSS2017(idxTrial,:),...
                         'DisplayName',...
                         txtName);
                    hold on;
    
                    
                case 'f1HNPreload'
                    %txtName = expTRSS2017.activeLengtheningData(idxTrial).seriesName;
                    %i0=strfind(txtName,'Exp.');
                    %txtName(1,i0:4)='Sim.';
                    seriesName = expTRSS2017Raw.trials{idxTrial};
                    txtName = ['Sim. ', seriesName];            
                    
                    plot(benchRecord.normFiberLength(:,idxTrial).*lceOptMdl,...
                         benchRecord.normFiberForce(:,idxTrial),...
                         '-','Color',lineColorsSimTRSS2017(idxTrial,:),...
                         'DisplayName',...
                         txtName);
                    hold on;
    
                case 'l1HNOffset'
                    %txtName = expTRSS2017.activeLengtheningData(idxTrial).seriesName;
                    %i0=strfind(txtName,'Exp.');
                    %txtName(1,i0:4)='Sim.';
                    seriesName = expTRSS2017Raw.trials{idxTrial};
                    txtName = ['Sim. ', seriesName];

                    plot(benchRecord.normFiberLength(:,idxTrial).*lceOptMdl,...
                         benchRecord.normFiberForce(:,idxTrial),...
                         '-','Color',lineColorsSimTRSS2017(idxTrial,:),...
                         'DisplayName',...
                         txtName);
                    hold on;

                case 'FpeQToF'
                    %txtName = expTRSS2017.activeLengtheningData(idxTrial).seriesName;
                    %i0=strfind(txtName,'Exp.');
                    %txtName(1,i0:4)='Sim.';
                    seriesName = expTRSS2017Raw.trials{idxTrial};
                    txtName = ['Sim. ', seriesName];            
                    
                    plot(benchRecord.normFiberLength(:,idxTrial).*lceOptMdl,...
                         benchRecord.normFiberForce(:,idxTrial),...
                         '-','Color',lineColorsSimTRSS2017(idxTrial,:),...
                         'DisplayName',...
                         txtName);
                    hold on;

                otherwise
                    assert(0,'Error: invalid optParams.name');
            end

    
        end
    end   
