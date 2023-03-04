function [success] = runLeonardJoumaaHerzog2010SimulationsDampedEquilibrium(...
                          nominalNormalizedFiberLength,...
                          passiveForceKeyPoints,...
                          timeSpan,...
                          lengthRampKeyPoints,...
                          stimulationKeyTimes,...
                          flag_useElasticTendon,...
                          flag_useFiberDamping,...
                          fiberDampingCoefficient,...
                          musculotendonProperties,...
                          sarcomereProperties,...
                          normMuscleCurves,...
                          outputFileEndingHill, ...
                          outputFolder,...
                          flag_simulateActiveStretch,...
                          flag_simulatePassiveStretch,...
                          flag_useOctave)
                        
success = 0;

disp('Running Damped-Equilibrium Leonard, Joumaa, and Herzog 2010 Simulations');

%figBasic = figure;

%%
% Solve for the initial state of the fiber
%%
idx=1;
benchRecord = [];
%for idxNormFiberLength = 1:1:length(normFiberLength)
%  for idxActivation = 1:1:length(nominalForce)
%    flag_setWarning=0;
%    for i=1:1:length(amplitudeMM)        
%      for j=1:1:length(bandwidthHz)

        assert(size(nominalNormalizedFiberLength,1)==1 && ...
               size(nominalNormalizedFiberLength,2)==1);


        fpe0 = calcBezierYFcnXDerivative(...
                passiveForceKeyPoints(1,1),...
                normMuscleCurves.fiberForceLengthCurve,...
                0);
        
        fpe1 = calcBezierYFcnXDerivative(...
                passiveForceKeyPoints(2,1),...
                normMuscleCurves.fiberForceLengthCurve,...
                0);
        
        scaleFPE = passiveForceKeyPoints(2,2)/fpe1;

        normMuscleCurves.fiberForceLengthCurve.ypts = ...
          normMuscleCurves.fiberForceLengthCurve.ypts.*scaleFPE; 
        normMuscleCurves.fiberForceLengthCurve.yEnd = ...
          normMuscleCurves.fiberForceLengthCurve.yEnd.*scaleFPE; 
        normMuscleCurves.fiberForceLengthCurve.dydxEnd = ...
          normMuscleCurves.fiberForceLengthCurve.dydxEnd.*scaleFPE; 
        normMuscleCurves.fiberForceLengthCurve.d2ydx2End = ...
          normMuscleCurves.fiberForceLengthCurve.d2ydx2End.*scaleFPE; 


        %%
        % Evaluate the model when it is passive and at its nominal length
        %%
                    



        alphaOpt  = musculotendonProperties.pennationAngle;
        lceOpt    = musculotendonProperties.optimalFiberLength;
        lce       = nominalNormalizedFiberLength*lceOpt;

        fibKin    = calcFixedWidthPennatedFiberKinematicsAlongTendon(...
                                       lce,0,lceOpt,alphaOpt);
        alphaNom     = fibKin.pennationAngle;            
        lceATNom     = fibKin.fiberLengthAlongTendon;
        
        rampStartTime   = lengthRampKeyPoints(1,1);
        rampEndTime     = lengthRampKeyPoints(2,1);            
        rampStartLength = lengthRampKeyPoints(1,2)*lceOpt;
        rampEndLength   = lengthRampKeyPoints(2,2)*lceOpt;
        

        %Here we assume all most all of the length change is borne by the
        %fiber. For Herzog & Leonard 2002 this is a good assumption because
        %the tendon is short and stiff and the experiments are done 
        %at fiber lengths that have low passive forces
        lceAT0 = lceATNom + rampStartLength; 

        fiberKinematics = calcFixedWidthPennatedFiberKinematics(...
                              lceAT0,0,lceOpt,alphaOpt);

        lce0   = fiberKinematics.fiberLength;
        dlce0  = fiberKinematics.fiberVelocity;
        alpha0 = fiberKinematics.pennationAngle;

        lceN0  = lce0/lceOpt;

        %Evaluate the passive fiber force
        fiso   = musculotendonProperties.fiso;
        fpeN0  = calcBezierYFcnXDerivative(lceN0, ...
                    normMuscleCurves.fiberForceLengthCurve  ,0)*scaleFPE;

        %Evaluate the passive path length
        ft0 = fpeN0*cos(alpha0);
        tendonForceLengthCurveInverse = ...
            createInverseCurve(normMuscleCurves.tendonForceLengthCurve);
        ltN0    = calcBezierYFcnXDerivative(ft0, tendonForceLengthCurveInverse,0);
        if(flag_useElasticTendon==0)
          ltN0 = 1;
        end            
        ltslk  = musculotendonProperties.tendonSlackLength;

        pathStartLength  = lceAT0+ltN0*ltslk;

        %%
        % Now evaluate the length of the fiber when the muscle is at
        % the beginning of the ramp (stretched by rampStart) and is activated
        %%
        falN0 = calcBezierYFcnXDerivative(lceN0, normMuscleCurves.activeForceLengthCurve  ,0);
        fvN0  = 1;        
        ft1    = (falN0*fvN0 + fpeN0)*cos(alpha0);
        ltN1   = calcBezierYFcnXDerivative(ft1,...
                   tendonForceLengthCurveInverse,0);
        if(flag_useElasticTendon==0)
          ltN1 = 1;
        end
        lceAT1  = pathStartLength - ltN1*ltslk;
        dlceAT1 = 0;

        fibKin = calcFixedWidthPennatedFiberKinematics(...
                                       lceAT1,dlceAT1,lceOpt,alphaOpt);
        alpha1   = fibKin.pennationAngle;
        lceN1   = fibKin.fiberLength/lceOpt;
        dlceN1  = fibKin.fiberVelocity...
                    /musculotendonProperties.maximumNormalizedFiberVelocity;

%           %%
%           %Evaluate the fiber length when its developing the desired
%           %force
%           %%
%           falN          = calcBezierYFcnXDerivative(lceN1,  normMuscleCurves.activeForceLengthCurve ,0);
%           DfalN_DlceN   = calcBezierYFcnXDerivative(lceN1,  normMuscleCurves.activeForceLengthCurve ,1);
% 
%           fvN            = calcBezierYFcnXDerivative(dlceN1, normMuscleCurves.fiberForceVelocityCurve,0);    
%           DfvN_DdlceN    = calcBezierYFcnXDerivative(dlceN1, normMuscleCurves.fiberForceVelocityCurve,1);
% 
%           fpeN        = calcBezierYFcnXDerivative(lceN1,  normMuscleCurves.fiberForceLengthCurve  ,0);
%           DfpeN_DlceN = calcBezierYFcnXDerivative(lceN1,  normMuscleCurves.fiberForceLengthCurve  ,1);
% 
%           fcnN   = ft1/cos(alpha);
%           activation(1,idxActivation) = max(0,(fcnN-fpeN)/(falN*fvN));
% 
%           if( fcnN < fpeN && flag_setWarning==0)              
%             fprintf('  Warning: \tdesired force (%1.3f) > passive force (%1.3f) at length (%1.3f)\n',...
%               fcnN, fpeN,lceN1);
%             flag_setWarning=1;
%           end


          %Now we have the path length offset


          %%
          % Setup the function handles
          %%

          activation=1;
          excitationSquareFcn = @(argT)calcStepFunction(argT,...
                                  stimulationKeyTimes(1,1),...
                                  stimulationKeyTimes(2,1),...
                                  activation);
                            
          excitationZeroFcn = @(argT)calcStepFunction(argT,...
                                  stimulationKeyTimes(1,1),...
                                  stimulationKeyTimes(2,1),...
                                  0);
                            
          activationFcn = @(argU,argA)calcFirstOrderActivationDerivative(...
                            argU,argA, sarcomereProperties.activationTimeConstant,...
                                       sarcomereProperties.deactivationTimeConstant,0);

          rampSlope = (rampEndLength-rampStartLength) ...
                     /(rampEndTime-rampStartTime);

          pathLengthRampFcn = @(argT)calcRampStateSharp(...
                                   argT,rampStartTime+0.1,rampEndTime+0.1,...
                                   pathStartLength,rampSlope);
                                 
          pathLengthStaticFcn = @(argT)calcRampStateSharp(...
                                   argT,rampStartTime,rampEndTime,...
                                   pathStartLength,0);

          %%
          % Set up the simulation structs
          %%
            minActivation=0;
            if(flag_useFiberDamping==0)
              minActivation=0.001;
            end

            dampedFiberElasticTendonConfig.npts          = round(timeSpan(1,2)-timeSpan(1,1))*100;
            dampedFiberElasticTendonConfig.relTol        = 1e-6;
            dampedFiberElasticTendonConfig.absTol        = 1e-6;
            dampedFiberElasticTendonConfig.minActivation = minActivation;
            dampedFiberElasticTendonConfig.color0        = [0,0,1].*0.5;
            dampedFiberElasticTendonConfig.color1        = [0,0,1];

            nStates=0;
            labelStates={''};
            if(flag_useElasticTendon==1)
              nStates=1;
              labelStates = {'$$\ell_{CE}$$'};
            end

            dampedFiberElasticTendonConfig.numberOfMuscleStates  = nStates;
            dampedFiberElasticTendonConfig.stateLabels = labelStates;
            dampedFiberElasticTendonConfig.name                  = '';
            dampedFiberElasticTendonConfig.initialState          = [];
            dampedFiberElasticTendonConfig.initialActivation     = 0;
            dampedFiberElasticTendonConfig.pathFcn               = [];
            dampedFiberElasticTendonConfig.excitationFcn         = [];
            dampedFiberElasticTendonConfig.activationFcn         = activationFcn; 
            dampedFiberElasticTendonConfig.tspan                 = timeSpan;  

            dampedFiberElasticTendonConfig.useFiberDamping  = flag_useFiberDamping;
            dampedFiberElasticTendonConfig.useElasticTendon = flag_useElasticTendon;
            dampedFiberElasticTendonConfig.damping          = fiberDampingCoefficient;
            dampedFiberElasticTendonConfig.iterMax          = 100;
            dampedFiberElasticTendonConfig.tol              = 1e-6;

            calcDampedFiberElasticTendonMuscleInfoFcn =...
                @(actState1,pathState2,mclState3)...
                calcMillard2012DampedEquilibriumMuscleInfo(  ...
                                            actState1,...
                                            pathState2, ... 
                                            mclState3,...                                                                           
                                            musculotendonProperties,...
                                            normMuscleCurves,...
                                            dampedFiberElasticTendonConfig);   

            calcDampedFiberElasticTendonInitialMuscleStateFcn = ...
                @(actState1,pathState2,calcMuscleInfo3, initConfig4) ...
                    calcInitialMuscleState(actState1,...
                                           pathState2,...
                                           musculotendonProperties,...
                                           calcMuscleInfo3,...
                                           initConfig4);


            dampedFiberElasticTendonConfig.numberOfMuscleStates = nStates;
            dampedFiberElasticTendonConfig.name = 'DFE';
            dampedFiberElasticTendonConfig.eventFcn = [];

            flag_appendEnergetics =0;
            numberOfSimulations = flag_simulateActiveStretch...
                                + flag_simulatePassiveStretch;            

            benchRecord = [];
            idx = 1;
            if(flag_simulateActiveStretch == 1)
              dampedFiberElasticTendonConfig.pathFcn       = pathLengthRampFcn;
              dampedFiberElasticTendonConfig.excitationFcn = excitationSquareFcn;
              benchRecord = runPrescribedLengthActivationSimulation(...
                             calcDampedFiberElasticTendonMuscleInfoFcn,...
                             calcDampedFiberElasticTendonInitialMuscleStateFcn,...
                             dampedFiberElasticTendonConfig,...
                             benchRecord,...
                             idx, numberOfSimulations,...
                             flag_appendEnergetics,...
                             flag_useOctave);   
              fprintf('%i / %i\n', idx, numberOfSimulations); 
              idx=idx+1;
            end
           
            if(flag_simulatePassiveStretch == 1)
              dampedFiberElasticTendonConfig.pathFcn               = pathLengthRampFcn;
              dampedFiberElasticTendonConfig.excitationFcn         = excitationZeroFcn;
              benchRecord = runPrescribedLengthActivationSimulation(...
                             calcDampedFiberElasticTendonMuscleInfoFcn,...
                             calcDampedFiberElasticTendonInitialMuscleStateFcn,...
                             dampedFiberElasticTendonConfig,...
                             benchRecord,...
                             idx, numberOfSimulations,...
                             flag_appendEnergetics,...
                             flag_useOctave); 
              fprintf('%i / %i\n', idx, numberOfSimulations); 
              idx=idx+1;
            end
                              
            nameModification = '';
            if(dampedFiberElasticTendonConfig.useElasticTendon == 1)
              nameModification = 'ElasticTendon';
            else
              nameModification = 'RigidTendon';          
            end

            modelName = 'benchRecordHill_';

          save([outputFolder,modelName,outputFileEndingHill,'.mat'],...
            'benchRecord',...
            'nominalNormalizedFiberLength',...
            'timeSpan',...
            'lengthRampKeyPoints',...
            'stimulationKeyTimes',...
            'flag_useElasticTendon',...
            'musculotendonProperties',...
            'sarcomereProperties',...
            'normMuscleCurves');
         
          
          
%       end
%     end
%   end
% end


success = 1;