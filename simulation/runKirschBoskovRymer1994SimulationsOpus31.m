function [success] = runKirschBoskovRymer1994SimulationsOpus31(...
                          inputFunctions,...                          
                          normFiberLength,...
                          nominalForce,...
                          amplitudeMM,...
                          bandwidthHz,... 
                          numberOfSimulations,...
                          flag_useElasticTendon,...
                          muscleArchitecture,...
                          sarcomereProperties,...
                          normMuscleCurves,...
                          outputFileEndingOpus31, ...
                          outputFolder,...
                          flag_useOctave)
success = 0;

disp('Running Opus 31  Kirsch, Boskov, and Rymer 1994  Simulations');

  idx=1;
  benchRecord = [];
  for idxNormFiberLength = 1:1:length(normFiberLength)
    for idxActivation = 1:1:length(nominalForce)  
      flag_setWarning=0;
      for i=1:1:length(amplitudeMM)        
        for j=1:1:length(bandwidthHz)
          
          
          idxWave = 0;
          for z=1:1:length(inputFunctions.amplitudeMM)
            if( abs(inputFunctions.amplitudeMM(z)-amplitudeMM(i))<sqrt(eps) ...
              && abs(inputFunctions.bandwidthHz(z)-bandwidthHz(j))<sqrt(eps))
              idxWave = z;
            end          
          end
      

            x0     = inputFunctions.x(1,idxWave);
            xdot0  = inputFunctions.xdot(1,idxWave);

            %Evaluate the fiber kinematics
            alphaOpt = muscleArchitecture.pennationAngle;
            lceOpt   = muscleArchitecture.optimalFiberLength;

            lce0 = normFiberLength(idxNormFiberLength,1)*lceOpt + x0;
            lceN0= lce0/lceOpt;

            dlce0= 0 + xdot0;
            dlceN0= dlce0 / muscleArchitecture.maximumNormalizedFiberVelocity;

            fibKin   = calcFixedWidthPennatedFiberKinematicsAlongTendon(...
                                           lce0,dlce0,lceOpt,alphaOpt);
            alpha = fibKin.pennationAngle;
            lceAT = fibKin.fiberLengthAlongTendon;

           
            %Evaluate the fiber force
            fiso   = muscleArchitecture.fiso;
            fpeN   = calcBezierYFcnXDerivative(lceN0,  normMuscleCurves.fiberForceLengthCurve  ,0);

            %Evaluate the baseline length: tendon is developing fpe
            %Evaluate the tendon length
            ft0 = fpeN*cos(alpha);
            tendonForceLengthCurveInverse = ...
                createInverseCurve(normMuscleCurves.tendonForceLengthCurve);
            ltN    = calcBezierYFcnXDerivative(ft0, tendonForceLengthCurveInverse,0);
            if(flag_useElasticTendon==0)
              ltN = 1;
            end            
            ltslk  = muscleArchitecture.tendonSlackLength;
            
            lp0  = lceAT+ltN*ltslk;

            %Now evaluate the length of the fiber when the muscle is
            %activated
            ft1    = nominalForce(1,idxActivation)/fiso;
            ltN    = calcBezierYFcnXDerivative(ft1, tendonForceLengthCurveInverse,0);
            if(flag_useElasticTendon==0)
              ltN = 1;
            end
            lceAT = lp0 - ltN*ltslk;

            fibKin   = calcFixedWidthPennatedFiberKinematics(...
                                           lceAT,0,lceOpt,alphaOpt);
            alpha = fibKin.pennationAngle;
            lceN1 = fibKin.fiberLength/lceOpt;
            dlceN1= fibKin.fiberVelocity/muscleArchitecture.maximumNormalizedFiberVelocity;
           
            %%
            %Evaluate the fiber length when its developing the desired
            %force
            %%
            falN          = calcBezierYFcnXDerivative(lceN1,  normMuscleCurves.activeForceLengthCurve ,0);
            DfalN_DlceN   = calcBezierYFcnXDerivative(lceN1,  normMuscleCurves.activeForceLengthCurve ,1);

            fvN            = calcBezierYFcnXDerivative(dlceN1, normMuscleCurves.fiberForceVelocityCurve,0);    
            DfvN_DdlceN    = calcBezierYFcnXDerivative(dlceN1, normMuscleCurves.fiberForceVelocityCurve,1);

            fpeN        = calcBezierYFcnXDerivative(lceN1,  normMuscleCurves.fiberForceLengthCurve  ,0);
            DfpeN_DlceN = calcBezierYFcnXDerivative(lceN1,  normMuscleCurves.fiberForceLengthCurve  ,1);
            
            fcnN   = ft1/cos(alpha);
            activation(1,idxActivation) = max( (fcnN-fpeN)/(falN*fvN), 0);
            
            if( fcnN < fpeN && flag_setWarning==0)              
              fprintf('  Warning: \tdesired force (%1.3f) > passive force (%1.3f) at length (%1.3f)\n',...
                fcnN, fpeN,lceN1);
              flag_setWarning=1;
            end       

          %%
          % Setup the function handles
          %%

          excitationFcn = @(argT)calcStepFunction(argT,0,max(inputFunctions.time)*2,activation(1,idxActivation));
          activationFcn = @(argU,argA)calcFirstOrderActivationDerivative(...
                            argU,argA, sarcomereProperties.activationTimeConstant,...
                                      sarcomereProperties.deactivationTimeConstant,0);
          pathLengthFcn = @(argT)calcOffsetWaveform(argT,lp0, inputFunctions.time, ...
                                          inputFunctions.x(:,idxWave), ...
                                          inputFunctions.xdot(:,idxWave));



          %%
          % Set up the simulation structs
          %%
            benchConfig.npts                  = inputFunctions.totalpoints;
            benchConfig.relTol                = 1e-6;
            benchConfig.absTol                = 1e-6;
            benchConfig.minActivation         = 0;
            benchConfig.color0                = [0,0,1].*0.5;
            benchConfig.color1                = [0,0,1];

            nStates = 0;
            labelStates = {''};
            if(flag_useElasticTendon==1)
              nStates = 5;
              labelStates= {'$$\ell_{CE}$$','$$\dot{\ell}_{a}$$',...
                 '$$\ell_{a}$$','$$\ell_1$$','$$\lambda$$'};%,'$$f_{e}^1$$', '$$f_{e}^2$$'};
            else
              nStates = 4;
              labelStates= {'$$\dot{\ell}_{a}$$',...
                 '$$\ell_{a}$$','$$\ell_1$$','$$\lambda$$'};%,'$$f_{e}^1$$', '$$f_{e}^2$$'};            
            end

            benchConfig.numberOfMuscleStates  = nStates;
            benchConfig.stateLabels  = labelStates;
            benchConfig.name                  = '';
            benchConfig.initialState          = [];
            benchConfig.initialActivation     = activation(1,idxActivation);
            benchConfig.pathFcn               = pathLengthFcn;
            benchConfig.excitationFcn         = excitationFcn;
            benchConfig.activationFcn         = activationFcn; 
            benchConfig.tspan                 = inputFunctions.time;  

            benchConfig.useFiberDamping  = 1;
            benchConfig.useElasticTendon = flag_useElasticTendon;
            benchConfig.damping          = 0.1;
            benchConfig.iterMax          = 100;
            benchConfig.tol              = 1e-6;

            loopTolerance = min(benchConfig.relTol,benchConfig.absTol)/100;
            
            modelConfig = struct( ...
              'iterMax'                 , 100             , ...
              'tol'                     , loopTolerance   , ... 
              'tolInit'                 , sqrt(eps)       , ...
              'minActivation'           , 0.0             , ...
              'useElasticTendon'        , flag_useElasticTendon , ...
              'initializeState'         , 0                );          

              modelConfig.initializeState =1;
              activationState0 = [0;benchConfig.excitationFcn(0)];
              pathState0 = benchConfig.pathFcn(0);
              muscleState0 = zeros(nStates,1);
              mtInfo =calcMillard2019MuscleInfoOpus31(activationState0,...
                                                      pathState0,...
                                                      muscleState0,...
                                                      muscleArchitecture,...
                                                      sarcomereProperties,...
                                                      normMuscleCurves,...
                                                      modelConfig);
              muscleState0 =mtInfo.state.value;
              modelConfig.initializeState =0;           

              benchConfig.numberOfMuscleStates = length(muscleState0);
              benchConfig.initialState         = muscleState0;

              benchConfig.minimumActivation    = 0;
              benchConfig.name = 'Opus31';
              benchConfig.eventFcn = [];            


              calcMillard2019MuscleInfoOpus31Fcn = ...
                   @(activationState1,pathState2,mclState3) ...
                   calcMillard2019MuscleInfoOpus31(   activationState1,...
                                                      pathState2,...
                                                      mclState3,...
                                                      muscleArchitecture,...
                                                      sarcomereProperties,...
                                                      normMuscleCurves,...
                                                      modelConfig);

              flag_appendEnergetics =0;

              benchRecord = runPrescribedLengthActivationSimulation(...
                                           calcMillard2019MuscleInfoOpus31Fcn,...
                                           [],...
                                           benchConfig,...
                                           benchRecord,...
                                           idx, numberOfSimulations,...
                                           flag_appendEnergetics,...
                                           flag_useOctave);

          nameModification = '';
          if(benchConfig.useElasticTendon == 1)
            nameModification = 'ElasticTendon';
          else
            nameModification = 'RigidTendon';          
          end


          save([outputFolder,'benchRecordOpus31_',nameModification,outputFileEndingOpus31,'.mat'],...
            'benchRecord',...
            'muscleArchitecture',...
            'sarcomereProperties',...
            'normMuscleCurves',...
            'inputFunctions',...                          
            'normFiberLength',...
            'nominalForce',...
            'amplitudeMM',...
            'bandwidthHz',... 
            'numberOfSimulations',...
            'flag_useElasticTendon');
          fprintf('%i / %i\n', idx, numberOfSimulations);
          idx=idx+1;
          pause(0.1);
        end
      end
    end
  end

success = 1;
                        
                        