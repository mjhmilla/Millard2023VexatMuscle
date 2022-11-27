function [success] = runKirschBoskovRymer1994SimulationsDampedEquilibrium( ...
                          inputFunctions,...                          
                          normFiberLength,...
                          nominalForce,...
                          amplitudeMM,...
                          bandwidthHz,... 
                          numberOfSimulations,...
                          flag_useElasticTendon,...
                          flag_useFiberDamping,...
                          scalePassiveForceLengthCurve,...
                          muscleArchitecture,...
                          sarcomereProperties,...
                          normMuscleCurves,...
                          outputFileEndingHill, ...
                          outputFolder,...
                          flag_useOctave)
                        
success = 0;

fiberDamping = 0.1;  
%figBasic = figure;

muscleArchitecture.scaleHillFpe = scalePassiveForceLengthCurve;

%%
% Solve for the initial state of the fiber
%%
idx=1;
benchRecord = [];
disp('Running Damped-Equilibrium Kirsch, Boskov, and Rymer 1994 Simulations');
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
          activation(1,idxActivation) = max(0,(fcnN-fpeN)/(falN*fvN));

          if( fcnN < fpeN && flag_setWarning==0)              
            fprintf('  Warning: \tdesired force (%1.3f) > passive force (%1.3f) at length (%1.3f)\n',...
              fcnN, fpeN,lceN1);
            flag_setWarning=1;
          end


          %Now we have the path length offset


          %%
          % Setup the function handles
          %%

          excitationFcn = @(argT)calcStepFunction(argT,-inf,...
                                      inf,...
                                      activation(1,idxActivation));
          activationFcn = @(argU,argA)calcFirstOrderActivationDerivative(...
                            argU,argA, sarcomereProperties.activationTimeConstant,...
                                      sarcomereProperties.deactivationTimeConstant,0);
          pathLengthFcn = @(argT)calcOffsetWaveform(argT,lp0, inputFunctions.time, ...
                                          inputFunctions.x(:,idxWave), ...
                                          inputFunctions.xdot(:,idxWave));



          %%
          % Set up the simulation structs
          %%
            minActivation=0;
            if(flag_useFiberDamping==0)
              minActivation=0.001;
            end

            dampedFiberElasticTendonConfig.npts                  = inputFunctions.totalpoints;
            dampedFiberElasticTendonConfig.relTol                = 1e-6;
            dampedFiberElasticTendonConfig.absTol                = 1e-6;
            dampedFiberElasticTendonConfig.minActivation         = minActivation;
            dampedFiberElasticTendonConfig.color0                = [0,0,1].*0.5;
            dampedFiberElasticTendonConfig.color1                = [0,0,1];

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
            dampedFiberElasticTendonConfig.initialActivation     = activation(1,idxActivation);
            dampedFiberElasticTendonConfig.pathFcn               = pathLengthFcn;
            dampedFiberElasticTendonConfig.excitationFcn         = excitationFcn;
            dampedFiberElasticTendonConfig.activationFcn         = activationFcn; 
            dampedFiberElasticTendonConfig.tspan                 = inputFunctions.time;  

            dampedFiberElasticTendonConfig.useFiberDamping  = flag_useFiberDamping;
            dampedFiberElasticTendonConfig.useElasticTendon = flag_useElasticTendon;
            dampedFiberElasticTendonConfig.damping          = 0.1;
            dampedFiberElasticTendonConfig.iterMax          = 100;
            dampedFiberElasticTendonConfig.tol              = 1e-6;

            calcDampedFiberElasticTendonMuscleInfoFcn =...
                @(actState1,pathState2,mclState3)...
                calcMillard2012DampedEquilibriumMuscleInfo(  ...
                                            actState1,...
                                            pathState2, ... 
                                            mclState3,...                                                                           
                                            muscleArchitecture,...
                                            normMuscleCurves,...
                                            dampedFiberElasticTendonConfig);   

            calcDampedFiberElasticTendonInitialMuscleStateFcn = ...
                @(actState1,pathState2,calcMuscleInfo3, initConfig4) ...
                    calcInitialMuscleState(actState1,...
                                           pathState2,...
                                           muscleArchitecture,...
                                           calcMuscleInfo3,...
                                           initConfig4);


            dampedFiberElasticTendonConfig.numberOfMuscleStates = nStates;
            dampedFiberElasticTendonConfig.name = 'DFE';
            dampedFiberElasticTendonConfig.eventFcn = [];

            flag_appendEnergetics =0;

            benchRecord = ...
                runPrescribedLengthActivationSimulation(...
                   calcDampedFiberElasticTendonMuscleInfoFcn,...
                   calcDampedFiberElasticTendonInitialMuscleStateFcn,...
                   dampedFiberElasticTendonConfig,...
                   benchRecord,...
                   idx, numberOfSimulations,...
                   flag_appendEnergetics,...
                   flag_useOctave);    

          nameModification = '';
          if(dampedFiberElasticTendonConfig.useElasticTendon == 1)
            nameModification = 'ElasticTendon';
          else
            nameModification = 'RigidTendon';          
          end


          modelName = 'benchRecordHill_';


          save([outputFolder,modelName,nameModification,outputFileEndingHill,'.mat'],...
            'benchRecord',...
            'inputFunctions',...                          
            'normFiberLength',...
            'nominalForce',...
            'amplitudeMM',...
            'bandwidthHz',... 
            'numberOfSimulations',...
            'flag_useElasticTendon',...
            'flag_useFiberDamping',...
            'scalePassiveForceLengthCurve',...
            'muscleArchitecture',...
            'sarcomereProperties',...
            'normMuscleCurves');
          fprintf('%i / %i\n', idx, numberOfSimulations);
          idx=idx+1;
          pause(0.1);
       end
     end
   end
 end


success = 1;