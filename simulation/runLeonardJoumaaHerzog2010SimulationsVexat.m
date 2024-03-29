%%
% SPDX-FileCopyrightText: 2023 Matthew Millard <millard.matthew@gmail.com>
%
% SPDX-License-Identifier: MIT
%
% If you use this code in your work please cite the pre-print of this paper
% or the most recent peer-reviewed version of this paper:
%
%    Matthew Millard, David W. Franklin, Walter Herzog. 
%    A three filament mechanistic model of musculotendon force and impedance. 
%    bioRxiv 2023.03.27.534347; doi: https://doi.org/10.1101/2023.03.27.534347 
%
%%

function [success] = runLeonardJoumaaHerzog2010SimulationsVexat(...
                          nominalNormalizedFiberLength,...
                          activeForceKeyPoints,...
                          passiveForceKeyPoints,...
                          timeSpan,...
                          lengthRampKeyPoints,...
                          stimulationKeyTimes,...
                          flag_useElasticTendon,...
                          musculotendonProperties,...
                          sarcomereProperties,...
                          normMuscleCurves,...
                          outputFileEndingVexat, ...
                          outputFolder,...
                          flag_simulateActiveStretch,...
                          flag_simulatePassiveStretch,...
                          flag_useOctave)
success = 0;

disp('Running Opus 31 Leonard, Joumaa, and Herzog 2010 Simulations');

  idx=1;
  benchRecord = [];
%  for idxNormFiberLength = 1:1:length(normFiberLength)
%    for idxActivation = 1:1:length(nominalForce)  
%      flag_setWarning=0;
%      for i=1:1:length(amplitudeMM)        
%        for j=1:1:length(bandwidthHz)
            assert(size(nominalNormalizedFiberLength,1)==1 && ...
                   size(nominalNormalizedFiberLength,2)==1);

            %%
            %Adjust the passive curves
            %%

            %Note:  the ECM and titin curves that Vexat uses have been
            %       fitted s.t. they sum together to match the 
            %       passive-force-length curve. Thus we can more
            %       conveniently evaluate the passive force produced by 
            %       Opus 31 by evaluting the passive curve.

    
            lTitinFixedHN = sarcomereProperties.ZLineToT12NormLengthAtOptimalFiberLength ...
                          + sarcomereProperties.IGDFixedNormLengthAtOptimalFiberLength;

            fTiN = 1;
            lerr=1;
            i=1;
            
            %Fit the titin forces to a point n between lopt and the passive
            %breaking point
            n = 0.75;
            passiveForceMidKeyPoint = zeros(1,2);
            passiveForceMidKeyPoint(1,1) = passiveForceKeyPoints(1,1) ...
              +n*(passiveForceKeyPoints(2,1)-passiveForceKeyPoints(1,1));
            passiveForceMidKeyPoint(1,2) = passiveForceKeyPoints(1,2) ...
              +n*(passiveForceKeyPoints(2,2)-passiveForceKeyPoints(1,2));


            while(abs(lerr) > 1e-3 && i < 100)

                l1N = calcBezierYFcnXDerivative(fTiN,...
                    normMuscleCurves.forceLengthProximalTitinInverseCurve,0);            
                l2N = calcBezierYFcnXDerivative(fTiN,...
                    normMuscleCurves.forceLengthDistalTitinInverseCurve,0);
    
                D_l1N_D_fN = calcBezierYFcnXDerivative(fTiN,...
                    normMuscleCurves.forceLengthProximalTitinInverseCurve,1);            
                D_l2N_D_fN = calcBezierYFcnXDerivative(fTiN,...
                    normMuscleCurves.forceLengthDistalTitinInverseCurve,1);
                
                lerr = 2*(l1N+l2N+lTitinFixedHN)-passiveForceMidKeyPoint(1,1);
    
                D_lerr_D_fN = 2*(D_l1N_D_fN+D_l2N_D_fN);
    
                dfTiN = -lerr/D_lerr_D_fN;
    
                fTiN = fTiN+dfTiN;
                i=i+1;
            end
            
           
            fTiNRatio = passiveForceMidKeyPoint(1,2)/(fTiN);
            
              
              %if(flag_fitTitin==1)
              %  sarcomereProperties.scaleECM = 0;
                
                %A simulation is required to get this scaling value 
                %(the number below).
              %  activeScaling = 1;%activeForceKeyPoints(2,2)/3.77761; 

              %  sarcomereProperties.scaleTitinDistal= fTiNRatio*activeScaling;
              %  sarcomereProperties.scaleTitinProximal = fTiNRatio/activeScaling;
              %else
                sarcomereProperties.scaleECM = 0;
                sarcomereProperties.scaleTitinProximal = 1;%fTiNRatio;
                sarcomereProperties.scaleTitinDistal   = 1;%fTiNRatio;
              %end
              

%             else
%               %Remove ECM until the desired target is met.
%               ecmFraction   = sarcomereProperties.extraCellularMatrixPassiveForceFraction;
%               titinFraction = 1-ecmFraction;
% 
%               sarcomereProperties.scaleECM = (fpe1-ecmFraction)/ecmFraction;              
%             end

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
                        normMuscleCurves.fiberForceLengthCurve  ,0)*fTiNRatio;

            %Evaluate the passive path length
            ft0 = fpeN0*cos(alpha0);
            ltN0 = 1;          
            if(flag_useElasticTendon==1)
                tendonForceLengthCurveInverse = ...
                    createInverseCurve(normMuscleCurves.tendonForceLengthCurve);
                ltN0    = calcBezierYFcnXDerivative(ft0, tendonForceLengthCurveInverse,0);                
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
            ltN1   = 1;
            if(flag_useElasticTendon==1)
                tendonForceLengthCurveInverse = ...
                    createInverseCurve(normMuscleCurves.tendonForceLengthCurve);
                ltN1    = calcBezierYFcnXDerivative(ft1, tendonForceLengthCurveInverse,0);
            end            
            lceAT1  = pathStartLength - ltN1*ltslk;
            dlceAT1 = 0;

            fibKin = calcFixedWidthPennatedFiberKinematics(...
                                           lceAT1,dlceAT1,lceOpt,alphaOpt);
            alpha1  = fibKin.pennationAngle;
            lceN1   = fibKin.fiberLength/lceOpt;
            dlceN1  = fibKin.fiberVelocity...
                        /musculotendonProperties.maximumNormalizedFiberVelocity;
           
            %%
            %  Evaluate the fiber length when its developing the desired force
            %
            %  N.B. Since Opus 31 has elastic cross bridges the relationship between
            %       fiber length and active force does not exactly follow the curve.
            %       However, since the cross bridges are so stiff I'm going to 
            %       ignore this fact when evalulating the activation necessary
            %       to initialize the model: the error should be small.
            %%

%             falN          = calcBezierYFcnXDerivative(lceN1,  ...
%                               normMuscleCurves.activeForceLengthCurve ,0);
%             DfalN_DlceN   = calcBezierYFcnXDerivative(lceN1,  ...
%                               normMuscleCurves.activeForceLengthCurve ,1);
% 
%             fvN           = calcBezierYFcnXDerivative(dlceN1, ...
%                               normMuscleCurves.fiberForceVelocityCurve,0);    
%             DfvN_DdlceN   = calcBezierYFcnXDerivative(dlceN1, ...
%                               normMuscleCurves.fiberForceVelocityCurve,1);
% 
%             fpeN          = calcBezierYFcnXDerivative(lceN1, ...
%                               normMuscleCurves.fiberForceLengthCurve  ,0);
%             DfpeN_DlceN   = calcBezierYFcnXDerivative(lceN1, ...
%                               normMuscleCurves.fiberForceLengthCurve  ,1);
%             
%             fcnN   = ft1/cos(alpha1);
%             activation = max( (fcnN-fpeN)/(falN*fvN), 0);
%             
%             assert(activation <= 1.0);
% 
%             if( fcnN < fpeN && flag_setWarning==0)              
%               fprintf('  Warning: \tdesired force (%1.3f) > passive force (%1.3f) at length (%1.3f)\n',...
%                 fcnN, fpeN,lceN1);
%               flag_setWarning=1;
%             end       

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
                                   argT,rampStartTime+0.1,rampEndTime+0.1,...
                                   pathStartLength,0);


          %%
          % Set up the simulation structs
          %%
            benchConfig.npts                  = round(timeSpan(1,2)-timeSpan(1,1))*100;
            benchConfig.relTol                = 1e-6;
            benchConfig.absTol                = 1e-6;
            benchConfig.minActivation         = 0;
            benchConfig.color0                = [0,0,1].*0.5;
            benchConfig.color1                = [0,0,1];

            nStates = 0;
            labelStates = {''};
            if(flag_useElasticTendon==1)
              nStates = 4;
              labelStates= {'$$\ell_{CE}$$','$$\dot{\ell}_{a}$$',...
                 '$$\ell_{a}$$','$$\ell_1$$'};%,'$$f_{e}^1$$', '$$f_{e}^2$$'};
            else
              nStates = 3;
              labelStates= {'$$\dot{\ell}_{a}$$',...
                 '$$\ell_{a}$$','$$\ell_1$$'};%,'$$f_{e}^1$$', '$$f_{e}^2$$'};            
            end

            benchConfig.numberOfMuscleStates  = nStates;
            benchConfig.stateLabels  = labelStates;
            benchConfig.name                  = '';
            benchConfig.initialState          = [];
            benchConfig.initialActivation     = 0;
            benchConfig.pathFcn               = [];
            benchConfig.excitationFcn         = [];
            benchConfig.activationFcn         = activationFcn; 
            benchConfig.tspan                 = timeSpan;  

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
              'initializeState'         , 0                     );          

            modelConfig.initializeState =1;
            activationState0 = [0;excitationSquareFcn(0)];
            pathState0 = pathLengthRampFcn(0);
            muscleState0 = zeros(nStates,1);
            mtInfo =calcMillard2023VexatMuscleInfo(activationState0,...
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
                 calcMillard2023VexatMuscleInfo(   activationState1,...
                                                    pathState2,...
                                                    mclState3,...
                                                    musculotendonProperties,...
                                                    sarcomereProperties,...
                                                    normMuscleCurves,...
                                                    modelConfig);

            flag_appendEnergetics =0;
            numberOfSimulations = flag_simulateActiveStretch...
                                + flag_simulatePassiveStretch;

           benchRecord = [];
           idx = 1;
           if(flag_simulateActiveStretch == 1)
             benchConfig.pathFcn               = pathLengthRampFcn;
             benchConfig.excitationFcn         = excitationSquareFcn;             
             benchRecord = runPrescribedLengthActivationSimulation(...
                                           calcMillard2023VexatMuscleInfoFcn,...
                                           [],...
                                           benchConfig,...
                                           benchRecord,...
                                           idx, numberOfSimulations,...
                                           flag_appendEnergetics,...
                                           flag_useOctave);
             fprintf('%i / %i\n', idx, numberOfSimulations);                            
             idx=idx+1;
           end


           if(flag_simulatePassiveStretch == 1)
             
              benchConfig.pathFcn               = pathLengthRampFcn;
              benchConfig.excitationFcn         = excitationZeroFcn;
              benchRecord = runPrescribedLengthActivationSimulation(...
                                           calcMillard2023VexatMuscleInfoFcn,...
                                           [],...
                                           benchConfig,...
                                           benchRecord,...
                                           idx, numberOfSimulations,...
                                           flag_appendEnergetics,...
                                           flag_useOctave);
              fprintf('%i / %i\n', idx, numberOfSimulations);                           
              idx=idx+1;
          end
                                       
         nameModification = '';
         if(benchConfig.useElasticTendon == 1)
           nameModification = 'ElasticTendon';
         else
           nameModification = 'RigidTendon';          
         end

         strTitinAdj = '_TiDefault';
         %if(flag_fitTitin==1)
         %  strTitinAdj='_TiAdj';           
         %end
         

         save([outputFolder,'benchRecordVexat_',outputFileEndingVexat,'.mat'],...
            'benchRecord',...
            'nominalNormalizedFiberLength',...
            'timeSpan',...
            'lengthRampKeyPoints',...
            'stimulationKeyTimes',...
            'flag_useElasticTendon',...
            'musculotendonProperties',...
            'sarcomereProperties',...
            'normMuscleCurves');
          
          %idx=idx+1;
          %pause(0.1);
%        end
%      end
%    end
%  end

success = 1;
                        
                        