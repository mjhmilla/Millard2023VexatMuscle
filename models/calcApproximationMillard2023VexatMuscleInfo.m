function [muscleStateDerivative, mtInfo] ...
  = calcApproximationMillard2023VexatMuscleInfo(...
        time,...
        muscleIntegratedState,...
        muscleStateInitial,...           
        activationState,...                                            
        pathFcn, ... 
        pathStateInitial,...
        muscleArchitecture,...
        sarcomereProperties,...
        normMuscleCurves,...
        modelConfig,...
        approximationMethod)
%%
% This function implements approximations of the model that
% are used to fit parts of the model to experimental data.
% These approximations absolutely strip the model down to the 
% bare minimum amount of computation just for the purposes
% of a specific fitting routine. This is not intended for general 
% use.
%%

muscleStateDerivative=zeros(size(muscleIntegratedState)).*nan;

switch approximationMethod
  case 'quasiStaticSkinnedFiberApproximation'
       pathState=pathFcn(time);

       %Translate laH along with the path
       %Assume dlaH matches the path velocity       
       muscleStateUpd = muscleStateInitial;
       muscleStateUpd(1) = pathState(1)*0.5;
       muscleStateUpd(2) = muscleStateInitial(2)...
            + (pathState(2)-pathStateInitial(2))*0.5;
       muscleStateUpd(3) = muscleIntegratedState(1);
       assert(length(muscleIntegratedState)==1,...
              'Error: expected only one state in muscleIntegratedState');

       %Approximate state
       dlaH   = muscleStateUpd(1);  
       laH    = muscleStateUpd(2);  
       l1H    = muscleStateUpd(3); 
      

      if(normMuscleCurves.useCalibratedCurves == 1)
        calcFvDer      = @(arg1, arg2)calcBezierYFcnXDerivative(arg1, ...
                          normMuscleCurves.fiberForceVelocityCalibratedCurve,...
                          arg2);
        calcFalDer     = @(arg1, arg2)calcBezierYFcnXDerivative(arg1, ...
                          normMuscleCurves.activeForceLengthCalibratedCurve,...
                          arg2);  
      else
          
        calcFvDer      = @(arg1, arg2)calcBezierYFcnXDerivative(arg1, ...
                          normMuscleCurves.fiberForceVelocityCurve,...
                          arg2);
        calcFalDer     = @(arg1, arg2)calcBezierYFcnXDerivative(arg1, ...
                        normMuscleCurves.activeForceLengthCurve,...
                        arg2);
      end  

      calcF1HDer       = @(arg1, arg2)calcBezierYFcnXDerivative(arg1, ...
                         normMuscleCurves.forceLengthProximalTitinCurve, ...
                            arg2);   

      calcF2HDer       = @(arg1, arg2)calcBezierYFcnXDerivative(arg1, ...
                         normMuscleCurves.forceLengthDistalTitinCurve, ...
                            arg2);      

      calcFeHDer  = @(arg1,arg2)calcBezierYFcnXDerivative(arg1, ...
                        normMuscleCurves.forceLengthECMHalfCurve, ...
                        arg2);      

      %Constants
      lceOpt=muscleArchitecture.optimalFiberLength;
      dlceMaxN=muscleArchitecture.maximumNormalizedFiberVelocity;
      lmHN    = sarcomereProperties.normMyosinHalfLength;
      forceVelocityCalibrationFactor = ...
          sarcomereProperties.forceVelocityCalibrationFactor;
      lTitinFixedHN = ...
            sarcomereProperties.ZLineToT12NormLengthAtOptimalFiberLength ...
          + sarcomereProperties.IGDFixedNormLengthAtOptimalFiberLength;
      
      scaleTitinDistal    = sarcomereProperties.scaleTitinDistal;
      scaleTitinProximal  = sarcomereProperties.scaleTitinProximal;
      f1HNPreload = 0;
      if(isfield(sarcomereProperties,'f1HNPreload'))
          f1HNPreload    = sarcomereProperties.f1HNPreload;
      end      
      
      lAHN               = sarcomereProperties.normActinLength;  

      ZLineToT12NormLengthAtOptimalFiberLength = ...
        sarcomereProperties.ZLineToT12NormLengthAtOptimalFiberLength;
      smoothStepFunctionRadius = ...
         sarcomereProperties.smoothStepFunctionRadius;
      lceHNTitinActinBondMin = ...
        sarcomereProperties.normLengthTitinActinBondMinimum*0.5;
      activationThresholdTitin = ...
        sarcomereProperties.activationThresholdTitin;

      betaTAaHN = sarcomereProperties.normMaxActiveTitinToActinDamping;
      betaTApHN = sarcomereProperties.normPassiveTitinToActinDamping;


      %Evaluate arguments to input functions
      vce = pathState(1);
      lce = pathState(2);
      dadt= activationState(1);
      a   = activationState(2);
      
      lceHN = lce*0.5/lceOpt;
      dlceH = vce*0.5;

      % XE forces
      laHN  = laH/lceOpt;
      lamN  = 2*(lmHN+laHN); 
      flN = calcFalDer(lamN,0);
      
      dlfNN = ((dlaH*2)/(lceOpt*dlceMaxN));
      fvN=calcFvDer(dlfNN*forceVelocityCalibrationFactor,0);

      % Prox. Titin spring      
      l1HN  = l1H/lceOpt;
      f1kHN = scaleTitinProximal*(calcF1HDer(l1HN,0)+f1HNPreload );      

      % Distal titin spring
      l2H =lce*0.5-(l1H+lTitinFixedHN*lceOpt);
      l2HN=l2H/lceOpt;
      f2kHN = scaleTitinDistal*calcF2HDer(l2HN,0);

      %Titin-actin bond
      kaTi = a/activationThresholdTitin;
      aTi = 1.-exp((-kaTi*kaTi));
      
      dTiA = (lAHN-(l1HN+ZLineToT12NormLengthAtOptimalFiberLength));
      kTiA = dTiA/smoothStepFunctionRadius;
      uTiA = 0.5+0.5*tanh(kTiA);       

      dTiLMin = lceHN-lceHNTitinActinBondMin;%lceHNZeroFpeN;
      kTiLMin = dTiLMin/smoothStepFunctionRadius;
      uTiLMin = 0.5+0.5*tanh(kTiLMin);

      beta1HNN = betaTApHN + betaTAaHN*aTi*uTiA*uTiLMin; 
      beta2HNN = 0;

      dl1HN = (f2kHN-f1kHN)/beta1HNN;
      dl1H  = dl1HN*lceOpt;
  
      dl2H  = dlceH - dl1H; 
      dl2HN = dl2H/lceOpt;
  
      f1dHN = beta1HNN*dl1HN;
      f2dHN = 0;   
  
      f1HN  = f1kHN + f1dHN;    
      f2HN  = f2kHN;         

      mtInfo.muscleDynamicsInfo.activation                = a;
      mtInfo.muscleDynamicsInfo.activationDerivative      = dadt;

      mtInfo.muscleLengthInfo.fiberPassiveForceLengthMultiplier = f2HN;
      mtInfo.muscleLengthInfo.fiberActiveForceLengthMultiplier  = flN;
      mtInfo.fiberVelocityInfo.fiberForceVelocityMultiplier     = fvN;

      mtInfo.muscleDynamicsInfo.activation                      = a;
      mtInfo.muscleDynamicsInfo.activationDerivative            = dadt;

      mtInfo.muscleDynamicsInfo.normFiberForce  = a*flN*fvN + f2HN;
      mtInfo.muscleDynamicsInfo.normFiberLength = lce/lceOpt;

      muscleStateDerivative = dl1H;
  otherwise
    assert(0,'Error: unrecognized approximationMethod');
end
