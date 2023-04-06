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

function [success] = runStandardForceLengthVelocitySimulationsDampedEquilibrium(...
                          normFiberLengthAtForceVelocitySample,...
                          flag_useElasticTendon,...
                          flag_useFiberDamping,...
                          fiberDampingCoefficient,...
                          musculotendonProperties,...
                          sarcomereProperties,...
                          normMuscleCurves,...
                          outputFileEndingHill, ...
                          outputFolder,...
                          flag_activeForceLengthSimulations,...
                          flag_passiveForceLengthSimulations,...
                          flag_forceVelocitySimulations,...
                          numberOfLengthSteps,...
                          numberOfVelocitySteps,...
                          maxShorteningVelocity,...
                          forceVelocityNormFiberHalfLength,...
                          flag_useOctave)
                        
success = 0;


nameModification = '';
if(flag_useElasticTendon == 1)
 nameModification = 'ElasticTendon';
else
 nameModification = 'RigidTendon';          
end

%Standard boiler plate setup
excitationMaxFcn = @(argT)calcStepFunction(argT,...
                              0.0,...
                              Inf,...
                              1.0);

excitationMaxWithStartDelayFcn = @(argT)calcStepFunction(argT,...
                              0.1,...
                              Inf,...
                              1.0);

excitationZeroFcn = @(argT)calcStepFunction(argT,...
                              0.,...
                              Inf,...
                              0);

activationFcn = @(argU,argA)calcFirstOrderActivationDerivative(...
                  argU,argA, sarcomereProperties.activationTimeConstant,...
                             sarcomereProperties.deactivationTimeConstant,0);

%%
%Fields that are updated per simulation
%%
benchConfig.npts                  = []; 
benchConfig.initialState          = [];
benchConfig.initialActivation     = 0;
benchConfig.pathFcn               = [];
benchConfig.excitationFcn         = [];
benchConfig.activationFcn         = activationFcn; 
benchConfig.tspan                 = [];  

%%
%Fields that stay fixed
%%
benchConfig.relTol                = 1e-6;
benchConfig.absTol                = 1e-6;
benchConfig.minActivation         = 0;
benchConfig.color0                = [0,0,1].*0.5;
benchConfig.color1                = [0,0,1];

nStates=0;
labelStates={''};
if(flag_useElasticTendon==1)
  nStates=1;
  labelStates = {'$$\ell_{CE}$$'};
end

minActivation=0;
if(flag_useFiberDamping==0)
  minActivation=0.001;
end

benchConfig.numberOfMuscleStates  = nStates;
benchConfig.stateLabels           = labelStates;
benchConfig.name                  = '';

benchConfig.useFiberDamping       = flag_useFiberDamping;
benchConfig.useElasticTendon      = flag_useElasticTendon;
benchConfig.damping               = fiberDampingCoefficient;
benchConfig.iterMax               = 100;
benchConfig.tol                   = 1e-6;
benchConfig.minActivation         = minActivation;

ltSlk     =  musculotendonProperties.tendonSlackLength;
fiso      =  musculotendonProperties.fiso;
lceOpt    =  musculotendonProperties.optimalFiberLength;
alphaOpt  =  musculotendonProperties.pennationAngle;

lpEnd = 0;
lpStart = 0;

if(flag_activeForceLengthSimulations==1)

  
  disp('Running Damped-Equilibrium Active Force Length Simulations');
  benchRecord = [];


  lceNormMin            = normMuscleCurves.activeForceLengthCurve.xEnd(1);
  lceNormMax            = normMuscleCurves.activeForceLengthCurve.xEnd(2);
  lceDelta = (lceNormMax-lceNormMin)/20;
  lceSteps              = numberOfLengthSteps;
  normFiberLengthVector = [(lceNormMin):(lceNormMax-lceNormMin)/(lceSteps-1):(lceNormMax)]';

  [val, idx] = min(abs(normFiberLengthVector-1));
  normFiberLengthVector(idx)=1; %The optimal fiber length occurs perfectly once.

  tendonForceLengthCurveInverse = ...
      createInverseCurve(normMuscleCurves.tendonForceLengthCurve);

  %Fixed-path length simulations
  lpStart = 0;
  lpEnd   = 0;
 



  %Active force-length simulations
  for z=1:1:length(normFiberLengthVector)
  
   
    z01     = (z-1)/(length(normFiberLengthVector)-1);
    tmax    = 2.5 + max(normFiberLengthVector)*1;
               
    lceN    = normFiberLengthVector(z,1);
    lce     = lceN*lceOpt;    
    fibKin  = calcFixedWidthPennatedFiberKinematicsAlongTendon(...
                  lce,0,lceOpt,alphaOpt);
    lceAT   = fibKin.fiberLengthAlongTendon;
    
    falN    = calcBezierYFcnXDerivative(lceN, normMuscleCurves.activeForceLengthCurve,0);
    fpeN    = calcBezierYFcnXDerivative(lceN, normMuscleCurves.fiberForceLengthCurve,0);
    fceATN  = (falN + fpeN)*cos(fibKin.pennationAngle);
    ltN     = calcBezierYFcnXDerivative(fceATN, tendonForceLengthCurveInverse,0);
    if(ltN<1)
        ltN=1;
    end
    if(flag_useElasticTendon==0)
      ltN=1;
    end
    
    lp = lceAT + ltSlk*ltN;

    if(z==1)
      lpStart=lp;
    end
    if(z==length(normFiberLengthVector))
      lpEnd=lp;
    end
    
    %Update benchConfig       
    benchConfig.npts                  = round(100*tmax); 
    benchConfig.initialActivation     = minActivation;
    benchConfig.pathFcn               = @(argT)calcConstantFunction(argT, lp);
    benchConfig.excitationFcn         = excitationMaxWithStartDelayFcn;
    benchConfig.activationFcn         = activationFcn; 
    benchConfig.tspan                 = [0, tmax]; 


    calcDampedFiberElasticTendonMuscleInfoFcn =...
                @(actState1,pathState2,mclState3)...
                calcMillard2012DampedEquilibriumMuscleInfo(  ...
                                    actState1,...
                                    pathState2, ... 
                                    mclState3,...                                                                           
                                    musculotendonProperties,...
                                    normMuscleCurves,...
                                    benchConfig);  
                                  
    calcDampedFiberElasticTendonInitialMuscleStateFcn = ...
                @(actState1,pathState2,calcMuscleInfo3, initConfig4) ...
                    calcInitialMuscleState(actState1,...
                                           pathState2,...
                                           musculotendonProperties,...
                                           calcMuscleInfo3,...
                                           initConfig4);
                                          
    benchConfig.numberOfMuscleStates = nStates;
    benchConfig.name = 'DFE';
    benchConfig.eventFcn = [];                                  
                                        

    flag_appendEnergetics=0;           
    benchRecord = runPrescribedLengthActivationSimulation(...
                     calcDampedFiberElasticTendonMuscleInfoFcn,...
                     calcDampedFiberElasticTendonInitialMuscleStateFcn,...
                     benchConfig,...
                     benchRecord,...
                     z, length(normFiberLengthVector),...
                     flag_appendEnergetics,...
                     flag_useOctave);
    fprintf('%i / %i\n', z, length(normFiberLengthVector));                                 

  end

  save([outputFolder,'activeForceLengthHill_',nameModification,outputFileEndingHill,'.mat'],...
      'benchRecord',...
      'flag_useElasticTendon',...
      'musculotendonProperties',...
      'normMuscleCurves');

  %Complete the passive-force-length experiments
end

if(flag_passiveForceLengthSimulations==1)

  disp('Running Damped-Equilibrium Passive Force Length Simulations');
  
  benchRecord = [];

  lceNormMin            = normMuscleCurves.activeForceLengthCurve.xEnd(1);
  lceNormMax            = normMuscleCurves.activeForceLengthCurve.xEnd(2);  
  
  lceMin      = lceNormMin*lceOpt;    
  fibKin      = calcFixedWidthPennatedFiberKinematicsAlongTendon(...
                lceMin,0,lceOpt,alphaOpt);
  lceATMin    = fibKin.fiberLengthAlongTendon;
    
  lceMax      = lceNormMax*lceOpt;    
  fibKin      = calcFixedWidthPennatedFiberKinematicsAlongTendon(...
                lceMax,0,lceOpt,alphaOpt);
  lceATMax    = fibKin.fiberLengthAlongTendon;  
  
  lpStart = lceATMin + ltSlk;
  lpEnd   = lceATMax + ltSlk;
  
  tRampStart  = 0.1;
  tRampEnd    = ((lpEnd-lpStart)/lceOpt)*10 + tRampStart;
  rampSlope   = (lpEnd-lpStart)/(tRampEnd-tRampStart);
  pathFcn     = @(argT)calcRampStateSharp(argT,tRampStart,tRampEnd,lpStart,rampSlope);

  benchConfig.npts               = round(100*tRampEnd); 
  benchConfig.initialActivation  = minActivation;
  benchConfig.pathFcn            = pathFcn;
  benchConfig.excitationFcn      = excitationZeroFcn;
  benchConfig.activationFcn      = activationFcn; 
  benchConfig.tspan              = [0, tRampEnd]; 


  calcDampedFiberElasticTendonMuscleInfoFcn =...
              @(actState1,pathState2,mclState3)...
              calcMillard2012DampedEquilibriumMuscleInfo(  ...
                                  actState1,...
                                  pathState2, ... 
                                  mclState3,...                                                                           
                                  musculotendonProperties,...
                                  normMuscleCurves,...
                                  benchConfig);  
                                
  calcDampedFiberElasticTendonInitialMuscleStateFcn = ...
              @(actState1,pathState2,calcMuscleInfo3, initConfig4) ...
                  calcInitialMuscleState(actState1,...
                                         pathState2,...
                                         musculotendonProperties,...
                                         calcMuscleInfo3,...
                                         initConfig4);
                                        
  benchConfig.numberOfMuscleStates = nStates;
  benchConfig.name = 'DFE';
  benchConfig.eventFcn = [];       

  flag_appendEnergetics=0;           
  benchRecord = runPrescribedLengthActivationSimulation(...
                   calcDampedFiberElasticTendonMuscleInfoFcn,...
                   calcDampedFiberElasticTendonInitialMuscleStateFcn,...
                   benchConfig,...
                   benchRecord,...
                   1, 1,...
                   flag_appendEnergetics,...
                   flag_useOctave);

  save([outputFolder,'passiveForceLengthHill_',...
        nameModification,outputFileEndingHill,'.mat'],...
      'benchRecord',...
      'flag_useElasticTendon',...
      'musculotendonProperties',...
      'normMuscleCurves');  

  benchRecord = [];
end

if(flag_forceVelocitySimulations==1)

  disp('Running Damped-Equilibrium Force Velocity Simulations');
  benchRecord = [];

  lceNdelta   = forceVelocityNormFiberHalfLength;%0.20;
  dlceNNSteps = numberOfVelocitySteps;
  dlceNNMax   = 1;
  normVelocitySeries = [0:(dlceNNMax/(dlceNNSteps-1)):dlceNNMax]';

  dlceNNMin = (dlceNNMax/(dlceNNSteps-1)); 

  tendonForceLengthCurveInverse = ...
      createInverseCurve(normMuscleCurves.tendonForceLengthCurve);

  
  %%
  %Solve for the starting, mid, and ending lengths of the musculotendon
  %%
  lceNCenter = normFiberLengthAtForceVelocitySample;
  
  lceN  = (lceNCenter+lceNdelta);        
  lce     = lceN*lceOpt;    
  fibKin  = calcFixedWidthPennatedFiberKinematicsAlongTendon(...
                lce,0,lceOpt,alphaOpt);
  lceAT   = fibKin.fiberLengthAlongTendon;
  falN    = calcBezierYFcnXDerivative(lceN, normMuscleCurves.activeForceLengthCurve,0);
  fpeN    = calcBezierYFcnXDerivative(lceN, normMuscleCurves.fiberForceLengthCurve,0);
  fceATN  = (falN + fpeN)*cos(fibKin.pennationAngle);
  ltN     = calcBezierYFcnXDerivative(fceATN, tendonForceLengthCurveInverse,0);
  if(flag_useElasticTendon==0)
      ltN=1;
  end
  if(ltN < 1)
      ltN=1;
  end
  lpLong  = lceAT + ltSlk*ltN;    

  lceN  = (lceNCenter);        
  lce     = lceN*lceOpt;    
  fibKin  = calcFixedWidthPennatedFiberKinematicsAlongTendon(...
                lce,0,lceOpt,alphaOpt);
  lceAT   = fibKin.fiberLengthAlongTendon;
  falN    = calcBezierYFcnXDerivative(lceN, normMuscleCurves.activeForceLengthCurve,0);
  fpeN    = calcBezierYFcnXDerivative(lceN, normMuscleCurves.fiberForceLengthCurve,0);
  fceATN  = (falN + fpeN)*cos(fibKin.pennationAngle);
  ltN     = calcBezierYFcnXDerivative(fceATN, tendonForceLengthCurveInverse,0);
  if(ltN < 1)
      ltN=1;
  end  
  if(flag_useElasticTendon==0)
      ltN=1;
  end
  
  lpMid  = lceAT + ltSlk*ltN;        
  
  lceN    = (lceNCenter-lceNdelta);
  lce     = lceN*lceOpt;    
  fibKin  = calcFixedWidthPennatedFiberKinematicsAlongTendon(...
                lce,0,lceOpt,alphaOpt);
  lceAT   = fibKin.fiberLengthAlongTendon;
  falN    = calcBezierYFcnXDerivative(lceN, normMuscleCurves.activeForceLengthCurve,0);
  fpeN    = calcBezierYFcnXDerivative(lceN, normMuscleCurves.fiberForceLengthCurve,0);
  fceATN  = (falN + fpeN)*cos(fibKin.pennationAngle);
  ltN     = calcBezierYFcnXDerivative(fceATN, tendonForceLengthCurveInverse,0);
  if(ltN < 1)
      ltN=1;
  end  
  if(flag_useElasticTendon==0)
      ltN=1;
  end
  
  lpShort = lceAT + ltSlk*ltN;    
  
  %%
  %Perform the simulations
  %%
  
  dlceATSlowest = dlceNNMin;
  tmaxLongest = abs(lpEnd-lpStart)/abs(dlceATSlowest) + 2;

  for k=1:1:2
    velDir = -1;

    if(k > 1)
      velDir = 1;
    end      
    
    for z=1:1:(dlceNNSteps)

      idx = (k-1)*length(normVelocitySeries) + z;

      %Construct the path function 
      pathFcn = [];
      tmax = 0;
      tRampDuration = 0;
      if(abs(normVelocitySeries(z,1)) > sqrt(eps))
        lpStart = lpLong;
        lpEnd   = lpShort;
        if(velDir > 0)
          lpStart   = lpShort;        
          lpEnd     = lpLong;
        end
        vceAtlceOpt = velDir*normVelocitySeries(z,1)*(lceOpt*maxShorteningVelocity);
        fibKin = calcFixedWidthPennatedFiberKinematicsAlongTendon(...
          lceOpt,vceAtlceOpt,lceOpt,alphaOpt);
        dlceAT = fibKin.fiberVelocityAlongTendon;       

        tRampStart    = 0.1;
        tRampDuration = (lpEnd-lpStart)/dlceAT;
        tRampEnd      = tRampStart + tRampDuration;
        tRampMid      = tRampStart + tRampDuration*0.5;
        rampSlope     = (lpEnd-lpStart)/tRampDuration;
        assert(tRampEnd > tRampStart);

        pathFcn = @(argT)calcRampStateSharp(argT, tRampStart, tRampEnd,...
                                                     lpStart, rampSlope);
        tmax = tRampEnd+0.1;
      else
        pathFcn = @(argT)calcConstantFunction(argT,lpMid);
        tmax = 0.2;
        tRampMid = tmax/2;
      end


      
      benchConfig.pathFcn = pathFcn;

      %Update benchConfig       
      benchConfig.npts                  = round(100*tmaxLongest); 
      benchConfig.initialActivation     = 1.0;
      benchConfig.excitationFcn         = excitationMaxFcn;
      benchConfig.activationFcn         = activationFcn; 
      benchConfig.tspan                 = [0, tmax]; 


      calcDampedFiberElasticTendonMuscleInfoFcn =...
                  @(actState1,pathState2,mclState3)...
                  calcMillard2012DampedEquilibriumMuscleInfo(  ...
                                      actState1,...
                                      pathState2, ... 
                                      mclState3,...                                                                           
                                      musculotendonProperties,...
                                      normMuscleCurves,...
                                      benchConfig);  
                                    
      calcDampedFiberElasticTendonInitialMuscleStateFcn = ...
                  @(actState1,pathState2,calcMuscleInfo3, initConfig4) ...
                      calcInitialMuscleState(actState1,...
                                             pathState2,...
                                             musculotendonProperties,...
                                             calcMuscleInfo3,...
                                             initConfig4);
                                            
      benchConfig.numberOfMuscleStates = nStates;
      benchConfig.name = 'DFE';
                                                                       
      lceThreshold = lceNCenter*lceOpt;
                  
      idxThresholdState   = 1;
      eventThresholdValue = lceThreshold; 
      eventDirection      = velDir;
      if(flag_useElasticTendon == 0 || abs(normVelocitySeries(z,1)) <= sqrt(eps))
        idxThresholdState   = -1;
        eventThresholdValue = tRampMid;
        eventDirection      = 1;
      end
                  
      idxThresholdState = idxThresholdState + 1; %activation state is pre-pended            
      
      benchConfig.eventFcn = @(argT,argY)eventStateThreshold(...
                                argT,argY,idxThresholdState,...
                                eventThresholdValue,eventDirection);



      flag_appendEnergetics=0;           
      benchRecord = runPrescribedLengthActivationSimulation(...
                      calcDampedFiberElasticTendonMuscleInfoFcn,...
                      calcDampedFiberElasticTendonInitialMuscleStateFcn,...
                      benchConfig,...
                      benchRecord,...
                      idx, length(normVelocitySeries)*2,...
                      flag_appendEnergetics,...
                      flag_useOctave);

      fprintf('%i / %i\n', idx, length(normVelocitySeries)*2);   
    end
  end

  save([outputFolder,'forceVelocityHill_',nameModification,outputFileEndingHill,'.mat'],...
    'benchRecord',...
    'flag_useElasticTendon',...
    'musculotendonProperties',...
    'normMuscleCurves');

end


success = 1;


success = 1;