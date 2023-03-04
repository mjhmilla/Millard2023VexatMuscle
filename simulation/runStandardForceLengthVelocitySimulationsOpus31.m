function [success] = runStandardForceLengthVelocitySimulationsOpus31(...
                          normFiberLengthAtForceVelocitySample,...
                          flag_useElasticTendon,...
                          musculotendonProperties,...
                          sarcomereProperties,...
                          normMuscleCurves,...
                          outputFileEndingOpus31, ...
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
benchConfig.relTol                = 1e-8;
benchConfig.absTol                = 1e-10;
benchConfig.minActivation         = 0;
benchConfig.color0                = [0,0,1].*0.5;
benchConfig.color1                = [0,0,1];

nStates = 0;
labelStates = {''};
if(flag_useElasticTendon==1)
  nStates = 4;
  labelStates= {'$$\ell_{CE}$$','$$\dot{\ell}_{a}$$','$$\ell_{a}$$','$$\ell_1$$'};%,'$$\lambda$$'};
else
  nStates = 3;
  labelStates= {'$$\dot{\ell}_{a}$$','$$\ell_{a}$$','$$\ell_1$$'};%,'$$\lambda$$'};        
end

benchConfig.numberOfMuscleStates  = nStates;
benchConfig.stateLabels           = labelStates;
benchConfig.name                  = '';

benchConfig.useElasticTendon      = flag_useElasticTendon;
benchConfig.iterMax               = 100;
benchConfig.tol                   = min(benchConfig.relTol,benchConfig.absTol);

loopTolerance = min(benchConfig.relTol,benchConfig.absTol)/100;

modelConfig = struct( ...
  'iterMax'                 , 100             , ...
  'tol'                     , loopTolerance   , ... 
  'tolInit'                 , sqrt(eps)       , ...
  'minActivation'           , 0.0             , ...
  'useElasticTendon'        , flag_useElasticTendon , ...
  'initializeState'         , 0                     );          

ltSlk     =  musculotendonProperties.tendonSlackLength;
fiso      =  musculotendonProperties.fiso;
lceOpt    =  musculotendonProperties.optimalFiberLength;
alphaOpt  =  musculotendonProperties.pennationAngle;



if(flag_activeForceLengthSimulations==1)

  
  disp('Running Opus 31 Active Force Length Simulations');
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
    if(ltN < 1)
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
    benchConfig.initialActivation     = excitationMaxWithStartDelayFcn(0);
    benchConfig.pathFcn               = @(argT)calcConstantFunction(argT, lp);
    benchConfig.excitationFcn         = excitationMaxWithStartDelayFcn;
    benchConfig.activationFcn         = activationFcn; 
    benchConfig.tspan                 = [0, tmax]; 

    %Initialize the model
    modelConfig.initializeState = 1;
    activationState0            = [0;benchConfig.excitationFcn(0)];
    pathState0                  = benchConfig.pathFcn(0);
    muscleState0                = zeros(nStates,1);

    mtInfo = calcMillard2019MuscleInfoOpus31(activationState0,...
                                            pathState0,...
                                            muscleState0,...
                                            musculotendonProperties,...
                                            sarcomereProperties,...
                                            normMuscleCurves,...
                                            modelConfig);

    muscleState0                     = mtInfo.state.value;

    modelConfig.initializeState      = 0;   
    benchConfig.initialState         = muscleState0;

    benchConfig.minimumActivation    = 0;
    benchConfig.name                 = 'Opus31';
    benchConfig.eventFcn             = [];           

    calcMillard2019MuscleInfoOpus31Fcn = ...
         @(activationState1,pathState2,mclState3) ...
         calcMillard2019MuscleInfoOpus31(   activationState1,...
                                            pathState2,...
                                            mclState3,...
                                            musculotendonProperties,...
                                            sarcomereProperties,...
                                            normMuscleCurves,...
                                            modelConfig);

    flag_appendEnergetics=0;           
    benchRecord = runPrescribedLengthActivationSimulation(...
                     calcMillard2019MuscleInfoOpus31Fcn,...
                     [],...
                     benchConfig,...
                     benchRecord,...
                     z, length(normFiberLengthVector),...
                     flag_appendEnergetics,...
                     flag_useOctave);
    fprintf('%i / %i\n', z, length(normFiberLengthVector));                                 

  end

  save([outputFolder,'activeForceLengthOpus31_',nameModification,outputFileEndingOpus31,'.mat'],...
      'benchRecord',...
      'flag_useElasticTendon',...
      'musculotendonProperties',...
      'sarcomereProperties',...
      'normMuscleCurves');
end

if(flag_passiveForceLengthSimulations==1)

  %Complete the passive-force-length experiments
  disp('Running Opus 31 Passive Force Length Simulations');
  
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
  benchConfig.initialActivation  = 0;
  benchConfig.pathFcn            = pathFcn;
  benchConfig.excitationFcn      = excitationZeroFcn;
  benchConfig.activationFcn      = activationFcn; 
  benchConfig.tspan              = [0, tRampEnd]; 

  %Initialize the model
  modelConfig.initializeState = 1;
  activationState0            = [0;benchConfig.excitationFcn(0)];
  pathState0                  = benchConfig.pathFcn(0);
  muscleState0                = zeros(nStates,1);

  mtInfo = calcMillard2019MuscleInfoOpus31(activationState0,...
                                          pathState0,...
                                          muscleState0,...
                                          musculotendonProperties,...
                                          sarcomereProperties,...
                                          normMuscleCurves,...
                                          modelConfig);

  muscleState0                     = mtInfo.state.value;

  modelConfig.initializeState      = 0;   
  benchConfig.initialState         = muscleState0;

  benchConfig.minimumActivation    = 0;
  benchConfig.name                 = 'Opus31';
  benchConfig.eventFcn             = [];           

  calcMillard2019MuscleInfoOpus31Fcn = ...
       @(activationState1,pathState2,mclState3) ...
       calcMillard2019MuscleInfoOpus31(   activationState1,...
                                          pathState2,...
                                          mclState3,...
                                          musculotendonProperties,...
                                          sarcomereProperties,...
                                          normMuscleCurves,...
                                          modelConfig);

  flag_appendEnergetics=0;           
  benchRecord = runPrescribedLengthActivationSimulation(...
                   calcMillard2019MuscleInfoOpus31Fcn,...
                   [],...
                   benchConfig,...
                   benchRecord,...
                   1, 1,...
                   flag_appendEnergetics,...
                   flag_useOctave);

  save([outputFolder,'passiveForceLengthOpus31_',...
        nameModification,outputFileEndingOpus31,'.mat'],...
      'benchRecord',...
      'flag_useElasticTendon',...
      'musculotendonProperties',...
      'sarcomereProperties',...
      'normMuscleCurves');  

  benchRecord = [];
end

if(flag_forceVelocitySimulations==1)

  disp('Running Opus 31 Force Velocity Simulations');
  benchRecord = [];

  lceNdelta   = forceVelocityNormFiberHalfLength;%...;0.20;
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
  if(ltN < 1)
      ltN=1;
  end
  if(flag_useElasticTendon==0)
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
  tmaxLongest = abs(lpLong-lpShort)/abs(dlceATSlowest) + 2;

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
        tmax = 2;
        tRampMid = tmax/2;
      end


      
      benchConfig.pathFcn = pathFcn;

      %Update benchConfig       
      benchConfig.npts                  = round(100*tmaxLongest); 
      benchConfig.initialActivation     = 1.0;
      benchConfig.excitationFcn         = excitationMaxFcn;
      benchConfig.activationFcn         = activationFcn; 
      benchConfig.tspan                 = [0, tmax]; 

      %Initialize the model
      modelConfig.initializeState = 1;
      activationState0            = [0;benchConfig.excitationFcn(0)];
      pathState0                  = benchConfig.pathFcn(0);
      muscleState0                = zeros(nStates,1);

      mtInfo = calcMillard2019MuscleInfoOpus31(activationState0,...
                                              pathState0,...
                                              muscleState0,...
                                              musculotendonProperties,...
                                              sarcomereProperties,...
                                              normMuscleCurves,...
                                              modelConfig);

      muscleState0                     = mtInfo.state.value;

      modelConfig.initializeState      = 0;   
      benchConfig.initialState         = muscleState0;

      benchConfig.minimumActivation    = 0;
      benchConfig.name                 = 'Opus31';


      lceThreshold = lceNCenter*lceOpt;
      
      idxThresholdState   = 1;
      eventThresholdValue = lceThreshold; 
      eventDirection      = velDir;

      %Fiber length is not a state for a rigid tendon model. Since the
      %middle of the ramp, by construction, will give us the desired 
      %fiber length we can trigger the event to occur at the middle of the
      %ramp.
      if(flag_useElasticTendon == 0 || abs(normVelocitySeries(z,1)) <= sqrt(eps)) %
        idxThresholdState   = -1;
        eventThresholdValue = tRampMid;
        eventDirection      = 1;
      end
                  
      idxThresholdState = idxThresholdState + 1; %activation state is pre-pended            
      
      benchConfig.eventFcn = @(argT,argY)eventStateThreshold(...
                                argT,argY,idxThresholdState,...
                                eventThresholdValue,eventDirection);

      calcMillard2019MuscleInfoOpus31Fcn = ...
           @(activationState1,pathState2,mclState3) ...
           calcMillard2019MuscleInfoOpus31(   activationState1,...
                                              pathState2,...
                                              mclState3,...
                                              musculotendonProperties,...
                                              sarcomereProperties,...
                                              normMuscleCurves,...
                                              modelConfig);

      flag_appendEnergetics=0;           
      benchRecord = runPrescribedLengthActivationSimulation(...
                       calcMillard2019MuscleInfoOpus31Fcn,...
                       [],...
                       benchConfig,...
                       benchRecord,...
                       idx, length(normVelocitySeries)*2,...
                       flag_appendEnergetics,...
                       flag_useOctave);

      fprintf('%i / %i\n', idx, length(normVelocitySeries)*2);   
    end
  end

  save([outputFolder,'forceVelocityOpus31_',nameModification,outputFileEndingOpus31,'.mat'],...
    'benchRecord',...
    'flag_useElasticTendon',...
    'musculotendonProperties',...
    'sarcomereProperties',...
    'normMuscleCurves');


end


success = 1;
                        
                        