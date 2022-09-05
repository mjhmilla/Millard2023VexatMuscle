function [muscleArchitectureUpd,sarcomerePropertiesUpd] = ...
  updateOpus31CrossBridgeParameters(nominalForceN,...
                                    nominalNormFiberLengthAtSlack,...
                                    flag_figureNumberToFitTo,...
                                    flag_frequencyToFitTo,...
                                    dataKBR1994Fig3Gain,...
                                    dataKBR1994Fig3Phase,...
                                    dataKBR1994Fig12K,...
                                    dataKBR1994Fig12D,...
                                    normTendonDampingConstant,...
                                    normTendonDampingLinear,...
                                    scaleSlidingTimeConstant,...
                                    scaleCrossBridgeCyclingDamping,...
                                    flag_useElasticTendon,...
                                    muscleArchitecture,...
                                    sarcomereProperties,...
                                    normMuscleCurves,...
                                    flag_usingOctave)
                                  
muscleArchitectureUpd  = muscleArchitecture;
sarcomerePropertiesUpd = sarcomereProperties;

%%
% Simple updates
%%
sarcomerePropertiesUpd.slidingTimeConstant ...
  = sarcomerePropertiesUpd.slidingTimeConstant*scaleSlidingTimeConstant;


sarcomerePropertiesUpd.normCrossBridgeCyclingDamping  ...
  = sarcomereProperties.normCrossBridgeCyclingDamping*scaleCrossBridgeCyclingDamping;     

muscleArchitectureUpd.normTendonDampingLinear   = normTendonDampingLinear;
muscleArchitectureUpd.normTendonDampingConstant = normTendonDampingConstant;

%%
% Evaluate the stiffness and damping of the data at the desired nominal
% force
%%
kmt = 0;
dmt = 0;

if(flag_figureNumberToFitTo==3)
  %Fit stiffness and damping values to the gain and phase curves
  
  k0 = 3000;
  d0 = 10;
  idx90 = 1;
  idx15 = 2;
  idxFig3Freq=0;
  freqHz = [];
  
  if(flag_frequencyToFitTo==90)
    idxFig3Freq = idx90;
    freqHz = [4:1:90]';
  elseif(flag_frequencyToFitTo==15)
    idxFig3Freq = idx15;
    freqHz = [4:1:20]';
  else 
    assert(0,'Error: flag_frequencyToFitTo must be 15 or 90');
  end

  
  gain = zeros(size(freqHz));
  phase= zeros(size(freqHz));
  
  
  %Interpolate the gain and phase data to get 1 set.

  for i=1:1:length(freqHz)
    gain(i,1) = interp1(dataKBR1994Fig3Gain(idxFig3Freq).x,...
                          dataKBR1994Fig3Gain(idxFig3Freq).y,...
                          freqHz(i,1),'linear');
    phase(i,1) = interp1(dataKBR1994Fig3Phase(idxFig3Freq).x,...
                          dataKBR1994Fig3Phase(idxFig3Freq).y,...
                          freqHz(i,1),'linear');                      
  end
  
  freqRads = freqHz .* (2*pi);
  gainNm   = gain .* 1000;
  phaseRads= phase .* (pi/180);

  %If fmincon is available, then constrain damping to be positive
  A0 = [0 -1];
  b0 = [ 0 ];

  for z=1:1:2
    
      x0 = [1,1];
      argScaling =[k0,d0];
      argScalingMin=1;
      argScaling(argScaling<argScalingMin) = argScalingMin;
      objScaling = 1;

      err0 = calcFrequencyDomainSquaredError(x0, ...
                                            freqRads,...
                                            gainNm,...
                                            phaseRads,...
                                            argScaling, ...
                                            objScaling);

      objScaling  = 1/max(abs(err0),sqrt(eps));

      errFcn0 = @(argX)calcFrequencyDomainSquaredError(argX, ...
                          freqRads,...
                          gainNm,...
                          phaseRads,...
                          argScaling,...
                          objScaling);

      paramOpt  = []; 
      fval      = [];
      exitFlag  = 0;
      options   = [];
      if(flag_usingOctave == 0)               
        options=optimoptions('fmincon','Display','none','Algorithm','sqp');
        [paramOpt, fval, exitFlag] = fmincon(errFcn0, x0,A0,b0,[],[],[],[],[],options);        
      else
        options=optimset('Display','none','Algorithm','sqp');%optimset('Display','none');
        [paramOpt, fval, exitFlag] = fminsearch(errFcn0, x0,options);                        
      end
      k0 = paramOpt(1)*argScaling(1);
      d0 = paramOpt(2)*argScaling(2);
          
  end 
  
  kmt = k0;
  dmt = d0;
  
  %dataKBR1994Fig3Gain
 %dataKBR1994Fig3Phase;
end

if(flag_figureNumberToFitTo==12)

  dataF = dataKBR1994Fig12K(1:1:end).x;
  dataK = dataKBR1994Fig12K(1:1:end).y;
  fitK = fit(dataF,dataK,'poly1');

  %Ignore the small constant offset
  kmt = (fitK.p1*nominalForceN + fitK.p2)*1000;


  dataF = dataKBR1994Fig12D(1:1:end).x;
  dataD = dataKBR1994Fig12D(1:1:end).y;
  fitD = fit(dataF,dataD,'poly1');

  %Ignore the small constant offset
  dmt = (fitD.p1*nominalForceN + fitD.p2)*1000;
end

if(flag_figureNumberToFitTo ~= 3 && flag_figureNumberToFitTo ~= 12)
    assert(0,'Error: flag_figureNumberToFitTo must be 3 or 12');
end

%%
%Given :
% data stiffness and damping data and the tendon model evaluate
%
%%

ltSlk   = muscleArchitecture.tendonSlackLength;
fiso    = muscleArchitecture.fiso;
lceOpt  = muscleArchitecture.optimalFiberLength;
alphaOpt= muscleArchitecture.pennationAngle;
vmax    = muscleArchitecture.maximumNormalizedFiberVelocity;
pennationAngleAtMinimumFiberLength = ...
          muscleArchitecture.pennationAngleAtMinimumFiberLength;

               
%%
% Solve for the length of the path such that the fiber reaches
% the desired length with an activation of 0.
%
% Note: this is identical to the code used to compute the baseline length
%       in the 'runPerturbationSimulationsOpus31' and 
%              'runPerturbationSimulationsHill' code
%%

x0     = 0;
xdot0  = 0;

%Evaluate the fiber kinematics
alphaOpt = muscleArchitecture.pennationAngle;
lceOpt   = muscleArchitecture.optimalFiberLength;

lce0 = nominalNormFiberLengthAtSlack*lceOpt + x0;
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


ft0 = fpeN/cos(alpha);
tendonForceLengthCurveInverse = ...
    createInverseCurve(normMuscleCurves.tendonForceLengthCurve);
ltN    = calcBezierYFcnXDerivative(ft0, tendonForceLengthCurveInverse,0);
if(flag_useElasticTendon==0)
  ltN = 1;
end            
ltslk  = muscleArchitecture.tendonSlackLength;

lp0  = lceAT+ltN*ltslk;
dlp0 = 0;

%Now evalute the lengths when the fiber is activated

ft1    = nominalForceN/fiso;
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
%Evaluate the fiber length when its developing the desired and estimate
%the activation.
%%
falN  = calcBezierYFcnXDerivative(lceN1, ...
          normMuscleCurves.activeForceLengthCurve ,0);
fvN   = calcBezierYFcnXDerivative(dlceN1, ...
          normMuscleCurves.fiberForceVelocityCurve,0);    
fpeN  = calcBezierYFcnXDerivative(lceN1,  ...
          normMuscleCurves.fiberForceLengthCurve  ,0);

fcnN   = ft1/cos(alpha);
activation = max(0,(fcnN-fpeN)/(falN*fvN));

if( fpeN > fcnN )
  fprintf('  Warning: \tdesired force (%1.3f) > passive force (%1.3f) at length (%1.3f)\n',...
    fcnN, fpeN,lceN1);
end



%%
%Iteratively vary activation until the desired tendon force is reached
%%
errLeft=0;
errRight=0;
errBest = 1e10;

actBest = activation;
actDelta = 0.1;
actLeft = 0;
actRight= 0;

modelConfig = struct( ...
  'iterMax'                 , 100                   , ...
  'tol'                     , 1e-12                 , ... 
  'tolInit'                 , sqrt(eps)*1000        , ...
  'minActivation'           , 0.0                   , ...
  'useElasticTendon'        , flag_useElasticTendon , ...
  'initializeState'         , 1                     );  

nStates = 0;
if(flag_useElasticTendon==1)
  nStates = 4;
else
  nStates = 3;
end

activationState0 = [0;actBest];
pathState0       = [dlp0, lp0];
muscleState0     = zeros(nStates,1);


mtInfoBest = [];

mtInfoBest =calcMillard2019MuscleInfoOpus31([0,actBest],...
                                      pathState0,...
                                      muscleState0,...
                                      muscleArchitecture,...
                                      sarcomereProperties,...
                                      normMuscleCurves,...
                                      modelConfig);
muscleState0 =mtInfoBest.state.value;
modelConfig.initializeState =0;  
mtInfoBest =calcMillard2019MuscleInfoOpus31([0,actBest],...
                                      pathState0,...
                                      muscleState0,...
                                      muscleArchitecture,...
                                      sarcomereProperties,...
                                      normMuscleCurves,...
                                      modelConfig);              

errBest = abs(mtInfoBest.muscleDynamicsInfo.tendonForce-nominalForceN);

for i=1:1:10

  actLeft = max(0,actBest-actDelta);
  modelConfig.initializeState =1; 
  mtInfoL =calcMillard2019MuscleInfoOpus31([0,actLeft],...
                                        pathState0,...
                                        muscleState0,...
                                        muscleArchitecture,...
                                        sarcomereProperties,...
                                        normMuscleCurves,...
                                        modelConfig);
  muscleState0 =mtInfoL.state.value;
  modelConfig.initializeState =0;  
  mtInfoL =calcMillard2019MuscleInfoOpus31([0,actLeft],...
                                        pathState0,...
                                        muscleState0,...
                                        muscleArchitecture,...
                                        sarcomereProperties,...
                                        normMuscleCurves,...
                                        modelConfig);                                      
  errL = abs(mtInfoL.muscleDynamicsInfo.tendonForce-nominalForceN);
  
  actRight = min(1,actBest+actDelta);
  modelConfig.initializeState =1; 
  mtInfoR =calcMillard2019MuscleInfoOpus31([0,actRight],...
                                        pathState0,...
                                        muscleState0,...
                                        muscleArchitecture,...
                                        sarcomereProperties,...
                                        normMuscleCurves,...
                                        modelConfig);
  muscleState0 =mtInfoR.state.value;
  modelConfig.initializeState =0;  
  mtInfoR =calcMillard2019MuscleInfoOpus31([0,actRight],...
                                        pathState0,...
                                        muscleState0,...
                                        muscleArchitecture,...
                                        sarcomereProperties,...
                                        normMuscleCurves,...
                                        modelConfig);                                        
  errR = abs(mtInfoR.muscleDynamicsInfo.tendonForce-nominalForceN);

  if(errL < errR && errL < errBest)
    errBest=errL;
    actBest = actLeft;
    mtInfoBest = mtInfoL;
  end
  if(errR < errL && errR < errBest)
    errBest = errR;
    actBest = actRight;
    mtInfoBest = mtInfoR;    
  end
  
  actDelta = actDelta*0.5;
  
end
  
  

%%
% Using the solution in mtInfo extract the stiffness and damping of the
% relevant components and solve for the new cross-bridge and damping 
% parameters
%%

kT    = mtInfoBest.muscleDynamicsInfo.tendonStiffness;
dT    = mtInfoBest.muscleDynamicsInfo.tendonDamping;
ke    = mtInfoBest.muscleDynamicsInfo.ecmStiffness;
de    = mtInfoBest.muscleDynamicsInfo.ecmDamping;
k1  = mtInfoBest.muscleDynamicsInfo.titin1Stiffness;
k2  = mtInfoBest.muscleDynamicsInfo.titin2Stiffness;
d1  = mtInfoBest.muscleDynamicsInfo.titin1Damping;
d2  = mtInfoBest.muscleDynamicsInfo.titin2Damping;

%The stiffness contribution of titin depends on the position of the
%N2A element and how activated the muscle is. In one extreme the muscle is 
%activated and thus the N2A element is effectively locked. Then titin's 
%stiffness contribution would be PEVK. In the other extreme the N2A element
%is free to move and titin's stiffness contribution is the stiffness of the
%igprox and pevk regions in series.
%
%Here we are assuming ignorance and
%saying that titin's contribution is about half of what it could be.

  
a = mtInfoBest.muscleDynamicsInfo.activation;
falN = mtInfoBest.muscleLengthInfo.fiberActiveForceLengthMultiplier;

%ktitinHigh = k2; 
%ktitinLow  = 1/((1/k1)+(1/k2));
%ktitin     = (1-a)*ktitinLow + a*ktitinHigh;
ktitin = k2;

%PEVK section has no damping
dtitin = 0;%1/((1/d1)+(1/d2));

assert(kT > kmt || flag_useElasticTendon == 0,...
    'The tendon must be stiffer than the target MT stiffness');
assert(dT > dmt || flag_useElasticTendon == 0,...
    'The tendon must be more damped than the target MT damping');

%Target fiber stiffness and damping along the tendon;
kfAT = 0;
dfAT = 0;
if(flag_useElasticTendon==1)
  kfAT = ( (1/kmt) - (1/kT) )^-1;
  dfAT = ( (1/dmt) - (1/dT) )^-1;
else
  kfAT = kmt;
  dfAT = dmt;
end


alpha = mtInfoBest.muscleLengthInfo.pennationAngle;
lce   = mtInfoBest.muscleLengthInfo.fiberLength;

fibKinDer = calcFixedWidthPennationPartialDerivatives(alpha,...
                                                      lce,...
                                                      lceOpt,...
                                                      alphaOpt);
Dalpha_Dlce   = fibKinDer.Dalpha_Dlce;

fT = mtInfoL.muscleDynamicsInfo.tendonForce;

% kfAT = d/dlce ( ff*cos(alpha) )(dlce/dlceAT)
%      = (kf*cos(alpha) - ff*sin(alpha)*Dalpha_Dlce)(1/cos(alpha))
% kf   = (kfAT cos(alpha) + ff*sin(alpha)*Dalpha_Dlce)



ff = fT/cos(alpha);
Dlce_DlceAT = 1/cos(alpha);

kf =  (kfAT/Dlce_DlceAT + ff*sin(alpha)*Dalpha_Dlce)/cos(alpha);

df = dfAT*cos(alpha)/(cos(alpha) - lce*sin(alpha)*Dalpha_Dlce);

%kf_dep = (kfAT+fT*sin(alpha)*Dalpha_Dlce)/cos(alpha);
%assert(0,'Check this');

%df_dep = dfAT/cos(alpha);  
%assert(0,'Re-derive this');


assert( kf > (ke+ktitin), ...
  ['Something is wrong: stiffness of ECM and Titin',...
   ' exceed desired active stiffness']);
assert( df > (de+dtitin), ...
  [ 'Something is wrong: damping ECM and Titin',...
    ' exceed desired active damping']);

  
kx = kf-(ktitin+ke);
dx = df-(dtitin+de);

kxMax = kx/(a*falN);
dxMax = dx/(a*falN);

kxMaxH = kxMax*2;
dxMaxH = dxMax*2;

kxMaxHN = kxMaxH/(fiso/lceOpt);
dxMaxHN = dxMaxH/(fiso/lceOpt);

sarcomerePropertiesUpd.normCrossBridgeDamping   = dxMaxHN;
sarcomerePropertiesUpd.normCrossBridgeStiffness = kxMaxHN;


                             