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
function mtInfo = calcMillard2023VexatMuscleInfo(...
                                            activationState,...                                            
                                            pathState, ... 
                                            muscleState,...                                                                                       
                                            muscleArchitecture,...
                                            sarcomereProperties,...
                                            normMuscleCurves,...
                                            modelConfig)
%%              
%                                 
%                                                                 *pennated here
%         la             lx             lm (fixed)
% |---------------->|<-------->|<---------------------->|        
% |              
% |================[|]===============|                  
%                   |     kx                
%                   |--|/\/\|--|                        |
%                   |     bx   |========================|
%                   |--[    ]--|                        |        kt
%                              |                        |   |--|/\/\|--|
% |     IgP       PEVK    IgD  |                        |---|    bt    |---
% |---|\/\/\|--|\/\/\/|--------|                        |   |--[    ]--|
% |===========[|]====================|                  |
% |---- l1 --->|                                        |
% |---------------------|\/\/\/\/\/\|-------------------| fecm
% |---------------------[           ]-------------------| becm
%
% |-------------------------1/2 lce --------------------|-------lt---------|  
%
%
% State vector : x = [lce,dla,la,l1], Stretch goal: e1]
%dlce    : fiber velocity                            (units: m/s)
% lce    : fiber length                              (units: m)
% la     : sliding length                            (units: m)
% lx     : cross bridge stretch                      (units: m)
% l1     : length of the titin segment between       (units: m)
%          the Z-line and the PEVK/distal IG border
% Stretch goal
% e1     : calcium + stretch induced enhancement     (integral of act*m/s)
%
% Ascii Variable Notation convention
%
% ce : fiber 
% t  : tendon 
% 1  : titin segment from the z-line to the PEVK/IG2 border
% 2  : titin segment from the PEVK/IG2 border to the m-line
% s  : actin myosin attachement
% x  : cross bridge
% a  : actin
% m  : m-line
% z  : z-line
% f  : length of the fiber excluding the cross bridge strain
% N  : normalized
% H  : half
%
% 
% lceH  = lce*0.5
% lceNH = lceH/lopt
% dlceH = dlce*0.5
% dlceNH= dlceH/(lopt*vmax)
%

useElasticTendon = modelConfig.useElasticTendon;
%assert(useElasticTendon==1,'Rigid tendon not yet coded');


if(useElasticTendon==1)
  assert(length(muscleState) == 4,'This is a 4 state muscle model!');
else
  assert(length(muscleState) == 3,'This is a 3 state muscle model!');
end
%%
%Create function handles for the curves to make the code easier to read                                                  
%%

%Active force length curve

     
%Force velocity curve
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

calcFeHDer  = @(arg1,arg2)calcBezierYFcnXDerivative(arg1, ...
                  normMuscleCurves.forceLengthECMHalfCurve, ...
                  arg2);



calcF1HDer       = @(arg1, arg2)calcBezierYFcnXDerivative(arg1, ...
                      normMuscleCurves.forceLengthProximalTitinCurve, ...
                      arg2);  
calcF2HDer       = @(arg1, arg2)calcBezierYFcnXDerivative(arg1, ...
                      normMuscleCurves.forceLengthDistalTitinCurve, ...
                      arg2); 

%Just to check to make sure the curves are consistent
if(normMuscleCurves.useTwoSidedTitinCurves==1)
    assert(normMuscleCurves.forceLengthProximalTitinCurve.xEnd(1,1) < 0);
else
    assert(normMuscleCurves.forceLengthProximalTitinCurve.xEnd(1,1) >= 0);
end


lceHNZeroFpeN = normMuscleCurves.fiberForceLengthCurve.xEnd(1,1)*0.5;
lceHNTitinActinBondMin= sarcomereProperties.normLengthTitinActinBondMinimum*0.5;

%Extract the slack length
l1HNZeroForce = 0;

if(normMuscleCurves.forceLengthProximalTitinCurve.ypts(1,1) > 0)
    %Single sided curve
    l1HNZeroForce = normMuscleCurves.forceLengthProximalTitinCurve.xpts(1,1);
else
    %Double sided curve
    for i=1:1:size(normMuscleCurves.forceLengthProximalTitinCurve.ypts,2)
        if((min(normMuscleCurves.forceLengthProximalTitinCurve.ypts(:,i)) <= 0) ...
         && (max(normMuscleCurves.forceLengthProximalTitinCurve.ypts(:,i)) > 0))
            if(normMuscleCurves.forceLengthProximalTitinCurve.ypts(1,i) > 0)
                l1HNZeroForce = normMuscleCurves.forceLengthProximalTitinCurve.xpts(1,i);
            else
                l1HNZeroForce = normMuscleCurves.forceLengthProximalTitinCurve.xpts(end,i);
            end
        end
    end
end
if(abs(l1HNZeroForce) < 0.01)
    here=1;
end

lTitinFixedHN = sarcomereProperties.ZLineToT12NormLengthAtOptimalFiberLength ...
              + sarcomereProperties.IGDFixedNormLengthAtOptimalFiberLength;



%Tendon force-length curve            
calcFtDer       = @(arg1, arg2)calcBezierYFcnXDerivative(arg1, ...
                  normMuscleCurves.tendonForceLengthCurve, ...
                  arg2); 
 
%Tendon damping-length curve            
calcDtDer       = @(arg1, arg2)calcBezierYFcnXDerivative(arg1, ...
                  normMuscleCurves.tendonStiffnessCurve, ...
                  arg2);                 
                
%Compressive force-lenght curve                
calcCpDer = @(arg1, arg2)calcBezierYFcnXDerivative(arg1, ...
                  normMuscleCurves.compressiveForceLengthCurve, ...
                  arg2); 

DftfcnN_DltN_max = normMuscleCurves.tendonForceLengthCurve.dydxEnd(1,2);

calcFpeDer = @(arg1, arg2)calcBezierYFcnXDerivative(arg1, ...
                  normMuscleCurves.fiberForceLengthCurve, ...
                  arg2);




%%
%Muscle Architectural Properties
%%

fiso      =   muscleArchitecture.fiso;
lceOpt    =   muscleArchitecture.optimalFiberLength;
alphaOpt  =   muscleArchitecture.pennationAngle;
ltSlk     =   muscleArchitecture.tendonSlackLength;
dlceMaxN  =   muscleArchitecture.maximumNormalizedFiberVelocity;
lceMin    =   muscleArchitecture.minimumFiberLength;
lceATMin  =   muscleArchitecture.minimumFiberLengthAlongTendon;
betaNum   =   muscleArchitecture.normNumericalDamping;
betafTN   =   muscleArchitecture.normTendonDampingLinear;
betaCTN   =   muscleArchitecture.normTendonDampingConstant;



%Scaled to be appropriate for a half sarcomere
scaleECM   = sarcomereProperties.scaleECM;
scaleTitinProximal  = sarcomereProperties.scaleTitinProximal;
scaleTitinDistal    = sarcomereProperties.scaleTitinDistal;

kAXHN      =   sarcomereProperties.normCrossBridgeStiffness;
betaAXHN   =   sarcomereProperties.normCrossBridgeDamping;

forceVelocityCalibrationFactor =   sarcomereProperties.forceVelocityCalibrationFactor;

betafEcmHN =   sarcomereProperties.normECMDamping;

betaTAaHN = sarcomereProperties.normMaxActiveTitinToActinDamping;
betaTApHN = sarcomereProperties.normPassiveTitinToActinDamping;

betaPevkPHN = sarcomereProperties.normPassivePevkDamping;
betaPevkAHN = sarcomereProperties.normActivePevkDamping;
activationThresholdTitin = ...
  sarcomereProperties.activationThresholdTitin;

slidingTimeConstant = sarcomereProperties.slidingTimeConstant;

normCrossBridgeCyclingDamping = sarcomereProperties.normCrossBridgeCyclingDamping;
lowActivationThreshold        = sarcomereProperties.lowActivationThreshold;
lowActivationGain             = sarcomereProperties.lowActivationGain;

normActinLength               = sarcomereProperties.normActinLength;

%%
%Sarcomere Architectural Properties
%%

%laN     = sarcomereProperties.normActinLength;  
lmHN    = sarcomereProperties.normMyosinHalfLength;



%la  = laN*lceOpt;
lmH = lmHN*lceOpt;


%Path state
dlp = pathState(1);
lp  = pathState(2);

%Activation state
%dadt  = activationState(1); this model does not require an activation state.

lce    = 0.;
dlaH   = 0.;  
laH    = 0.;  
l1H    = 0.;
lambda = 0.;

% I've removed the calcium-based force enhancement states for the two 
% titin segments from the model for now: this makes at most 5% difference
% in the stiffness of these two segments, and I am running out of time.
% This will have to wait for the next model.

if(useElasticTendon == 1)
  lce    = muscleState(1);
  dlaH   = muscleState(2);  
  laH    = muscleState(3);  
  l1H    = muscleState(4);
  lambda = 0;%muscleState(5);
else
  dlaH   = muscleState(1);  
  laH    = muscleState(2);  
  l1H    = muscleState(3);
  lambda = 0;%muscleState(4);
end  

%Future: update state to be the fiber length along the tendon.
% if(lce < muscleArchitecture.minimumFiberLength)
%   lce=muscleArchitecture.minimumFiberLength;
% end

                                            
dadt   = activationState(1);
a      = activationState(2);

%assert(a >= 0);


if( (lce >= muscleArchitecture.minimumFiberLength ...
           || modelConfig.initializeState == 1 ) == 0 )
  here=1;
end

if( ~(lce >= muscleArchitecture.minimumFiberLength ...
              || modelConfig.initializeState == 1 ...
              || useElasticTendon == 0) )
here=1;

end
%  assert( (lce >= muscleArchitecture.minimumFiberLength ...
%                || modelConfig.initializeState == 1 ...
%                || useElasticTendon == 0));

if(sum(isnan(muscleState))>0)
   disp('muscleState has a NaN entry'); 
end


sqrtEps = sqrt(eps);



%Notation:  In_Out
%   Dli_DliN: li in, liN out.
%   This is backwards in the sense of the chain rule but it
%   is a bit easier to read since the correct ordering will
%   yield neighbors that are similar
%
%   lceN = lce * lce_lceH;
%           |      |
%           similar neighbors, correct
%
%   lceN = lce * lceH_lce;
%           |      |
%           dissimilar similar neighbors, incorrect.
%

lce_lceH      = 0.5;
lce_lceHN     = 0.5/lceOpt;
dlce_dlceN    = 1.0/lceOpt;
dlce_dlceNN  = 1/(lceOpt*dlceMaxN);
                     
lce_lceN    = 1/lceOpt;
lt_ltN      = 1/ltSlk;
lk_lkN      = 1/lceOpt;
li_liN      = 1/lceOpt; % Here 'i' denotes titin.
li_liHN     = lce_lceHN;           

lceH_lce     = 1/lce_lceH    ; 
lceHN_lce    = 1/lce_lceHN   ; 
dlceN_dlce   = 1/dlce_dlceN  ; 
dlceNN_dlce  = 1/dlce_dlceNN ; 

lceN_lce    = 1/lce_lceN;
ltN_lt      = 1/lt_ltN  ;
lkN_lk      = 1/lk_lkN  ;
liN_li      = 1/li_liN  ;
liHN_li     = lceHN_lce;
                   
%%
%Set up utility structs
%%
modelCurves = struct( ...
        'calcFalDer' , calcFalDer ,...
        'calcFvDer'  , calcFvDer   ,...  
        'calcFeHDer' , calcFeHDer  ,...        
        'calcF1HDer' , calcF1HDer  ,...
        'calcF2HDer' , calcF2HDer  ,...
        'calcFtDer'  , calcFtDer   ,...
        'calcDtDer'  , calcDtDer   ,...
        'calcCpDer'  , calcCpDer   ,...
        'calcFpeDer' , calcFpeDer);




 
      
modelConstants = struct( ...
    'fiso'          , fiso             , ...
    'lceOpt'        , lceOpt           , ...    
    'alphaOpt'      , alphaOpt         , ...
    'ltSlk'         , ltSlk            , ...
    'dlceMaxN'      , dlceMaxN         , ...
    'lceMin'        , lceMin           , ...
    'lceHNZeroFpeN' , lceHNZeroFpeN    , ...
    'lceHNTitinActinBondMin', lceHNTitinActinBondMin,...
    'kAXHN'         , kAXHN            , ...
    'betaAXHN'      , betaAXHN         , ...
    'normCrossBridgeCyclingDamping'  , normCrossBridgeCyclingDamping, ...
    'lowActivationThreshold', lowActivationThreshold,...
    'lowActivationGain',lowActivationGain,...
      'betaNum'     , betaNum          , ...                          
      'betafTN'      , betafTN           , ...    
      'betaCTN'      , betaCTN           , ...          
    'betafEcmHN'      , betafEcmHN         , ...                    
    'betaTAaHN'    , betaTAaHN       , ...
    'betaTApHN'    , betaTApHN       , ...  
    'betaPevkPHN'   , betaPevkPHN      , ...
    'betaPevkAHN'   , betaPevkAHN      , ...  
    'activationThresholdTitin', activationThresholdTitin,...
    'forceVelocityCalibrationFactor'  , forceVelocityCalibrationFactor,...
    'tau'       , slidingTimeConstant, ...
    'lTitinFixedHN' , lTitinFixedHN    , ...
    'ZLineToT12NormLengthAtOptimalFiberLength',...
      sarcomereProperties.ZLineToT12NormLengthAtOptimalFiberLength,...
    'lmHN'          , lmHN             , ...
    'lmH'           , lmH              , ...    
    'lAHN',normActinLength,...
      'lce_lceH'    , lce_lceH         , ... 
      'lceH_lce'    , lceH_lce         , ...       
      'lce_lceHN'   , lce_lceHN        , ...    
      'dlce_dlceN'  , dlce_dlceN       , ...                            
      'dlce_dlceNN',  dlce_dlceNN      , ...                                             
      'lce_lceN'    , lce_lceN         , ...
      'lceN_lce'    , lceN_lce         , ...
    'lt_ltN'        , lt_ltN             , ...                                
    'lk_lkN'        , lk_lkN             , ...                                
    'li_liN'        , li_liN             , ...                                             
    'ltN_lt'        , ltN_lt             , ...                                
    'lkN_lk'        , lkN_lk             , ...                                
    'liN_li'        , liN_li             , ...
    'scaleECM',  scaleECM,...
    'scaleTitinProximal', scaleTitinProximal,...
    'scaleTitinDistal', scaleTitinDistal,...
    'DftfcnN_DltN_max',DftfcnN_DltN_max,...
    'smoothStepFunctionRadius',...
      sarcomereProperties.smoothStepFunctionRadius,...
    'titinModelType', sarcomereProperties.titinModelType,...
    'initialization', modelConfig.initializeState,...
    'iterMax',modelConfig.iterMax);



modelCachedValues = struct(...
    'dlt'            , NaN       , ...
    'lt'             , NaN       , ...
    'dltN'           , NaN       , ...    
    'ltN'            , NaN       , ...
    'et'             , NaN       , ...
      'lp'             , lp        , ...
      'dlp'            , dlp       , ...
      'dw_lp'          , NaN      , ...
    'lce'            , lce       , ...
    'lceAT'          , NaN       , ...
    'lceN'           , NaN       , ...
    'lceH'           , NaN       , ...      
    'lceHN'          , NaN       , ...
    'dlce'           , NaN       , ...
    'dlceH'          , NaN       , ...
    'dlceHN'         , NaN       , ...    
    'dlceAT'         , NaN       , ...    
      'laH'            , laH       , ...
      'dlaH'           , dlaH      , ...
      'ddlaH'          , NaN       , ...
      'laHN'           , NaN        , ...
      'dlaHN'          , NaN        , ...
      'ddlaHN'         , NaN        , ...    
    'l1H'           , l1H        , ...
    'l1HN'          , NaN        , ...
    'dl1H'          , NaN        , ...
    'dl1HN'         , NaN        , ...
    'dw_n2aH'        , NaN        , ...
      'l2H'           , NaN      , ...
      'l2HN'          , NaN      , ...
      'dl2H'          , NaN      , ...
      'dl2HN'         , NaN      , ...
    'lxH'           , NaN        , ...
    'lxHN'          , NaN        , ...        
    'dlxH'          , NaN        , ...
    'dlxHN'         , NaN        , ...    
    'dw_xkH'         , NaN        , ...
    'dw_xdH'         , NaN        , ...
    'dw_xaH'         , NaN        , ...
      'alpha'         , NaN       , ...
      'cosAlpha'      , NaN       , ...
      'sinAlpha'      , NaN       , ...
      'dalpha'        , NaN       , ...
    'a'             , a       , ...
    'dadt'          , dadt    , ...
    'lfN'           , NaN     , ...
    'dlfNN'         , NaN     , ...
    'flN'           , NaN     , ...    
    'fvN'           , NaN     , ...
      'f1H'           , NaN       , ...
      'f1HN'          , NaN       , ...
      'f1kHN'         , NaN       , ...
      'f1dHN'         , NaN       , ...      
      'k1HNN'         , NaN       , ...
      'beta1HNN'      , NaN       , ...
      'dw_1kH'         , NaN       , ...
    'f2H'           , NaN       , ...
    'f2HN'          , NaN       , ...
    'f2kHN'         , NaN       , ...
    'f2dHN'         , NaN       , ...      
    'k2HNN'         , NaN       , ...
    'beta2HNN'      , NaN       , ...
    'dw_2kH'         , NaN       , ...
      'fxHN'        , NaN       , ... 
      'kxHNN'           , NaN       , ...
      'betaxHNN'        , NaN       , ...    
    'fceN'            , NaN       , ...   
    'fTfcnN'          , NaN       , ...  
    'DfTfcnN_DltN'    , NaN       , ...    
    'fTkN'            , NaN       , ...
    'fTdN'            , NaN       , ...
    'fTN'             , NaN       , ...
    'kTNN'            , NaN       , ...  
      'fCpN'          , NaN       , ...
      'kCpN'          , NaN       , ...
    'betaTNN'         , NaN       , ...
    'dw_Tk'           , NaN       , ...
    'dw_Td'           , NaN       , ...     
      'fEcmfcnHN'         , NaN       , ...
      'DfEcmfcnHN_DlceHN' , NaN       , ...      
      'fEcmkHN'           , NaN       , ...
      'fEcmdHN'           , NaN       , ...
      'fEcmHN'            , NaN       , ...
      'kEcmHNN'           , NaN       , ...
      'betaEcmHNN'        , NaN       , ...
      'dw_EcmkH'           , NaN       , ...
      'dw_EcmdH'           , NaN       , ... 
      'dw_Cp'             , NaN       , ...
    'Dalpha_Dlce'     , NaN       , ...
    'Ddalpha_Ddlce'   , NaN       , ...
    'Dlt_Dlce'        , NaN       , ...
    'DltN_DlceN'      , NaN       , ...
    'D_f2HN_D_dlce'       , NaN, ...   
    'D_fxHN_D_dlce'       , NaN, ...
    'D_fEcmHN_D_dlce'     , NaN, ...
    'D_fTN_D_dlce'        , NaN, ...
      'DfEcmkHN_DlceHN'   , NaN, ...
      'DfEcmdHN_DlceHN'   , NaN, ...
      'DfEcmHN_DlceHN'    , NaN, ...
      'DfCpN_DlceN'       , NaN, ...
    'DfTkN_DltN'      , NaN ,...
    'DfTdN_DltN'      , NaN ,...
    'DfTN_DltN'       , NaN ,...
    'dT', NaN,...
    'dV', NaN,...
    'dW', NaN,...
        'tau',NaN,...
        'lambda',lambda,...
        'dlambda',NaN,...
        'ddlaHN_HillError',NaN,...
        'ddlaHN_Damping',NaN,...
        'ddlaHN_Tracking',NaN);

%Careful: 
% :Quantities below are in units of Joules
% :For the whole muscle
%dw_Tk
%dw_Tbeta
%dw_ECMk
%dw_ECMbeta
%dw_IG1
%dw_PEVK
%dw_N2A
%dw_Xk
%dw_Xbeta
%dw_Xcycle

%dwdt_XBridgeDamping
%dwdt_XBridgeCycling


if(modelConfig.initializeState==1)

  %assert(abs(dlp) < sqrt(eps)); %General case is not yet handled.

  errF    = 0;
  errFJac = 0;
  errI    = zeros(1+useElasticTendon,1);
  errIJac = zeros(1+useElasticTendon,1+useElasticTendon);
  argsBest = zeros(2,1);
  flag_initializeAtRest     = 0;

  if(  flag_initializeAtRest == 0)
      flag_updatePositionLevel                =0;
      flag_evaluateJacobian                   =0; 
      flag_evaluateDerivatives                =0; 
      flag_updateModelCache                   =0;
      flag_evaluateInitializationFunctions    =1;
      flag_useArgs                            =1;
    
    
      vars=zeros(2,1);

      numMaxBisections = modelConfig.iterMax;
    
      %Solve for an initial state that leads to dlx->0 and an acceleration
      %of zero.

      for i=1:1:2
    
          delta=0;
          switch i
              case 1
                  vars(1,1)=lp*0.5;
                  delta=vars(1,1)*0.5;
              case 2
                  x = vars(1,1);
                  y = lceOpt*sin(alphaOpt);
                  vars(2,1)=0.5*sqrt(x*x + y*y);
                  delta=vars(2,1)*0.5;                  
              otherwise
                  assert(0);
          end
          
          flag_evaluateInitializationFunctions=i;
          [errF, errFJac, errI, errIJac, modelCache] ...
                = updateMillard2023VexatCache(...
                                      vars                                     ,... 
                                      modelCachedValues                        ,...
                                      modelCurves                             ,...
                                      modelConstants                          ,...
                                      sarcomereProperties                     ,...
                                      useElasticTendon                        ,...
                                      flag_updatePositionLevel                ,...
                                      flag_evaluateJacobian                   ,... 
                                      flag_evaluateDerivatives                ,... 
                                      flag_updateModelCache                   ,...
                                      flag_evaluateInitializationFunctions    ,...
                                      flag_useArgs);
          switch i 
              case 1
                errBest = abs(errI(1,1));
                if(useElasticTendon==0)
                    vars(i,1)=modelCache.lceAT;
                end
                varsBest=vars;
                argsBest(i,1)=vars(i,1);
              case 2
                errBest = abs(errI(2,1));
                varsBest=vars;                  
                argsBest(i,1)=vars(i,1);
                
              otherwise
                  assert(0);
          end
          

        if(i==1 && useElasticTendon==0)
            numBisections = 0;
        else
            numBisections=numMaxBisections;
        end

        for j=1:1:numBisections
    
          stepSign = 0;
    
          errL=0;
          errR=0;
          varsL=0;
          varsR=0;
          for k=1:1:2
            switch k 
              case 1 
                stepSign=-1;
              case 2
                stepSign= 1;
              otherwise
                assert(0);
            end
            vars = varsBest;
            vars(i,1)=vars(i,1)+stepSign*delta; 
      
            if(vars(1,1) < lceATMin)
              vars(1,1) = lceATMin;
            end
            
            flag_evaluateInitializationFunctions=i;

            [errF, errFJac, errI, errIJac, modelCache] ...
              = updateMillard2023VexatCache(...
                                    vars                                   ,... 
                                    modelCachedValues                       ,...
                                    modelCurves                             ,...
                                    modelConstants                          ,...
                                    sarcomereProperties                     ,...
                                    useElasticTendon                        ,...
                                    flag_updatePositionLevel                ,...
                                    flag_evaluateJacobian                   ,... 
                                    flag_evaluateDerivatives                ,... 
                                    flag_updateModelCache                   ,...
                                    flag_evaluateInitializationFunctions    ,...
                                    flag_useArgs);
            switch k
              case 1
                varsL = vars;
                errL  = abs(errI(i,1));
              case 2
                varsR = vars;
                errR  = abs(errI(i,1));          
              otherwise 
                assert(0)
            end
          end

          if(errL < errBest && errL <= errR )                                       
            errBest   = errL;
            varsBest  = varsL;
            argsBest(i,1)=varsL(i,1);

          end
          if(errR < errBest && errR < errL )                                       
            errBest   = errR;
            varsBest  = varsR;
            argsBest(i,1)=varsR(i,1);
          end
          delta = delta*0.5;
        end 
        here=1;
      end

    %Update the output of the initialization values in the cache
    flag_updatePositionLevel                =0;
    flag_evaluateJacobian                   =0; 
    flag_evaluateDerivatives                =0; 
    flag_updateModelCache                   =0;
    flag_evaluateInitializationFunctions    =1;
    flag_useArgs                            =1;
    errRecord = zeros(useElasticTendon+1,1);
    for i=1:1:2
        flag_evaluateInitializationFunctions=i;
    
        [errF, errFJac, errI, errIJac, modelCache] ...
          = updateMillard2023VexatCache(...
                                argsBest                                ,... 
                                modelCachedValues                       ,...
                                modelCurves                             ,...
                                modelConstants                          ,...
                                sarcomereProperties                     ,...
                                useElasticTendon                        ,...
                                flag_updatePositionLevel                ,...
                                flag_evaluateJacobian                   ,... 
                                flag_evaluateDerivatives                ,... 
                                flag_updateModelCache                   ,...
                                flag_evaluateInitializationFunctions    ,...
                                flag_useArgs);
        errRecord(i,1)=errI(i,1);
        modelCachedValues=modelCache;
    end

    flag_updatePositionLevel                =1;
    flag_evaluateJacobian                   =1; 
    flag_evaluateDerivatives                =1; 
    flag_updateModelCache                   =1;
    flag_evaluateInitializationFunctions    =0;
    flag_useArgs                            =0;


    [errF,errFJac,errI,errIJac,modelCachedValuesUpd] = ...
      updateMillard2023VexatCache( vars,...
                                  modelCachedValues,...
                                  modelCurves,...
                                  modelConstants,...
                                  sarcomereProperties,...
                                  useElasticTendon,...
                                  flag_updatePositionLevel,...
                                  flag_evaluateJacobian, ...
                                  flag_evaluateDerivatives, ...
                                  flag_updateModelCache,...
                                  flag_evaluateInitializationFunctions,...
                                  flag_useArgs); 

    modelCachedValues = modelCachedValuesUpd;  
  end



  if(flag_initializeAtRest == 1)
    lceAT  = max(lp-ltSlk*1.01,lceATMin+lceOpt*0.01);
    fibKin = calcFixedWidthPennatedFiberKinematics(lceAT,0,lceOpt,alphaOpt);
    lce    = fibKin.fiberLength;%,lceMin+lceOpt*0.01);  
    lceN = lce*lce_lceN;
    laHN = 0.5*lceN - lmHN;  
    laH  = laHN*lceN_lce;
    dlaH = 0;

    
    l1HN = l1HNZeroForce;
    l1H  = l1HN * liN_li;
      
    flag_updatePositionLevel                = 1;
    flag_evaluateJacobian                   = 1;
    flag_evaluateDerivatives                = 1;
    flag_updateModelCache                   = 1;
    flag_evaluateInitializationFunctions    = 1;
    flag_useArgs                            = 0;

    vars = [];
    if(useElasticTendon==1)
      vars = 0.;
    end
      
    modelCachedValues.lce  = lce;
    modelCachedValues.laH  = laH;
    modelCachedValues.dlaH = dlaH;
    modelCachedValues.l1H  = l1H;
    modelCachedValues.lambda= 0;

    [errF,errFJac,errI,errIJac,modelCachedValuesUpd] = ...
      updateMillard2023VexatCache( vars,...
                                  modelCachedValues,...
                                  modelCurves,...
                                  modelConstants,...
                                  sarcomereProperties,...
                                  useElasticTendon,...
                                  flag_updatePositionLevel,...
                                  flag_evaluateJacobian, ...
                                  flag_evaluateDerivatives, ...
                                  flag_updateModelCache,...
                                  flag_evaluateInitializationFunctions,...
                                  flag_useArgs);  
    
                                
    iter=1;
    iterMax = 100;
    tolInit = 1e-8;
    
    while(iter < iterMax && errI'*errI > tolInit)

      modelCachedValues.lce  = lce;
      modelCachedValues.laH  = laH;
      modelCachedValues.dlaH = 0;
      modelCachedValues.l1H  = l1H;
      modelCachedValues.lambda= 0;


      [errF,errFJac,errI,errIJac,modelCachedValuesUpd] = ...
        updateMillard2023VexatCache( vars,...
                                    modelCachedValues,...
                                    modelCurves,...
                                    modelConstants,...
                                    sarcomereProperties,...
                                    useElasticTendon,...
                                    flag_updatePositionLevel,...
                                    flag_evaluateJacobian, ...
                                    flag_evaluateDerivatives, ...
                                    flag_updateModelCache,...
                                    flag_evaluateInitializationFunctions,...
                                    flag_useArgs);      

      %Numerically build the jacobian
      errL = errI;
      errR = errI;
      errJacNum = errIJac;
      delta = sqrt(eps);
      delta_lce = 0;
      delta_laH = 0;
      delta_l1H = 0;

        for i=1:1:length(errI)
          
          delta_lce = 0;
          delta_laH = 0;
          delta_l1H = 0;

          if(useElasticTendon==1)
            switch(i)
              case 1
                delta_lce = delta;
              case 2
                delta_laH = delta;        
              case 3
                delta_l1H = delta;
            end
          else
            switch(i)
              case 1
                delta_laH = delta;        
              case 2
                delta_l1H = delta;
            end
          end

          modelCachedValues.lce  = lce + delta_lce;
          modelCachedValues.laH  = laH + delta_laH;
          modelCachedValues.dlaH = 0;
          modelCachedValues.l1H  = l1H + delta_l1H;

          [errF,errFJac,errIR,errJacEmpty,modelCachedValuesUpd] = ...
            updateMillard2023VexatCache( vars,...
                                        modelCachedValues,...
                                        modelCurves,...
                                        modelConstants,...
                                        sarcomereProperties,...
                                        useElasticTendon,...
                                        flag_updatePositionLevel,...
                                        flag_evaluateJacobian, ...
                                        flag_evaluateDerivatives, ...
                                        flag_updateModelCache,...
                                        flag_evaluateInitializationFunctions,...
                                        flag_useArgs);      

          modelCachedValues.lce  = lce - delta_lce;
          modelCachedValues.laH  = laH - delta_laH;
          modelCachedValues.dlaH = 0;
          modelCachedValues.l1H  = l1H - delta_l1H;

          [errF,errFJac,errIL,errJacEmpty,modelCachedValuesUpd] = ...
            updateMillard2023VexatCache( vars,...
                                        modelCachedValues,...
                                        modelCurves,...
                                        modelConstants,...
                                        sarcomereProperties,...
                                        useElasticTendon,...
                                        flag_updatePositionLevel,...
                                        flag_evaluateJacobian, ...
                                        flag_evaluateDerivatives, ...
                                        flag_updateModelCache,...
                                        flag_evaluateInitializationFunctions,...
                                        flag_useArgs); 
          errIJac(:,i) = (errIR-errIL)./(2*delta);
        end
         
        deltaVar = -pinv(errIJac)*errI;
         
        if(useElasticTendon == 1)
          if(lce + deltaVar(1) < lceMin)
            stepLength = (lceMin-lce)/deltaVar(1);
            stepLength = stepLength*0.5;
            deltaVar = deltaVar.*stepLength;
          end
          lce = lce + deltaVar(1);
          laH = laH + deltaVar(2);
          l1H = l1H + deltaVar(3);   
        else
          laH = laH + deltaVar(1);
          l1H = l1H + deltaVar(2);           
        end
               
        iter=iter+1;
    end

    assert(iter <= iterMax);
    assert(errI'*errI <= tolInit,...
        ['Error: calcMillard2023VexatMuscleInfo: initialization error ',...
          sprintf('%1.2e',errI'*errI),' exceeds tolerance ', ...
          sprintf('%1.2e',tolInit)]);
    
    flag_updatePositionLevel                = 1;
    flag_evaluateJacobian                   = 1;
    flag_evaluateDerivatives                = 1;
    flag_updateModelCache                   = 1;
    
    modelCachedValues.lce  = lce;
    modelCachedValues.laH  = laH;
    modelCachedValues.dlaH = 0;
    modelCachedValues.l1H  = l1H;

    [errF,errFJac,errI,errIJac,modelCachedValuesUpd] = ...
      updateMillard2023VexatCache( vars,...
                                  modelCachedValues,...
                                  modelCurves,...
                                  modelConstants,...
                                  sarcomereProperties,...
                                  useElasticTendon,...
                                  flag_updatePositionLevel,...
                                  flag_evaluateJacobian, ...
                                  flag_evaluateDerivatives, ...
                                  flag_updateModelCache,...
                                  flag_evaluateInitializationFunctions,...
                                  flag_useArgs); 

    modelCachedValues = modelCachedValuesUpd;  
  end
  
elseif(useElasticTendon ==0 && modelConfig.initializeState==0) 

  flag_updatePositionLevel                = 1;
  flag_evaluateJacobian                   = 1;
  flag_evaluateDerivatives                = 1;
  flag_updateModelCache                   = 1;
  flag_evaluateInitializationFunctions    = 0;
  flag_useArgs                            = 0;
  vars = [];
  
  [errF,errFJac,errI,errIJac,modelCachedValuesUpd] = ...
    updateMillard2023VexatCache( vars,...
                                modelCachedValues,...
                                modelCurves,...
                                modelConstants,...
                                sarcomereProperties,...
                                useElasticTendon,...
                                flag_updatePositionLevel,...
                                flag_evaluateJacobian, ...
                                flag_evaluateDerivatives, ...
                                flag_updateModelCache,...
                                flag_evaluateInitializationFunctions,...
                                flag_useArgs);

  modelCachedValues = modelCachedValuesUpd;
  
elseif(useElasticTendon ==1 && modelConfig.initializeState==0) 
      
    vars  = 0;

    flag_updatePositionLevel              = 1;
    flag_evaluateJacobian                 = 1;
    flag_evaluateDerivatives              = 1;
    flag_updateModelCache                 = 1;
    flag_evaluateInitializationFunctions  = 0;
    flag_useArgs                          = 0;
    
    errF    = 0;
    errFJac = 0;
    errI = zeros(3,1);
    errIJac = zeros(3,3);    
    
    [errF,errFJac,errI,errIJac,modelCachedValuesUpd] = ...
        updateMillard2023VexatCache( vars,modelCachedValues,...
                                    modelCurves,...
                                    modelConstants,...
                                    sarcomereProperties,...
                                    useElasticTendon,...
                                    flag_updatePositionLevel,...
                                    flag_evaluateJacobian, ...
                                    flag_evaluateDerivatives, ...
                                    flag_updateModelCache, ...
                                    flag_evaluateInitializationFunctions,...
                                    flag_useArgs);

    errNorm = sqrt(errF'*errF);
    assert(errNorm <= modelConfig.tol,'Error tolerance not met');

    modelCachedValues = modelCachedValuesUpd;

    %
    % Now that dlce is being analytically evaluated the iterative code
    % below is no longer needed. For now I'm keeping it in place just
    % incase I update the model later ... and need to once again double
    % check if I have made a dumb mistake.
    %
    flag_iterativelySolveForCEVelocity = 0;
    if(flag_iterativelySolveForCEVelocity == 1)    

        flag_updatePositionLevel = 0;
        iterMax = modelConfig.iterMax;
        tol     = modelConfig.tol;
        iter    = 0;
        err     = 0;
        errNorm = tol*100;
    
        errJac        = zeros(length(vars),length(vars));
        errJacNum     = zeros(length(vars),length(vars));
    
        h = sqrt(eps); %purturbation for building the error
                       %Jacobian using finite differences.
    
        tmpCache              = modelCachedValues;
    
        eqVars       = zeros(4,1);
        eqVarsL      = zeros(4,1);
        eqVarsR      = zeros(4,1);
        eqVarsJac    = zeros(4,1);
        eqVarsNumJac = zeros(4,1);
    
        errL = zeros(1,1);
        errR = zeros(1,1);
    
        flag_calcNumJac = 0;
        %This flag, and the while loop below, is only in place to 
        %numerically check that the analytic Jacobians that I'm 
        %evaluating are equal to a numerical derivative (within tolerance)
        
        if(a > 0.01)
          here=1;
        end
        delta=1;
        
        while( (errNorm > tol) && iter < iterMax)
          flag_updatePositionLevel                = 0;
          flag_evaluateJacobian                   = 1;
          flag_evaluateDerivatives                = 0;
          flag_updateModelCache                   = 1;
          flag_useArgs                            = 0;
          
          if(a > 0.01 && abs(dlp) > 0.01 && useElasticTendon==1)
            here=1;
          end
    
          [errF,errFJac,errI,errIJac,modelCachedValuesUpd] = ...
            updateMillard2023VexatCache(...
                                        vars,...
                                        modelCachedValues,...
                                        modelCurves,...
                                        modelConstants,...
                                        sarcomereProperties,...
                                        useElasticTendon,...
                                        flag_updatePositionLevel,...
                                        flag_evaluateJacobian, ...
                                        flag_evaluateDerivatives, ...
                                        flag_updateModelCache, ...
                                        flag_evaluateInitializationFunctions,...
                                        flag_useArgs);
    
          vars = modelCachedValuesUpd.dlce;
          errNorm = sqrt(errF'*errF);
    
            eqVars(1,1) = modelCachedValuesUpd.f2HN;
            eqVars(2,1) = modelCachedValuesUpd.fEcmHN;
            eqVars(3,1) = modelCachedValuesUpd.fxHN;
            eqVars(4,1) = modelCachedValuesUpd.fTN;      
          
          if(flag_calcNumJac==1)
            eqVars(1,1) = modelCachedValuesUpd.f2HN;
            eqVars(2,1) = modelCachedValuesUpd.fEcmHN;
            eqVars(3,1) = modelCachedValuesUpd.fxHN;
            eqVars(4,1) = modelCachedValuesUpd.fTN;      
    
            eqVarsJac(1,1) = modelCachedValuesUpd.D_f2HN_D_dlce;
            eqVarsJac(2,1) = modelCachedValuesUpd.D_fEcmHN_D_dlce;
            eqVarsJac(3,1) = modelCachedValuesUpd.D_fxHN_D_dlce;
            eqVarsJac(4,1) = modelCachedValuesUpd.D_fTN_D_dlce;      
    
            eqVarsL = zeros(4,1);
            eqVarsR = zeros(4,1);
                    
            flag_updatePositionLevel                = 0;
            flag_evaluateJacobian                   = 1;
            flag_evaluateDerivatives                = 0;
            flag_updateModelCache                   = 1;
            flag_useArgs                            = 1;
            for i=1:1:length(vars)
              for(j=1:1:2)
                if(j==1)
                  varsL     = vars;
                  varsL(i)  = varsL(i)-h;            
                  tmpCache = modelCachedValues;
                  [errL,errJacTmp,errI,errIJac,tmpCache] = ...
                    updateMillard2023VexatCache(...
                                        varsL,...
                                        modelCachedValues,...
                                        modelCurves,...
                                        modelConstants,...
                                        sarcomereProperties,...
                                        useElasticTendon,...
                                        flag_updatePositionLevel,...
                                        flag_evaluateJacobian, ...
                                        flag_evaluateDerivatives, ...
                                        flag_updateModelCache,...
                                        flag_evaluateInitializationFunctions,...
                                        flag_useArgs);
    
                  eqVarsL(1,i) = tmpCache.f2HN;
                  eqVarsL(2,i) = tmpCache.fEcmHN;
                  eqVarsL(3,i) = tmpCache.fxHN;
                  eqVarsL(4,i) = tmpCache.fTN;
                else
                  varsR     = vars;
                  varsR(i)  = varsR(i)+h;
                  tmpCache = modelCachedValues;
                  [errR,errJacTmp,errI, errIJac,tmpCache] = ...
                    updateMillard2023VexatCache(...
                                        varsR,...
                                        modelCachedValues,...
                                        modelCurves,...
                                        modelConstants,...
                                        sarcomereProperties,...
                                        useElasticTendon,...
                                        flag_updatePositionLevel,...
                                        flag_evaluateJacobian, ...
                                        flag_evaluateDerivatives, ...
                                        flag_updateModelCache, ...
                                        flag_evaluateInitializationFunctions,...
                                        flag_useArgs);
    
                  eqVarsR(1,i) = tmpCache.f2HN;
                  eqVarsR(2,i) = tmpCache.fEcmHN;
                  eqVarsR(3,i) = tmpCache.fxHN;
                  eqVarsR(4,i) = tmpCache.fTN;
                end                  
              end
              eqVarsNumJac(:,i) = (eqVarsR(:,i)-eqVarsL(:,i))./(2*h);
              errJacNum(:,i) = (errR-errL)./(2*h);
            end
    
            eqVarsJacErr = eqVarsJac - eqVarsNumJac;
            eqVarsJacErrNorm = sqrt(sum(sum(eqVarsJacErr.^2)));
    
            errJacErr = errFJac - errJacNum;
            errJacErrNorm = sqrt(sum(sum(errJacErr.^2)));
            here=1;
          end
    
          %Debugging
          delta = -errFJac\errF;                    
          vars  = vars+delta;
          iter=iter+1;
          
          flag_debugging=0;
          if(iter > 30)
            here=1; 
            flag_debugging=1;
          end
          
          
          if(flag_debugging==1)
            fprintf('%1.2e\t%1.2e\n',...
                    errF,errFJac);
            here=1;
          end
          %if(errJacErrNorm > 0.1)
          %  here=1;
          %end
            
          
        end
    
        
        
        
        assert(errNorm <= tol,'Error tolerance not met');
    
          flag_updatePositionLevel                = 0;
          flag_evaluateJacobian                   = 1;
          flag_evaluateDerivatives                = 1;
          flag_updateModelCache                   = 1;
          flag_useArgs                            = 0;
          [errF,errFJac,errI,errIJac,modelCachedValuesUpd] = ...
            updateMillard2023VexatCache( vars,...
                                        modelCachedValues,...
                                        modelCurves,...
                                        modelConstants,...
                                        sarcomereProperties,...
                                        useElasticTendon,...
                                        flag_updatePositionLevel,...
                                        flag_evaluateJacobian, ...
                                        flag_evaluateDerivatives, ...
                                        flag_updateModelCache,...
                                        flag_evaluateInitializationFunctions,...
                                        flag_useArgs);
    
          modelCachedValues = modelCachedValuesUpd;

    end
else
    assert(0, 'Error: somehow an impossible(!) case has been reached');
end

%%
% Update the state derivative
dlce    = modelCachedValues.dlce;

%%==============================================================================
%%==============================================================================

%%
% Calculate the stiffness & damping of the tendon
%%
ktN     = modelCachedValues.kTNN*(1/ltN_lt);
kt      = fiso*ktN;

dtN     = modelCachedValues.betaTNN*(1/ltN_lt);
dt      = fiso*dtN;


%%
% Calculate the stiffness of titin segments 1 and 2
%%
k1HN  = modelCachedValues.k1HNN*(1/liN_li);
k1H   = fiso*k1HN; %N/m
k1N   = k1HN*0.5;
k1    = k1H*0.5;

k2HN  = modelCachedValues.k2HNN*(1/liN_li);
k2H   = fiso*k2HN; %N/m
k2N   = k2HN*0.5;
k2    = k2H*0.5;

%%
% Calculate the damping of titin segments 1 and 2
%   Note: beta1HNN and beta2HNN set accordingly to the titin model selected
%         in updateMillard2023VexatCache
%%

d1HN = modelCachedValues.beta1HNN*(1/liN_li);
d1H  = fiso*d1HN; %N/(m/s)
d1N  = d1HN*0.5;
d1   = d1H*0.5;

d2HN = modelCachedValues.beta2HNN*(1/liN_li);
d2H  = fiso*d2HN; %N/(m/s)
d2N  = d2HN*0.5;
d2   = d2H*0.5;

%%
% ECM stiffness & damping
%%
keHN = modelCachedValues.kEcmHNN*(1/lceN_lce);
keH  = fiso*keHN;
keN  = keHN*0.5;
ke   = keH*0.5;

deHN = modelCachedValues.betaEcmHNN*(1/lceN_lce);
deH  = fiso*deHN;
deN  = deHN*0.5;
de   = deH*0.5;

%%
% Compressive element stiffness & damping
%%
kCpN = modelCachedValues.DfCpN_DlceN*(1/lceN_lce);
kCp  = fiso*kCpN;

%Philosophy point:
%  The bond between titin and actin is not really a physical damper. Its 
%  probably electrostatic and has a high stiffness and zero damping. However
%  the way the bond slides along actin is reministent of damping - much like the
%  modelled bond between the XE element and actin.
%
%  And thus the damping contribution that we evaluate here, while mathematically
%  correct, should be ignored because this damping only arises because of the
%  simplified model we are using to model the titin-actin bond. Here I leave
%  the de value unchanged to be mathematically consistent. However, I do not
%  expect to see any evidence of damping from the titin-actin bond.
%
%  By leaving this de value as is, during the Herzog & Leonard 2002 simulations
%  there is a curious spike in the damping coefficient when the muscle is de
%  activated.

%%
% Calculate the stiffness & damping of the cross-bridges
%%
kxHN = modelCachedValues.kxHNN*(1/lceN_lce);
kxH = kxHN*fiso;
kxN = 0.5*kxHN;
kx  = 0.5*kxH;

dxHN = modelCachedValues.betaxHNN*(1/lceN_lce);
dxH  = dxHN*fiso;
dxN  = 0.5*dxHN;
dx   = 0.5*dxH;

%%
% Calculate the stiffness and damping of the fiber
%%
z1  = complex(k1,d1); %z: impedance
z2  = complex(k2,d2);
a12 = (1/z1)+(1/z2) ; %a: admittance (reciprocal of impedance)
z12 = 1/a12;

zEcm = complex(ke,de);

zx = complex(kx,dx);

zCp = complex(kCp,0);
zfCp= complex(modelCachedValues.fCpN*fiso,0);

Dalpha_Dlce = modelCachedValues.Dalpha_Dlce;
cosAlpha    = modelCachedValues.cosAlpha;
sinAlpha    = modelCachedValues.sinAlpha;

%Resolve the compressive element into the direction of the fiber
zf = z12+zEcm+zx ...
    - zCp/cosAlpha ...
    - ( (-zfCp/(cosAlpha*cosAlpha)) * (-sinAlpha*Dalpha_Dlce) );

kf = real(zf);
df = imag(zf);

%%
% Resolve it along the direction of hte tendon
%%

fceN  = modelCachedValues.fceN;

kfAT = kf*cosAlpha - fiso*fceN*sinAlpha*Dalpha_Dlce; 
dfAT = df*cosAlpha;      

%%
% Calculate the impedance of the whole musculotendon
%%

km = 0;
dm = 0;

if(useElasticTendon)
  zf = complex(kfAT,dfAT);
  zt = complex(kt,dt);
  am = (1/zf)+(1/zt);
  zm = 1/am;
  km = real(zm);
  dm = imag(zm);
  
else
  km = kfAT;
  dm = dfAT;
end


% fxAT = fx*cos(alpha)
% kxAT = (d_fx/d_l)*cos(alpha) - fx*sin(alpha)*Dalpha_Dlce
% dxAT = (d_fx/d_dlce) * cos(alpha)

% Dalpha_Dlce = modelCachedValues.Dalpha_Dlce;
% cosAlpha    = modelCachedValues.cosAlpha;
% sinAlpha    = modelCachedValues.sinAlpha;

%%
%Calculate the impedance of the contractile element
%%


%fce   = fceN*fiso;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4. Populate the muscle length and velocity information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Position Information
%%

mtInfo.muscleLengthInfo.tendonLength            = modelCachedValues.lt;        %length        m
mtInfo.muscleLengthInfo.normTendonLength        = modelCachedValues.ltN;       %length/length m/m        
mtInfo.muscleLengthInfo.tendonStrain            = modelCachedValues.et;        %length/length m/m        

mtInfo.muscleLengthInfo.normFiberSlidingLength  = modelCachedValues.laH;
mtInfo.muscleLengthInfo.fiberSlidingLength      = modelCachedValues.laHN*lceOpt;

mtInfo.muscleLengthInfo.fiberLength             = modelCachedValues.lce;       %length        m  
mtInfo.muscleLengthInfo.fiberLengthAlongTendon  = modelCachedValues.lceAT;     %length        m
mtInfo.muscleLengthInfo.normFiberLength         = modelCachedValues.lceN;      %length/length m/m   

mtInfo.muscleLengthInfo.pennationAngle          = modelCachedValues.alpha;     %angle         1/s    
mtInfo.muscleLengthInfo.cosPennationAngle       = modelCachedValues.cosAlpha;  %NA            NA         
mtInfo.muscleLengthInfo.sinPennationAngle       = modelCachedValues.sinAlpha;  %NA            NA         	

mtInfo.muscleLengthInfo.titin1Length            = modelCachedValues.l1H;        %length         m  
mtInfo.muscleLengthInfo.normTitin1Length        = modelCachedValues.l1HN;       %length/length  m/m  
mtInfo.muscleLengthInfo.titin2Length            = modelCachedValues.l2H;        %length         m  
mtInfo.muscleLengthInfo.normTitin2Length        = modelCachedValues.l2HN;       %length/length  m/m  

mtInfo.muscleLengthInfo.crossBridgeLength     = modelCachedValues.lxH;
mtInfo.muscleLengthInfo.normCrossBridgeLength = modelCachedValues.lxHN;



mtInfo.muscleLengthInfo.isClamped               = 0;%isFiberClamped;% || isCrossBridgeClamped;



mtInfo.muscleLengthInfo.fiberPassiveForceLengthMultiplier = modelCachedValues.fEcmfcnHN;
mtInfo.muscleLengthInfo.fiberActiveForceLengthMultiplier = modelCachedValues.flN;

%%
%Velocity information
%%

mtInfo.fiberVelocityInfo.tendonVelocity               = modelCachedValues.dlt;     %length/time           m/s
mtInfo.fiberVelocityInfo.normTendonVelocity           = modelCachedValues.dltN;    %(length/time)/length (m/s)/m

mtInfo.muscleLengthInfo.normFiberSlidingVelocity  = modelCachedValues.dlfNN;
mtInfo.muscleLengthInfo.fiberSlidingVelocity      = modelCachedValues.dlfNN*lceOpt*dlceMaxN;

mtInfo.fiberVelocityInfo.fiberVelocity                = modelCachedValues.dlce;    %length/time           m/s
mtInfo.fiberVelocityInfo.fiberVelocityAlongTendon     = modelCachedValues.dlceAT;  %length/time           m/s
mtInfo.fiberVelocityInfo.normFiberVelocity            = modelCachedValues.dlce/(lceOpt*dlceMaxN);   %(length/time)/length (m/s)/m
mtInfo.fiberVelocityInfo.pennationAngularVelocity     = modelCachedValues.dalpha;  %angle/time            rad/s
mtInfo.fiberVelocityInfo.fiberForceVelocityMultiplier = modelCachedValues.fvN;                       %force/force

mtInfo.fiberVelocityInfo.titin1Velocity            = modelCachedValues.dl1H;        %length/time           m/s  
mtInfo.fiberVelocityInfo.normTitin1Velocity        = modelCachedValues.dl1HN;       %(length/time)/length (m/s)/m  
mtInfo.fiberVelocityInfo.titin2Velocity            = modelCachedValues.dl2H;        %length/time           m/s  
mtInfo.fiberVelocityInfo.normTitin2Velocity        = modelCachedValues.dl2HN;       %(length/time)/length (m/s)/m  

mtInfo.muscleLengthInfo.crossBridgeVelocity     = modelCachedValues.dlxH;
mtInfo.muscleLengthInfo.normCrossBridgeVelocity = modelCachedValues.dlxHN;


   
%%
%7. Evaluate the potential energy stored in the musculotendon
%%
pt  = 0;
pi1 = 0;
pi2 = 0;
pc  = 0;

% Leaving this for later: this code requires a (time consuming) update
% and this is not needed to assess the physical correctness of the model

%if(isempty(normMuscleCurves.tendonForceLengthCurve.integral) == 0)
%    if(useElasticTendon == 1)
%       ptN = calcFtDer(modelCachedValues.ltN, -1); 
       
%       tendonStrainAtOneNormForce = ...
%           normMuscleCurves.tendonForceLengthCurve.integral.xScaling;
%       tendonStretchAtOneNormForce = tendonStrainAtOneNormForce*ltSlk;
       
       %To get this scaling we're getting calculating the energy in
       %a square that is 100% strain by 1 maximum isometric force
       %      first in units of Nm 
       %      second in normalized units
       %And then taking the ratio to obtain the scaling
%       tendonEnergyScaling = tendonStretchAtOneNormForce*fiso ...
%                           / tendonStrainAtOneNormForce*1;
%
%       pt = ptN*tendonEnergyScaling;

%    end
%end



pe = pt+pi1+pi2+pc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%8. Populate the dynamics and energy information structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mtInfo.muscleDynamicsInfo.activation                = a;
mtInfo.muscleDynamicsInfo.activationDerivative      = dadt;

mtInfo.muscleDynamicsInfo.lambda                    = modelCachedValues.lambda;
mtInfo.muscleDynamicsInfo.dlambda                   = modelCachedValues.dlambda;

mtInfo.muscleDynamicsInfo.tendonForce               = modelCachedValues.fTN*fiso;     % force                N
mtInfo.muscleDynamicsInfo.normTendonForce           = modelCachedValues.fTN;          % force/force          N/N
mtInfo.muscleDynamicsInfo.tendonStiffness           = kt;                             % force/length         N/m
mtInfo.muscleDynamicsInfo.tendonDamping             = dt;

mtInfo.muscleDynamicsInfo.fiberForce                = fceN*fiso;                             % force                N
mtInfo.muscleDynamicsInfo.fiberForceAlongTendon     = fceN*fiso*modelCachedValues.cosAlpha;  % force                N
mtInfo.muscleDynamicsInfo.normFiberForce            = fceN;                                 % force/force          N/N


%The active force acting between the myosin & actin
mtInfo.muscleDynamicsInfo.activeFiberForce          = modelCachedValues.fxHN*fiso;   % force                N
mtInfo.muscleDynamicsInfo.normActiveFiberForce      = modelCachedValues.fxHN;        % force                N

mtInfo.muscleDynamicsInfo.passiveFiberForce         = (modelCachedValues.f2HN+modelCachedValues.fEcmHN)*fiso;  % force                N
mtInfo.muscleDynamicsInfo.normPassiveFiberForce     = (modelCachedValues.f2HN+modelCachedValues.fEcmHN);  % force                N

%This is not so simple any more thanks to the titin segment:
%  no damping      : easy - titin is two springs in series
%  infinite damping: easy - just the PEVK region stretches
%  between         : this is an LTI system with a frequency dependence.
mtInfo.muscleDynamicsInfo.crossBridgeStiffness      = kx;  % force/length         N/m
mtInfo.muscleDynamicsInfo.crossBridgeDamping        = dx;  % force/length         N/m

%mtInfo.muscleDynamicsInfo.muscleStiffness           = NaN;                % force/length         N/m
%mtInfo.muscleDynamicsInfo.normMuscleStiffness       = NaN;  % force/length         N/m

%Doing this later.
mtInfo.muscleDynamicsInfo.fiberActivePower          = 0;   % force*velocity       W
mtInfo.muscleDynamicsInfo.fiberPassivePower         = 0;    % force*velocity       W
mtInfo.muscleDynamicsInfo.tendonPower               = 0;    % force*velocity       W
mtInfo.muscleDynamicsInfo.musclePower               = 0;  % force*velocity       W
mtInfo.muscleDynamicsInfo.fiberParallelElementPower = 0;

mtInfo.muscleDynamicsInfo.ecmStiffness = ke;
mtInfo.muscleDynamicsInfo.ecmDamping   = de;
mtInfo.muscleDynamicsInfo.titin1Stiffness = k1;
mtInfo.muscleDynamicsInfo.titin2Stiffness = k2;
mtInfo.muscleDynamicsInfo.titin1Damping = d1;
mtInfo.muscleDynamicsInfo.titin2Damping = d2;


mtInfo.muscleDynamicsInfo.titin1Force                = modelCachedValues.f1HN*fiso;  %length/time           m/s  
mtInfo.muscleDynamicsInfo.normTitin1Force            = modelCachedValues.f1HN;       %(length/time)/length (m/s)/m  
mtInfo.muscleDynamicsInfo.titin2Force                = modelCachedValues.f2HN*fiso;  %length/time           m/s  
mtInfo.muscleDynamicsInfo.normTitin2Force            = modelCachedValues.f2HN;       %(length/time)/length (m/s)/m  

mtInfo.muscleDynamicsInfo.titin1SpringForce          = modelCachedValues.f1HN*fiso;  %length/time           m/s  
mtInfo.muscleDynamicsInfo.normTitin1SpringForce      = modelCachedValues.f1HN;       %(length/time)/length (m/s)/m  
mtInfo.muscleDynamicsInfo.titin2SpringForce          = modelCachedValues.f2HN*fiso;  %length/time           m/s  
mtInfo.muscleDynamicsInfo.normTitin2SpringForce      = modelCachedValues.f2HN;       %(length/time)/length (m/s)/m  

mtInfo.muscleDynamicsInfo.titin1DampingForce          = (modelCachedValues.f2HN-modelCachedValues.f1HN)*fiso;  %length/time           m/s  
mtInfo.muscleDynamicsInfo.normTitin1DampingForce      = (modelCachedValues.f2HN-modelCachedValues.f1HN);       %(length/time)/length (m/s)/m  
mtInfo.muscleDynamicsInfo.titin2DampingForce          = 0;  %length/time           m/s  
mtInfo.muscleDynamicsInfo.normTitin2DampingForce      = 0;       %(length/time)/length (m/s)/m  

mtInfo.muscleDynamicsInfo.crossBridgeForce            = modelCachedValues.fxHN*fiso;
mtInfo.muscleDynamicsInfo.normCrossBridgeForce        = modelCachedValues.fxHN;

mtInfo.muscleDynamicsInfo.crossBridgeSpringForce      = modelCachedValues.kxHNN*modelCachedValues.lxHN*fiso;
mtInfo.muscleDynamicsInfo.normCrossBridgeSpringForce  = modelCachedValues.kxHNN*modelCachedValues.lxHN;

mtInfo.muscleDynamicsInfo.crossBridgeDampingForce     = modelCachedValues.betaxHNN*modelCachedValues.dlxHN*fiso;
mtInfo.muscleDynamicsInfo.normCrossBridgeDampingForce = modelCachedValues.betaxHNN*modelCachedValues.dlxHN;

%Work: doing this later
%mtInfo.muscleDynamicsInfo.crossBridgeForce          = modelCachedValues.faN;
mtInfo.muscleDynamicsInfo.dampingForces             =  0;
mtInfo.muscleDynamicsInfo.dampingPower              =  0;
mtInfo.muscleDynamicsInfo.boundaryPower             =  modelCachedValues.fTN*fiso*modelCachedValues.dlp;

%
mtInfo.muscleDynamicsInfo.fiberStiffness      = kf;
mtInfo.muscleDynamicsInfo.fiberDamping        = df;
mtInfo.muscleDynamicsInfo.normFiberStiffness  = kf/(fiso/lceOpt);
mtInfo.muscleDynamicsInfo.normFiberDamping    = df/(fiso/lceOpt);

mtInfo.muscleDynamicsInfo.musculotendonStiffness      = km;
mtInfo.muscleDynamicsInfo.musculotendonDamping        = dm;
mtInfo.muscleDynamicsInfo.normMusculotendonStiffness  = km/(fiso/lceOpt);
mtInfo.muscleDynamicsInfo.normMusculotendonDamping    = dm/(fiso/lceOpt);


%Energy: doing this later
mtInfo.musclePotentialEnergyInfo.fiberPotentialEnergy  = pi1+pi2+pc;    %force*distance    J   
mtInfo.musclePotentialEnergyInfo.tendonPotentialEnergy = pt;            %force*distance    J   
mtInfo.musclePotentialEnergyInfo.musclePotentialEnergy = pi1+pi2+pc+pt; % + pbr;%force*distance    J (Nm)

mtInfo.muscleDynamicsInfo.dT = modelCachedValues.dT;
mtInfo.muscleDynamicsInfo.dW = modelCachedValues.dW;
mtInfo.muscleDynamicsInfo.dV = modelCachedValues.dV;


mtInfo.state.value      = NaN;
mtInfo.state.derivative = NaN;


if useElasticTendon == 1
  mtInfo.state.value      = [modelCachedValues.lce;...
                             modelCachedValues.dlaH; ...
                             modelCachedValues.laH;...
                             modelCachedValues.l1H];%...
                             %modelCachedValues.lambda];


  mtInfo.state.derivative = [modelCachedValues.dlce;...
                             modelCachedValues.ddlaH; ...
                             modelCachedValues.dlaH;...
                             modelCachedValues.dl1H];%...
                             %modelCachedValues.dlambda];
else
  mtInfo.state.value      = [modelCachedValues.dlaH; ...
                             modelCachedValues.laH;...
                             modelCachedValues.l1H];%...
                             %modelCachedValues.lambda];


  mtInfo.state.derivative = [modelCachedValues.ddlaH; ...
                             modelCachedValues.dlaH;...
                             modelCachedValues.dl1H];%...
                             %modelCachedValues.dlambda];
end


mtInfo.extra = [modelCachedValues.l1HN,...               
                modelCachedValues.f1kHN,...
                modelCachedValues.f1dHN,...
                modelCachedValues.f1HN,...
                modelCachedValues.l2HN,...
                modelCachedValues.f2kHN,...
                modelCachedValues.f2dHN,...
                modelCachedValues.f2HN];
mtInfo.extraLabels={'l1HN','f1kN','f1dN','f1N',...
                    'l2HN','f2kN','f2dN','f2N'};

if(sum(isnan(mtInfo.state.derivative)) > 0)
   here = 1; 
end
