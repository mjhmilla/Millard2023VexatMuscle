%%
% 2nd order model.
%
%%
function mtInfo = calcMillard2019MuscleInfoOpus31(...
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
%                                                       |   |--|/\/\|--|
% |     IgP       PEVK    IgD                           |---|    bt    |---
% |---|\/\/\|--|\/\/\/|--------|------------------------|   |--[    ]--|
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


%Titin proximal force-length curve: Z-line to N2A
if(normMuscleCurves.useTitinCurvesWithRigidIgDSegment==1)  
    calcF1HDer       = @(arg1, arg2)calcBezierYFcnXDerivative(arg1, ...
                          normMuscleCurves.forceLengthIgpCurveB, ...
                          arg2);  
    calcF2HDer       = @(arg1, arg2)calcBezierYFcnXDerivative(arg1, ...
                          normMuscleCurves.forceLengthPevkCurveB, ...
                          arg2);  
else
    calcF1HDer       = @(arg1, arg2)calcBezierYFcnXDerivative(arg1, ...
                          normMuscleCurves.forceLengthIgpCurve, ...
                          arg2);  
    calcF2HDer       = @(arg1, arg2)calcBezierYFcnXDerivative(arg1, ...
                          normMuscleCurves.forceLengthPevkIgdCurve, ...
                          arg2);  
end                    
                    
%Normalized length of the rigid IG2 section
%disp('Replace this with a fixed parameter so that these curves are not required')
%lIG2HN  = 0.5 ...
%          - (normMuscleCurves.titinPEVK2ForceLengthHalfCurve.xpts(1,1)) ... 
%          - (normMuscleCurves.titinZN2AForceLengthHalfCurve.xpts(1,1));
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
                
                
DftfcnN_DltN_max = normMuscleCurves.tendonForceLengthCurve.dydxEnd(1,2);


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
scalePEVK  = sarcomereProperties.scalePEVK;
scaleIGP   = sarcomereProperties.scaleIGP;

kAXHN      =   sarcomereProperties.normCrossBridgeStiffness;
betaAXHN   =   sarcomereProperties.normCrossBridgeDamping;

forceVelocityCalibrationFactor =   sarcomereProperties.forceVelocityCalibrationFactor;

betafEcmHN =   sarcomereProperties.normECMDamping;

betaN2AaHN = sarcomereProperties.normMaxActiveTitinToActinDamping;
betaN2ApHN = sarcomereProperties.normPassiveTitinToActinDamping;

slidingTimeConstant = sarcomereProperties.slidingTimeConstant;

normCrossBridgeCyclingDamping = sarcomereProperties.normCrossBridgeCyclingDamping;
lowActivationThreshold        = sarcomereProperties.lowActivationThreshold;

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
 assert( (lce >= muscleArchitecture.minimumFiberLength ...
               || modelConfig.initializeState == 1 ...
               || useElasticTendon == 0));

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
        'calcFalDer', calcFalDer ,...
        'calcFvDer'  , calcFvDer   ,...  
        'calcFeHDer' , calcFeHDer  ,...        
        'calcF1HDer' , calcF1HDer  ,...
        'calcF2HDer' , calcF2HDer  ,...
        'calcFtDer'  , calcFtDer   ,...
        'calcDtDer'  , calcDtDer );




 
      
modelConstants = struct( ...
    'fiso'          , fiso             , ...
    'lceOpt'        , lceOpt           , ...    
    'alphaOpt'      , alphaOpt         , ...
    'ltSlk'         , ltSlk            , ...
    'dlceMaxN'      , dlceMaxN         , ...
    'lceMin'        , lceMin           , ...
    'kAXHN'         , kAXHN            , ...
    'betaAXHN'      , betaAXHN         , ...
    'normCrossBridgeCyclingDamping'  , normCrossBridgeCyclingDamping, ...
    'lowActivationThreshold', lowActivationThreshold,...
      'betaNum'     , betaNum          , ...                          
      'betafTN'      , betafTN           , ...    
      'betaCTN'      , betaCTN           , ...          
    'betafEcmHN'      , betafEcmHN         , ...                    
    'betaN2AaHN'    , betaN2AaHN       , ...
    'betaN2ApHN'    , betaN2ApHN       , ...    
    'forceVelocityCalibrationFactor'  , forceVelocityCalibrationFactor,...
    'tau'       , slidingTimeConstant, ...
    'lTitinFixedHN' , lTitinFixedHN    , ...
    'lmHN'          , lmHN             , ...
    'lmH'           , lmH              , ...    
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
    'scalePEVK', scalePEVK,...
    'scaleIGP', scaleIGP,...
    'DftfcnN_DltN_max',DftfcnN_DltN_max);



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
      'k1HNN'         , NaN       , ...
      'beta1HNN'      , NaN       , ...
      'dw_1kH'         , NaN       , ...
    'f2H'           , NaN       , ...
    'f2HN'          , NaN       , ...
    'k2HNN'         , NaN       , ...
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

  assert(abs(dlp) < sqrt(eps)); %General case is not yet handled.

  errF    = 0;
  errFJac = 0;
  errI    = zeros(2+useElasticTendon,1);
  errIJac = zeros(2+useElasticTendon,2+useElasticTendon);
  


  lceAT  = max(lp-ltSlk,lceATMin);
  fibKin = calcFixedWidthPennatedFiberKinematics(lceAT,0,lceOpt,alphaOpt);
  lce    = max(fibKin.fiberLength,lceMin);  
  lceN = lce*lce_lceN;
  laHN = 0.5*lceN - lmHN;  
  laH  = laHN*lceN_lce;
  dlaH = 0;

  
  l1HN = normMuscleCurves.forceLengthIgpCurve.xEnd(1)*lceN;
  l1H = l1HN * liN_li;
    
  flag_updatePositionLevel                = 1;
  flag_evaluateJacobian                   = 1;
  flag_evaluateDerivatives                = 1;
  flag_updateModelCache                   = 1;
  flag_evaluateInitializationFunctions    = 1;
  
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
    calcEquilibriumErrorOpus31( vars,...
                                modelCachedValues,...
                                modelCurves,...
                                modelConstants,...
                                sarcomereProperties,...
                                useElasticTendon,...
                                flag_updatePositionLevel,...
                                flag_evaluateJacobian, ...
                                flag_evaluateDerivatives, ...
                                flag_updateModelCache,...
                                flag_evaluateInitializationFunctions);  
  
                              
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
      calcEquilibriumErrorOpus31( vars,...
                                  modelCachedValues,...
                                  modelCurves,...
                                  modelConstants,...
                                  sarcomereProperties,...
                                  useElasticTendon,...
                                  flag_updatePositionLevel,...
                                  flag_evaluateJacobian, ...
                                  flag_evaluateDerivatives, ...
                                  flag_updateModelCache,...
                                  flag_evaluateInitializationFunctions);      

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
          calcEquilibriumErrorOpus31( vars,...
                                      modelCachedValues,...
                                      modelCurves,...
                                      modelConstants,...
                                      sarcomereProperties,...
                                      useElasticTendon,...
                                      flag_updatePositionLevel,...
                                      flag_evaluateJacobian, ...
                                      flag_evaluateDerivatives, ...
                                      flag_updateModelCache,...
                                      flag_evaluateInitializationFunctions);      

        modelCachedValues.lce  = lce - delta_lce;
        modelCachedValues.laH  = laH - delta_laH;
        modelCachedValues.dlaH = 0;
        modelCachedValues.l1H  = l1H - delta_l1H;

        [errF,errFJac,errIL,errJacEmpty,modelCachedValuesUpd] = ...
          calcEquilibriumErrorOpus31( vars,...
                                      modelCachedValues,...
                                      modelCurves,...
                                      modelConstants,...
                                      sarcomereProperties,...
                                      useElasticTendon,...
                                      flag_updatePositionLevel,...
                                      flag_evaluateJacobian, ...
                                      flag_evaluateDerivatives, ...
                                      flag_updateModelCache,...
                                      flag_evaluateInitializationFunctions); 
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
  assert(errI'*errI <= tolInit);
  
  flag_updatePositionLevel                = 1;
  flag_evaluateJacobian                   = 1;
  flag_evaluateDerivatives                = 1;
  flag_updateModelCache                   = 1;
  
  modelCachedValues.lce  = lce;
  modelCachedValues.laH  = laH;
  modelCachedValues.dlaH = 0;
  modelCachedValues.l1H  = l1H;

  [errF,errFJac,errI,errIJac,modelCachedValuesUpd] = ...
    calcEquilibriumErrorOpus31( vars,...
                                modelCachedValues,...
                                modelCurves,...
                                modelConstants,...
                                sarcomereProperties,...
                                useElasticTendon,...
                                flag_updatePositionLevel,...
                                flag_evaluateJacobian, ...
                                flag_evaluateDerivatives, ...
                                flag_updateModelCache,...
                                flag_evaluateInitializationFunctions); 

  modelCachedValues = modelCachedValuesUpd;  
  
elseif(useElasticTendon ==0 && modelConfig.initializeState==0) 

  flag_updatePositionLevel                = 1;
  flag_evaluateJacobian                   = 1;
  flag_evaluateDerivatives                = 1;
  flag_updateModelCache                   = 1;
  flag_evaluateInitializationFunctions    = 0;
  vars = [];
  
  [errF,errFJac,errI,errIJac,modelCachedValuesUpd] = ...
    calcEquilibriumErrorOpus31( vars,...
                                modelCachedValues,...
                                modelCurves,...
                                modelConstants,...
                                sarcomereProperties,...
                                useElasticTendon,...
                                flag_updatePositionLevel,...
                                flag_evaluateJacobian, ...
                                flag_evaluateDerivatives, ...
                                flag_updateModelCache,...
                                flag_evaluateInitializationFunctions);

  modelCachedValues = modelCachedValuesUpd;
  
elseif(useElasticTendon ==1 && modelConfig.initializeState==0) 
    
    
    iterMax = modelConfig.iterMax;
    tol     = modelConfig.tol;
    iter    = 0;
    err     = 0;
    errNorm = tol*100;
  
    vars  = 0;

    flag_updatePositionLevel              = 1;
    flag_evaluateJacobian                 = 0;
    flag_evaluateDerivatives              = 0;
    flag_updateModelCache                 = 0;
    flag_evaluateInitializationFunctions  = 0;
    
    errF    = 0;
    errFJac = 0;
    errI = zeros(3,1);
    errIJac = zeros(3,3);
    
    [errF,errFJac,errI,errIJac,modelCachedValuesUpd] = ...
        calcEquilibriumErrorOpus31( vars,modelCachedValues,...
                                    modelCurves,...
                                    modelConstants,...
                                    sarcomereProperties,...
                                    useElasticTendon,...
                                    flag_updatePositionLevel,...
                                    flag_evaluateJacobian, ...
                                    flag_evaluateDerivatives, ...
                                    flag_updateModelCache, ...
                                    flag_evaluateInitializationFunctions);
    flag_updatePositionLevel = 0;
    modelCachedValues = modelCachedValuesUpd;
    
% 
%     stepSize = 0.5*lceOpt*dlceMaxN;
%     stepSign = [1 0 1 -1 0 -1 ; ...
%                 0 1 1 0 -1 -1];
%          
%     errNormBest = sqrt(errF'*errF);
%     errNorm     = errNormBest;
%     varsBest    = vars;
%     
%     for i=1:1:10
%         for j=1:1:size(stepSign,2)
%             vars = varsBest + stepSign(:,j).*stepSize;
%             [err,errJac,modelCachedValuesUpd] = ...
%                 calcEquilibriumErrorOpus31(...
%                                     vars,modelCachedValues,...
%                                     modelCurves,...
%                                     modelConstants,...
%                                     useElasticTendon,...
%                                     sarcomereProperties,...
%                                     flag_updatePositionLevel,...
%                                     flag_evaluateJacobian, ...
%                                     flag_evaluateDerivatives, ...
%                                     flag_updateModelCache); 
%             errNorm = sqrt(err'*err);
%             if(errNorm < errNormBest)
%                varsBest = vars; 
%                errNormBest=errNorm;
%             end
%         end
%         stepSize = stepSize*0.5;
%     end
%     
%     
%     vars = varsBest;


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
    
    if(a > 0.01)
      here=1;
    end
    delta=1;
    
    while( (errNorm > tol) && iter < iterMax)
      flag_updatePositionLevel                = 0;
      flag_evaluateJacobian                   = 1;
      flag_evaluateDerivatives                = 0;
      flag_updateModelCache                   = 1;
      
      [errF,errFJac,errI,errIJac,modelCachedValuesUpd] = ...
        calcEquilibriumErrorOpus31(...
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
                                    flag_evaluateInitializationFunctions);

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
        for i=1:1:length(vars)
          for(j=1:1:2)
            if(j==1)
              varsL     = vars;
              varsL(i)  = varsL(i)-h;            
              tmpCache = modelCachedValues;
              [errL,errJacTmp,errI,errIJac,tmpCache] = ...
                calcEquilibriumErrorOpus31(...
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
                                    flag_evaluateInitializationFunctions);

              eqVarsL(1,i) = tmpCache.f2HN;
              eqVarsL(2,i) = tmpCache.fEcmHN;
              eqVarsL(3,i) = tmpCache.fxHN;
              eqVarsL(4,i) = tmpCache.fTN;
            else
              varsR     = vars;
              varsR(i)  = varsR(i)+h;
              tmpCache = modelCachedValues;
              [errR,errJacTmp,errI, errIJac,tmpCache] = ...
                calcEquilibriumErrorOpus31(...
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
                                    flag_evaluateInitializationFunctions);

              eqVarsR(1,i) = tmpCache.f2HN;
              eqVarsR(2,i) = tmpCache.fEcmHN;
              eqVarsR(3,i) = tmpCache.fxHN;
              eqVarsR(4,i) = tmpCache.fTN;
            end      
            eqVarsNumJac(:,i) = (eqVarsR(:,i)-eqVarsL(:,i))./(2*h);
          end
          errJacNum(:,i) = (errR-errL)./(2*h);
        end

        eqVarsJacErr = eqVarsJac - eqVarsNumJac;
        eqVarsJacErrNorm = sqrt(sum(sum(eqVarsJacErr.^2)));

        errJacErr = errFJac - errJacNum;
        errJacErrNorm = sqrt(sum(sum(errJacErr.^2)));
        
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
      [errF,errFJac,errI,errIJac,modelCachedValuesUpd] = ...
        calcEquilibriumErrorOpus31( vars,...
                                    modelCachedValues,...
                                    modelCurves,...
                                    modelConstants,...
                                    sarcomereProperties,...
                                    useElasticTendon,...
                                    flag_updatePositionLevel,...
                                    flag_evaluateJacobian, ...
                                    flag_evaluateDerivatives, ...
                                    flag_updateModelCache,...
                                    flag_evaluateInitializationFunctions);

      modelCachedValues = modelCachedValuesUpd;
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
% Calculate the stiffness & damping of titin segment 1
%%
k1HN  = modelCachedValues.k1HNN*(1/liN_li);
k1H   = fiso*k1HN; %N/m
k1N   = k1HN*0.5;
k1    = k1H*0.5;

d1HN           = modelCachedValues.beta1HNN*(1/liN_li);
d1H            = fiso*d1HN; %N/(m/s)

d1N = d1HN*0.5;
d1  = d1H*0.5;

%%
% Calculate the stiffness & damping of titin segment 2
%%
k2HN  = modelCachedValues.k2HNN*(1/liN_li);
k2H   = fiso*k2HN; %N/m
k2N   = k2HN*0.5;
k2    = k2H*0.5;

d2HN  = 0;
d2H   = 0; %N/(m/s)
d2N   = 0;
d2    = 0; %N/(m/s)

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

zf = z12+zEcm+zx;

kf = real(zf);
df = imag(zf);

%%
% Resolve it along the fiber
%%

Dalpha_Dlce = modelCachedValues.Dalpha_Dlce;
cosAlpha    = modelCachedValues.cosAlpha;
sinAlpha    = modelCachedValues.sinAlpha;

fceN  = modelCachedValues.fceN;

kfAT = kf*cosAlpha - fiso*fceN*sinAlpha*Dalpha_Dlce; 
dfAT = df*cosAlpha;      

%%
% Calculate the impedance of the whole muscle
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

mtInfo.extra = [modelCachedValues.tau,...
                modelCachedValues.ddlaHN_HillError,...
                modelCachedValues.ddlaHN_Damping,...
                modelCachedValues.ddlaHN_Tracking];
mtInfo.extraLabels={'$$\tau$$','Hill Tracking','Cycle Damping','CE Tracking'};

if(sum(isnan(mtInfo.state.derivative)) > 0)
   here = 1; 
end
