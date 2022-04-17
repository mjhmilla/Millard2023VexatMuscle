function [errF, errFJac, errI, errIJac, modelCache] ...
          = calcEquilibriumErrorOpus31(...
                                args                                    ,... 
                                modelCache                              ,...
                                modelCurves                             ,...
                                modelConstants                          ,...
                                sarcomereProperties                     ,...
                                flag_useElasticTendon                   ,...
                                flag_updatePositionLevel                ,...
                                flag_evaluateJacobian                   ,... 
                                flag_evaluateDerivatives                ,... 
                                flag_updateModelCache                   ,...
                                flag_evaluateInitializationFunctions)

%%              
%                                 
%                                                                 *pennated here
%         la             lx             lm (fixed)
% |---------------->|<-------->|<---------------------->|        
% |              
% |================[|]===============|                     assuming
%                   |     kx                               symmetry
%                   |--|/\/\|--|                        |   
%                   |     bx   |========================|   
%                   |--[    ]--|                        |        kt
%                                                       |   |--|/\/\|--|
% |     IgP       PEVK    IgD                           |---|    bt    |---
% |---|\/\/\|--|\/\/\/|--------|------------------------|   |--[    ]--|
% |===========[|]====================|                  |   |   |
% |---- l1 --->|                                        |   |   |
% |---------------------|\/\/\/\/\/\|-------------------|   |---| fecm
% |---------------------[           ]-------------------|   |---| becm
%
% |-------------------------1/2 lce --------------------|1/2 lce|-lt---------|  
%
%
% State vector : x = [lce,dla,la,l1] stretch goal: e1
%dlce    : fiber velocity                            (units: m/s)
% lce    : fiber length                              (units: m)
% la     : sliding length                            (units: m)
% lx     : cross bridge stretch                      (units: m)
% l1     : length of the titin segment between       (units: m)
%          the Z-line and the PEVK/distal IG border
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
% 
% N  : normalized
% H  : half
%


if(flag_useElasticTendon == 1)
  assert(length(args)==1);
else
  assert(length(args)==0);
end

errF        = 0;
errFJac     = 0;

errI      = zeros(2+flag_useElasticTendon,1);
errIJac   = zeros(2+flag_useElasticTendon,2+flag_useElasticTendon);



%Curves
calcFalDer  = modelCurves.calcFalDer;
%calcFpeDer  = modelCurves.calcFpeDer;
calcFvDer   = modelCurves.calcFvDer  ;

calcFeHDer  = modelCurves.calcFeHDer ;
calcF1HDer  = modelCurves.calcF1HDer ;
calcF2HDer  = modelCurves.calcF2HDer ;
calcFtDer   = modelCurves.calcFtDer  ;
calcDtDer   = modelCurves.calcDtDer  ;

%Constants

fiso        = modelConstants.fiso        ;
lceOpt      = modelConstants.lceOpt      ;
alphaOpt    = modelConstants.alphaOpt    ;
ltSlk       = modelConstants.ltSlk       ;

dlceMaxN    = modelConstants.dlceMaxN    ;
lceMin      = modelConstants.lceMin      ;

kAXHN        = modelConstants.kAXHN        ;
betaAXHN     = modelConstants.betaAXHN     ;

betaCXHN  = modelConstants.normCrossBridgeCyclingDamping   ;

betaNum     = modelConstants.betaNum     ;
betafTN     = modelConstants.betafTN      ;
betaCTN     = modelConstants.betaCTN; 
betafEcmHN   = modelConstants.betafEcmHN    ;

betaN2AaHN  = modelConstants.betaN2AaHN  ;
betaN2ApHN  = modelConstants.betaN2ApHN  ;

forceVelocityCalibrationFactor = modelConstants.forceVelocityCalibrationFactor;
tau                = modelConstants.tau       ;
lowActivationThreshold = modelConstants.lowActivationThreshold;


lTitinFixedHN = modelConstants.lTitinFixedHN      ;
lmHN          = modelConstants.lmHN        ;
lmH           = modelConstants.lmH         ;

lce_lceH    = modelConstants.lce_lceH    ;
lce_lceHN   = modelConstants.lce_lceHN   ;
dlce_dlceN  = modelConstants.dlce_dlceN  ;
dlce_dlceNN = modelConstants.dlce_dlceNN ;
lce_lceN    = modelConstants.lce_lceN    ;
lceN_lce    = modelConstants.lceN_lce    ;
lceH_lce    = modelConstants.lceH_lce    ;


lt_ltN      = modelConstants.lt_ltN      ;
li_liN      = modelConstants.li_liN      ;

ltN_lt      = modelConstants.ltN_lt      ;
liN_li      = modelConstants.liN_li      ;

scaleECM    = modelConstants.scaleECM;
scalePEVK   = modelConstants.scalePEVK;
scaleIGP    = modelConstants.scaleIGP;

DftfcnN_DltN_max = modelConstants.DftfcnN_DltN_max;


%%
%Cached Quantities 
%%


%Kinematic quantities
dlt            = modelCache.dlt           ;
lt             = modelCache.lt            ;
dltN           = modelCache.dltN          ;
ltN            = modelCache.ltN           ;
et             = modelCache.et            ;

lp             = modelCache.lp            ;
dlp            = modelCache.dlp           ;

lce            = modelCache.lce           ;
lceAT          = modelCache.lceAT         ;
lceN           = modelCache.lceN          ;
lceH           = modelCache.lceH          ;
lceHN          = modelCache.lceHN         ;

dlce           = modelCache.dlce          ;
dlceH          = modelCache.dlceH         ;
dlceHN         = modelCache.dlceHN        ;
dlceAT         = modelCache.dlceAT        ;

laH             = modelCache.laH          ;
dlaH            = modelCache.dlaH         ;
ddlaH           = modelCache.ddlaH        ;
laHN            = modelCache.laHN         ; 
dlaHN           = modelCache.dlaHN        ;
ddlaHN          = modelCache.ddlaHN       ;    

lambda         = modelCache.lambda       ;
dlambda        = modelCache.dlambda      ;

l1H            = modelCache.l1H           ;
l1HN           = modelCache.l1HN          ;
dl1H           = modelCache.dl1H          ;
dl1HN          = modelCache.dl1HN         ;

l2H            = modelCache.l2H           ;
l2HN           = modelCache.l2HN          ;
dl2H           = modelCache.dl2H          ;
dl2HN          = modelCache.dl2HN         ;

lxH            = modelCache.lxH           ;
lxHN           = modelCache.lxHN          ;
dlxH           = modelCache.dlxH          ;
dlxHN          = modelCache.dlxHN         ;

alpha          = modelCache.alpha         ;
cosAlpha       = modelCache.cosAlpha      ;
sinAlpha       = modelCache.sinAlpha      ;
dalpha         = modelCache.dalpha        ;

%Force quantities
a              = modelCache.a             ;
dadt           = modelCache.dadt          ;
lfN            = modelCache.lfN           ; %Fiber length excluding x-bridge strain
dlfNN          = modelCache.dlfNN         ; %Fiber velocity excluding x-bridge strain rate
flN            = modelCache.flN           ; 
fvN            = modelCache.fvN           ;


f1H            = modelCache.f1H              ;
f1HN           = modelCache.f1HN             ;
k1HNN          = modelCache.k1HNN;
beta1HNN       = modelCache.beta1HNN;

f2H            = modelCache.f2H              ;
f2HN           = modelCache.f2HN             ;
k2HNN          = modelCache.k2HNN;

fxHN           = modelCache.fxHN             ;
kxHNN          = modelCache.kxHNN             ;
betaxHNN       = modelCache.betaxHNN          ;

fTfcnN         = modelCache.fTfcnN           ;
DfTfcnN_DltN   = modelCache.DfTfcnN_DltN     ;
fTkN           = modelCache.fTkN             ;
fTdN           = modelCache.fTdN             ;
fTN            = modelCache.fTN              ;
kTNN           = modelCache.kTNN             ;
betaTNN        = modelCache.betaTNN          ;

fEcmfcnHN         = modelCache.fEcmfcnHN          ;
DfEcmfcnHN_DlceHN = modelCache.DfEcmfcnHN_DlceHN  ;
fEcmkHN           = modelCache.fEcmkHN            ;
fEcmdHN           = modelCache.fEcmdHN            ;
fEcmHN            = modelCache.fEcmHN             ;
kEcmHNN           = modelCache.kEcmHNN            ;
betaEcmHNN        = modelCache.betaEcmHNN         ;

D_f2HN_D_dlce      = modelCache.D_f2HN_D_dlce       ;
D_fxHN_D_dlce      = modelCache.D_fxHN_D_dlce       ;
D_fEcmHN_D_dlce    = modelCache.D_fEcmHN_D_dlce       ;
D_fTN_D_dlce       = modelCache.D_fTN_D_dlce       ;



%This gets a flag because these quantities only need to be updated
%once, and some of the quantities involve evaluating a Bezier curve
%or a trancendential function
if(flag_updatePositionLevel == 1)

  if(flag_useElasticTendon == 0)
    lceAT = lp - ltSlk;
    dlceAT= dlp;
    fibKin = calcFixedWidthPennatedFiberKinematics(lceAT,...
                                        dlceAT,...
                                        lceOpt,...
                                        alphaOpt);
    lce = fibKin.fiberLength;
    modelCache.lce = lce;
  end

  lceN   = lce*lce_lceN;  % Norm. fiber length  
  lceH   = lce*lce_lceH;  % Half fiber length
  lceHN  = lce*lce_lceHN; % Half norm. fiber length

if(lceN > 2.0)
    here=1;
end

  lxH   = lceH - (lmH+laH);
  lxHN  = lxH*lce_lceN;

  l1HN  = l1H*li_liN;    % Norm. length titin segment: z_line-PEVK/IG2 border



  l2H   = lceH  - (l1H+lTitinFixedHN*lceOpt);     % Length IG2 (PEVK/IG2) - m-line
  l2HN  = l2H*lce_lceN;    % Norm. length IG2


  laHN  = laH*lce_lceN;


%   if(lce < lceOpt*sin(alphaOpt)+1e-5)
%     here=1;
%   end
  
  fibKin  = calcFixedWidthPennatedFiberKinematicsAlongTendon(...
                  lce,NaN,lceOpt,alphaOpt);

  lceAT  = fibKin.fiberLengthAlongTendon   ;
  alpha  = fibKin.pennationAngle           ;

  cosAlpha = cos(alpha);
  sinAlpha = sin(alpha);

  lt    = lp - lceAT;
  ltN   = lt*lt_ltN;
  et    = ltN-1.0;     %Cauchy/engineering strain of the tendon

  lamH      = (lmH+laH);
  lamN      = (2*lamH)*lce_lceN;
  
  flN       = calcFalDer(lamN,0);
  kxHNN     = a*flN*kAXHN;
  betaxHNN  = a*flN*betaAXHN;
  f1HN      = scaleIGP*calcF1HDer(l1HN,0);
  f2HN      = scalePEVK*calcF2HDer(l2HN,0);
  fEcmfcnHN = scaleECM*calcFeHDer(lceHN,0);

  if(flag_useElasticTendon == 1)
    fTfcnN    = calcFtDer(ltN,0);
    DfTfcnN_DltN = calcFtDer(ltN,1);
  else
    fTfcnN    = NaN;
  end
                                
  modelCache.lceN     = lceN ;
  modelCache.lceH     = lceH ;
  modelCache.lceHN    = lceHN; 

  modelCache.lxH  = lxH ;
  modelCache.lxHN = lxHN;

  modelCache.l1HN = l1HN;
  modelCache.l2H  = l2H ;
  modelCache.l2HN = l2HN; 

  modelCache.laHN = laHN;

  modelCache.lceAT    = lceAT   ;
  modelCache.alpha    = alpha   ;
  modelCache.cosAlpha = cosAlpha;
  modelCache.sinAlpha = sinAlpha;

  modelCache.lt       = lt      ;
  modelCache.ltN      = ltN     ;
  modelCache.et       = et      ;   

  modelCache.lfN = lfN;

  modelCache.flN        = flN     ; 
  modelCache.kxHNN      = kxHNN    ; 
  modelCache.betaxHNN   = betaxHNN ; 
  modelCache.f1HN       = f1HN    ; 
  modelCache.f2HN       = f2HN    ; 
  modelCache.fEcmfcnHN  = fEcmfcnHN  ; 
  modelCache.fTfcnN     = fTfcnN  ; 
  modelCache.DfTfcnN_DltN = DfTfcnN_DltN;
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Kinematics
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%Fiber
%

%%------------------------------------------------------------------------------
%Pennated fiber and tendon
%%------------------------------------------------------------------------------
if(flag_useElasticTendon == 1)
  dlce    = args(1);
  fibKin  = calcFixedWidthPennatedFiberKinematicsAlongTendon(...
                  lce,dlce,lceOpt,alphaOpt);
  dlceAT = fibKin.fiberVelocityAlongTendon ;
  dalpha = fibKin.pennationAngularVelocity ;  
  dlt   = dlp - dlceAT;    

else
  dlceAT= dlp;
  fibKin = calcFixedWidthPennatedFiberKinematics(lceAT,...
                                      dlceAT,...
                                      lceOpt,...
                                      alphaOpt);
  dlce   = fibKin.fiberVelocity;
  dalpha = fibKin.pennationAngularVelocity;    
  dlt = 0;
end

dlceH   = dlce*lce_lceH;
dlceHN  = dlceH*dlce_dlceN;  
dltN    = dlt*lt_ltN;                  

%%------------------------------------------------------------------------------
% Sliding velocity
%%------------------------------------------------------------------------------

dlaHN  = dlaH*dlce_dlceN;
dlaNN  = dlaH*dlce_dlceNN*lceH_lce;
dlfNN  = dlaNN;
dlceNN = dlce*dlce_dlceNN;
%fvN    = calcFvDer(dlceNN,0);

fvN=calcFvDer(dlfNN*forceVelocityCalibrationFactor,0);

%%------------------------------------------------------------------------------
%Titin segments 
%   1: lumped T12+IG1                    : z-line to the N2A epitope
%   2: PEVK section + compliance of IG2  : N2A epitope to the PEVK/IG2 border
%%------------------------------------------------------------------------------

%Note 1: This implicitly assumes that the N2A epitope is never longer than
%        the length of an actin filament. This is probably quite a safe assumption
%        as the forces required to stretch the IG1 section (with a normalized  
%        resting length of 0.064) to the full length of actin would be enormous 
%        - greater than 10 maximum isometric forces.
%
%Note 2: this damping coefficient is the same in lengthening as it is
%        in shortening. This might not be the case! Something interesting
%        to test with an experiment
beta1HNN = betaN2ApHN + betaN2AaHN*a; 

dl1HN = (f2HN-f1HN)/beta1HNN;
dl1H  = dl1HN*lceN_lce;

dl2H  = dlceH - dl1H; 
dl2HN = dl2H*lce_lceN;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Net forces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ECM
fEcmkHN     = fEcmfcnHN;
betaEcmHNN  = betaNum + betafEcmHN*fEcmkHN;
fEcmdHN     = betaEcmHNN*dlceHN;
fEcmHN      = fEcmkHN + fEcmdHN;

%Cross bridges
dlxH  = dlceH - (dlaH);
dlxHN = dlxH*lce_lceN;

%N.B. This function assumes that the duty cycle does not change with 
%     velocity. There is a good chance that it does.
fxHN  = kxHNN*lxHN + betaxHNN*dlxHN;


lambda              =   0;
ddlaHN_HillError    =   ((fxHN  - a*flN*(fvN))/(tau));
ddlaHN_Damping      = - betaCXHN*dlaNN;

ka                  = (a/lowActivationThreshold);
ddlaHN_Tracking     =   exp(-ka*ka)*(lxHN + dlxHN);
ddlaHN = ddlaHN_HillError + ddlaHN_Damping + ddlaHN_Tracking;

%An alternate update rule for cross-bridge cycling. Note that since
%the zetax term is greater than the physical damping term in the model, this
%has the property of positive feedback: the more active lengthening takes place
%the more resistant to the movement of ls becomes. This has some positive
%aspects: simulations of Herzog & Leonard 2002 are improved, and the time
%constant used can be much bigger. This has two points which give me pause: 
%first, is there evidence of a positive feedback mechanism?; second, such a 
%mechanism cannot continue without restraint, so what happens during a really 
%long active stretch?.

%m      = 0.01;
%wnx    = sqrt(kxHNN/m);%10*2*pi;
%zetax  = sqrt(a*flN)*betaAXHN/(sqrt(kAXHN/m)*2);

%ddlaHN = -(a*flN*fvN/m + (2*zetax*wnx)*dlaHN) + (2*zetax*wnx)*dlxHN + wnx*wnx*lxHN...
%           + ddlaHN_Tracking; 

ddlaH  = ddlaHN*lceN_lce;



%Tendon force
if(flag_useElasticTendon == 1)
  fTkN         = fTfcnN;
  betaTNN      = calcDtDer(ltN,0)*betafTN + DftfcnN_DltN_max*betaCTN;
  fTdN         = betaTNN*dltN;
  fTN          = fTkN + fTdN;
else
  fTfcnN =   (fxHN + f2HN + fEcmHN)*cosAlpha;
  fTkN    = fTfcnN;
  betaTNN = 0;
  fTdN    = 0; %It's rigid: its not stretching.
  fTN     = fTkN + fTdN; 
end


%Sign convention
% + forces that pull to the right
% - forces that pull to the left
%Equlilbrium Eqn 1: applied to fiber-tendon junction
%                   assuming both sarcomere halves are in equilibrium
errF    = -(fxHN + f2HN + fEcmHN)*cosAlpha + fTN;
errFJac = NaN;


if(flag_updateModelCache == 1)

    %Kinematic quantities

    modelCache.dlce        = dlce       ; 
    modelCache.dlceH       = dlceH     ;
    modelCache.dlceHN      = dlceHN     ;
    modelCache.dlceAT      = dlceAT     ;

    modelCache.dlaHN       = dlaHN;
    modelCache.dlaNN       = dlaNN;

    modelCache.dlt         = dlt        ; 
    modelCache.dltN        = dltN       ; 

    modelCache.dl1H        = dl1H       ; 
    modelCache.dl1HN       = dl1HN      ; 
    
    modelCache.beta1HNN    = beta1HNN;
    
    modelCache.dl2H        = dl2H       ; 
    modelCache.dl2HN       = dl2HN      ; 

    modelCache.dlxH        = dlxH       ; 
    modelCache.dlxHN       = dlxHN      ; 

            
    modelCache.dalpha      = dalpha     ; 


    
    %Force/Titin related quantities

    modelCache.fEcmfcnHN    = fEcmfcnHN     ;
    modelCache.betaEcmHNN   = betaEcmHNN  ;
    modelCache.fEcmkHN      = fEcmkHN        ;          
    modelCache.fEcmdHN      = fEcmdHN        ; 
    modelCache.fEcmHN       = fEcmHN         ;     
    
    modelCache.fTfcnN   = fTfcnN     ;
    modelCache.betaTNN  = betaTNN   ;
    modelCache.fTkN     = fTkN       ; 
    modelCache.fTdN     = fTdN       ; 
    modelCache.fTN      = fTN        ; 
    
    modelCache.fxHN     = fxHN;
    modelCache.ddlaHN   = ddlaHN;
    modelCache.ddlaH    = ddlaH;

    modelCache.tau      = tau;
    modelCache.lambda   = lambda;
    modelCache.dlambda  = dlambda;
    modelCache.ddlaHN_HillError = ddlaHN_HillError;
    modelCache.ddlaHN_Damping   = ddlaHN_Damping;
    modelCache.ddlaHN_Tracking  = ddlaHN_Tracking;
    
    modelCache.fvN      = fvN;
    modelCache.dlfNN    = dlfNN;
    
    modelCache.fceN     = (fxHN + f2HN + fEcmHN);

    %Power quantities

    modelCache.dw_lp    =   fiso*fTN*dlp;               % W: work    
    modelCache.dw_n2aH  = - fiso*(f2HN-f1HN)*dl1H;      % W: work
    modelCache.dw_xkH   =  (fiso*kxHNN*lxHN)*dlxH;      % V: potential
    modelCache.dw_xdH   = -(fiso*betaxHNN*dlxHN)*dlxH;  % W: work
    modelCache.dw_xaH   = -(fiso*fxHN)*dlaH;            % W: work
    modelCache.dw_1kH   =  (fiso*f1HN)*dl1H;            % V: potential
    %modelCache.dw_1dH   = 0;                            % element has no damping
    modelCache.dw_2kH   =  (fiso*f2HN)*dl2H;            % V: potential
    %modelCache.dw_2dH   = 0;                            % element has no damping
    modelCache.dw_Tk    =  (fiso*fTkN)*dlt;             % V: potential
    modelCache.dw_Td    = -(fiso*fTdN)*dlt;             % W: work
    modelCache.dw_EcmkH =  (fiso*fEcmkHN)*dlceH;        % V: potential
    modelCache.dw_EcmdH = -(fiso*fEcmdHN)*dlceH;        % W: work

    modelCache.dT = 0;
    modelCache.dV = 2*modelCache.dw_xkH ...
                    + 2*modelCache.dw_1kH ...
                    + 2*modelCache.dw_2kH ...
                    + 2*modelCache.dw_EcmkH ...
                    + modelCache.dw_Tk;

    modelCache.dW =   modelCache.dw_lp ...
                      + 2*modelCache.dw_n2aH ...
                      + 2*modelCache.dw_xdH ...
                      + 2*modelCache.dw_xaH ...
                      + 2*modelCache.dw_EcmdH ...
                      + modelCache.dw_Td;


end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Derivatives
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fibKinDer = calcFixedWidthPennationPartialDerivatives(alpha,...
                                                      lce,...
                                                      lceOpt,...
                                                      alphaOpt);
Dalpha_Dlce   = fibKinDer.Dalpha_Dlce;
Ddalpha_Ddlce = fibKinDer.Ddalpha_Ddlce;


%err: -(fxHN + f2HN + feHN)*cosAlpha + ftN
%
% To evaluate the Jacobian we need
%
% DfxHN_Ddlce  = DfxHN_DlxHN * DlxHN_DlceN * DlceN_Dlce 
%              + DfxHN_DdlxHN * DdlxHN_DdlceN * DdlceN_Ddlce 
% Df2HN_Ddlce  = Df2HN_Dl2HN * Dl2HN_DlceN * DlceN_Dlce
% DfeHN_Ddlce  = DfeHN_DlceHN * DlceHN_Dlce
%               +DfeHN_DdlceHN * DdlceHN_Ddlce
% Ddalpha_Ddlce= from the pennation model
% DftN_Ddlce   = DftN_DltN*DltN_DlceN*DlceN_Dlce
%                DftN_DdltN*DdltN_DdlceN*DdlceN_Ddlce
%

D_f2HN_D_l2HN = 0;
if(flag_updateModelCache==1 && flag_evaluateDerivatives || flag_evaluateJacobian)
  D_f2HN_D_l2HN     = scalePEVK*calcF2HDer( l2HN,1 );
end


if(flag_updateModelCache==1 && flag_evaluateDerivatives)

  %Evaluate the necessary curve derivatives
  D_f1HN_D_l1HN       = scaleIGP*calcF1HDer(l1HN,1);
  D_fEcmfcnHN_D_lceHN = scaleECM*calcFeHDer(lceHN,1 );

  D_fTfcnN_D_ltN = 0;
  if(flag_useElasticTendon)
    D_fTfcnN_D_ltN      = calcFtDer(   ltN,1 );
  else
    D_fTfcnN_D_ltN = Inf;
  end
  modelCache.D_f2HN_D_l2HN        = D_f2HN_D_l2HN;
  modelCache.D_fEcmfcnHN_D_lceHN  = D_fEcmfcnHN_D_lceHN;
  modelCache.D_fTfcnN_D_ltN       = D_fTfcnN_D_ltN;

  %Evaluate the stiffness and damping coefficients of each element
  k1HNN  = D_f1HN_D_l1HN;
  k2HNN  = D_f2HN_D_l2HN; 

  DfEcmkHN_DlceHN     = D_fEcmfcnHN_D_lceHN;   
  DfEcmdHN_DlceHN     = (betafEcmHN*DfEcmkHN_DlceHN)*dlceHN;  
  DfEcmHN_DlceHN      = DfEcmkHN_DlceHN + DfEcmdHN_DlceHN;  
  kEcmHNN = DfEcmkHN_DlceHN;

  %Cross bridges: coefficents has already been computed
  % kXHN
  % betaXHN

  %The sliding element is strictly an active element. Since its force 
  %does not vary with the instantaneous length and velocity and velocity
  %of the contractile element, but instead with la and dla, it has no
  %stiffness and damping


  DfTkN_DltN      = D_fTfcnN_D_ltN;
  %DfTdN_DltN      = D_fTfcnN_D_ltN*betafTN*dltN;
  
  %betaTNN       = calcDtDer(ltN,0)*betafTN + betaCTN;
  DfTdN_DltN     = ( calcDtDer(ltN,1)*betafTN )*dltN; 
  
  DfTN_DltN       = DfTdN_DltN + DfTkN_DltN;
  kTNN      = DfTkN_DltN;
  %betaTNN  already been evaluated

  %Derivatives: for now reporting a minimal amount.
  modelCache.Dalpha_Dlce        = Dalpha_Dlce   ;
  modelCache.Ddalpha_Ddlce      = Ddalpha_Ddlce ;

  %lt = lp - lce*cos(alpha)
  %Dlt_Dlce = -cos(alpha) + lce*sin(alpha)*Dalpha_Dlce  
  Dlt_Dlce    = 0;
  DltN_DlceN  = 0;  
  if(flag_useElasticTendon == 1)
    Dlt_Dlce    = -cosAlpha + lce*sinAlpha*Dalpha_Dlce;
    DltN_DlceN  = Dlt_Dlce*lt_ltN*lceN_lce;
  else
    Dlt_Dlce   = 0;
    DltN_DlceN = 0;
  end
  modelCache.Dlt_Dlce         = Dlt_Dlce;
  modelCache.DltN_DlceN       = DltN_DlceN;

  modelCache.k1HNN             = k1HNN;
  modelCache.k2HNN             = k2HNN;

  modelCache.DfEcmfcnHN_DlceHN   = D_fEcmfcnHN_D_lceHN ;
  modelCache.DfEcmkHN_DlceHN     = DfEcmkHN_DlceHN   ;
  modelCache.DfEcmdHN_DlceHN     = DfEcmdHN_DlceHN   ;
  modelCache.DfEcmHN_DlceHN      = DfEcmHN_DlceHN    ;
  modelCache.kEcmHNN             = kEcmHNN;
  modelCache.betaEcmHNN          = betaEcmHNN;

  modelCache.kxHNN    = kxHNN;
  modelCache.betaxHNN = betaxHNN;

  modelCache.DfTfcnN_DltN     = DfTfcnN_DltN  ;
  modelCache.DfTkN_DltN       = DfTkN_DltN    ;
  modelCache.DfTdN_DltN       = DfTdN_DltN    ;
  modelCache.DfTN_DltN        = DfTN_DltN     ;
  modelCache.kTNN             = kTNN;
  modelCache.betaTNN          = betaTNN;  


end


if(flag_evaluateJacobian == 1)

  %d/dt lt          = d/dt lp - d/dt lceAt
  %     vt          = vp - dlceAT
  %   Ddlt_DdlceAT  = -1
  %     lceAT       = lce*cos(alpha)
  %    dlceAT       = dlce*cos(alpha) - lce*sin(alpha)*dalpha
  %  DdlceAT_Ddlce  =    1*cos(alpha) - lce*sin(alpha)*Ddalpha_Ddlce
  Ddlt_DdlceAT  = 0;
  DdlceAT_Ddlce = 0;
  if(flag_useElasticTendon == 1)
    Ddlt_DdlceAT  = -1;
  else
    Ddlt_DdlceAT  = 0;
  end
  DdlceAT_Ddlce =  1*cosAlpha - lce*sinAlpha*Ddalpha_Ddlce;

  D_dltN_D_dlce = lt_ltN*(Ddlt_DdlceAT*DdlceAT_Ddlce);

  D_dlceH_D_dlce    = 0.5; 
  D_dlceN_D_dlce    = dlce_dlceN;
  D_dlceHN_D_dlce   = D_dlceH_D_dlce*D_dlceN_D_dlce; 

  %%----------------------------------------------------------------------------
  %Cross-bridges
  %%----------------------------------------------------------------------------
  %
  % fxHN          = kx*lxHN + betaXHN*dlxHN
  % D_fxHN_D_dlce =         + betaXHN*D_dlxHN_D_dlce
  % lxHN  = lceHN - (lmHN + laHN)
  % dlxHN = dlceHN - (0  + dlaHN)
  % D_dlxHN_D_dlce = D_dlceHN_D_dlce
  D_dlxHN_D_dlce = D_dlceHN_D_dlce;
  D_fxHN_D_dlce = betaxHNN*D_dlxHN_D_dlce;


  %%----------------------------------------------------------------------------
  %Titin segment between the N2A epitope and the M-line
  %%----------------------------------------------------------------------------
  %
  % f2HN          = scalePEVK*calcF2HDer(l2HN,0)
  % D_f2HN_D_l2HN = scalePEVK*calcF2HDer(l2HN,1)
  % 
  % l2HN           = lceHN - l1N
  % dl2HN          = dlceHN - dl1N
  % D_dl2HN_D_dlce = D_dlceHN_D_dlce
  %
  %D_dl2HN_D_dlce = D_dlceHN_D_dlce;
  %D_f2HN_D_dlce  = D_f2HN_D_l2HN*D_dl2HN_D_dlce;
  D_f2HN_D_dlce = 0; %This is a spring: there is no damping
  
  %%----------------------------------------------------------------------------
  % ECM Force Velocity Derivatives
  %%----------------------------------------------------------------------------
  %fekHN       = fefcnHN;
  %betaEcmHN   = betaNum + betafEcmHN*fekHN;
  %fedHN       = betaEcmHN*dlceHN;
  %feHN        = fekHN + fedHN;
  %
  % D_feHN_D_dlce = D_fekHN_D_dlce + D_fedHN_D_dlce
  % D_fekHN_D_dlce = 0
  % D_fedHN_D_dlce = betaEcmHN*D_dlceHN_D_dlce

  D_fEcmHN_D_dlce = betaEcmHNN*D_dlceHN_D_dlce;

  %%----------------------------------------------------------------------------
  %Tendon Velocity Derivatives
  %%----------------------------------------------------------------------------
  % ftkN         = ftfcnN;
  % betaTNN      = calcDtDer(ltN,0)*betafTN + betaCTN;
  % fTdN         = betaTNN*dltN;
  % fTN          = fTkN + fTdN;
  %
  % D_ftN_D_dlce = D_ftkN_D_dlt*D_dlt_D_dlce->0 + D_ftdN_D_dlt*D_dlt_D_dlce
  % D_ftN_D_dlt  = betaTNN*D_dltN_D_dlt
  %
  %

  D_fTN_D_dlce = betaTNN*D_dltN_D_dlce;  

  
  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Error Derivatives
  %%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Sign convention
  % + forces that pull to the right
  % - forces that pull to the left  
  %Equlilbrium Eqn 2: (x-bridge + Titin 2) = Tendon
  %errF   = -( fxN + f2HN + feHN + betaNum*dlceH)*cosAlpha + ftN;
  %errJac = -( D_fxN_D_dlce + D_f2HN_D_dlce + D_feHN_D_dlce )*cosAlpha 
  %         +  D_ftN_D_dlce

  errFJac  = - (D_fxHN_D_dlce + D_f2HN_D_dlce + D_fEcmHN_D_dlce)*cosAlpha ...
              + D_fTN_D_dlce;

      

  if(flag_updateModelCache==1)

    modelCache.D_fxHN_D_dlce    = D_fxHN_D_dlce ; 
    modelCache.D_f2HN_D_dlce    = D_f2HN_D_dlce ;                              
    modelCache.D_fEcmHN_D_dlce  = D_fEcmHN_D_dlce ;                                 
    modelCache.D_fTN_D_dlce     = D_fTN_D_dlce;                                                                
    
  end
 

end

if(flag_evaluateInitializationFunctions==1)
  %fpeN = calcFpeDer(lceN,0);
  
  errI    = zeros(2+flag_useElasticTendon,1);
  errIJac = zeros(2+flag_useElasticTendon,2+flag_useElasticTendon); %Numerical for now

  if(flag_useElasticTendon == 1)
    errI(1,1) = -(fxHN + f2HN + fEcmHN)*cosAlpha + fTN; 
    %When f1HN=f2HN the passive force developed by fEcm + f2 = fpe;  
    errI(2,1) = f1HN-f2HN;
    errI(3,1) = ddlaHN;
  else
    errI(1,1) = f1HN-f2HN;
    errI(2,1) = ddlaHN;    
  end


end
