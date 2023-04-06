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
function [errF, errFJac, errI, errIJac, modelCache] ...
          = updateMillard2023VexatCache(...
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
                                flag_evaluateInitializationFunctions    ,...
                                flag_useArgs)

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
% |===========[|]====================|                  |   |   
% |---- l1 --->|                                        |   |   
% |---------------------|\/\/\/\/\/\|-------------------|---| fecm
% |---------------------[           ]-------------------|---| becm
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
calcCpDer   = modelCurves.calcCpDer  ;
calcFpeDer  = modelCurves.calcFpeDer ;
%Constants

fiso        = modelConstants.fiso        ;
lceOpt      = modelConstants.lceOpt      ;
alphaOpt    = modelConstants.alphaOpt    ;
ltSlk       = modelConstants.ltSlk       ;

lceHNZeroFpeN= modelConstants.lceHNZeroFpeN;

dlceMaxN    = modelConstants.dlceMaxN    ;
lceMin      = modelConstants.lceMin      ;

kAXHN        = modelConstants.kAXHN        ;
betaAXHN     = modelConstants.betaAXHN     ;

betaCXHN  = modelConstants.normCrossBridgeCyclingDamping   ;

betaNum     = modelConstants.betaNum     ;
betafTN     = modelConstants.betafTN      ;
betaCTN     = modelConstants.betaCTN; 
betafEcmHN   = modelConstants.betafEcmHN    ;

betaTAaHN  = modelConstants.betaTAaHN  ;
betaTApHN  = modelConstants.betaTApHN  ;

forceVelocityCalibrationFactor = modelConstants.forceVelocityCalibrationFactor;
tau                     = modelConstants.tau       ;
lowActivationThreshold  = modelConstants.lowActivationThreshold;
lowActivationGain       = modelConstants.lowActivationGain;

lTitinFixedHN = modelConstants.lTitinFixedHN      ;
lmHN          = modelConstants.lmHN        ;
lmH           = modelConstants.lmH         ;
ZLineToT12NormLengthAtOptimalFiberLength = ...
  modelConstants.ZLineToT12NormLengthAtOptimalFiberLength;

smoothStepFunctionRadius = ... 
  modelConstants.smoothStepFunctionRadius;  

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
scaleTitinDistal      = modelConstants.scaleTitinDistal;
scaleTitinProximal    = modelConstants.scaleTitinProximal;

lAHN = modelConstants.lAHN; %Normalized actin length

betaPevkPHN = modelConstants.betaPevkPHN;
betaPevkAHN = modelConstants.betaPevkAHN; 
activationThresholdTitin = ...
  modelConstants.activationThresholdTitin;

titinModelType=modelConstants.titinModelType;

initialization = modelConstants.initialization;

DftfcnN_DltN_max = modelConstants.DftfcnN_DltN_max;


%%
% Initialization Routines
%%

if(flag_evaluateInitializationFunctions > 0)
  %fpeN = calcFpeDer(lceN,0);
  
  %errI    = zeros(2+flag_useElasticTendon,1);
  %errIJac = zeros(2+flag_useElasticTendon,2+flag_useElasticTendon); %Numerical for now
  errI = zeros(flag_useElasticTendon+1,1).*nan;  

  if(flag_useElasticTendon == 1)

    %Unknowns: lce, dlaH, laH, l1H
    
    % lce: Static elastic solution between CE and the tendon 
    %
    %1. Assuming that the X-bridge forces are equal to the Hill model
    %2. By using the fact that the fpe curve is the sum of the ecm and titin
    %   curves when the titin curve is unbound from actin.
    assert(length(args)==2);


    switch flag_evaluateInitializationFunctions
        case 1

            lceAT     = args(1,1);

            %Positions: CE length, pennation angle, and tendon length
            fiberKinematics = calcFixedWidthPennatedFiberKinematics(lceAT,0,lceOpt,alphaOpt);
        
            lce   = fiberKinematics.fiberLength;
            lceH  = lce*lce_lceH;
            lceN  = lce*lce_lceN;
            lceHN = lceH*lce_lceN;
            alpha = fiberKinematics.pennationAngle;
      
            sinAlpha  = sin(alpha);
            cosAlpha  = cos(alpha);
            lceATN    = lce*cosAlpha*lce_lceN;

            lp        = modelCache.lp;   
            dlp       = modelCache.dlp;
            lt        = lp-lce*cosAlpha;
            ltN       = lt*lt_ltN;
        
            %Length Dependent Forces
            fEcmfcnHN   = scaleECM*calcFeHDer(lceHN,0);
            fEcmkHN     = fEcmfcnHN;   

%            fpeN=calcFpeDer(lceN,0);
            
            fCpN=calcCpDer(lceATN,0);

            fTfcnN      = calcFtDer(ltN,0);            
            fTkN        = fTfcnN;
            betaTNN     = calcDtDer(ltN,0)*betafTN + DftfcnN_DltN_max*betaCTN + betaNum;


            %errI(1,1) = (modelCache.a*calcFalDer(lceN,0)+calcFpeDer(lceN,0))*cosAlpha ...
            %            - calcFtDer(ltN,0) - calcCpDer(lceATN,0);
            %-(fxHN + f2HN + fEcmHN)*cosAlpha + fCpN + fTN;
            %
            % fxHN  = kxHNN*lxHN + betaxHNN*dlxHN;

            %Approximate for now: lx=0, dlx=0;
            a     = modelCache.a;

            %Iterate over dlfNN until a value is found s.t. dlx=0
            solnBest = struct('dlceAT',dlp*0.5,'lce',0,'laH',0,'dlaH',0,'err',-1);
            solnL    = struct('dlceAT',0,'lce',0,'laH',0,'dlaH',0,'err',-1);
            solnR    = struct('dlceAT',0,'lce',0,'laH',0,'dlaH',0,'err',-1);

            %This loop iterates over the length of the titin segments

            numMaxBisections=modelConstants.iterMax;
            lTitinFreeH = (lceHN-lTitinFixedHN)*lceN_lce;
            argBest = 0.5;
            errBest =-Inf;
            argL    = 0;
            errL    = Inf;
            argR    = 0;
            errR    = Inf;

            for i=1:1:numMaxBisections
              stepDir = 2;

              if(i==1)
                delta = 0;
                stepDir=1;
              end 
              if(i==2)
                delta=0.25;                
              end

              errL=inf;
              errR=inf;

              for j=1:1:stepDir
                stepSign=0;
                switch j 
                  case 1 
                    stepSign=1;
                  case 2
                    stepSign=-1;
                  otherwise
                    assert(0);
                end
                arg = (argBest+stepSign*delta);
                l1H       = lTitinFreeH*arg;
                l2H       = lce*lce_lceH - (l1H+lTitinFixedHN*lceOpt);
                l1HN      = l1H*li_liN;
                l2HN      = l2H*li_liN;
                f1kHN     = scaleTitinProximal*calcF1HDer(l1HN,0);
                f2kHN     = scaleTitinDistal*calcF2HDer(l2HN,0);    
            
                errSoln = abs(f2kHN-f1kHN);
                if(errBest<0)
                    errBest=errSoln;
                elseif(errSoln < errBest && stepSign < 0)
                    errL=errSoln;
                    argL=arg;
                elseif(errSoln < errBest && stepSign > 0)
                    errR=errSoln;
                    argR=arg;                    
                end
              end
              if(errL<errBest && errL<errR)
                errBest=errL;
                argBest=argL;
              end
              if(errR<errBest && errR<errL)
                errBest=errR;
                argBest=argR;
              end
              delta=delta*0.5;
            end

            arg       = argBest;
            l1H       = lTitinFreeH*arg;
            l2H       = lce*lce_lceH - (l1H+lTitinFixedHN*lceOpt);
            l1HN      = l1H*li_liN;
            l2HN      = l2H*li_liN;
            f1kHN     = scaleTitinProximal*calcF1HDer(l1HN,0);
            f2kHN     = scaleTitinDistal*calcF2HDer(l2HN,0);  


            %This loop iterates over dlfNN in an effort to find a 
            %velocity at which the model is moving at a constant
            %velocity: when the strain rate of the XE is zero,
            %and the acceleration (ddlaHN) is zero.
            for i=1:1:numMaxBisections
              stepDir = 2;

              if(i==1)
                delta = 0;
                stepDir=1;
              end 
              if(i==2)
                delta=dlp*0.25;                
              end

              solnL.err=inf;
              solnR.err=inf;

              for j=1:1:stepDir
                stepSign=0;
                switch j 
                  case 1 
                    stepSign=1;
                  case 2
                    stepSign=-1;
                  otherwise
                    assert(0);
                end

                dlceAT = solnBest.dlceAT + stepSign*delta;

                fibKin = calcFixedWidthPennatedFiberKinematics(...
                            lceAT,dlceAT,lceOpt,alphaOpt);

                dlce = fibKin.fiberVelocity;
                dlt  = dlp-dlceAT;

                %Normalized velocities
                dlceH   = dlce*lce_lceH;
                dlceHN  = dlceH*dlce_dlceN;  
                dltN    = dlt*lt_ltN;                  
                dlceNN  = dlce*dlce_dlceNN;
                
                %Cross bridge strain rates
                %I want a solution in which the cross-bridge strain rate is zero
                dlxH  = 0;
                dlxHN = dlxH*lce_lceN;                
                dlaH  = dlceH - dlxH;                

                %Related velocities
                dlaHN  = dlaH*dlce_dlceN;
                dlaNN  = dlaH*dlce_dlceNN*lceH_lce;
                dlfNN  = dlaNN;

                %Solve for the XE strain as if the acceleration equation
                %evalates to zero.
                %
                %  ddlaHN  = ((kxHNN*lxHN + betaxHNN*dlxHN)-a*flN*fvN)/tau
                %             -betaCXHN*dlaNN
                %             +lowActivationGain*exp(-ka*ka)*(lxHN + dlxHN)   
                %          = 0
                %  since we've set dlx=0, and have a candidate values for lceAT
                %  and dlceAT the only unknown in this equation is lx. Technically
                %  this function is nonlinear in lx since flN(lamN(lmH,laH))
                %
                %    lxHN*(kxHNN/tau + lowActivationGain*exp(-ka*ka)) = 
                %      -( betaxHNN*dlxHN - a*flN*fvN)/tau
                %      + betaCXHN*dlaNN
                %      -( lowActivationGain*exp(-ka*ka) )*dlxHN
                %
                fvN=calcFvDer(dlfNN*forceVelocityCalibrationFactor,0);

                lxHN      = fvN/kAXHN;

                %Fixed point interation: this converges very quickly
                for z=1:1:2
                    lxH       = lxHN*lceN_lce;
    
                    laH       = lceH-lmH-lxH;
                    lamH      = (lmH+laH);
                    lamN      = (2*lamH)*lce_lceN;
                    
                    flN              = calcFalDer(lamN,0); 
                    kxHNN            = a*flN*kAXHN;
                    betaxHNN         = a*flN*betaAXHN;               
                     
                    ka                  = (a/lowActivationThreshold);
    
                    t0 =  (kxHNN/tau + lowActivationGain*exp(-ka*ka));
                    t1 =  -(betaxHNN*dlxHN- (a*flN*fvN))/tau ...
                          + betaCXHN*dlaNN ...
                          -( lowActivationGain*exp(-ka*ka) )*dlxHN;
                    lxHN = t1/t0;                    
                end

                %Check


                %Solve for flN and related quantities
                lxH       = lxHN*lceN_lce;
                laH       = lceH-lmH-lxH;
                lamH      = (lmH+laH);
                lamN      = (2*lamH)*lce_lceN; 

                flN       = calcFalDer(lamN,0);
                kxHNN     = a*flN*kAXHN;
                betaxHNN  = a*flN*betaAXHN;

                A = -(kxHNN*lxHN ...
                        - (betaxHNN*dlaH/lceOpt) ...
                        + f2kHN ...
                        + fEcmfcnHN ...
                     )*cosAlpha ...
                     + fCpN;
              
                B = (fTkN + betaTNN*(dlp/ltSlk));
              
                C = - ((betaxHNN*0.5/lceOpt) + (betaNum + betafEcmHN*fEcmkHN)*(0.5/lceOpt))*cosAlpha ...
                    - (betaTNN/(cosAlpha*ltSlk));
                              
                errSoln =  abs(-C*dlce-(A+B));
                %MM 2022/11/8
                % I've also tried an error that is the absolute value
                % of ddlaHN, and the sum of the equilibrium error (above)
                % and the absolute value of ddlaHN: these lead to worse
                % solutions

                if(solnBest.err < 0)
                  solnBest.err    = errSoln;
                  solnBest.dlceAT = dlceAT;
                  solnBest.lce    = lce;
                  solnBest.dlaH   = dlaH;
                  solnBest.laH    = laH;
                else
                  if(errSoln<solnBest.err && stepSign > 0)
                    solnR.err    = errSoln;
                    solnR.dlceAT = dlceAT;
                    solnR.lce    = lce;
                    solnR.dlaH   = dlaH;
                    solnR.laH    = laH;                                       
                  elseif(errSoln<solnBest.err && stepSign < 0)
                    solnL.err    = errSoln;
                    solnL.dlceAT = dlceAT;
                    solnL.lce    = lce;
                    solnL.dlaH   = dlaH;
                    solnL.laH    = laH;
                  end
                end

              end

              if(solnL.err < solnBest.err && solnL.err < solnR.err)
                solnBest=solnL;
              end
              if(solnR.err < solnBest.err && solnR.err < solnL.err)
                solnBest=solnR;
              end
              delta=delta*0.5;
            end


            errI(1,1) = solnBest.err;
            modelCache.dlceAT = solnBest.dlceAT;
            modelCache.lce    = solnBest.lce;
            modelCache.laH    = solnBest.laH;
            modelCache.dlaH   = solnBest.dlaH;

            return;
        case 2
            % l1H:  Static elastic solution between the two titin segments
            lceAT           = args(1,1);
            fiberKinematics = calcFixedWidthPennatedFiberKinematics(lceAT,0,lceOpt,alphaOpt);        
            lce             = fiberKinematics.fiberLength;

            l1H       = args(2,1);
            l2H       = lce*lce_lceH - (l1H+lTitinFixedHN*lceOpt);
            l1HN      = l1H*li_liN;
            l2HN      = l2H*li_liN;
            f1kHN     = calcF1HDer(l1HN,0);
            f2kHN     = calcF2HDer(l2HN,0);    
        
            errI(2,1) = f2kHN-f1kHN;
            modelCache.l1H=l1H;
            return;
        otherwise
            assert(0);
    end

  else
%%
% Rigid Tendon Initialization
%%
    assert(length(args)==2);

    switch flag_evaluateInitializationFunctions
        case 1
            lp        = modelCache.lp;   
            dlp       = modelCache.dlp;
            lceAT     = lp-ltSlk;
            dlceAT    = dlp;
            fibKin = calcFixedWidthPennatedFiberKinematics(...
                lceAT,dlceAT,lceOpt,alphaOpt);
            lce    = fibKin.fiberLength;
            lceH   = lce*lce_lceH;
            dlce    = fibKin.fiberVelocity;
            if(lce <= lceMin)
                lce=lceMin;
                fibKin=calcFixedWidthPennatedFiberKinematicsAlongTendon(...
                                        lce,...
                                        0,...
                                        lceOpt,...
                                        alphaOpt);
                lceAT = fibKin.fiberLengthAlongTendon;
                dlceAT = 0;
                dlce   = 0;
            end
            dlceH = dlce*lce_lceH;
            %Assume the XE has a strain rate of zero
            dlxH  = 0;
            dlxHN = dlxH*lce_lceN;
            dlaH  = dlceH-dlxH;


            %As with the elastic tendon, assume that the 
            %XE forces perfectly balance the Hill forces in the 
            %acceleration equation
            dlaHN  = dlaH*dlce_dlceN;
            dlaNN  = dlaH*dlce_dlceNN*lceH_lce;
            dlfNN  = dlaNN;
            

            fvN     =calcFvDer(dlfNN*forceVelocityCalibrationFactor,0);
            lxHN    = fvN/kAXHN;
            a       = modelCache.a;
            %Fixed point interation: this converges very quickly
            for z=1:1:2
                lxH       = lxHN*lceN_lce;

                laH       = lceH-lmH-lxH;
                lamH      = (lmH+laH);
                lamN      = (2*lamH)*lce_lceN;
                
                flN              = calcFalDer(lamN,0); 
                kxHNN            = a*flN*kAXHN;
                betaxHNN         = a*flN*betaAXHN;               
                 
                ka                  = (a/lowActivationThreshold);

                t0 =  (kxHNN/tau + lowActivationGain*exp(-ka*ka));
                t1 =  -(betaxHNN*dlxHN- (a*flN*fvN))/tau ...
                      + betaCXHN*dlaNN ...
                      -( lowActivationGain*exp(-ka*ka) )*dlxHN;
                lxHN = t1/t0;                    
            end

            errI(1,1) = 0;
            modelCache.dlceAT = dlceAT;
            modelCache.lce    = lce;
            modelCache.lceAT  = lceAT;
            modelCache.laH    = laH;
            modelCache.dlaH   = dlaH;
            
        case 2 
            % l1H:  Static elastic solution between the two titin segments
            lceAT           = args(1,1);
            fiberKinematics = calcFixedWidthPennatedFiberKinematics(lceAT,0,lceOpt,alphaOpt);        
            lce             = fiberKinematics.fiberLength;
    
            l1H       = args(2,1);
            l2H       = lce*lce_lceH - (l1H+lTitinFixedHN*lceOpt);
            l1HN      = l1H*li_liN;
            l2HN      = l2H*li_liN;
            f1kHN     = scaleTitinProximal * calcF1HDer(l1HN,0);
            f2kHN     = scaleTitinDistal   * calcF2HDer(l2HN,0);    
        
            errI(2,1) = f2kHN-f1kHN;
            modelCache.l1H=l1H;        
            return;
    end



    %errI(1,1) = f1kHN-f2kHN; %This error term works for titin models 0 and 1
    %errI(2,1) = ddlaHN;    
  end

  
end

%%
%Cached Quantities 
%%
if(flag_evaluateInitializationFunctions == 0)
    
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
    
    
    f1H            = modelCache.f1H           ;
    f1HN           = modelCache.f1HN          ;
    f1kHN          = modelCache.f1kHN         ;
    f1dHN          = modelCache.f1dHN         ;
    k1HNN          = modelCache.k1HNN         ;
    beta1HNN       = modelCache.beta1HNN      ;
    
    f2H            = modelCache.f2H           ;
    f2HN           = modelCache.f2HN          ;
    f2kHN          = modelCache.f2kHN         ;
    f2dHN          = modelCache.f2dHN         ;
    k2HNN          = modelCache.k2HNN         ;
    beta2HNN       = modelCache.beta2HNN      ;
    
    fxHN           = modelCache.fxHN          ;
    kxHNN          = modelCache.kxHNN         ;
    betaxHNN       = modelCache.betaxHNN      ;
    
    fTfcnN         = modelCache.fTfcnN        ;
    DfTfcnN_DltN   = modelCache.DfTfcnN_DltN  ;
    fTkN           = modelCache.fTkN          ;
    fTdN           = modelCache.fTdN          ;
    fTN            = modelCache.fTN           ;
    kTNN           = modelCache.kTNN          ;
    betaTNN        = modelCache.betaTNN       ;
    
    fEcmfcnHN         = modelCache.fEcmfcnHN          ;
    DfEcmfcnHN_DlceHN = modelCache.DfEcmfcnHN_DlceHN  ;
    fEcmkHN           = modelCache.fEcmkHN            ;
    fEcmdHN           = modelCache.fEcmdHN            ;
    fEcmHN            = modelCache.fEcmHN             ;
    kEcmHNN           = modelCache.kEcmHNN            ;
    betaEcmHNN        = modelCache.betaEcmHNN         ;
    
    fCpN              = modelCache.fCpN;
    kCpN              = modelCache.kCpN;
    
    
    D_f2HN_D_dlce      = modelCache.D_f2HN_D_dlce       ;
    D_fxHN_D_dlce      = modelCache.D_fxHN_D_dlce       ;
    D_fEcmHN_D_dlce    = modelCache.D_fEcmHN_D_dlce       ;
    D_fTN_D_dlce       = modelCache.D_fTN_D_dlce       ;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Position-level quantities
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
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
    
    
      lxH   = lceH - (lmH+laH);
      lxHN  = lxH*lce_lceN;
    
      l1HN  = l1H*li_liN;    % Norm. length titin segment: z_line-PEVK/IG2 border
    
    
      l2H   = lceH  - (l1H+lTitinFixedHN*lceOpt);     % Length IG2 (PEVK/IG2) - m-line
      l2HN  = l2H*li_liN;    % Norm. length IG2
    
    
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
    
      f1kHN     = scaleTitinProximal * calcF1HDer(l1HN,0);
      f2kHN     = scaleTitinDistal   * calcF2HDer(l2HN,0);
      k1HNN     = scaleTitinProximal * calcF1HDer(l1HN,1);
      k2HNN     = scaleTitinDistal   * calcF2HDer(l2HN,1);
    
      fEcmfcnHN = scaleECM*calcFeHDer(lceHN,0);
    
      if(flag_useElasticTendon == 1)
        fTfcnN    = calcFtDer(ltN,0);
        DfTfcnN_DltN = calcFtDer(ltN,1);
      else
        fTfcnN    = NaN;
      end
          
      lceATN=lceAT*lce_lceN;
      fCpN      = calcCpDer(lceATN,0);
      kCpN      = calcCpDer(lceATN,1);   
    
      if(fCpN > 1e-6 && flag_evaluateInitializationFunctions==0)
          here=1;
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
    
      modelCache.f1kHN      = f1kHN    ; 
      modelCache.f2kHN      = f2kHN    ; 
      modelCache.k1HNN      = k1HNN    ; 
      modelCache.k2HNN      = k2HNN    ; 
    
      modelCache.fEcmfcnHN  = fEcmfcnHN  ; 
      modelCache.fTfcnN     = fTfcnN  ; 
      modelCache.DfTfcnN_DltN = DfTfcnN_DltN;

      modelCache.fCpN       = fCpN;
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Velocities
    %
    %   If the tendon is rigid, all element positions and velocities can be
    %   evaluated
    %
    %   If the tendon is elastic, the CE velocity has to be solved for by
    %   evaluating the force inbalance between the CE and the tendon
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%------------------------------------------------------------------------------
    % Sliding velocity between XE and actin 
    %   Important:  Due to the viscoelasticity of the XE this is not the same as 
    %               the velocity between actin and myosin.
    %%------------------------------------------------------------------------------
    
    dlaHN  = dlaH*dlce_dlceN;
    dlaNN  = dlaH*dlce_dlceNN*lceH_lce;
    dlfNN  = dlaNN;
    
    fvN=calcFvDer(dlfNN*forceVelocityCalibrationFactor,0);
    
    
    %%------------------------------------------------------------------------------
    % Velocities: CE and CE dependent 
    %%------------------------------------------------------------------------------
    
    %ECM - position level forces
    fEcmkHN     = fEcmfcnHN;
    betaEcmHNN  = betafEcmHN*fEcmkHN + betaNum;
    
    %Tendon force - position level forces
    if(flag_useElasticTendon == 1)
      fTkN         = fTfcnN;
      betaTNN      = calcDtDer(ltN,0)*betafTN + DftfcnN_DltN_max*betaCTN + betaNum;
      fTdN         = nan; %Depends on the CE velocity: cannot yet solve for this
      fTN          = nan; %Depends on the CE velocity: cannot yet solve for this
    
    else
      fTfcnN  = nan; %Depends on the CE velocity: cannot yet solve for this
      fTkN    = nan; %Depends on the CE velocity: cannot yet solve for this
      betaTNN = 0;
      fTdN    = nan; %Depends on the CE velocity: cannot yet solve for this
      fTN     = nan; %Depends on the CE velocity: cannot yet solve for this
    end
    
    %Solve for dlce
    if(flag_useElasticTendon == 1)
    
      assert(titinModelType==0,['Error: the method that explicitly solves ',...
                    'for the CE velocity has not been generalized: the ',...
                    'current formulation only works with the titin-actin',...
                    'interaction.']);  
    
      %Contractile element - tendon force balance
      %  dlce is the only unknown in the equation fce*cosAlpha-ft = 0
    
      %%
      % Assumed system equations
      %   dlceAT     = dlce/cosAlpha; 
      %
      %   dltN       = (dlp-dlce/cosAlpha)/ltSlk;
      %
      %   dlxHN      = (dlce*0.5-dlaH)/lceOpt;
      %
      %   fce0      = (fxHN + f2HN + fEcmHN)*cosAlpha - fCpN;
      %   fce1      = (kxHNN*lxHN + betaxHNN*dlxHN0 
      %                + f2HN 
      %                + (fEcmfcnHN + (betaNum + betafEcmHN*fEcmfcnHN)*dlceHN))*cosAlpha
      %                - fCpN
      %
      % Where fce1 is the expanded version of fce0
      %
      %   ft0  = fTkN + betaTNN*(dlp - dlce/cosAlpha)/ltSlk;
      %
      % and ft0 is the expanded version of ftN
      %
      % This results in an equilibrum equation of
      %
      % errF0 = -(kxHNN*lxHN ...
      %            - (betaxHNN*dlaH/lceOpt) ...
      %            + f2HN ...
      %            + fEcmfcnHN ...
      %         )*cosAlpha ...  
      %         + fCpN
      %         + (fTkN + betaTNN*(dlp/ltSlk)) ...
      %         - dlce*( (betaxHNN*0.5/lceOpt) + (betaNum + betafEcmHN*fEcmfcnHN)*(0.5/lceOpt))*cosAlpha ...  
      %         - dlce*(betaTNN/(cosAlpha*ltSlk))
      %         = 0
      %
      % which can be solved for dlce since it is the only unknown. 
      % 
      % Since the damping model for the ECM and the tendo both include betaNum, which is a small
      % value above zero, there is guaranteed to be a finite solution for dlce.
      %%
    
    
      A = -(kxHNN*lxHN ...
              - (betaxHNN*dlaH/lceOpt) ...
              + f2kHN ...
              + fEcmfcnHN ...
           )*cosAlpha ...
           + fCpN;
    
      B = (fTkN + betaTNN*(dlp/ltSlk));
    
      C = - ((betaxHNN*0.5/lceOpt) ...
          +  (betaNum + betafEcmHN*fEcmkHN)*(0.5/lceOpt) )*cosAlpha ...
          - (betaTNN/(cosAlpha*ltSlk));
    
      dlce = -(A+B)/C;
    
    else
    
      %Rigid tendon kinematic model
    
      dlceAT= dlp;
      fibKin = calcFixedWidthPennatedFiberKinematics(lceAT,...
                                          dlceAT,...
                                          lceOpt,...
                                          alphaOpt);
      dlce   = fibKin.fiberVelocity;
      
    end
    
    
    
    %%------------------------------------------------------------------------------
    % Update velocities dependent on the CE velocity
    %%------------------------------------------------------------------------------
    
    %Pennation angle and tendon
    fibKin  = calcFixedWidthPennatedFiberKinematicsAlongTendon(...
                    lce,dlce,lceOpt,alphaOpt);
    
    dlceAT = fibKin.fiberVelocityAlongTendon ;
    dalpha = fibKin.pennationAngularVelocity ;  
    dlt   = dlp - dlceAT;    
    
    %Normalized velocities
    dlceH   = dlce*lce_lceH;
    dlceHN  = dlceH*dlce_dlceN;  
    dltN    = dlt*lt_ltN;                  
    dlceNN  = dlce*dlce_dlceNN;
    
    %Cross bridge strain rates
    dlxH  = dlceH - (dlaH);
    dlxHN = dlxH*lce_lceN;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Forces
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%------------------------------------------------------------------------------
    % Update forces dependent on the CE velocity
    %%------------------------------------------------------------------------------
    
    
    
    %%------------------------------------------------------------------------------
    %Titin segments 
    %   1: lumped T12+IG1+Q(PEVK)                : z-line to the titin/actin bonding location
    %   2: (1-Q)PEVK section + compliance of IG2 : N2A epitope to the PEVK/IG2 border
    %%------------------------------------------------------------------------------
    
    %Note 1: this damping coefficient is the same in lengthening as it is
    %        in shortening. This might not be the case! Something interesting
    %        to test with an experiment
    
    beta1HNN = 0; 
    beta2HNN = 0;
    
    dl1HN = 0;
    dl1H  = 0;
    
    dl2H  = 0; 
    dl2HN = 0;
    
    %
    %
    %From Fukutani and Herzog it is known that the active properties of titin
    %appear to be saturated at very low levels of calcium concentration. Here
    %we use a function that smoothly saturates to 1 as activation becomes 
    %greater than activationThresholdTitin.
    %
    % Fukutani A, Herzog W. Residual Force Enhancement Is Preserved for Conditions 
    % of Reduced Contractile Force. Medicine and Science in Sports and Exercise. 
    % 2018 Jun 1;50(6):1186-91.
    %
    
    kaTi = a/activationThresholdTitin;
    aTi = 1.-exp((-kaTi*kaTi));
    
    switch titinModelType
      case 0
        %Sticky-Spring active titin model:
        %   Damping acts between actin and some point on the PEVK segment
        %
        %|     IgP       PEVK    IgD                           |
        %|---|\/\/\|--|\/\/\/|--------|------------------------|
        %|===========[b]====================|                  |
        %|---- l1 --->|                    lAHN                |   
        % 
    
        %uActin is a step function that is 1 provided the titin-actin attachement
        %point is still within actin, otherwise it goes to zero. Because this is
        %in normalized coordinates and its a step function that goes to zero when
        %the titin-actin bond slides off of titin I'm directly putting in this
        %
    
        dTiA = (lAHN-(l1HN+ZLineToT12NormLengthAtOptimalFiberLength));
        kTiA = dTiA/smoothStepFunctionRadius;
        uTiA = 0.5+0.5*tanh(kTiA);
    
        if(uTiA <= 0.1)
            here=1;
        end
    
        %
        % From Hisey et al. it is known that force enhancement does not
        % take place for act streches that remain less than ~ 1 optimal
        % fiber length. I'm going to assume here that this has nothing
        % to do with the optimal fiber length per se, but that this 
        % minimal length is due to the prox. Ig segment having a minimum
        % slack length beyond which it does not shorten.
        %
        % Hisey B, Leonard TR, Herzog W. Does residual force enhancement increase 
        % with increasing stretch magnitudes?. Journal of biomechanics. 
        % 2009 Jul 22;42(10):1488-92.
        %
    
        dTiLMin = lceHN-lceHNZeroFpeN;
        kTiLMin = dTiLMin/smoothStepFunctionRadius;
        uTiLMin = 0.5+0.5*tanh(kTiLMin);
        % To break beta1HNN down:
        %
        %   beta1HNN = betaTApHN + betaTAaHN*aTi*uTia
        %
        %   beta1HNN : titin-actin damping at the bond site
        %   betaTApHN: small passive damping coefficient that never goes to zero
        %   betaTAaHN: maximum active titin-actin damping coefficient
        %   aTi      : value between 0-1 that indicates if the bond is active
        %   uTia     : a value that is 1 provided the bond overlaps with actin,
        %              otherwise its zero.    
        %
        % Note: the active titin-actin bond is disabled for initialization
        beta1HNN = betaTApHN + betaTAaHN*aTi*uTiA*uTiLMin*(1-initialization); 
        beta2HNN = 0;
    
        dl1HN = (f2kHN-f1kHN)/beta1HNN;
        dl1H  = dl1HN*lceN_lce;
    
        dl2H  = dlceH - dl1H; 
        dl2HN = dl2H*lce_lceN;
    
        f1dHN = beta1HNN*dl1HN;
        f2dHN = 0;   
    
        f1HN  = f1kHN + f1dHN;    
        f2HN  = f2kHN;   
    
        assert(abs(f2dHN) < eps && abs(beta2HNN) < eps && abs(f2HN-f2kHN) < eps*10,...
          ['Error: Method used to directly solve dlce',...
            ' for an elastic tendon is no longer valid']);       
    
      case 1 
    
        %Stiff-spring model:
        %  Damping acts on the PEVK element
        %
        %|     IgP       PEVK    IgD                           |
        %|---|\/\/\|--|\/\/\/\|--------|------------------------|
        %|            |       |
        %|            |--[ ]--|
        %|                b
        %|---- l1 --->|                                        |   
        %
        % This is a very weird model: PEVK becomes very viscous when exposed
        % to cycling crossbriges within actin-myosin overlap, but does not lose
        % this viscous state until Ca2+ levels drop. Even when it goes outside
        % of the reach of actin.
        %
        % This model was made after simulating Leonard, Joumaa, and Herzog
        % and finding that a titin-actin bond would slip off well before
        % the lengths that Leonard, Joumaa, and Herzog measured. If one assumes
        % that there is no experimental error, and all assumptions are correct,
        % then a titin-actin bond couldn't support enhanced forces at these
        % long lengths. And so I made a model that could.
        %
        % However, reality is not so simple. There is ample indirect 
        % and direct evidence that actin and titin interact. The experiments
        % are also challenging to undertake and some room for error should be
        % given - especially now that I've worked with the Aurora machines
        % myself. Finally, there are a number of plausible options within
        % the sarcomere that could explain the long lengths: it's possible
        % according to Herzog that some of the distal Ig segment becomes 
        % separated from myosin. This would free up a lot of extra titin
        % to stretch.
        %
    
        beta1HNN = 0;
        beta2HNN = betaPevkPHN + betaPevkAHN*(1-exp(-aTi)); 
    
        dl2HN = (f1kHN-f2kHN)/beta2HNN;
        dl2H  = dl2HN*lceN_lce;
    
        dl1H  = dlceH - dl2H; 
        dl1HN = dl1H*lce_lceN;
    
        f1dHN = 0;              %beta1HNN*dl1HN;
        f2dHN = beta2HNN*dl2HN;
        f1HN  = f1kHN;          % + f1dHN;    
        f2HN  = f2kHN + f2dHN;
      otherwise
        assert(0,'titinModelType must be 0 (sticky-spring) or 1 (stiff spring');
    end
    
    %ECM
    fEcmdHN     = betaEcmHNN*dlceHN;
    fEcmHN      = fEcmkHN + fEcmdHN;
    
    %N.B. This function assumes that the duty cycle does not change with 
    %     velocity. There is a good chance that it does.
    fxHN  = kxHNN*lxHN + betaxHNN*dlxHN;
    
    %Tendon force
    if(flag_useElasticTendon == 1)
      fTdN         = betaTNN*dltN;
      fTN          = fTkN + fTdN;
    else
      fTdN    = 0; %It's rigid: its not stretching.
      fTN     = (fxHN + f2HN + fEcmHN)*cosAlpha -fCpN; 
    end
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Acceleration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    lambda              =   0;
    ddlaHN_HillError    =   ((fxHN - a*flN*(fvN))/(tau));
    ddlaHN_Damping      = - betaCXHN*dlaNN;
    
    ka                  = (a/lowActivationThreshold);
    ddlaHN_Tracking     =   lowActivationGain*exp(-ka*ka)*(lxHN + dlxHN);
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
    
    
    %Sign convention
    % + forces that pull to the right
    % - forces that pull to the left
    %Equlilbrium Eqn 1: applied to fiber-tendon junction
    %                   assuming both sarcomere halves are in equilibrium
    errF    = -(fxHN + f2HN + fEcmHN)*cosAlpha + fCpN + fTN;
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
            
        modelCache.dl2H        = dl2H       ; 
        modelCache.dl2HN       = dl2HN      ; 
    
        modelCache.dlxH        = dlxH       ; 
        modelCache.dlxHN       = dlxHN      ; 
                
        modelCache.dalpha      = dalpha     ; 
    
        
        %Force/Titin related quantities
    
        modelCache.f1dHN    = f1dHN;    
        modelCache.beta1HNN = beta1HNN;
        modelCache.f1HN     = f1HN;
        modelCache.f1H      = f1HN*fiso;
    
        modelCache.f2dHN    = f2dHN;    
        modelCache.beta2HNN = beta2HNN;
        modelCache.f2HN     = f2HN;
        modelCache.f2H      = f2HN*fiso;
    
    
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
        
        modelCache.fceN     = (fxHN + f2HN + fEcmHN) - (fCpN/cosAlpha);
    
        %Power quantities
    
        modelCache.dw_lp    =   fiso*fTN*dlp;               % W: work    
    
        switch titinModelType
          case 0
            modelCache.dw_12AH  = - fiso*(f2HN-f1HN)*dl1H;      % W: work
            modelCache.dw_1kH   =  (fiso*f1HN)*dl1H;            % V: potential
            modelCache.dw_2kH   =  (fiso*f2HN)*dl2H;            % V: potential
    
          case 1
    
            modelCache.dw_12AH  = - fiso*(f1HN-f2HN)*dl2H;      % W: work
            modelCache.dw_1kH   =  (fiso*f1HN)*dl1H ;           % V: potential
            modelCache.dw_2kH   =  (fiso*f2HN)*dl2H;            % V: potential 
    
        otherwise
          assert(0,'titinModelType must be 0 (sticky-spring) or 1 (stiff spring');
        end
    
        modelCache.dw_xkH   =  (fiso*kxHNN*lxHN)*dlxH;      % V: potential
        modelCache.dw_xdH   = -(fiso*betaxHNN*dlxHN)*dlxH;  % W: work
        modelCache.dw_xaH   = -(fiso*fxHN)*dlaH;            % W: work
        %modelCache.dw_1dH   = 0;                            % element has no damping
        %modelCache.dw_2dH   = 0;                            % element has no damping
        modelCache.dw_Tk    =  (fiso*fTkN)*dlt;             % V: potential
        modelCache.dw_Td    = -(fiso*fTdN)*dlt;             % W: work
        modelCache.dw_EcmkH =  (fiso*fEcmkHN)*dlceH;        % V: potential
        modelCache.dw_EcmdH = -(fiso*fEcmdHN)*dlceH;        % W: work
        modelCache.dw_Cp    =  (fiso*fCpN)*dlce;            % V: potential
    
        modelCache.dT = 0;
        modelCache.dV = 2*modelCache.dw_xkH ...
                        + 2*modelCache.dw_1kH ...
                        + 2*modelCache.dw_2kH ...
                        + 2*modelCache.dw_EcmkH ...
                        + modelCache.dw_Tk;
    
        modelCache.dW =   modelCache.dw_lp ...
                          + 2*modelCache.dw_12AH ...
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
    
    
    %err: -(fxHN + f2HN + feHN)*cosAlpha + ftN + fCp
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
    
    
    
    if(flag_updateModelCache==1 && flag_evaluateDerivatives)
    
    
      D_fEcmfcnHN_D_lceHN = scaleECM*calcFeHDer(lceHN,1 );
    
      D_fTfcnN_D_ltN = 0;
      if(flag_useElasticTendon)
        D_fTfcnN_D_ltN      = calcFtDer(   ltN,1 );
      else
        D_fTfcnN_D_ltN = Inf;
      end
    
      modelCache.D_f1HN_D_l1HN        = k1HNN;
      modelCache.D_f2HN_D_l2HN        = k2HNN;
      modelCache.D_fEcmfcnHN_D_lceHN  = D_fEcmfcnHN_D_lceHN;
      modelCache.D_fTfcnN_D_ltN       = D_fTfcnN_D_ltN;
    
      %Evaluate the stiffness and damping coefficients of each element
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
    
      modelCache.DfCpN_DlceN = kCpN;
    end
    
    
    if(flag_evaluateJacobian == 1)
    
      %d/dt lt          = d/dt lp - d/dt lceAt
      %     vt          = vp - dlceAT
      %   Ddlt_DdlceAT  = -1
      %
      %Deriving dlceAT
      %
      %     lceAT       = lce*cos(alpha)
      %    dlceAT       = dlce*cos(alpha) - lce*sin(alpha)*dalpha
      %
      %Deriving dalpha (as in d/dt alpha)
      %
      %    lce*sin(alpha) = h
      %   dlce*sin(alpha) + lce*cos(alpha)*dalpha = 0
      %   dalpha = -dlce*sin(alpha)/( lce*cos(alpha))
      %
      %Substituting into dlceAT
      %dlceAT       = dlce*cos(alpha) 
      %             - lce*sin(alpha)*( -dlce*sin(alpha)/( lce*cos(alpha)) )
      %
      %dlceAT       =  dlce*cos^2(alpha)/( cos(alpha)) )
      %              + dlce*sin^2(alpha)/( cos(alpha)) )
      %
      %dlceAT       =  dlce*(cos^2(alpha)+ sin^2(alpha))/( cos(alpha)) )
      %
      %dlceAT       = dlce/cos(alpha)
      %
      %  DdlceAT_Ddlce  =    1/cos(alpha)
      %
      Ddlt_DdlceAT  = 0;
      DdlceAT_Ddlce = 0;
      if(flag_useElasticTendon == 1)
        Ddlt_DdlceAT  = -1;
      else
        Ddlt_DdlceAT  = 0;
      end
      DdlceAT_Ddlce =  1/cosAlpha;
    
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
      %Distal Titin Segment 
      %%----------------------------------------------------------------------------
    
      %  f2kHN = calcF2HDer(l2HN,0)
      %
      % Stiff spring model
      %  f2dHN = (betaPevkP + (1-exp(a/aThresh))*betaPevkA)*dl2HN 
      %
      % Sticky spring model
      %  f2dHN = 0 
      %
      %  f2HN           = scaleTitinDistal*(f2kHN + f2dHN)
      %  D_f2HN_D_dl2HN = scaleTitinDistal*(
      %                     D_f2kHN_D_dl2HN
      %                    +D_f2dHN_D_dl2HN)
      %
      %  D_f2kHN_D_dl2HN = d/dl2HN (calcF2HDer(l2HN,0)) = 0
      %  D_f2dHN_D_dl2HN = d/dl2HN ((betaPevkP + (1-exp(a/aThresh))*betaPevkA))*dl2HN 
      %                  = ((betaPevkP + (1-exp(a/aThresh))*betaPevkA))
      %                  = beta                 
      %
      %  l2HN           = lceHN - l1N - (LIg)
      %  dl2HN_dt       = dlceHN - dl1N
      %  D_dl2HN_D_dlce = D_dlceHN_D_dlce
      %
      %  D_f2HN_D_dlce  = D_f2HN_D_l2HN*D_dl2HN_D_dlce;
      D_f2HN_D_dlce  = 0;  
    
    
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
end


