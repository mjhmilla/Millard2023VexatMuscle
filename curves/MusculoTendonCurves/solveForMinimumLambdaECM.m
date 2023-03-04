function [lambdaECM, forceLengthPevkIgdCurve, forceLengthIgpCurve] = ...
  solveForMinimumLambdaECM(fiberForceLengthData, ...
      normStiffnessHalfPevkIgd, loptPevkIgd, ...
      normHalfStiffnessIgp, loptIgp, ...
      normLengthZToT12, normLengthHalfMyosin,...
      normPassiveFiberLengthAtOneNormForce,flag_useOctave)
    
lambdaECM = 0.5;
delta = 0.5;

%This is a bit nasty: if lambdaECM is too low, then it is possible that
%it is not possible to construct a monotonic igp curve because its slope
%might be too gradual to meet the demands.
    
%normStiffnessHalfPevkIgd          = kNormPEVKIGD; 
%loptPevkIgd                       = lPEVK+lIGD;
%normHalfStiffnessIgp              = kNormIGPeRatio; 
%loptIgp                           = lIGP;
%normLengthZToT12                  = lT12;
%normLengthHalfMyosin              = lmH;
%normFiberLengthAtOnePassiveForce  = lceFpeFiso;

iterMax = 10;

[forceLengthPevkIgdCurve, forceLengthIgpCurve] = ...
    constructTitinCurves(lambdaECM, ...
                    normStiffnessHalfPevkIgd, loptPevkIgd, ...
                    normHalfStiffnessIgp, loptIgp, ...
                    normLengthZToT12, normLengthHalfMyosin,...
                    normPassiveFiberLengthAtOneNormForce,... 
                    flag_useOctave);

forceLengthPevkIgdInverseCurve  = createInverseCurve(forceLengthPevkIgdCurve); 
forceLengthIgpInverseCurve      = createInverseCurve(forceLengthIgpCurve);

forceError = calcTitinCurveError( forceLengthPevkIgdInverseCurve, ...
                                  forceLengthIgpInverseCurve, ...
                                  fiberForceLengthData,...
                                  normLengthZToT12, normLengthHalfMyosin);                                

err = max(forceError);                                
                                
for i=1:1:iterMax
  %Left
  [forceLengthPevkIgdCurve, forceLengthIgpCurve] = ...
    constructTitinCurves(lambdaECM-delta, ...
                    normStiffnessHalfPevkIgd, loptPevkIgd, ...
                    normHalfStiffnessIgp, loptIgp, ...
                    normLengthZToT12, normLengthHalfMyosin,...
                    normPassiveFiberLengthAtOneNormForce,... 
                    flag_useOctave);

  %Check for monotonicity in the constructed curves: this will be broken
  %if lambdaECM is set to an infeasble value given the terminal stiffness of
  %the Pevk-IgD and IgP regions. 
  dxPevkIgd = diff(forceLengthPevkIgdCurve.xpts);
  dyPevkIgd = diff(forceLengthPevkIgdCurve.ypts);
  dxIgp     = diff(forceLengthIgpCurve.xpts);
  dyIgp     = diff(forceLengthIgpCurve.ypts);
  
  flag_curvesInvalid = 0;
  if(min(dxPevkIgd) < 0 || min(dyPevkIgd) < 0 || min(min(dxIgp)) < 0 || min(min(dyIgp)) < 0)
    flag_curvesInvalid=1;
  end

  errLeft = Inf;
  if(flag_curvesInvalid == 0)
    forceLengthPevkIgdInverseCurve = createInverseCurve(forceLengthPevkIgdCurve); 
    forceLengthIgpInverseCurve = createInverseCurve(forceLengthIgpCurve);

    forceErrorLeft = calcTitinCurveError( forceLengthPevkIgdInverseCurve, ...
                                      forceLengthIgpInverseCurve, ...
                                      fiberForceLengthData,...
                                      normLengthZToT12, normLengthHalfMyosin);                                
    errLeft = max(forceErrorLeft);
  end
  %Right
  [forceLengthPevkIgdCurve, forceLengthIgpCurve] = ...
    constructTitinCurves(lambdaECM+delta, ...
                    normStiffnessHalfPevkIgd, loptPevkIgd, ...
                    normHalfStiffnessIgp, loptIgp, ...
                    normLengthZToT12, normLengthHalfMyosin,...
                    normPassiveFiberLengthAtOneNormForce,... 
                    flag_useOctave);

  dxPevkIgd = diff(forceLengthPevkIgdCurve.xpts);
  dyPevkIgd = diff(forceLengthPevkIgdCurve.xpts);
  dxIgp     = diff(forceLengthIgpCurve.xpts);
  dyIgp     = diff(forceLengthIgpCurve.xpts);
  
  flag_curvesInvalid = 0;
  if(min(dxPevkIgd) < 0 || min(dyPevkIgd) < 0 || min(min(dxIgp)) < 0 || min(min(dyIgp)) < 0)
    flag_curvesInvalid=1;
  end

  errRight = Inf;
  if(flag_curvesInvalid == 0)
    forceLengthPevkIgdInverseCurve = createInverseCurve(forceLengthPevkIgdCurve); 
    forceLengthIgpInverseCurve = createInverseCurve(forceLengthIgpCurve);

    forceErrorRight = calcTitinCurveError( forceLengthPevkIgdInverseCurve, ...
                                      forceLengthIgpInverseCurve, ...
                                      fiberForceLengthData,...
                                      normLengthZToT12, normLengthHalfMyosin);                                

    errRight=max(forceErrorRight);
  end
  
  if( errLeft < 0 &&  abs(errLeft) < abs(err) || ((err > 0) && (errLeft < err)) ) 
    err = errLeft;
    lambdaECM = lambdaECM-delta;
  end
  if( errRight < 0 &&  abs(errRight) < abs(err) || ((err > 0) && (errRight < err)) )
    err = errRight;
    lambdaECM = lambdaECM+delta;
  end
  delta = delta*0.5;                                
end

[forceLengthPevkIgdCurve, forceLengthIgpCurve] = ...
  constructTitinCurves(lambdaECM, ...
                  normStiffnessHalfPevkIgd, loptPevkIgd, ...
                  normHalfStiffnessIgp, loptIgp, ...
                  normLengthZToT12, normLengthHalfMyosin,...
                  normPassiveFiberLengthAtOneNormForce,... 
                  flag_useOctave);





