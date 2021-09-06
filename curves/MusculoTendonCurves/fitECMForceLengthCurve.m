function ecmCurve = fitECMForceLengthCurve(...
              fiberForceLengthData,...
              fiberForceLengthCurve, fiberForceLengthInverseCurve, ...
              forceLengthPevkIgdCurve,forceLengthPevkIgdInverseCurve,...
              forceLengthIgpCurve, forceLengthIgpInverseCurve,...
              normLengthZToT12, normLengthHalfMyosin,...
              flag_useOctave)
            
ecmCurve = [];

%Get the end point positions
x0 = fiberForceLengthData(1,1);
x0Ti = x0*0.5 - (normLengthZToT12 + normLengthHalfMyosin);

%Solve for lengths and forces of passive pevk-igd & igp sections  ...
[x0PevkIgd, x0Igp, y0Ti] = calcSeriesSpringStretch(x0Ti, ...
                 forceLengthPevkIgdCurve,forceLengthPevkIgdInverseCurve,...
                 forceLengthIgpCurve    ,forceLengthIgpInverseCurve);               

%Solve for the ECM force
y0Ecm = fiberForceLengthData(1,2) - y0Ti;
y0  = max(0,y0Ecm);

dydx0fpe     = calcBezierYFcnXDerivative(       x0, fiberForceLengthCurve,  1);
dydx0PevkIgD = calcBezierYFcnXDerivative(x0PevkIgd, forceLengthPevkIgdCurve,1);        
dydx0IgP     = calcBezierYFcnXDerivative(    x0Igp, forceLengthIgpCurve,    1);
dydx0Ti      = 1/((1/dydx0IgP) + (1/dydx0PevkIgD));

dydxEcm = dydx0fpe-dydx0Ti;
dydx0   = max(0,dydxEcm);
if(y0 < sqrt(eps))
  dydx0 = 0;
end

x1   = calcBezierYFcnXDerivative(1, fiberForceLengthInverseCurve,0);
x1Ti = x1*0.5 - (normLengthZToT12 + normLengthHalfMyosin);

%Solve for lengths and forces of passive pevk-igd & igp sections  ...
[x1PevkIgd, x1Igp, y1Ti] = calcSeriesSpringStretch(x1Ti, ...
                 forceLengthPevkIgdCurve,forceLengthPevkIgdInverseCurve,...
                 forceLengthIgpCurve    ,forceLengthIgpInverseCurve);               


y1 =  1-y1Ti;

dydx1fpe     = calcBezierYFcnXDerivative(x1, fiberForceLengthCurve,1);
dydx1PevkIgD = calcBezierYFcnXDerivative(x1PevkIgd, forceLengthPevkIgdCurve,1);
dydx1IgP     = calcBezierYFcnXDerivative(x1Igp, forceLengthIgpCurve,1);
dydx1Ti      = 1/((1/dydx1IgP) + (1/dydx1PevkIgD));

dydx1 = max(dydx1fpe,dydx1Ti);


%Now fit the curviness of this ECM section to come as close to the
%experimental data as possible

c       = 0.5;
iterMax = 10;

%We are making a half curve
x0=x0*0.5;
x1=x1*0.5;
dydx1 = dydx1*2;
dydx0 = dydx0*2;


p0 =  calcQuinticBezierCornerControlPoints(x0, y0, dydx0, 0, ...
                                           x1, y1, dydx1, 0, c);

%Create the curve structure
ecmCurve.xpts    = p0(:,1);
ecmCurve.ypts    = p0(:,2);
ecmCurve.xEnd         = [x0, x1];
ecmCurve.yEnd         = [y0, y1];
ecmCurve.dydxEnd      = [dydx0, dydx1];
ecmCurve.d2ydx2End    = [0, 0];
ecmCurve.integral = [];  

yErrTotal = 0;
ykErr = 0;
ykErrV = zeros(size(fiberForceLengthData,1),1);
for k=1:1:size(fiberForceLengthData,1)
  xk = fiberForceLengthData(k,1);
  xkTi = xk*0.5 - (normLengthZToT12 + normLengthHalfMyosin);
  %Solve for these lengths
  [xkPevkIgd, xkIgp, ykTi] = calcSeriesSpringStretch(xkTi, ...
                 forceLengthPevkIgdCurve,forceLengthPevkIgdInverseCurve,...
                 forceLengthIgpCurve    ,forceLengthIgpInverseCurve);               
  
  
  ykEcm     = calcBezierYFcnXDerivative(xk*0.5, ecmCurve,0);    
  ykErr = (ykEcm + ykTi) - fiberForceLengthData(k,2);
  ykErrV(k,1) = ykErr;
  yErrTotal = yErrTotal + ykErr*ykErr;  
end

delta = dydx1*0.5;

err = yErrTotal;
ykErrRV = zeros(size(fiberForceLengthData,1),1);
ykErrLV = zeros(size(fiberForceLengthData,1),1);
%Now use the bisection method to find the c of best fit
for i=1:1:iterMax
  
  %Form the left hand bound on c
  dydx1L = dydx1-delta;
  
  p0 =  calcQuinticBezierCornerControlPoints(x0, y0, dydx0, 0, ...
                                             x1, y1, dydx1L, 0, c);
  %Create the curve structure
  ecmCurve.xpts    = p0(:,1);
  ecmCurve.ypts    = p0(:,2);
  ecmCurve.dydxEnd      = [dydx0, dydx1L];

  %Evaluate the error 
  yErrTotal = 0;
  ykErr = 0;
  for k=1:1:size(fiberForceLengthData,1)
    xk = fiberForceLengthData(k,1);
    xkTi = xk*0.5 - (normLengthZToT12 + normLengthHalfMyosin);
    %Solve for these lengths
    [xkPevkIgd, xkIgp, ykTi] = calcSeriesSpringStretch(xkTi, ...
                   forceLengthPevkIgdCurve,forceLengthPevkIgdInverseCurve,...
                   forceLengthIgpCurve    ,forceLengthIgpInverseCurve);               


    ykEcm     = calcBezierYFcnXDerivative(xk*0.5, ecmCurve,0);    
    ykErr = (ykEcm + ykTi) - fiberForceLengthData(k,2);
    ykErrLV(k,1) = ykErr;
    yErrTotal = yErrTotal + ykErr*ykErr;  
  end
  errL = yErrTotal;
    
  %Form and evaluate the right hand bound
  dydx1R = dydx1+delta;
  p0 =  calcQuinticBezierCornerControlPoints(x0, y0, dydx0, 0, ...
                                             x1, y1, dydx1R, 0, c);

  %Create the curve structure
  ecmCurve.xpts    = p0(:,1);
  ecmCurve.ypts    = p0(:,2);
  ecmCurve.dydxEnd      = [dydx0, dydx1R];

  %Evaluate the error 
  yErrTotal = 0;
  ykErr = 0;
  for k=1:1:size(fiberForceLengthData,1)
    xk = fiberForceLengthData(k,1);
    xkTi = xk*0.5 - (normLengthZToT12 + normLengthHalfMyosin);
    %Solve for these lengths
    [xkPevkIgd, xkIgp, ykTi] = calcSeriesSpringStretch(xkTi, ...
                   forceLengthPevkIgdCurve,forceLengthPevkIgdInverseCurve,...
                   forceLengthIgpCurve    ,forceLengthIgpInverseCurve);               

    ykEcm     = calcBezierYFcnXDerivative(xk*0.5, ecmCurve,0);    
    ykErr = (ykEcm + ykTi) - fiberForceLengthData(k,2);
    ykErrRV(k,1) = ykErr;
    yErrTotal = yErrTotal + ykErr*ykErr;  
  end
  errR = yErrTotal;
      
  if(abs(errL) < abs(err))
    err = errL;
    dydx1 = dydx1L;
  end
  if(abs(errR) < abs(err))
    err = errR;
    dydx1 = dydx1R;    
  end
  
  delta=delta*0.5;
  
end


p0 =  calcQuinticBezierCornerControlPoints(x0, y0, dydx0, 0, ...
                                           x1, y1, dydx1, 0, c);

%Create the curve structure
ecmCurve.xpts    = p0(:,1);
ecmCurve.ypts    = p0(:,2);
ecmCurve.dydxEnd      = [dydx0, dydx1];

ecmCurve.name = 'ecm';


