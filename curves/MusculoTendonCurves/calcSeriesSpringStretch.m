function [x1, x2, y12] = calcSeriesSpringStretch(x12, y1Curve,y1CurveInv, y2Curve,y2CurveInv)

%calcBezierYFcnXDerivative(
y12 = 0.5*(calcBezierYFcnXDerivative(x12, y1Curve, 0) ...
         + calcBezierYFcnXDerivative(x12, y2Curve,0));
delta = y12;

x1 = calcBezierYFcnXDerivative(y12, y1CurveInv,0);
x2 = calcBezierYFcnXDerivative(y12, y2CurveInv,0);
err = (x1+x2)-x12;

%Use the Bisection method to get close to the root
for i=1:1:5
  y12L = y12-delta;
  x1L = calcBezierYFcnXDerivative(y12L, y1CurveInv,0);
  x2L = calcBezierYFcnXDerivative(y12L, y2CurveInv,0);
  errL = (x1L+x2L)-x12;
  
  y12R = y12+delta;
  x1R = calcBezierYFcnXDerivative(y12R,y1CurveInv,0);
  x2R = calcBezierYFcnXDerivative(y12R,y2CurveInv,0);
  errR = (x1R+x2R)-x12;
  
  if( abs(errL) < abs(err) )
    y12 = y12L;
    err = errL;
  end

  if( abs(errR) < abs(err) )
    y12 = y12R;
    err = errR;
  end

  delta = 0.5*delta;
end

iter=0;
iterMax=100;
tol=1e-12;


%Polish the root up using Newton's method
while(abs(err) > tol && iter < iterMax)
  
  x1 = calcBezierYFcnXDerivative(y12,y1CurveInv,0);
  x2 = calcBezierYFcnXDerivative(y12,y2CurveInv,0);
  err = (x1+x2)-x12;

  dx1 = calcBezierYFcnXDerivative(y12,y1CurveInv,1);
  dx2 = calcBezierYFcnXDerivative(y12,y2CurveInv,1);
  derr = (dx1+dx2);
  
  delta = -err/derr;
  y12=y12+delta;
  
  iter=iter+1;
end

assert(abs(err) <= tol);