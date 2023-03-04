function errVal = calcThirdDerivativeError( argPts, indexLeftSpline, curveParamsDefault, errScaling)

indexRightSpline=indexLeftSpline+1;

curveParams = curveParamsDefault;
curveParams.ypts(4,indexLeftSpline)  = argPts(1,1);
curveParams.ypts(3,indexRightSpline) = argPts(2,1);

nrow = size(curveParams.xpts,1);
ncol = size(curveParams.xpts,2);
xmin = curveParams.xEnd(1,1);
xmax = curveParams.xEnd(1,2);



%
% Evaluate the 3rd derivative at the end point of the left spline
%
u = 1.0;
col = indexLeftSpline;

  xV  = curveParams.xpts(:,col);
  x1V = diff(xV) .*(nrow-1);
  x2V = diff(x1V).*(nrow-2);
  x3V = diff(x2V).*(nrow-3);

  yV  =  curveParams.ypts(:,col);
  y1V =  diff(yV) .*(nrow-1);
  y2V =  diff(y1V).*(nrow-2);
  y3V =  diff(y2V).*(nrow-3);

  x1  = calc1DBezierCurveValue(u, x1V);
  y1  = calc1DBezierCurveValue(u, y1V);

  x2  = calc1DBezierCurveValue(u, x2V);
  y2  = calc1DBezierCurveValue(u, y2V);

  x3  = calc1DBezierCurveValue(u, x3V);
  y3  = calc1DBezierCurveValue(u, y3V);

  t1 = 1 / x1;
  t3 = x1*x1;
  t4 = 1 / t3;
  t11 = x2*x2;
  t14 = y1 * t4;

  %val = (  (x1*x1*y3-2*x1*x2*y2+(2*x2*x2-x1*x2)*y1)...
  %        /(x1*x1*x1);
  %val =(    (x2*y2+x1*y3-x3*y1-x2*y2)/(x1*x1) ...
  %              - 2*(x1*y2-x2*y1)*x2/(x1*x1*x1)  )/x1;
  d3ydx3L = ((y3*t1 - 2*y2*t4*x2 ...
        + 2*y1/t3/x1 * t11 - t14 * x3) * t1 ...
      - (y2*t1 - t14*x2)*t4*x2) * t1;      


%
% Evaluate the 3rd derivative at the end point of the left spline
%
u = 0.0;
col = indexRightSpline;

  xV  = curveParams.xpts(:,col);
  x1V = diff(xV) .*(nrow-1);
  x2V = diff(x1V).*(nrow-2);
  x3V = diff(x2V).*(nrow-3);

  yV  =  curveParams.ypts(:,col);
  y1V =  diff(yV) .*(nrow-1);
  y2V =  diff(y1V).*(nrow-2);
  y3V =  diff(y2V).*(nrow-3);

  x1  = calc1DBezierCurveValue(u, x1V);
  y1  = calc1DBezierCurveValue(u, y1V);

  x2  = calc1DBezierCurveValue(u, x2V);
  y2  = calc1DBezierCurveValue(u, y2V);

  x3  = calc1DBezierCurveValue(u, x3V);
  y3  = calc1DBezierCurveValue(u, y3V);

  t1 = 1 / x1;
  t3 = x1*x1;
  t4 = 1 / t3;
  t11 = x2*x2;
  t14 = y1 * t4;

  %val = (  (x1*x1*y3-2*x1*x2*y2+(2*x2*x2-x1*x2)*y1)...
  %        /(x1*x1*x1);
  %val =(    (x2*y2+x1*y3-x3*y1-x2*y2)/(x1*x1) ...
  %              - 2*(x1*y2-x2*y1)*x2/(x1*x1*x1)  )/x1;
  d3ydx3R = ((y3*t1 - 2*y2*t4*x2 ...
        + 2*y1/t3/x1 * t11 - t14 * x3) * t1 ...
      - (y2*t1 - t14*x2)*t4*x2) * t1;      
    
flag_debug=0;
if(flag_debug==1)
  d3ydx3LTest = calcBezierYFcnXDerivative( ...
    curveParams.xpts(end,indexLeftSpline)-1e-6,...
    curveParams,...
    3);
  d3ydx3RTest = calcBezierYFcnXDerivative( ...
    curveParams.xpts(1,indexRightSpline)+1e-6,...
    curveParams,...
    3);  
  
  relErrL = abs(d3ydx3LTest-d3ydx3L)/abs(d3ydx3LTest+d3ydx3L);
  relErrR = abs(d3ydx3RTest-d3ydx3R)/abs(d3ydx3RTest+d3ydx3R); 
  assert(relErrL < 1e-3);
  assert(relErrR < 1e-3);  
end
    
    
errVal = (d3ydx3L-d3ydx3R)*(d3ydx3L-d3ydx3R)*errScaling;

% 
%
%
%
% dydx   = dydu/dxdu
% d2ydx2 =      d2ydu2/dxdu
%         -dydu*d2xdu2/(dxdu*dxdu)
% d3ydx3 =        d3ydu3/dxdu
%       -2*d2ydu2*d2xdu2/(dxdu*dxdu)
%         -dydu*d3xdu3/(dxdu*dxdu)
%       +2*dydu*d2xdu2/(dxdu*dxdu*dxdu)
%
%