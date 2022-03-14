function invCurve = createInverseCurve(curve)

invCurve.xpts = curve.ypts;
invCurve.ypts = curve.xpts;

invCurve.xEnd = curve.yEnd;
invCurve.yEnd = curve.xEnd;

%assert(abs(curve.dydxEnd(1)) > eps && abs(curve.dydxEnd(2)), ...
%       'End-point slopes are flat: this curve is not invertible');

invCurve.dydxEnd = zeros(1,2);
invCurve.d2ydx2End=zeros(1,2);
for k=1:1:2
    if(isinf(curve.dydxEnd(k))==1)
       invCurve.dydxEnd(k) = 0; 
    elseif( curve.dydxEnd(k)==0)
       invCurve.dydxEnd(k) = Inf;     
    elseif( isnan(curve.dydxEnd(k)) == 1)
       invCurve.dydxEnd(k)=NaN;
    else
       delta = eps - (k-1)*2*eps;
       invCurve.dydxEnd(k) = calcBezierYFcnXDerivative(invCurve.xEnd(k)+delta,...
                                                       invCurve,...
                                                       1);
       invCurve.d2ydx2End(k) = calcBezierYFcnXDerivative(invCurve.xEnd(k)+delta,...
                                                       invCurve,...
                                                       2);
                                                   
    end
    
end

invCurve.integral = [];

invCurve.name = sprintf('%s.%s',curve.name,'inverse');

