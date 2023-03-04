function curve = createProximalIgCurve(x0, y0 ,dydx0, ... 
                                       x1, y1, dydx1,...       
                                       x2, y2, dydx2, c, ...
                                       computeIntegral, curveName, ...
                                       flag_usingOctave)
% x1,y1,dydx1,...                                      

c = scaleCurviness(c);

p0 = calcQuinticBezierCornerControlPoints(x0, y0, dydx0, 0, ...
                                          x1, y1, dydx1, 0,c);                                        
p1 = calcQuinticBezierCornerControlPoints(x1, y1, dydx1, 0, ...
                                          x2, y2, dydx2, 0,c);                                      

%p0 = calcQuinticBezierCornerControlPoints(x0, y0, dydx0, 0, ...
%                                          x2, y2, dydx2, 0,c);


%Create the curve structure
curve.xpts    = [p0(:,1),p1(:,1)];
curve.ypts    = [p0(:,2),p1(:,2)];

curve.xEnd         = [x0,x2];
curve.yEnd         = [y0,y2];
curve.dydxEnd      = [dydx0,dydx2];
curve.d2ydx2End    = [0, 0];

curve.integral = [];

curve.name = curveName;

if(computeIntegral == 1)
    xScaling = (x2-x0);
    curve.integral = ...
        createCurveIntegralStructure(curve, ...
                                     1000,...
                                     1e-12,...
                                     xScaling, flag_usingOctave,1);    
end                       
                                        
                                        