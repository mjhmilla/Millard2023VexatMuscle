function ramp = createC2ContinuousRamp(x0,y0,dydx0,x1,y1,dydx1,...
                                       computeIntegral,curveName,flag_usingOctave)


c = scaleCurviness(0.5);

%y1 = y0 + 0.5*(dydx1+dydx0)*(x1-x0);

p0 = calcQuinticBezierCornerControlPoints(x0, y0, dydx0, 0, ...
                                          x1, y1, dydx1, 0,c);

%Create the curve structure
ramp.xpts    = p0(:,1);
ramp.ypts    = p0(:,2);

ramp.xEnd         = [x0,x1];
ramp.yEnd         = [y0,y1];
ramp.dydxEnd      = [dydx0,dydx1];
ramp.d2ydx2End    = [0, 0];

ramp.integral = [];

ramp.name = curveName;

if(computeIntegral == 1)
    xScaling = (x1-x0);
    ramp.integral = ...
        createCurveIntegralStructure(ramp, ...
                                     1000,...
                                     1e-12,...
                                     xScaling, flag_usingOctave,1);    
end