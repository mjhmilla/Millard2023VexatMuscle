function curve = createCompressiveForceLengthCurve2022(...
                    fiberLengthAlongTendonAtOneNormForce,...
                    fiberLengthAlongTendonAtZeroNormForce,...
                    curviness, muscleName, flag_usingOctave)

xLeft = fiberLengthAlongTendonAtOneNormForce;
yLeft = 1;

xRight = fiberLengthAlongTendonAtZeroNormForce;
yRight = 0;

dydxLeft  = 5.0*(yRight-yLeft)/(xRight-xLeft);
dydxRight = 0;

p01 = calcQuinticBezierCornerControlPoints(...
         xLeft,  yLeft,  dydxLeft, 0,...
        xRight, yRight, dydxRight, 0, curviness);


curve.xpts = p01(:,1);
curve.ypts = p01(:,2);

curve.xEnd = [xLeft,xRight];
curve.yEnd = [yLeft,yRight];
curve.dydxEnd = [dydxLeft,dydxRight];
curve.d2ydx2End = [0,0];
curve.integral = [];

curve.name = sprintf('%s.%s',muscleName,'compressiveForceLengthCurve');