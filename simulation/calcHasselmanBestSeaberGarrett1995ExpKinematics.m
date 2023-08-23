function expKinematics =...
    calcHasselmanBestSeaberGarrett1995ExpKinematics(...
        musculotendonProperties, ...
        normMuscleCurves)

lceOpt      = musculotendonProperties.optimalFiberLength;
alphaOpt    = musculotendonProperties.pennationAngle;
ltSlk       = musculotendonProperties.tendonSlackLength;
etOne       = musculotendonProperties.tendonStrainAtOneNormForce;
fiso        = musculotendonProperties.fiso;

fNStart   = 0.050*9.81/fiso; %Starting point was the passive length at 50g load   

lceNStart = calcBezierFcnXGivenY(fNStart,normMuscleCurves.fiberForceLengthCurve,0.8);
ltNStart  = calcBezierFcnXGivenY(fNStart,normMuscleCurves.tendonForceLengthCurve,1.01);
fiberKin = calcFixedWidthPennatedFiberKinematicsAlongTendon(lceNStart,0,lceOpt,alphaOpt);
alpha = fiberKin.pennationAngle;
lpStart = lceNStart*cos(alpha)*lceOpt+ltNStart*ltSlk;



expKinematics.lceNStart= lceNStart;
expKinematics.lpStart  = lpStart;

expKinematics.lpOpt    = lceOpt*cos(alphaOpt) + ltSlk*(1+etOne); 