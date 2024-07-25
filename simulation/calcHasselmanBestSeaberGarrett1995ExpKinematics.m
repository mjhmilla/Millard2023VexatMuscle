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