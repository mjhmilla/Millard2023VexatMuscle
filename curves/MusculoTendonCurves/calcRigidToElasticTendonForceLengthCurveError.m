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

function errV = calcRigidToElasticTendonForceLengthCurveError( params, data, scaling, ...
    fixedParams, elasticTendonReferenceModel, flag_useOctave)

xshift = params(1)./scaling;
xwidth = params(2)./scaling;
kLow   = fixedParams(1,1);
kNum   = fixedParams(1,2);

normLengthZero = xshift;
normLengthToe  = xwidth + xshift;
fToe  = 1;
kZero = fixedParams(1,3);

kToe  = kNum/(normLengthToe-normLengthZero);
curviness= fixedParams(1,4);

computeIntegral = 0;



fiberForceLengthCurve = createFiberForceLengthCurve2021(normLengthZero,...
                    normLengthToe,...
                    fToe,...
                    kZero,...
                    kLow,...
                    kToe,...
                    curviness,...
                    0,...
                    'fitted',...
                    flag_useOctave);
                                                
errV = zeros(size(data,1),1);

lopt    = elasticTendonReferenceModel.musculotendon.optimalFiberLength;
penOpt  = elasticTendonReferenceModel.musculotendon.pennationAngle;
ltSlk   = elasticTendonReferenceModel.musculotendon.tendonSlackLength;

for(i=1:1:size(data,1))
   lceN_ET = data(i,1);
   fpeN_ET = calcBezierYFcnXDerivative(data(i,1), ...
       elasticTendonReferenceModel.curves.fiberForceLengthCurve,0);

   fiberKin = calcFixedWidthPennatedFiberKinematicsAlongTendon(lceN_ET,0,lopt,penOpt);
   alpha = fiberKin.pennationAngle;

   ftN_ET = fpeN_ET*cos(alpha);

   ltN_ET = calcBezierYFcnXDerivative(ftN_ET, ...
            elasticTendonReferenceModel.curves.tendonForceLengthInverseCurve, 0);
   lp = (lopt * lceN_ET)*cos(alpha) + ltN_ET*ltSlk;

   lceAT_RT = lp - ltSlk;
   
   fiberKinRT = calcFixedWidthPennatedFiberKinematics(lceAT_RT,0,lopt,penOpt);
   lce_RT = fiberKinRT.fiberLength;

   lceN_RT = lce_RT/lopt;

   fpeN_RT = calcBezierYFcnXDerivative(lceN_RT, fiberForceLengthCurve,0);

   errV(i) = (fpeN_RT - fpeN_ET)*scaling;
   
end
