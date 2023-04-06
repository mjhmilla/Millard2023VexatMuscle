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

function errV = calcFittedFiberForceLengthCurveError( params, data, scaling, ...
    fixedParams,flag_useOctave)

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

for(i=1:1:size(data,1))
   errV(i) = (calcBezierYFcnXDerivative(data(i,1),fiberForceLengthCurve,0)...
             - data(i,2))*scaling;
   
end
