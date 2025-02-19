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
kNum   = params(3)./scaling;

kLow        = fixedParams(1,1);
kZero       = fixedParams(1,2);
%kNum        = fixedParams(1,3);
curviness   = fixedParams(1,3);

normLengthZero = xshift;
normLengthToe  = xwidth + xshift;
fToe  = 1;


kToe  = kNum/(normLengthToe-normLengthZero);

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

dataMinForce= 0.0;
count = 0;
for i=1:1:length(data)
    if(data(i,2) >= dataMinForce)
        count=count+1;
    end
end

errV = zeros(count,1);



j=1;
for i=1:1:size(data,1)
    if(data(i,2)>=dataMinForce)        
       errV(j) = (calcBezierYFcnXDerivative(data(i,1),fiberForceLengthCurve,0)...
                 - data(i,2))*scaling;
       j=j+1;
    end   
end
