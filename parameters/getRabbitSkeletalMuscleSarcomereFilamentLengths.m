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
function [activeForceLengthKeyPoints, ...
          halfMyosinBareLength, ...
          halfMyosinLength,...
          zLineLength,...
          actinLength] = getRabbitSkeletalMuscleSarcomereFilamentLengths()

%From Higuchi et al.
actinLength       = 1.12;
myosinLength      = 1.63;
zLineLength       = 0.07;
myosinBareLength  = 0.16;

halfMyosinBareLength = myosinBareLength*0.5;
halfMyosinLength     = myosinLength*0.5;

lasc  = max(actinLength,myosinLength)    + 2*zLineLength;    
loptA = 2*actinLength  -myosinBareLength + 2*zLineLength;
loptB = 2*actinLength  +myosinBareLength + 2*zLineLength;
lmax  = 2*actinLength  +myosinLength     + 2*zLineLength;

geoRabbit = [1.27,lasc, loptA, loptB,  lmax]; 

activeForceLengthKeyPoints = geoRabbit;

flag_debug=1;
if(flag_debug==1)
    halfMyosinBareLength =     ( activeForceLengthKeyPoints(1,4) ...
                               - activeForceLengthKeyPoints(1,3) )*0.25;
    
    assert(abs(halfMyosinBareLengthTest-halfMyosinBareLength) < sqrt(eps));

    halfMyosinLength      = 0.5*(activeForceLengthKeyPoints(1,5) ...
                                -activeForceLengthKeyPoints(1,4)) ...
                                +halfMyosinBareLength;  

    assert(abs(halfMyosinLengthTest-halfMyosinLength) < sqrt(eps));
                       
    actinLength           = 0.5*( (activeForceLengthKeyPoints(1,5))...
                            -2*halfMyosinLength ...
                            -2*zLineLength); 
    assert(abs(actinLengthTest-actinLength) < sqrt(eps));
end
here=1;