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
function path = oscillatingSmoothStepFunction(t, length, deltaLength, ...
                                            omega, timeStepStart, timeStepEnd,...
                                            stepHeight)


y  = length+deltaLength*sin(omega*t);
dy = (deltaLength*omega)*cos(omega*t);

if(t > timeStepStart && t < timeStepEnd)
    dt = 1/(timeStepEnd-timeStepStart);
    delta = (t-timeStepStart)*dt;
    delta2=delta*delta;
    delta3=delta2*delta;
    s = stepHeight*(3*delta2-2*delta3);

    ds = stepHeight*(6*delta-6*delta2)*dt;

    y=y+s;
    dy=dy+ds;

end

if(t >= timeStepEnd)
    y=y+stepHeight;
end

path = [y;dy];