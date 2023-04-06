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

function sinusoidState = calcSinusoidState(t, tStart, tEnd, ...
    lengthStart, sinusoidFrequencyInRadiansPerSecond, sinusoidAmplitude)
%
% This function implements:
%
% if t < tStart
%   l(t) = lengthStart
% else 
%   l(t) = lengthStart + sinusoidAmplitude*sin((t-tStart)*sinusoidFrequencyInRadiansPerSecond)
% @param t: time
% @param tStart: when the ramp starts
% @param tEnd: when the ramp begins
% @param lengthStart: the starting rest length
% @param sinusoidFrequencyInRadiansPerSecond
% @param sinusoidAmplitude: amplitude of the sinusoid
% @returns sinusoidState 2x1 vector:
%          sinusoidState(1) = dl(t)/dt 
%          sinusoidState(2) =  l(t)
%
%%

sinusoidState = zeros(2,1);

dldt = NaN;
l    = NaN;

if(t < tStart)
  l    = lengthStart;
  dldt = 0;
elseif(t >= tStart && t <= tEnd)
  l    = lengthStart + sinusoidAmplitude*sin(sinusoidFrequencyInRadiansPerSecond*(t-tStart));
  dldt = (sinusoidAmplitude*sinusoidFrequencyInRadiansPerSecond) ...
            *cos(sinusoidFrequencyInRadiansPerSecond*(t-tStart));
else
  l    = lengthStart + sinusoidAmplitude*sin(sinusoidFrequencyInRadiansPerSecond*(tEnd-tStart));
  dldt = 0;
end

sinusoidState(1) = dldt;
sinusoidState(2) = l;
