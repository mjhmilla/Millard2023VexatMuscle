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

function rampState = calcRampStateSharp(t, tStart, tEnd, lengthStart, rampSlope)
%tStart,tEnd,lengthStart,rampSlope)
%%
% Calculates the value of a ramp function that begins with an intial length
% and extends following a shape like this
%
%                   tEnd
%                   /------------ lengthStart + rampSlope*(tEnd-tStart)
%                  / rampSlope
%   --------------/
%   lengthStart   tstart
%
%
% @param t: time
% @param tStart: when the ramp starts
% @param tEnd: when the ramp begins
% @param lengthStart: the starting rest length
% @param rampSlope: the slope of the ramp
%
% @returns rampState 2x1 vector:
%          rampState(1) = dl(t)/dt 
%          rampState(2) =  l(t)
%
%%

rampState = zeros(2,1);

dldt = NaN;
l    = NaN;

if(t < tStart)
  l    = lengthStart;
  dldt = 0;
elseif(t >= tStart && t <= tEnd)
  l    = lengthStart + rampSlope*(t-tStart);
  dldt = rampSlope;
else
  l    = lengthStart + rampSlope*(tEnd-tStart);
  dldt = 0;
end





rampState(1) = dldt;
rampState(2) = l;
