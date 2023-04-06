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

function [uC, isClampedToLB,isClampedToUB]= clampWithinInterval(u,umin,umax)
%%
% This function clamps u to be within the interval [0,1]
%
% @param u
% @return uC, version of u that is guaranteed to be within [0,1]
%%
uC = u;
isClampedToLB=0;
isClampedToUB=0;
if(u<=umin)
    uC=umin;
    isClampedToLB=1;
end
if(u>=umax)
    uC=umax;
    isClampedToUB=1;
end

