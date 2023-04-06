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

function uC = clampU(u)
%%
% This is used to ensure that the argument for the Bezier curve's Bernstein
% polynomial is limited to [0,1] and not a bit outside
%
% @param u
% @returns uC, version of u that is guaranteed to be within [0,1]
%%
uC = u;
if(u<0.0)
    uC=0;
end
if(u>1.0)
    uC=1;
end

