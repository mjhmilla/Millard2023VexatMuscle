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

function fiberState = clampFiberState(lce,dlce,lceMin)
%%
% This function returns the state of the fiber, and if necessary, 
% clamps it to a lower bound defined by lceMin
%
% @param lce: length of the fiber (m)
% @param dlce: fiber lengthening velocity (m/s)
% @param lceMin: the minimum allowable length of the fiber
%
% @return fiberState a structure with fields of
%           .lce : fiber length
%           .dlce: fiber velocity
%           .isClamped: 0 - not clamped
%                       1 - clamped
%%
isClamped = 0;
if(lce < lceMin || (lce == lceMin && dlce < 0))            
    %Clamp the fiber length along the tendon
    lce  = lceMin;
    dlce = 0;    
    isClamped = 1;
end

fiberState.lce     = lce;
fiberState.dlce    = dlce;
fiberState.isClamped = isClamped; 