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

function [position,isterminal,direction] = ...
  eventStateThreshold(t,y, idxState, stateValue, stateDirection)

position   = -1;
isterminal = 0;
direction  = 1;

if(idxState <= 0)
  position   = t-stateValue;
  direction  = stateDirection;
  isterminal = 0;
else
  position   = y(idxState)-stateValue;
  isterminal = 0;
  direction  = stateDirection;    
end

if position*stateDirection >= 0
  here=1;
end
