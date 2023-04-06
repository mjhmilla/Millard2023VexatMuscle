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

function clampedActivation = clampActivation(activation, activationMin, activationMax)
%%
% This function clamps the activation to be within the domain of 
% activationMin and activationMax
%
% @param activation: the value of the activation of the muscle
% @param activationMin: the minimum allowed activation value
% @param activationMax: the maximum allowed activation value.
%
% @return clampedActivation
%%
clampedActivation = activation;

if activation < activationMin
    clampedActivation = activationMin;
end

if activation > activationMax
   clampedActivation = activationMax; 
end
