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

function waveState = calcOffsetWaveform(t, lengthOffset, timeSeries,lengthSeries,lengthDotSeries)

waveState = zeros(2,1);

waveState(1) = interp1(timeSeries,lengthDotSeries,t);
waveState(2) = interp1(timeSeries,lengthSeries,t)+lengthOffset;
