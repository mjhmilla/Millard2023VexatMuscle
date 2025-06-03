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

function rampState = calcLinearlyInterpolatedState(t, ...
                        timeSeries,lengthSeries, flag_returnDerivative)

rampState=[];
if(flag_returnDerivative==1)
    rampState = zeros(2,1);
else
    rampState = 0;
end


l = interp1(timeSeries,lengthSeries,t,...
            'linear','extrap');

velocitySeries = diff(lengthSeries)./diff(timeSeries);

dldt = 0;
for i=2:1:length(timeSeries)
    if( t >= timeSeries(i-1,1) && t < timeSeries(i,1) )
        dldt = velocitySeries(i-1,1);
    end
end



if(flag_returnDerivative==1)
    rampState(1) = dldt;
    rampState(2) = l;
else
    rampState = l;
end