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
function sse = calcFrequencyDomainSquaredError(x, dataFreqRadians, dataGain,...
   dataPhaseRadians, argScaling, objScaling, gainScaling, phaseScaling)

sse = 0;
assert(length(x)==2)

k    = x(1)*argScaling(1);
beta = x(2)*argScaling(2);

for i=1:1:length(dataFreqRadians)
  
  %modelResponse = k + (beta*complex(0,-1)*dataFreqRadians(i,1));      
  modelResponse = calcFrequencyModelResponse(k,beta,dataFreqRadians(i,:));
  gainError = (abs(modelResponse) - dataGain(i,1))/(dataGain(end,1));
  phaseError= (angle(modelResponse) - dataPhaseRadians(i,1))/(dataPhaseRadians(end,1));
  
  errorSq = (gainError*gainError)*gainScaling ...
           +(phaseError*phaseError)*phaseScaling;
  
  sse = sse + errorSq/(length(dataFreqRadians));
end

sse = sse.*objScaling;
