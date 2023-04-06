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

function  idx =  getSignalIndex(amplitudeMM, bandwidthHz, inputFunctions)

idx = 0;

for z=1:1:length(inputFunctions.amplitudeMM)
  if( abs(inputFunctions.amplitudeMM(z)-amplitudeMM<sqrt(eps)) ...
    && abs(inputFunctions.bandwidthHz(z)-bandwidthHz)<sqrt(eps) )
    idx = z;
  end          
end  