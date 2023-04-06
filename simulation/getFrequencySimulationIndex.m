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

function idxSim = getFrequencySimulationIndex(amplitudeMM,bandwidthHz, ...
                  nominalForceN,normFiberLength, freqSimData)

idxSim = 0;
tol = 1e-6;
for m=1:1:size(freqSimData.force,2)     
  if( abs(freqSimData.amplitudeMM(1,m)     - amplitudeMM   ) <= tol && ...
      abs(freqSimData.bandwidthHz(1,m)     - bandwidthHz   ) <= tol && ...
      abs(freqSimData.nominalForceDesired(1,m)    - nominalForceN ) <= tol && ...
      abs(freqSimData.normFiberLength(1,m) - normFiberLength   ) <= tol)
    if(idxSim == 0)
      idxSim = m;
    else
      assert(0); %Error condition: there should not be 2 simulations with 
                 %the same configuration
    end
  end
end                