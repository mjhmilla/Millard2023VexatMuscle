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

function idx = getIndexIntoVectors(entryA, entryB, vectorA, vectorB)
idx = 0;
tol = 1e-6;
i = 1;

while(idx == 0 && i <= length(vectorA))
  if( abs( vectorA(i)-entryA ) <= tol && abs( vectorB(i)-entryB ) <= tol)
    idx = i;
  end
  i=i+1;
end
