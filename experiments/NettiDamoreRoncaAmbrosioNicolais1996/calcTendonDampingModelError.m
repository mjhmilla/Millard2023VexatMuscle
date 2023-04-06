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

function errSq = calcTendonDampingModelError(arg, stiffnessData, dampingData,errSqNorm)

errSq=0;

beta0 = arg(1,1);
beta1 = arg(2,1);


for z=1:1:length(dampingData)
  beta = beta0 + beta1.*mean(stiffnessData(z).y);
  err = beta - mean(dampingData(z).y);
  errSq = errSq + err*err;
end

errSq = errSq./errSqNorm;