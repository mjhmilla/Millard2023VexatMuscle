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
function dy = stateDerivativeDwellDependency(t, y, pathFcn,...
    omega, zeta,timeConstant)

if(t > 0.1)
    here=1;
end

path = pathFcn(t);
l  = path(1);
dl = path(2);

x  = y(1);
dx = y(2);

f = (dl-dx)/timeConstant;

ddx = f - (2*zeta*omega)*dx -(omega*omega)*x;

dy = [dx;ddx];

%t0 = dl/velocityStick;
%s  = exp(-t0*t0);
%dy = dl*(1-s) - s*(lsrs/timeConstant);



