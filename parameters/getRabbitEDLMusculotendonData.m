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

function [activeForceLengthData,...
          passiveForceLengthData] = getRabbitEDLMusculotendonData(...
                                      musculotendonPropertiesExp,...                              
                                      musculotendonProperties,...
                                      normPlateauShift,...
                                      useElasticTendon,...
                                      shiftPassiveForceLengthCurveEDLInM,...
                                      projectFolders,...
                                      flag_useOctave)



%Evaluate the passive force length curve up to a normalized force of 0.2
%the maximum passive force measured
k1 = 0.001;
k2 = 0.82;
lpec0 = 13.1;

fpecMax = 0.2*musculotendonPropertiesExp.fiso;
dlMax = log((fpecMax/k1)+1)/k2;

dlce = ([0:0.1:1]').*dlMax;
fpe  = k1*(exp(k2.*dlce)-1);
lce  = lpec0 + dlce + shiftPassiveForceLengthCurveEDLInM*1000;

lce = lce ./1000; %Convert to meters
lceN = lce./musculotendonPropertiesExp.optimalFiberLength;
fpeN = fpe./musculotendonPropertiesExp.fiso;

passiveForceLengthData = [lceN, fpeN];

%Evaluate the active force length curve, which is just linear
lceFLPts = [0.39; 0.82; 1; 1.21;2.03].*musculotendonPropertiesExp.optimalFiberLengthSiebert;
fceFLPts = [   0; 0.91; 1; 1   ;   0].*musculotendonPropertiesExp.fiso;

lceFLmid = lceFLPts(1:end-1,1)+diff(lceFLPts).*0.5;

lceFL = [lceFLPts; lceFLmid ];
lceFL = sort(lceFL);
fceFL = interp1(lceFLPts, fceFLPts, lceFL);

lceFLN = lceFL./musculotendonPropertiesExp.optimalFiberLength;
fceFLN = fceFL./musculotendonPropertiesExp.fiso;

activeForceLengthData = [lceFLN, fceFLN];

flag_debug=1;
if(flag_debug==1)
    figDebug=figure;
    plot(activeForceLengthData(:,1),activeForceLengthData(:,2),'-k');
    hold on;
    plot(passiveForceLengthData(:,1),passiveForceLengthData(:,2),'-b');
    hold on;
    xlabel('$$\ell^{M} / \ell^M_o$$');
    ylabel('$$f^{M} / f^M_o$$');
    here=1;
end




