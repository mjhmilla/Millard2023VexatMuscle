%%
% SPDX-FileCopyrightText: 2023 Matthew Millard <mfreqillard.matthew@gmail.com>
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
function [stiffness,damping,exitFlag] = ...
    fitSpringDamperToGainAndPhaseProfiles(freq,gain,phase, ...
                                            flag_usingOctave)


k0=1;
d0=1;
if(flag_usingOctave == 0)
    f0 = fit(freq,...
             gain,...
             'poly1');
    
    k0 = f0.p2;
    d0 = f0.p1;
end



%If fmincon is available, then constrain damping to be positive
A0 = [0 -1];
b0 = [ 0];

lb = [];
ub = [];

x0 = [1;1];
argScaling =[k0;d0];
argScalingMin=1;
argScaling(argScaling<argScalingMin) = argScalingMin;
objScaling = 1;


err0 = calcFrequencyDomainSquaredError(x0, ...
        freq,...
        gain,...
        phase,...
        argScaling,objScaling);

objScaling  = 1/max(sum(err0),1);

errFcn = @(argX)calcFrequencyDomainSquaredError(argX, ...
        freq,...
        gain,...
        phase,...
        argScaling,objScaling);

errStarting = errFcn(x0);

paramOpt  = []; 
fval      = [];
exitFlag  = 0;
options   = [];

if(flag_usingOctave == 0)               
    options.Algorithm = 'levenberg-marquardt';
    options.Display='off';
    [paramOpt, fval,residual,exitflag, output]=...
        lsqnonlin(errFcn,x0,[],[],options);

else
    options=optimset('Display','off','Algorithm','sqp');
    [paramOpt, fval, exitFlag] = fminsearch(errFcn0, x0,options);                        
end

stiffness = paramOpt(1)*argScaling(1);
damping = paramOpt(2)*argScaling(2);
