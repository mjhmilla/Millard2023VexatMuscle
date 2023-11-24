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


k0=100;
d0=1;
if(flag_usingOctave == 0)
f0 = fit(freq,...
         gain,...
         'poly1');

k0 = f0.p2;
d0 = f0.p1;
end



%If fmincon is available, then constrain damping to be positive
A0 = [];%[-1 0; 0 -1];
b0 = [];%[ 0 ; 0];

lb = [0;0];
ub = [];

x0 = [1,1];
argScaling =[k0,d0];
argScalingMin=1;
argScaling(argScaling<argScalingMin) = argScalingMin;
objScaling = 1;

gainScaling=1;
phaseScaling=1;

err0 = calcFrequencyDomainSquaredError(x0, ...
        freq,...
        gain,...
        phase,...
        argScaling,objScaling,gainScaling,phaseScaling);

objScaling  = 1/max(abs(err0),sqrt(eps));

errFcn0 = @(argX)calcFrequencyDomainSquaredError(argX, ...
        freq,...
        gain,...
        phase,...
        argScaling,objScaling,gainScaling,phaseScaling);


paramOpt  = []; 
fval      = [];
exitFlag  = 0;
options   = [];
if(flag_usingOctave == 0)               
options=optimoptions('fmincon','Display','none','Algorithm','sqp');%optimset('Display','none');
[paramOpt, fval, exitFlag] = fmincon(errFcn0, x0,[],[],[],[],lb,ub,[],options);        
else
options=optimset('Display','none','Algorithm','sqp');%optimset('Display','none');
[paramOpt, fval, exitFlag] = fminsearch(errFcn0, x0,options);                        
end
k0 = paramOpt(1)*argScaling(1);
d0 = paramOpt(2)*argScaling(2);

%Now that we've got a decent initial solution, use it to
%adjust the scaling to the inputs and run the optimization 
%routine again.

argScaling = [k0,d0];
argScaling(argScaling<argScalingMin) = argScalingMin;

x1 = [0,0];
x1(1,1) = k0/argScaling(1);
x1(1,2) = d0/argScaling(2);
A1 = [0 -1];
b1 = [ 0 ];


err1 = calcFrequencyDomainSquaredError(x1, ...
        freq,...
        gain,...
        phase,...
        argScaling,1,gainScaling,phaseScaling);

objScaling = 1/max(sqrt(eps),abs(err1));

errFcn1 = @(argX)calcFrequencyDomainSquaredError(argX, ...
        freq,...
        gain,...
        phase,...
        argScaling,objScaling,gainScaling,phaseScaling);     

%options=optimset('Display','none');
%[paramOpt, fval, exitFlag] = fminsearch(errFcn, [k1, d1],options); 

if(flag_usingOctave==0)       
[paramOpt, fval, exitFlag] = fmincon(errFcn1, x1,[],[],[],[],lb,ub,[],options);        
else
[paramOpt, fval, exitFlag] = fminsearch(errFcn1, x1,options);
end


stiffness = paramOpt(1)*argScaling(1);
damping   = paramOpt(2)*argScaling(2);