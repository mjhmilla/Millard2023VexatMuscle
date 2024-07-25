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
function stiffness = calcActinMyosinLoadPathStiffness(...
                      numAttachedCrossbridges, ...
                      lengthMLineToMyosinMeanAttachmentPoint, ...
                      sarcomereLength,...
                      sarcomereProperties)
% k: stiffness
% l: length
% x: crossbridge
% m: myosin
% a: actin
kx = sarcomereProperties.crossbridgeStiffness;
ka = sarcomereProperties.actinStiffness;
km = sarcomereProperties.myosinStiffness;
lm = sarcomereProperties.myosinLength;
la = sarcomereProperties.actinLength;

% N: number
xN = numAttachedCrossbridges;

%Mean filament attachment point
%A : attachment
lmA= lengthMLineToMyosinMeanAttachmentPoint; 
laA = (0.5*sarcomereLength - lmA);
assert(laA <= la, 'Error: actin is too short for the desired attachment point');

%n: normalized
nlmA = lmA/lm;
nlaA = laA/la;

%Stiffness of the loaded filament length:
% Z-line to the actin's point of attachement
% M-line to myosin's point of attachment
%
% Note: stiffness varies inversely with length
kaA = ka./nlmA;
kmA = km./nlaA;

%Total stiffness of all components.
% To avoid the combinatorial problem of assigning crossbridges to
% specific actin filaments, I'm going to assume that crossbridges can
% be fractionally divided among all 6 actin filaments. This assumption
% will overestimate the stiffness in the case that all crossbridges are 
% attached to 1 actin filament, or are more heterogeneously distributed
% in general. Because actin filaments are much stiffer than crossbridges
% the error caused by this assumption is low. For example, if 3
% crossbridges attach to 3 actin filaments the stiffness is 0.6582-3.4469
% pN/nm. If instead the 3 crossbridges attach to 6 actin filaments (by
% allowing half a crossbridge to join to each actin filament) the stiffness
% is 0.6591-3.4634 pN/nm.

kx1 = (numAttachedCrossbridges/6).*kx;
ka1 = kaA;

%c: compliance, the reciprocal of stiffness
cx1 = 1./kx1;
ca1 = 1./ka1;
cax1= cx1+ca1;
kax1 = 1./cax1;

kax = 6.*kax1;

%Total load path stiffness
cax = 1./kax;
cm = 1./kmA;
cp = cax+cm;
kp = 1./cp;

stiffness=kp;

