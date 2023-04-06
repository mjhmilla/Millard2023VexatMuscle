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

function [forceLengthPevkIgdCurve, forceLengthIgpCurve] = ...
  constructTitinCurves(lambdaECM, ...
                        normStiffnessHalfPevkIgd, loptPevkIgd, ...
                        normHalfStiffnessIgp, loptIgp, ...
                        normLengthZToT12, normLengthHalfFixedIgd,...
                        normFiberLengthAtOnePassiveForce,... 
                        flag_useOctave)


x0    = -loptPevkIgd*0.025;
y0    = sqrt(eps);
dydx0 = sqrt(eps);

dydx1 = normStiffnessHalfPevkIgd;        
x1    = loptPevkIgd*0.025;
y1    = dydx1*(x1);

computeIntegral = 0;
forceLengthPevkIgdCurve = createC2ContinuousRamp(x0,y0,dydx0,x1,y1,dydx1,...
                                computeIntegral,'PevkIgd',flag_useOctave);

forceLengthPevkIgdInverseCurve = createInverseCurve(forceLengthPevkIgdCurve);

x0    = 0;
y0    = sqrt(eps);
dydx0 = sqrt(eps);

lPevkIgDRef = calcBezierYFcnXDerivative((1-lambdaECM),...
                    forceLengthPevkIgdInverseCurve,0);

lTiFlex = normFiberLengthAtOnePassiveForce*0.5 ...
          - normLengthHalfFixedIgd - normLengthZToT12;

lIgPRef = lTiFlex-lPevkIgDRef;

x2    = lIgPRef;   
dydx2 = normHalfStiffnessIgp;   
y2    = (1-lambdaECM);

x12 = x2 - y2/dydx2;


x1    = x0 + 0.5*(x12-x0);   
dydx1 = normHalfStiffnessIgp*0.005;
y1 = dydx1*(x1-x0);





computeIntegral=0;
curveName = sprintf('Igp');  
c=0.75;
forceLengthIgpCurve = createProximalIgCurve(x0,y0,dydx0,...
                                            x1,y1,dydx1,...
                                            x2,y2,dydx2, c,...                                                
                            computeIntegral,curveName,flag_useOctave);


