

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

function [activeForceLengthCurve, ...
  activeForceLengthCurveAnnotationPoints]= ...
  createSiebert2015FiberActiveForceLengthCurve(l1N, l2N, l3N, l4N, fcN,...
                                      normCrossbridgeStiffness,...
                                      curviness, ...
                                      flag_compensateForCrossbridgeStiffness,...
                                      flag_enableNumericallyNonZeroGradient, ...
                                      smallNumericallyNonZeroNumber,...
                                      computeIntegral, ...
                                      muscleName, ...
                                      flag_usingOctave)




activeForceLengthCurve = [];
activeForceLengthCurve.name = sprintf('%s.%s',muscleName,'activeForceLengthCurveSiebert2015');




%%
%Check inputs
%%
%rootEPS = eps^0.5;
tol = 1e-6;

%Check signs
assert(curviness >= 0 && curviness <= 1);

%%
%Create the curve
%%
c = scaleCurviness(curviness);

% Ascii force-length curve
%                  2 3 3
%             1  2 c_p_c           p : control point
%             c_ p       \         c : corner
%           1/             \ 4
%           p                p
%      0  0/                   \4  5
%------p--c                     c--p---------
%       ~0.5       1.0         ~1.6
%



dxpN = 0.5*(l3N-1); %Shift to account for the fact that lopt is 
                    %defined differently for Siebert 2015 (left most point 
                    % of plateau) and the active force length curves here
                    %(center of plateau).


l1N = l1N-dxpN;
l2N = l2N-dxpN;
loptSN = 1-dxpN; %Siebert et al.'s normalized optimal fiber length 
                 %expressed using our definition here
l3N = l3N-dxpN;
l4N = l4N-dxpN;


shoulder = 0.05*(l4N - l1N);



p0x    = l1N - shoulder;
p0y    = 0;
p0DyDx = 0;

if(flag_enableNumericallyNonZeroGradient==0)
  p0y    = smallNumericallyNonZeroNumber;
  p0DyDx = smallNumericallyNonZeroNumber/10.;    
end    

p1x    = l1N + 0.5*(l2N-l1N);
p1y    = fcN*0.5;
p1DyDx = 1.1*(fcN)/(l2N-l1N);

p2x    = l2N + 0.5*(loptSN-l2N);
p2y    = fcN + 0.5*(1-fcN);
p2DyDx = (1-fcN)/(loptSN-l2N);

p3x = 1;
p3y = 1;
p3DyDx=0;

p4x     = l3N + 0.5*(l4N-l3N);
p4y     = 1 + 0.5*(0-1);
p4DyDx  = -1.1/(l4N-l3N);

p5x     = l4N + 0.5*shoulder;
p5y     = 0;
p5DyDx  = 0;

if(flag_enableNumericallyNonZeroGradient==0)
  p5y    = smallNumericallyNonZeroNumber;
  p5DyDx =-smallNumericallyNonZeroNumber/10.;    
end   


if(flag_compensateForCrossbridgeStiffness==0)
    assert( abs(p3x - 1.0) < 1e-6);
    assert( abs(p3y - 1.0) < 1e-6);
    assert( abs(p3DyDx )   < 1e-6);
end
 

activeForceLengthCurveAnnotationPoints = struct('x',[],'y',[]);
activeForceLengthCurveAnnotationPoints.x ...
                    = [l1N;...
                      l2N;...
                      (1.0-(l3N-1.0));...
                      1.0;...
                      l3N;...
                      l4N];
activeForceLengthCurveAnnotationPoints.y ...
                    = [0;...
                      fcN;...
                      1;...
                      1;...
                      1;...
                      0];                    

if(flag_compensateForCrossbridgeStiffness==1)
    %Since the stiffness is f*Kx, the strain introduced by a half
    %cross bridge under an active force f is f/(f*Kx) or 1/Kx. For the 
    %entire sarcomere its 2/Kx.
    p1x = p1x - 2/(normCrossbridgeStiffness);
    p2x = p2x - 2/(normCrossbridgeStiffness);
    p3x = p3x - 2/(normCrossbridgeStiffness);
    p4x = p4x - 2/(normCrossbridgeStiffness);    
end

%Compute the locations of the control points
b0 = calcQuinticBezierCornerControlPoints( p0x, p0y,  p0DyDx, 0, ...
                                           p1x, p1y,  p1DyDx, 0, c);
b1 = calcQuinticBezierCornerControlPoints( p1x, p1y,  p1DyDx, 0, ...
                                           p2x, p2y,  p2DyDx, 0, c);
b2 = calcQuinticBezierCornerControlPoints( p2x, p2y,  p2DyDx, 0, ...
                                           p3x, p3y,  p3DyDx, 0, c);
b3 = calcQuinticBezierCornerControlPoints( p3x, p3y,  p3DyDx, 0, ...
                                           p4x, p4y,  p4DyDx, 0, c);
b4 = calcQuinticBezierCornerControlPoints( p4x, p4y,  p4DyDx, 0, ...
                                           p5x, p5y,  p5DyDx, 0, c);
                                         
   
xpts = [b0(:,1) b1(:,1) b2(:,1) b3(:,1) b4(:,1)];
ypts = [b0(:,2) b1(:,2) b2(:,2) b3(:,2) b4(:,2)];

%Create the curve structure
activeForceLengthCurve.xpts    = xpts;
activeForceLengthCurve.ypts    = ypts;

activeForceLengthCurve.xEnd         = [p0x,     p5x];
activeForceLengthCurve.yEnd         = [p0y,     p5y];
activeForceLengthCurve.dydxEnd      = [p0DyDx,  p5DyDx];
activeForceLengthCurve.d2ydx2End    = [0, 0];
activeForceLengthCurve.integral = [];



        