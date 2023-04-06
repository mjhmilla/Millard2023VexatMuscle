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
  createFiberActiveForceLengthCurve(  normMyosinLength, ...
                                      normMyosinBareLength, ...
                                      normActinLength, ...
                                      normZLineThickness,...
                                      normSarcomereLengthZeroForce,...
                                      normCrossbridgeStiffness,...
                                      curviness, ...
                                      shiftLengthActiveForceLengthCurveDescendingCurve,...
                                      flag_compensateForCrossbridgeStiffness,...
                                      flag_enableNumericallyNonZeroGradient, ...
                                      smallNumericallyNonZeroNumber,...
                                      computeIntegral, ...
                                      muscleName, ...
                                      flag_usingOctave)




activeForceLengthCurve = [];
activeForceLengthCurve.name = sprintf('%s.%s',muscleName,'activeForceLengthCurve');




%%
%Check inputs
%%
%rootEPS = eps^0.5;
tol = 1e-6;

%Check signs
assert(normMyosinLength     >  0);
assert(normMyosinBareLength >= 0);
assert(normActinLength      >  0);
assert(normZLineThickness   >  0);
assert( normSarcomereLengthZeroForce > 0);
assert(curviness >= 0 && curviness <= 1);

assert( normActinLength > 0.5*normMyosinLength,... 
       ['Error: this curve construction function assumes that actin',...
        ' is longer than half the length of myosin']);

%Check normalization
assert(abs(normActinLength*2 + normZLineThickness*2 + normMyosinBareLength - 1.0) <= tol,...
       ['Sarcomere geometry incorrectly normalized: ',...
        '2*(normActinLength+normZLineThickness)+ normMyosinBareLength==1'] )
assert(normSarcomereLengthZeroForce < 1.0);
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

dx45 = shiftLengthActiveForceLengthCurveDescendingCurve;

c0y    = 0;
c4y    = 0;

p0DyDx = 0;
p5DyDx = 0;

% if(flag_enableNumericallyNonZeroGradient == 1)
%   c0y    =  10*smallNumericallyNonZeroNumber;
%   c4y    =  10*smallNumericallyNonZeroNumber;
% 
%   p0DyDx =  smallNumericallyNonZeroNumber;  
%   p5DyDx = -smallNumericallyNonZeroNumber;
% end


% Corner 0: zero force happens ... for some (largely) unknown reason.
%
%         : Originally it was thought that the myosin ends are developing 
%           compressive forces as they contact the z-line. There exists data
%           showing that the myosin tips can pass through the z-line (mentioned
%           by Christoph Rode in conversation).
%
%         : Christoph Rode also mentioned that when myosin protrudes into its
%           neighbors sarcomere that the cross-bridges can still attach to
%           its neighbors actin filaments. This would reduce the available
%           sites for attachment and attentuate the tension developed. In this
%           case a tension of zero is reached as the titin segments get 
%           stretched. This cannot (to me) explain the rapid drop in force: 
%           titin is too compliant for this I think.
%
%         : Another idea is that perhaps the cross-bridges are not able to 
%           reach the actins in a short sarcomere because it has expanded too
%           much (due to the constant volume property) which is an idea put 
%           forward by Robert Rockenfeller & Michael Gunther. This does explain
%           the rapid decrease in force with length.
%
%    |--------------------------------
%    |--------------------------------|
%    |                                |
% ||=================----===============||
%    |                                |
%    |--------------------------------|
%     --------------------------------|
% 
c0x = normSarcomereLengthZeroForce;
%c0y : already set


% Corner 1: myosin tips touching the z-line
%  
% |--------------------------------
% |      --------------------------------|
% |                                      |
% ||=================----===============||
% |                                      |
% |--------------------------------      |
%        --------------------------------|
%
c1x = normMyosinLength + 2*normZLineThickness;


%At the point where the myosin ends just touch the Z-line evaluate the 
%length of the section where the actins from both ends of the sarcomere
%overlap: here assume cross-bridges have a 50% chance of attaching to the
%correct actin element. The 50% that attach to the wrong actin generate
%no force because they would be compressing actin. Since actin is a short
%stiff fiber (but not a rod) it cannout support compression and this section
%at the end would go slack. The 50% that attach to the correct actin generate
%force.
normShallowPlateauInterference = 2*(normActinLength ...
                                   - 0.5*normMyosinLength ...
                                   - 0.5*normMyosinBareLength);

normShallowPlateauOverlap      = 2*(normMyosinLength - normActinLength);

normMaxOverlap                 = normMyosinLength - normMyosinBareLength;

%With half of the cross-bridges pulling in one direction and the other half 
%pulling in the opposite direction the interference section contributes no force

netNormInterferenceTension     = 0.5; 

c1y = ( normShallowPlateauInterference*netNormInterferenceTension ...
      + normShallowPlateauOverlap )/normMaxOverlap;

% if(flag_compensateForCrossbridgeStiffness==1)
%   c1x = c1x - c1y/(c1y*normCrossbridgeStiffness);
% end


% Corner 3: maximum overlap, short end of the plateau
%  
% |--------------------------------
% |                            --------------------------------|
% |                                                            |
% |          |=================----===============|            |
% |                                                            |
% |--------------------------------                            |
%                              --------------------------------|
%
c2x = 1.0 - normMyosinBareLength; %The
c2y = 1.0;

% if(flag_compensateForCrossbridgeStiffness==1)
%   c2x = c2x - c2y/(c2y*normCrossbridgeStiffness);
% end

% Corner 4: maximum overlap, long end of the plateau
%  
% |--------------------------------
% |                                    --------------------------------|
% |                                                                    |
% |              |=================----===============|                |
% |                                                                    |
% |--------------------------------                                    |
%                                      --------------------------------|
%
c3x = 1.0 + normMyosinBareLength; %The
c3y = 1.0;

% if(flag_compensateForCrossbridgeStiffness==1)
%   c3x = c3x - c3y/(c3y*normCrossbridgeStiffness);
% end

% Corner 3: Overlap is lost
%
% |--------------------------------
% |                                                                     --------------------------------|
% |                                                                                                     |
% |                                |=================----===============|                               |
% |                                                                                                     |
% |--------------------------------                                                                     |
%                                                                       --------------------------------|

c4x = 2*normZLineThickness + 2*normActinLength + normMyosinLength;
%c4y : already set 



%%
% Calculate the control point locations and slopes
%%
c0c1x = (c1x-c0x);
c1c2x = (c2x-c1x);
c2c3x = (c3x-c2x);
c3c4x = (c4x-c3x);

c0c1y = (c1y-c0y);
c1c2y = (c2y-c1y);
c2c3y = (c3y-c2y);
c3c4y = (c4y-c3y);


shoulder = 0.125*(c4x - c0x);

p0x   = c0x        - 0.5*shoulder + dx45;
p0y   = 0;%c0y - p0DyDx*(0.5*shoulder);

if(flag_enableNumericallyNonZeroGradient==0)
  p0y    = smallNumericallyNonZeroNumber;
  p0DyDx = smallNumericallyNonZeroNumber/10.;    
end    

p1x     = c0x + 0.5*c0c1x;  
p1y     = c0y + 0.5*c0c1y;
p1DyDx  = (c0c1y)/(c0c1x);

p2x     = c1x + 0.5*c1c2x;
p2y     = c1y + 0.5*c1c2y;
p2DyDx  = (c1c2y)/(c1c2x);

p3x     = c2x + 0.5*c2c3x;
p3y     = c2y + 0.5*c2c3y;
p3DyDx  = (c2c3y)/(c2c3x);

if(flag_compensateForCrossbridgeStiffness==0)
    assert( abs(p3x - 1.0) < 1e-6);
    assert( abs(p3y - 1.0) < 1e-6);
    assert( abs(p3DyDx )   < 1e-6);
end

p3x     = 1.;
p3y     = 1.;
p3DyDx  = 0.;

p4x     = c3x + 0.5*c3c4x + dx45;
p4y     = c3y + 0.5*c3c4y;
p4DyDx  = (c3c4y)/(c3c4x);


p5x     = c4x + 0.5*shoulder + dx45;
p5y     = c4y + p5DyDx*(0.5*shoulder);
if(flag_enableNumericallyNonZeroGradient==0)
  p5y    = smallNumericallyNonZeroNumber;
  p5DyDx =-smallNumericallyNonZeroNumber/10.;    
end    

activeForceLengthCurveAnnotationPoints = struct('x',[],'y',[]);
activeForceLengthCurveAnnotationPoints.x ...
                    = [c0x;...
                      c1x;...
                      c2x;...
                      c3x;...
                      c4x];
activeForceLengthCurveAnnotationPoints.y ...
                    = [c0y;...
                      c1y;...
                      c2y;...
                      c3y;...
                      c4y];                    

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



        