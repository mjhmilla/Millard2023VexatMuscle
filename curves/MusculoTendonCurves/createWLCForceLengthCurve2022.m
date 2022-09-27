% -------------------------------------------------------------------------- %
%                OpenSim:  SmoothSegmentedFunctionFactory.cpp                %
% -------------------------------------------------------------------------- %
% The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  %
% See http:%opensim.stanford.edu and the NOTICE file for more information.  %
% OpenSim is developed at Stanford University and supported by the US        %
% National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    %
% through the Warrior Web program.                                           %
%                                                                            %
% Copyright (c) 2005-2012 Stanford University and the Authors                %
% Author(s): Matthew Millard                                                 %
%                                                                            %
% Licensed under the Apache License, Version 2.0 (the "License"); you may    %
% not use this file except in compliance with the License. You may obtain a  %
% copy of the License at http:%www.apache.org/licenses/LICENSE-2.0.         %
%                                                                            %
% Unless required by applicable law or agreed to in writing, software        %
% distributed under the License is distributed on an "AS IS" BASIS,          %
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   %
% See the License for the specific language governing permissions and        %
% limitations under the License.                                             %
% -------------------------------------------------------------------------- %
%
% Derivative work
% Date      : March 2022
% Authors(s): Millard
% Updates   : ported to code to Matlab
%
% If you use this code in your work please cite this paper
%
%  Millard, M., Uchida, T., Seth, A., & Delp, S. L. (2013). 
%    Flexing computational muscle: modeling and simulation of 
%    musculotendon dynamics. Journal of biomechanical engineering, 
%    135(2), 021005.
%%
function fiberForceLengthCurve = createWLCForceLengthCurve2022(...
               normLengthZero, normLengthToe, normForceToe,...
               normLengthContour, normForceFailure,...
               kZero, kLow, kToe, curviness,...
               computeIntegral, flag_useWLCTitinModel, muscleName,...
               flag_usingOctave)
%%

%This function will generate a C2 continuous curve that fits a fiber's 
%tensile force length curve.
%
%@param normLengthZero The normalized fiber length at which the fiber 
%      begins to develop force.
%
%@param normLengthToe The normalized fiber length at which the fiber 
%     transitions from the toe region, to a linear region.
%
%@param normForceToe The normalized force developed at a length of 
%                    normLengthToe
%
%@param normLengthContour: the normalized contour length (as defined by a 
%           worm-like chain model)
% 
%@param normFailureForce: the normalized force at which the muscle is
%       expected to mechanically fail. At the present time I'm using a
%       value of 5.14 as this is the average failure force that appears in
%       Leonard, Joumaa, Herzog 2010. 
%
%@param kZero   The normalized stiffness (or slope) of the curve 
%			  at normLengthZero
%
%@param kLow   The normalized stiffness (or slope) of the fiber curve 
%			  close to the location where the force-length curve 
%			  approaches a normalized force of 0. This is usually 
%			  chosen to be a small, but non-zero fraction of kToe 
%			  (kLow = 0.025 kToe is typical).
%
%@param kToe   The normalized stiffness (or slope) of the fiber curve 
%			  at a length of normLengthToe 
%
%@param curviness    The dimensionless 'curviness' parameter that 
%					can vary between 0 (a line) to 1 (a smooth, but 
%					sharply bent elbow)
%
%@param computeIntegral  If this is true, the integral for this curve
%						is numerically calculated and splined. If false, 
%						this integral is not computed, and a call to 
%						.calcIntegral will throw an exception
%
% @param curveName The name of the muscle this curve applies to. This 
%				  curve name should have the name of the muscle and the
%				  curve in it (e.g. "bicep_fiberForceLengthCurve") 
%				  so that if this curve ever causes an exception, a 
%				  userfriendly error message can be displayed to the
%				  end user to help them debug their model.
%
%@throws exception unless the following conditions are met
%	-normLengthToe > normLengthZero            
%	-kToe > 1/(normLengthToe-normLengthZero)
%	-0 < kLow < kToe
%	-0 <= curviness <= 1
%
%@return fiberForceLengthCurve 
%       A structure that the function calcNormalizedMuscleCurveDerivative 
%       can use to evaluate the active force length curve value
%       or up to the 3rd derivative.
%
%<B>Example:</B>
%@code
%	 normLengthToe      = 1.6;
%	 normLengthZero     = 1.0;
%	 kToe      = 4.0/(normLengthToe-normLengthZero);
%	 kNearZero = 0.025*kToe
%	 c         = 0.5;
%
%	 fiberFLCurve = createFiberForceLengthCurve(normLengthZero, normLengthToe,
%								  kLow, kToe, c, true,"test");
%@endcode
%%
fiberForceLengthCurve = [];
fiberForceLengthCurve.name = sprintf('%s.%s',muscleName,'fiberForceLengthCurve');

%%
%Check the input arguments
%%
assert( normLengthToe > normLengthZero , ...
    sprintf('%s: The following must hold: normLengthToe  > normLengthZero',fiberForceLengthCurve.name));

assert( kToe > (normForceToe/(normLengthToe-normLengthZero)) , ...
    sprintf('%s: kToe must be greater than 1/(normLengthToe-normLengthZero) (%f)',...
    fiberForceLengthCurve.name, (1.0/(normLengthToe-normLengthZero))));

assert(kLow > 0.0 && kLow < normForceToe/(normLengthToe-normLengthZero),...
    sprintf('%s: kLow must be greater than 0 and less than or equal to 1',...
    fiberForceLengthCurve.name));

assert( (curviness>=0 && curviness <= 1),...      
    sprintf('%s: curviness must be between 0.0 and 1.0',...
            fiberForceLengthCurve.name));
assert(normLengthContour > normLengthToe, ...
       sprintf('%s: normContourLength must be greater than normLengthToe',...
            fiberForceLengthCurve.name))

%%
%Translate the user parameters to quintic Bezier points
%%
c = scaleCurviness(curviness);
xZero = normLengthZero;
yZero = kZero*xZero;

xIso = normLengthToe;
yIso = normForceToe;

deltaX = min(0.2*(normForceToe/kToe), 0.2*(xIso-xZero));

xLow     = xZero + deltaX;
xfoot    = xZero + 0.5*(xLow-xZero);
yfoot    = yZero;
yLow     = yfoot + kLow*(xLow-xfoot);


p01 = calcQuinticBezierCornerControlPoints(xZero, yZero,kZero, 0, ...
                                            xIso, yIso, kToe, 0,c);

if(flag_useWLCTitinModel==0)
    xpts = [p01(:,1)];
    ypts = [p01(:,2)];
        
    %Create the curve structure
    fiberForceLengthCurve.xpts    = xpts;
    fiberForceLengthCurve.ypts    = ypts;
    
    fiberForceLengthCurve.xEnd         = [xZero, xIso];
    fiberForceLengthCurve.yEnd         = [yZero, yIso];
    fiberForceLengthCurve.dydxEnd      = [kZero, kToe];
    fiberForceLengthCurve.d2ydx2End    = [0, 0];
    
    fiberForceLengthCurve.integral = [];

else
    xo=xIso-(yIso/kToe);
    
    d = 1;
    dSq = d*d;
    
    yWlc = calcWormLikeChainModelDer(xIso-xo,normLengthContour, d,dSq,[0,0]);
    
    d = yIso/yWlc;
    dSq=d*d;
    
    xFailure = calcWormLikeChainModelInvDer(normForceFailure,normLengthContour,...
                    d,dSq,[0,0]);
    xFailure = xFailure + xo;
    
    yFailure = calcWormLikeChainModelDer(xFailure-xo, ...
                    normLengthContour,d,dSq,[0,0]);
    
    kFailure = calcWormLikeChainModelDer(xFailure-xo, ...
                    normLengthContour,d,dSq,[1,0]);

    cWlc = 0.5;
    p12 = calcQuinticBezierCornerControlPoints(xIso, yIso,kToe, 0, ...
                               xFailure, yFailure, kFailure, 0,cWlc);
    
    xpts = [p01(:,1),p12(:,1)];
    ypts = [p01(:,2),p12(:,2)];
    
    xEnd = xFailure;
    yEnd = yFailure;
    kEnd = kFailure;
    
    
    %Create the curve structure
    fiberForceLengthCurve.xpts    = xpts;
    fiberForceLengthCurve.ypts    = ypts;
    
    fiberForceLengthCurve.xEnd         = [xZero, xEnd];
    fiberForceLengthCurve.yEnd         = [yZero, yEnd];
    fiberForceLengthCurve.dydxEnd      = [kZero, kEnd];
    fiberForceLengthCurve.d2ydx2End    = [0, 0];
    
    fiberForceLengthCurve.integral = [];
    
    fiberForceLengthCurveL=fiberForceLengthCurve;
    fiberForceLengthCurveR=fiberForceLengthCurve;
    %Iterate over the second curviness parameter to get a better fit with the 
    %WLC model
    xWlc = xIso + 0.5*(xFailure-xIso);
    delta = 0.25;
    fWlc = calcWormLikeChainModelDer(xWlc-xo, ...
                    normLengthContour,d,dSq,[0,0]);
    errBest = abs(calcBezierYFcnXDerivative(xWlc,fiberForceLengthCurve,0)...
                   -fWlc);
    
    for i=1:1:10
        cL = cWlc-delta;
        p12 = calcQuinticBezierCornerControlPoints(xIso, yIso,kToe, 0, ...
                               xFailure, yFailure, kFailure, 0,cL);
        fiberForceLengthCurveL.xpts = [p01(:,1),p12(:,1)];
        fiberForceLengthCurveL.ypts = [p01(:,2),p12(:,2)];
    
        errL = abs(calcBezierYFcnXDerivative(xWlc,fiberForceLengthCurveL,0)...
                   -fWlc);
        
        cR = cWlc+delta;
        p12 = calcQuinticBezierCornerControlPoints(xIso, yIso,kToe, 0, ...
                               xFailure, yFailure, kFailure, 0,cR);
        fiberForceLengthCurveR.xpts = [p01(:,1),p12(:,1)];
        fiberForceLengthCurveR.ypts = [p01(:,2),p12(:,2)];
    
        errR = abs(calcBezierYFcnXDerivative(xWlc,fiberForceLengthCurveR,0)...
                   -fWlc);
    
        if(errL < errR && errL < errBest)
            cWlc = cL;
            errBest=errL;
        end
        if(errR < errL && errR < errBest)
            cWlc = cR;
            errBest=errR;
        end
    
        delta=delta*0.5;
    
    end
    
    p12 = calcQuinticBezierCornerControlPoints(xIso, yIso,kToe, 0, ...
                           xFailure, yFailure, kFailure, 0,cWlc);
    fiberForceLengthCurve.xpts = [p01(:,1),p12(:,1)];
    fiberForceLengthCurve.ypts = [p01(:,2),p12(:,2)];

end

if(computeIntegral == 1)
    xScaling = normLengthToe;
    fiberForceLengthCurve.integral = ...
        createCurveIntegralStructure(fiberForceLengthCurve, ...
                                     1000,...
                                     1e-12,...
                                     xScaling, flag_usingOctave,1);    
end
                                   