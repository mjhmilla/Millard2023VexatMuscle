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
% Licensed under the Apache License, Version 2.0 (the 'License'); you may    %
% not use this file except in compliance with the License. You may obtain a  %
% copy of the License at http:%www.apache.org/licenses/LICENSE-2.0.         %
%                                                                            %
% Unless required by applicable law or agreed to in writing, software        %
% distributed under the License is distributed on an 'AS IS' BASIS,          %
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   %
% See the License for the specific language governing permissions and        %
% limitations under the License.                                             %
% -------------------------------------------------------------------------- %
%
% Derivative work
% Date      : March 2015
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
function tendonStiffnessCurve = ...
          createTendonStiffnessCurve2021(tendonForceLengthCurve,...
                                        curviness,...
                                        muscleName, flag_usingOctave)
%%

%Will generate a C2 continous (continuous to the second derivative) 
%curve in a MuscleFunctionObject object that fits the derivative of the
%tendon force lenth curve
%
%
%%


xV0 = zeros(1,size(tendonForceLengthCurve.xpts,2)+1);

for i=1:1:size(tendonForceLengthCurve.xpts,2)
  xV0(1,i)=tendonForceLengthCurve.xpts(1,i);
end
xV0(1,end)=tendonForceLengthCurve.xpts(end,end);



xV=xV0;

yVal = zeros(size(xV));
dydxVal = zeros(size(xV));

for i=1:1:length(xV)
  yVal(1,i)    = calcBezierYFcnXDerivative(xV(1,i),tendonForceLengthCurve,1);
  dydxVal(1,i) = calcBezierYFcnXDerivative(xV(1,i),tendonForceLengthCurve,2);  
end

    
         
xpts = zeros(6,length(xV)-1);          
ypts = zeros(6,length(xV)-1);          


for i=1:1:size(xpts,2)
  c = scaleCurviness(curviness); 
  p01 = calcQuinticBezierCornerControlPoints(...
            xV(1,i), yVal(1,i), dydxVal(1,i), 0,...
            xV(1,i+1), yVal(1,i+1), dydxVal(1,i+1), 0, c);
  xpts(:,i)=p01(:,1);
  ypts(:,i)=p01(:,2);
end


tendonStiffnessCurve.xpts    = xpts;
tendonStiffnessCurve.ypts    = ypts;

tendonStiffnessCurve.xEnd         = [xV(1,1), xV(1,end)];
tendonStiffnessCurve.yEnd         = [yVal(1,1), yVal(1,end)];
tendonStiffnessCurve.dydxEnd      = [dydxVal(1,1), dydxVal(1,end)];
tendonStiffnessCurve.d2ydx2End    = [0, 0];

tendonStiffnessCurve.integral = [];




