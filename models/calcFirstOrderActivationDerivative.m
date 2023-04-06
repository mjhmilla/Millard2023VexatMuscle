%%
% SPDX-FileCopyrightText: 2015 Matthew Millard <millard.matthew@gmail.com>
%
% SPDX-License-Identifier: Apache-2.0
%
%%

%/*-------------------------------------------------------------------------- *
%*              OpenSim:  FirstOrderMuscleActivationDynamics.cpp              *
%* -------------------------------------------------------------------------- *
%* The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
%* See http://opensim.stanford.edu and the NOTICE file for more information.  *
%* OpenSim is developed at Stanford University and supported by the US        *
%* National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
%* through the Warrior Web program.                                           *
%*                                                                            *
%* Copyright (c) 2005-2013 Stanford University and the Authors                *
%* Author(s): Thomas Uchida, Matthew Millard, Ajay Seth                       *
%*                                                                            *
%* Licensed under the Apache License, Version 2.0 (the "License"); you may    *
%* not use this file except in compliance with the License. You may obtain a  *
%* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
%*                                                                            *
%* Unless required by applicable law or agreed to in writing, software        *
%* distributed under the License is distributed on an "AS IS" BASIS,          *
%* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
%* See the License for the specific language governing permissions and        *
%* limitations under the License.                                             *
%* -------------------------------------------------------------------------- */
%
% Derivative work
% Date      : March 2015
% Authors(s): Millard
% Update    : ported to Matlab.
%
function [dactivation, a] = calcFirstOrderActivationDerivative(excitation,...
                                                          activation,...
                                              activationTimeConstant,...
                                            deActivationTimeConstant,...
                                                   minimumActivation)

if(activation < 0)
  here=1;
end
%assert(activation >= 0);
                                                 
[u,uIsClampedToLB,uIsClampedToUB] = clampWithinInterval(excitation,minimumActivation,1);
[a,aIsClampedToLB,aIsClampedToUB] = clampWithinInterval(activation,minimumActivation,1);
tau = NaN;
    
if(u >= a)
    tau = activationTimeConstant * (0.5 + 1.5*activation);
else
    tau = deActivationTimeConstant / (0.5 + 1.5*activation);
end

dactivation = (u-a)/tau;

if(aIsClampedToLB && dactivation < 0)
  dactivation=0;
end
if(aIsClampedToUB && dactivation > 0)
  dactivation=0;
end
%if( activation < 1e-3 && dactivation < 0)
%  here=1;
%end
