%%
% SPDX-FileCopyrightText: 2015 Matthew Millard <millard.matthew@gmail.com>
%
% SPDX-License-Identifier: Apache-2.0
%
%%
% -------------------------------------------------------------------------- %
%                OpenSim:  SmoothSegmentedFunctionFactory.cpp                %
% -------------------------------------------------------------------------- %
% The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  %
% See http://opensim.stanford.edu and the NOTICE file for more information.  %
% OpenSim is developed at Stanford University and supported by the US        %
% National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    %
% through the Warrior Web program.                                           %
%                                                                            %
% Copyright (c) 2005-2012 Stanford University and the Authors                %
% Author(s): Matthew Millard                                                 %
%                                                                            %
% Licensed under the Apache License, Version 2.0 (the "License"); you may    %
% not use this file except in compliance with the License. You may obtain a  %
% copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         %
%                                                                            %
% Unless required by applicable law or agreed to in writing, software        %
% distributed under the License is distributed on an "AS IS" BASIS,          %
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
function c = scaleCurviness(curviness)
%%
%Scales a [0,1] curviness value to be within (0,1).
%This function is widely used to create the Bezier corner elements 
%used throughout the normalized muscle curves, and so while this is a very
%simple function, it is important that it be consistently used throughout.
%
% @param curviness: a [0,1] number
% @returns c      : a (0,1) number
%%
c = 0.1 + 0.8*curviness;
