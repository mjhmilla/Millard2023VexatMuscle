%%
% SPDX-FileCopyrightText: 2015 Matthew Millard <millard.matthew@gmail.com>
%
% SPDX-License-Identifier: Apache-2.0
%
%%
%/* -------------------------------------------------------------------------- *
% *                           singleMuscleBench.cpp                            *
% * -------------------------------------------------------------------------- *
% * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
% * See http://opensim.stanford.edu and the NOTICE file for more information.  *
% * OpenSim is developed at Stanford University and supported by the US        *
% * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
% * through the Warrior Web program.                                           *
% *                                                                            *
% * Copyright (c) 2005-2012 Stanford University and the Authors                *
% * Author(s): Matthew Millard                                                 *
% *                                                                            *
% * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
% * not use this file except in compliance with the License. You may obtain a  *
% * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
% *                                                                            *
% * Unless required by applicable law or agreed to in writing, software        *
% * distributed under the License is distributed on an "AS IS" BASIS,          *
% * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
% * See the License for the specific language governing permissions and        *
% * limitations under the License.                                             *
% * -------------------------------------------------------------------------- */
%
% Derivative work
% Date      : March 2015
% Authors(s): Millard
% Updates   : Ported to code to Matlab which included some major rework
%
% If you use this code in your work please cite this paper
%
%  Millard, M., Uchida, T., Seth, A., & Delp, S. L. (2013). 
%    Flexing computational muscle: modeling and simulation of 
%    musculotendon dynamics. Journal of biomechanical engineering, 
%    135(2), 021005.
%%
function benchRecord = ...
    runPrescribedLengthActivationSimulation(calcMuscleInfoFcn,...   
                                         calcInitialMuscleStateFcn,...
                                         benchConfig,...
                                         benchRecord,...
                                         simulationIndex,...
                                         numberOfSimulations,...
                                         flag_appendEnergetics,...
                                         flag_useOctave)
%%
% This function will run a single prescribed-length and prescribed-activation
% muscle simulation and put its results in the simulationIndex^th column of
% the matrices in benchRecord. 
%
%
% @param calcMuscleInfoFcn: a function handle to a muscle function
%                      like (calcMillard2012DampedEquilibriumMuscleInfo.m) 
%                      and takes arguments of activation, pathState, and 
%                      muscle state
%
% @param calcInitialMuscleStateFcn: a function handle to a function that
%                                 takes arguments of activation, path
%                                 state, calcMuscleInfoFcn, and
%                                 a struct initConfig (see
%                                 calcInitialMuscleState.m for details)
%                                 and returns a muscle state that satisfies
%                                 the equilibrium equation
%
% @param benchConfig: a structure that configures this constant activation
% sinusoidal stretch benchmark simulation:
%
%   benchConfig.npts                 : number of points to evaluate the
%                                      model at during a 1 second constant 
%                                      activation sinusoidal stretch 
%                                      simulation
%
%   benchConfig.numberOfMuscleStates : 0 for rigid tendon, 
%                                      1 for elastic tendon
%
%   benchConfig.minimumActivation    : 0 unless the classic elastic tendon
%                                      model is being used, then it should 
%                                      be > 0, and probably 0.05
%
%   benchConfig.name                 : short musclemodel name e.g.
%                                      'soleusRigidTendon'   
%
% @param benchRecord : see return struct for details.
%
% @param simulationIndex: the current simulation index - this is just used
%                         for indexing the correct column for the matrices 
%                         in benchRecord.
%
% @param numberOfSimulations: total number of simulations in this series.
%                             This variable is used for sizing the number of
%                             columns in the matricies that appear in 
%                             benchRecord
%
% @param figBasicInfo   : handle to an empty figure
%
% @param figEnergyInfo  : handle to an empty figure
%
% @param figPowerInfo   : handle to an empty figure
%
% @return benchRecord: A structure containing a set of matricies, where
%                      each column represents the data of 1 
%                      constant-activation sinusoidal-stretch simulation.
% 
%     benchRecord.activation                  
%     benchRecord.cpuTime                     
%     benchRecord.normFiberForceAlongTendon   
%     benchRecord.normFiberLength              
%     benchRecord.pennationAngle              
%     benchRecord.normFiberVelocity           
%     benchRecord.pennationAngVelocity        
%     benchRecord.fiberStiffnessAlongTendon   
%     benchRecord.tendonStiffness  
%     benchRecord.muscleStiffness             
%     benchRecord.pathLength                  
%     benchRecord.pathVelocity                
% 
%     benchRecord.dSystemEnergyLessWork       
%     benchRecord.systemEnergyLessWork        
%     benchRecord.tendonPotentialEnergy       
%     benchRecord.fiberPotentialEnergy        
%     benchRecord.fiberActiveWork             
%     benchRecord.dampingWork                 
%     benchRecord.boundaryWork                
% 
%     benchRecord.tendonPower                 
%     benchRecord.fiberParallelElementPower   
%     benchRecord.fiberActivePower            
%     benchRecord.dampingPower                
%     benchRecord.boundaryPower               
%
%%

%%
%
% Initialize the output structure if needed
%
%%
npts = benchConfig.npts;

if(isempty(benchRecord)==1)
  benchRecord.time                        = zeros(npts,numberOfSimulations);
  benchRecord.activation                  = zeros(npts,numberOfSimulations);
  benchRecord.activationDot               = zeros(npts,numberOfSimulations);  
  benchRecord.cpuTime                     = zeros(numberOfSimulations);
  benchRecord.normFiberForceAlongTendon   = zeros(npts,numberOfSimulations);
  benchRecord.normFiberForce              = zeros(npts,numberOfSimulations);
  benchRecord.activeFiberForceAlongTendon = zeros(npts,numberOfSimulations);
  benchRecord.activeFiberForce            = zeros(npts,numberOfSimulations);   
  benchRecord.normActiveFiberForce        = zeros(npts,numberOfSimulations); 
  benchRecord.normPassiveFiberForce        = zeros(npts,numberOfSimulations);   
  benchRecord.fiberForceAlongTendon       = zeros(npts,numberOfSimulations);
  benchRecord.tendonForce                 = zeros(npts,numberOfSimulations);
  benchRecord.fiberForce                  = zeros(npts,numberOfSimulations);  
  benchRecord.normFiberLength             = zeros(npts,numberOfSimulations);    
  benchRecord.pennationAngle              = zeros(npts,numberOfSimulations);
  benchRecord.normFiberVelocity           = zeros(npts,numberOfSimulations);
  benchRecord.pennationAngularVelocity    = zeros(npts,numberOfSimulations);
  %benchRecord.fiberStiffnessAlongTendon   = zeros(npts,numberOfSimulations);
  %benchRecord.fiberDampingAlongTendon     = zeros(npts,numberOfSimulations);  
  benchRecord.tendonStiffness             = zeros(npts,numberOfSimulations);
  benchRecord.tendonDamping               = zeros(npts,numberOfSimulations);
  
  benchRecord.crossBridgeStiffness        = zeros(npts,numberOfSimulations);
  benchRecord.crossBridgeDamping          = zeros(npts,numberOfSimulations);   
  benchRecord.pathLength                  = zeros(npts,numberOfSimulations);
  benchRecord.pathVelocity                = zeros(npts,numberOfSimulations);

  benchRecord.dSystemEnergyLessWork       = zeros(npts,numberOfSimulations);
  benchRecord.systemEnergyLessWork        = zeros(npts,numberOfSimulations);
  benchRecord.tendonPotentialEnergy       = zeros(npts,numberOfSimulations);
  benchRecord.fiberPotentialEnergy        = zeros(npts,numberOfSimulations);
  benchRecord.fiberActiveWork             = zeros(npts,numberOfSimulations);
  benchRecord.dampingWork                 = zeros(npts,numberOfSimulations);
  benchRecord.boundaryWork                = zeros(npts,numberOfSimulations);

  benchRecord.tendonPower                 = zeros(npts,numberOfSimulations);
  benchRecord.fiberParallelElementPower   = zeros(npts,numberOfSimulations);
  benchRecord.fiberActivePower            = zeros(npts,numberOfSimulations);
  benchRecord.dampingPower                = zeros(npts,numberOfSimulations);
  benchRecord.boundaryPower               = zeros(npts,numberOfSimulations);

%   benchRecord.bristleLength               = zeros(npts,numberOfSimulations);
%   benchRecord.normBristleLength           = zeros(npts,numberOfSimulations);
%   benchRecord.slidingLength               = zeros(npts,numberOfSimulations);
%   benchRecord.normBristleSpringForce      = zeros(npts,numberOfSimulations);
%   benchRecord.normBristleSpringForceNoFal = zeros(npts,numberOfSimulations);
%   benchRecord.normBristleDampingForce     = zeros(npts,numberOfSimulations);
%   benchRecord.normBristleStiffness        = zeros(npts,numberOfSimulations);
%   benchRecord.normBristleDamping          = zeros(npts,numberOfSimulations);
%   benchRecord.bristleBlendingVariable     = zeros(npts,numberOfSimulations);
%   
%   benchRecord.crossBridgeForce           = zeros(npts,numberOfSimulations);
%   benchRecord.slowForceVelocityForce     = zeros(npts,numberOfSimulations);
%   benchRecord.fastForceVelocityForce     = zeros(npts,numberOfSimulations);

  benchRecord.musculotendonStiffness  = zeros(npts,numberOfSimulations);
  benchRecord.musculotendonDamping    = zeros(npts,numberOfSimulations);
  benchRecord.fiberStiffness  = zeros(npts,numberOfSimulations);
  benchRecord.fiberDamping    = zeros(npts,numberOfSimulations);


  benchRecord.fiberActiveForceLengthMultiplier = zeros(npts,numberOfSimulations);
  benchRecord.fiberPassiveForceLengthMultiplier= zeros(npts,numberOfSimulations);

  benchRecord.eventTime             = zeros(npts,numberOfSimulations);
  benchRecord.eventTendonForce      = zeros(npts,numberOfSimulations);
  benchRecord.eventNormFiberForce   = zeros(npts,numberOfSimulations);
  benchRecord.eventNormActiveFiberForce   = zeros(npts,numberOfSimulations);  
  benchRecord.eventNormFiberForce   = zeros(npts,numberOfSimulations);
  benchRecord.eventNormFiberLength  = zeros(npts,numberOfSimulations);
  benchRecord.eventPathLength       = zeros(npts,numberOfSimulations);  
  benchRecord.eventNormFiberVelocity= zeros(npts,numberOfSimulations);
   

  benchRecord.extra =  zeros(  npts,...
                             numberOfSimulations,...
                             10 );
  benchRecord.extraLabels = {};
  

  benchRecord.state = zeros(  npts,...
                              numberOfSimulations,...
                              benchConfig.numberOfMuscleStates);
  benchRecord.dstate = zeros(  npts,...
                              numberOfSimulations,...
                              benchConfig.numberOfMuscleStates);                            
end

%%
%
% Set up function handles
%
%%

pathFcn       = benchConfig.pathFcn;
excitationFcn = benchConfig.excitationFcn;
activationFcn = benchConfig.activationFcn;
              
dfcn = @(argt,argState)...
        calcPrescribedMusculotendonStateDerivativeWrapper(...
                          argt,...
                          argState,...                                                      
                          pathFcn,...
                          excitationFcn,...
                          activationFcn,...
                          calcMuscleInfoFcn,...
                          flag_appendEnergetics);


%%
%
% Initialize the muscle if necessary
%
%%
muscleState0 = [];

if( isempty(benchConfig.initialState)==0)
    muscleState0 =  benchConfig.initialState;
end

if(benchConfig.numberOfMuscleStates ~= 0 ...
        && isempty(benchConfig.initialState)==1)  

    initConfig.iterMax = 100;
    initConfig.tol     = 1e-8;
    initConfig.useStaticFiberSolution = 0;
    e0 = excitationFcn(0);
    a0 = e0; %Assume that the activation is such that d/dt a = 0;

    initSoln = calcInitialMuscleStateFcn([0,a0],...
                                      pathFcn(0),...                                          
                                      calcMuscleInfoFcn,...
                                      initConfig);

    if(initSoln.converged == 0 && initSoln.isClamped == 0)
        %%
        %If we're here then we've been unlucky enough to be in the
        %negative stiffness region of the active force length curve and
        %we will have to settle for an initial condition where the
        %velocity of the fiber is 0 ... which often leads to force
        %transients at the beginning of the simulation.
        %%
        initConfig.useStaticFiberSolution = 1;
        initSoln = calcInitialMuscleStateFcn([0,a0],...
                                      pathFcn(0),...
                                      calcMuscleInfoFcn,...
                                      initConfig);

        assert(initSoln.converged == 1 || initSoln.isClamped == 1,...
               'Failed to bring the muscle to a valid initial solution');
    end

    muscleState0 = initSoln.muscleState(:);       
end

initialState = [benchConfig.initialActivation; muscleState0];


%%
%
% Initialize working variables
%
%%
muscleStateV  = [];

npts    = benchConfig.npts;
tmin    = min(benchConfig.tspan(1));
tmax    = max(benchConfig.tspan(2));
tV = benchConfig.tspan;

if(length(tV) == 2)
  tV      = [tmin:((tmax-tmin)/(npts-1)):tmax];
end

assert(length(tV) == npts);

t0        = 0;
cpuTime   = 0;
xe        = [];
ye        = [];
cpuTime   = [];



%%
%
% Integrate 
%
%%

if flag_useOctave
   
  lsode_options('relative tolerance', benchConfig.relTol);
  lsode_options('absolute tolerance', benchConfig.absTol);
   
  dfcnOctave = @(argX,argT)dfcn(argT,argX);

  x0 = [];
  if(flag_appendEnergetics==1)
    x0 = [initialState(:);0;0;0];
  else
    x0 = initialState(:);
  end
  
  t0 = clock();          
  ye = lsode(dfcnOctave,x0,tV);                                                   
  xe = tV;
  cpuTime = etime(clock(),t0);
else
  options = [];
  xe = [];
  ye = [];
  eventTime = [];
  eventState = [];
  eventIndex = [];
  
  x0 = [];
  if(flag_appendEnergetics==1)
    x0 = [initialState(:);0;0;0];
  else
    x0 = initialState(:);
  end
  
  
  if(isempty(benchConfig.eventFcn) == 1)                 
    options = odeset('RelTol',benchConfig.relTol,...
                     'AbsTol',benchConfig.absTol,...
                     'NormControl','on',...
                     'MaxStep',0.001,...
                     'Stats','off');
    t0 = tic;
    [xe, ye] = ode15s(dfcn,tV,x0,options);                                                   
    cpuTime = toc(t0);  
  else
    options = odeset('RelTol',benchConfig.relTol,...
                     'AbsTol',benchConfig.absTol,...
                     'Stats','off',...
                     'NormControl','on',...
                     'MaxStep',0.001,...                     
                     'Events',benchConfig.eventFcn);    
    
    t0 = tic;
    [xe, ye, eventTime, eventState, eventIndex] = ...
      ode15s(dfcn,tV,x0,options);                                                   
    cpuTime = toc(t0);  
    
    
  end
end

activation        = zeros(size(xe));
workOfBoundary    = zeros(size(xe));
workOfActiveFiber = zeros(size(xe));
workOfDamping     = zeros(size(xe));

activation = ye(:,1);

n = benchConfig.numberOfMuscleStates;
if(n >= 1)            
    muscleStateV = ye(:,2:1:(n+1));
end

energeticsT    = [];
energeticsV    = [];
energeticsWNeg = [];

if(flag_appendEnergetics ==1)
  energeticsT    = ye(:,n+2);
  energeticsV    = ye(:,n+3);
  energeticsWNeg = ye(:,n+4);
else
  energeticsT    = zeros(length(ye(:,1)),1);
  energeticsV    = zeros(length(ye(:,1)),1);
  energeticsWNeg = zeros(length(ye(:,1)),1);  
end

if(benchConfig.numberOfMuscleStates == 0) 
    t0 = tic;
end

%%
%
% Post-process the simulation
%
%%

idxSim = simulationIndex;

if(isempty(benchConfig.eventFcn)==0)
  
  
  for j=1:1:length(eventTime)
    benchRecord.eventTime(j,idxSim) = eventTime(j,1);

    excitationEvent = excitationFcn(eventTime(j,1));
    activationEvent = eventState(j,1);
    
    dactivationEvent = activationFcn(excitationEvent, activationEvent);
    activationStateEvent = [dactivationEvent,activationEvent];
    pathStateEvent       = pathFcn(eventTime(j,1));

    
    mtInfo = calcMuscleInfoFcn(activationStateEvent,...                                   
                               pathStateEvent,...
                               eventState(j,2:1:(n+1)));
    
    benchRecord.eventNormActiveFiberForce(j,idxSim) = ...
      mtInfo.muscleDynamicsInfo.normActiveFiberForce;    
    
    benchRecord.eventNormFiberLength(j,idxSim) = ...
      mtInfo.muscleLengthInfo.normFiberLength;    
    
    benchRecord.eventPathLength(j,idxSim) = ...
      pathStateEvent(2,1); 
        
    benchRecord.eventNormFiberVelocity(j,idxSim) = ...
      mtInfo.fiberVelocityInfo.normFiberVelocity;  
    
    if(mtInfo.fiberVelocityInfo.normFiberVelocity > 0)
        here=1;
    end

    benchRecord.eventTendonForce(j,idxSim) = ...
      mtInfo.muscleDynamicsInfo.tendonForce;
    
            
    benchRecord.eventNormFiberForce(j,idxSim) = ...
      mtInfo.muscleDynamicsInfo.normFiberForce;
    
    
  end
end


T0V0 = 0;
for j=1:1:length(tV)
    
    if(tV(j) >= 0.3)
       here=1; 
    end
    
    muscleState     = zeros(benchConfig.numberOfMuscleStates,1);
    dmuscleState     = zeros(benchConfig.numberOfMuscleStates,1);
    
     if(benchConfig.numberOfMuscleStates ~= 0)
        if( j > size(muscleStateV,1))
            here=1;
        else
            muscleState = muscleStateV(j,:)'; 
        end
     end
    
    
    dactivation = activationFcn(excitationFcn(tV(j)), activation(j));
    activationState = [dactivation,activation(j)];
    pathState       = pathFcn(tV(j));

    dlp = pathState(1);
    lp  = pathState(2);
    
    mtInfo = calcMuscleInfoFcn(activationState,...                                   
                               pathState,...
                               muscleState);
                             

    if(isempty(mtInfo.extra)==0)
        numberOfExtra = length(mtInfo.extra);
        benchRecord.extra(j,idxSim,1:numberOfExtra)=mtInfo.extra';
        benchRecord.extraLabels = mtInfo.extraLabels;
    end

                             
    benchRecord.state(j,idxSim,:)=mtInfo.state.value';
    benchRecord.dstate(j,idxSim,:)=mtInfo.state.derivative';

    tendonForce = mtInfo.muscleDynamicsInfo.tendonForce;
           
    lceN      = mtInfo.muscleLengthInfo.normFiberLength;
    alpha     = mtInfo.muscleLengthInfo.pennationAngle;
    dlceN     = mtInfo.fiberVelocityInfo.normFiberVelocity;
    dalpha    = mtInfo.fiberVelocityInfo.pennationAngularVelocity;
            
    falN      = mtInfo.muscleLengthInfo.fiberActiveForceLengthMultiplier;
    fpeN      = mtInfo.muscleDynamicsInfo.normPassiveFiberForce;
    
    faeN      = mtInfo.muscleDynamicsInfo.normActiveFiberForce;
    fceNAT    = mtInfo.muscleDynamicsInfo.normFiberForce*cos(alpha);
    fceN      = mtInfo.muscleDynamicsInfo.normFiberForce;
    fceAT     = mtInfo.muscleDynamicsInfo.fiberForce*cos(alpha);
    fce       = mtInfo.muscleDynamicsInfo.fiberForce;    
    %kFiberAT  = mtInfo.muscleDynamicsInfo.fiberStiffnessAlongTendon;
    %dFiberAT  = mtInfo.muscleDynamicsInfo.fiberDampingAlongTendon;
    kx    = mtInfo.muscleDynamicsInfo.crossBridgeStiffness;
    dx    = mtInfo.muscleDynamicsInfo.crossBridgeDamping;    
    kTendon   = mtInfo.muscleDynamicsInfo.tendonStiffness;
    dTendon   = mtInfo.muscleDynamicsInfo.tendonDamping;
    
    km = mtInfo.muscleDynamicsInfo.musculotendonStiffness;
    dm = mtInfo.muscleDynamicsInfo.musculotendonDamping;
    kf = mtInfo.muscleDynamicsInfo.fiberStiffness;
    df = mtInfo.muscleDynamicsInfo.fiberDamping;
    
    %kMuscle   = mtInfo.muscleDynamicsInfo.normMuscleStiffness;

    %%
    %force and kinematic information
    %%
    benchRecord.activation(j,idxSim)                  = activation(j);
    benchRecord.activationDot(j,idxSim)               = dactivation;
    benchRecord.tendonForce(j,idxSim)                 = tendonForce;
    benchRecord.normFiberForceAlongTendon(j,idxSim)   = fceNAT;
    benchRecord.normFiberForce(j,idxSim)              = fceN;
    benchRecord.normActiveFiberForce(j,idxSim)        = faeN;
    benchRecord.fiberForceAlongTendon(j,idxSim)       = fceAT;
    benchRecord.fiberForce(j,idxSim)                  = fce;    
    benchRecord.tendonForce(j,idxSim)                 = tendonForce;
    benchRecord.normFiberLength(j,idxSim)             = lceN;    
    benchRecord.pennationAngle(j,idxSim)              = alpha;
    benchRecord.normFiberVelocity(j,idxSim)           = dlceN;
    benchRecord.pennationAngularVelocity(j,idxSim)    = dalpha;
    
    %benchRecord.fiberStiffnessAlongTendon(j,idxSim)   = kFiberAT;
    benchRecord.crossBridgeStiffness(j,idxSim)        = kx;    
    %benchRecord.fiberDampingAlongTendon(j,idxSim)     = dFiberAT;
    benchRecord.crossBridgeDamping(j,idxSim)          = dx;   
    
    benchRecord.tendonStiffness(j,idxSim)             = kTendon;
    benchRecord.tendonDamping(j,idxSim)               = dTendon;
    
    benchRecord.pathLength(j,idxSim)                  = lp;
    benchRecord.pathVelocity(j,idxSim)                = dlp;
    
    benchRecord.fiberActiveForceLengthMultiplier(j,idxSim) = falN;
    benchRecord.fiberPassiveForceLengthMultiplier(j,idxSim) = fpeN;
    benchRecord.normPassiveFiberForce(j,idxSim) = fpeN;
    
    %if(benchConfig.numberOfMuscleStates >= 4)
    %  benchRecord.crossBridgeForce(j,idxSim)       = ...
    %       mtInfo.muscleDynamicsInfo.crossBridgeForce;
    %  benchRecord.slowForceVelocityForce(j,idxSim)     = ...
    %       mtInfo.muscleDynamicsInfo.slowForceVelocityForce;
    %  benchRecord.fastForceVelocityForce(j,idxSim)     = ...
    %       mtInfo.muscleDynamicsInfo.fastForceVelocityForce; 
    %end
    
    %%
    %Energetics & power information
    %%
    
    fiberVelocity    = mtInfo.fiberVelocityInfo.fiberVelocity;
    activeFiberForce = mtInfo.muscleDynamicsInfo.activeFiberForce;
    tendonForce      = mtInfo.muscleDynamicsInfo.tendonForce;
    dampingForce     = mtInfo.muscleDynamicsInfo.dampingForces;

    fae = mtInfo.muscleDynamicsInfo.activeFiberForce;
    benchRecord.activeFiberForce(j,idxSim)             = fae;
    benchRecord.activeFiberForceAlongTendon(j,idxSim)  = fae*cos(alpha);
    
    benchRecord.musculotendonStiffness(j,idxSim)  = km;
    benchRecord.musculotendonDamping(j,idxSim)    = dm;
    
    benchRecord.fiberStiffness(j,idxSim)  = kf;
    benchRecord.fiberDamping(j,idxSim)    = df;
    
    benchRecord.tendonPower(j,idxSim)               = mtInfo.muscleDynamicsInfo.tendonPower;
    benchRecord.fiberParallelElementPower(j,idxSim) = mtInfo.muscleDynamicsInfo.fiberParallelElementPower;
    benchRecord.fiberActivePower(j,idxSim)          = mtInfo.muscleDynamicsInfo.fiberActivePower;
    benchRecord.dampingPower(j,idxSim)              = mtInfo.muscleDynamicsInfo.dampingPower;
    benchRecord.boundaryPower(j,idxSim)             = mtInfo.muscleDynamicsInfo.boundaryPower;

    boundaryPower       = mtInfo.muscleDynamicsInfo.boundaryPower;
    activeFiberPower    = mtInfo.muscleDynamicsInfo.fiberActivePower;
    dampingPower        = mtInfo.muscleDynamicsInfo.dampingPower;

    tendonPower               = mtInfo.muscleDynamicsInfo.tendonPower;
    fiberParallelElementPower = mtInfo.muscleDynamicsInfo.fiberParallelElementPower;
    
    %dTpVmW = - tendonPower...
    %         - fiberParallelElementPower ...
    %         - fiberBristleSpringPower ...
    %       - activeFiberPower ...
    %       - dampingPower ...
    %       - boundaryPower;
    
    TpVmW  = NaN;
    dTpVmW = NaN;
    if(flag_appendEnergetics == 1)
      assert(0,'Update!');
    end 
    
%     if(j==1)
%         T0V0 = mtInfo.musclePotentialEnergyInfo.fiberPotentialEnergy ...
%             +  mtInfo.musclePotentialEnergyInfo.tendonPotentialEnergy ...
%             + fiberBristleSpringPotentialEnergy;
%         
%         
%     end

    %TpVmW = energeticsT(j,idxSim)+energeticsV(j,idxSim)+energeticsWNeg(j,idxSim);
    
%     TpVmW   =  mtInfo.musclePotentialEnergyInfo.fiberPotentialEnergy ...
%              + mtInfo.musclePotentialEnergyInfo.tendonPotentialEnergy ...
%              + fiberBristleSpringPotentialEnergy ...
%              - workOfActiveFiber(j) ...
%              - workOfDamping(j) ...
%              - workOfBoundary(j) ...
%              - T0V0;
            
    benchRecord.systemEnergyLessWork(j,idxSim)        = TpVmW;
    benchRecord.dSystemEnergyLessWork(j,idxSim)       = dTpVmW;
    benchRecord.tendonPotentialEnergy(j,idxSim)       = ...
        mtInfo.musclePotentialEnergyInfo.tendonPotentialEnergy;
    benchRecord.fiberPotentialEnergy(j,idxSim)        = ...
        mtInfo.musclePotentialEnergyInfo.fiberPotentialEnergy;
    benchRecord.fiberActiveWork(j,idxSim)             = NaN;%workOfActiveFiber(j);
    benchRecord.dampingWork(j,idxSim)                 = NaN;%workOfDamping(j);
    benchRecord.boundaryWork(j,idxSim)                = NaN;%workOfBoundary(j);
    
    
    
end       
benchRecord.time(:,idxSim) = tV;

if(benchConfig.numberOfMuscleStates == 0)
    cpuTime = toc(t0);
end

benchRecord.cpuTime(idxSim) = cpuTime;


                                         
                                     