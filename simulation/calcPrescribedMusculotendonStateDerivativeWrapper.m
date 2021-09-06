function dState = ...
    calcPrescribedMusculotendonStateDerivativeWrapper(t,...
                                                      state,...                                                      
                                                      calcPrescribedPathFcn,...
                                                      calcPrescribedExcitationFcn,...
                                                      calcActivationDot,...
                                                      calcMuscleInfoFcn,...
                                                      flag_appendEnergetics)
%%
% This is a wrapper ultimately to create a derivative function that takes t
% and muscle state as arguments and returns the muscle state derivative,
% and the various powers of interest, as the muscle undergoes a constant
% activation sinusoidal stretch benchmark simulation
%
% @param t: time in seconds
%
% @param state: the state vector for this benchmark simulation. This vector
%               contains:
%
%               [muscleState;
%                boundaryPower;
%                activeFiberPower;
%                dampingPower]
% 
%                The last three entries are required to numerically
%                evaluate T + V - W = const to ensure that the model is
%                conservative.
%
% @param calcPrescribedPathFcn : a handle to a function that given time t
%                            produces a 2x1 vector containing the velocity 
%                            and length of the path the muscle lies on.
%
% @param calcActivationDot: a handle to a function that given time t
%                            produces a 1x1 scalar of the activation of the
%                            muscle.
%
% @param calcMuscleInfoFcn: a function handle to a muscle function
%                      like (calcMillard2012DampedEquilibriumMuscleInfo.m) 
%                      and takes arguments of activation, pathState, and 
%                      muscle state
%
% @returns dState: the first time derivative of the state vector
%
%               [muscleState;
%                boundaryPower;
%                activeFiberPower;
%                dampingPower]
%
%
%%
% 
% disp(t);
% 
% if( t <= 0.9)
%    here=1;
% end
% 
% if(t >= 1.5 && t < 1.6)
%    here=1;
% end


pathState       = calcPrescribedPathFcn(t);

stateLength     = length(state);
muscleState     = [];
if(flag_appendEnergetics==1)
  muscleState     = state(2:(stateLength-3));  
else
  muscleState     = state(2:end);
end


activation      = state(1);

if(activation < 0)
  here=1;
end

excitation       = calcPrescribedExcitationFcn(t);
[activationDot, activationClamped] = calcActivationDot(excitation,activation);
activationState = [activationDot;activationClamped]; 

if(excitation > 0.01)
   here=1; 
end

mtInfo = calcMuscleInfoFcn(activationState,...
                           pathState,...
                           muscleState);      

              
                       
dlp = pathState(1);                       


boundaryPower       = mtInfo.muscleDynamicsInfo.boundaryPower;
activeFiberPower    = mtInfo.muscleDynamicsInfo.fiberActivePower;
dampingPower        = mtInfo.muscleDynamicsInfo.dampingPower;



dState = [];

if(flag_appendEnergetics == 1)
    dT = mtInfo.muscleDynamicsInfo.dT;
    dV = mtInfo.muscleDynamicsInfo.dV;
    dW = mtInfo.muscleDynamicsInfo.dW;
  
    if(sum(isnan(mtInfo.state.derivative)) >0)
      dState = [ activationState(1);...
                 dT; dV; -dW];          
    else
      dState = [ activationState(1);...
                 mtInfo.state.derivative; ...
                 dT; dV; -dW];          
    end
    
           
else
    if(sum(isnan(mtInfo.state.derivative)) >0)
      dState = [ activationState(1)];          
    else
      dState = [ activationState(1);...
                 mtInfo.state.derivative];          
    end
end

