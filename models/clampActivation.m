function clampedActivation = clampActivation(activation, activationMin, activationMax)
%%
% This function clamps the activation to be within the domain of 
% activationMin and activationMax
%
% @param activation: the value of the activation of the muscle
% @param activationMin: the minimum allowed activation value
% @param activationMax: the maximum allowed activation value.
%
% @return clampedActivation
%%
clampedActivation = activation;

if activation < activationMin
    clampedActivation = activationMin;
end

if activation > activationMax
   clampedActivation = activationMax; 
end
