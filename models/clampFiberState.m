function fiberState = clampFiberState(lce,dlce,lceMin)
%%
% This function returns the state of the fiber, and if necessary, 
% clamps it to a lower bound defined by lceMin
%
% @param lce: length of the fiber (m)
% @param dlce: fiber lengthening velocity (m/s)
% @param lceMin: the minimum allowable length of the fiber
%
% @return fiberState a structure with fields of
%           .lce : fiber length
%           .dlce: fiber velocity
%           .isClamped: 0 - not clamped
%                       1 - clamped
%%
isClamped = 0;
if(lce < lceMin || (lce == lceMin && dlce < 0))            
    %Clamp the fiber length along the tendon
    lce  = lceMin;
    dlce = 0;    
    isClamped = 1;
end

fiberState.lce     = lce;
fiberState.dlce    = dlce;
fiberState.isClamped = isClamped; 