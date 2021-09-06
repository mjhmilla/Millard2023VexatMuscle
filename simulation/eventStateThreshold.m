function [position,isterminal,direction] = ...
  eventStateThreshold(t,y, idxState, stateValue, stateDirection)

position   = -1;
isterminal = 0;
direction  = 1;

if(idxState <= 0)
  position   = t-stateValue;
  direction  = stateDirection;
  isterminal = 0;
else
  position   = y(idxState)-stateValue;
  isterminal = 0;
  direction  = stateDirection;    
end

if position*stateDirection >= 0
  here=1;
end
