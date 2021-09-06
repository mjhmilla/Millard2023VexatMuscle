function fiberState = clampFiberStateAlongTendon(lceAT,dlceAT,lceATMin)
%%
%
% @param lceAT: fiber length along the tendon
% @param dlceAT: d/dt lceAt - fiber velocity along the tendon
% @param lceATMin: minimum fiber length permitted along the tendon
%
% @returns fiberState, a struct with fields of
%
%       .lceAT : length of the fiber along the tendon
%       .dlceAT: fiber velocity along the tendon
%       .isClamped: 0 if the fiber value has not been clamped
%                   1 if the fiber velocity has been clamped
%%
isClamped = 0;
if(lceAT <= lceATMin || (lceAT == lceATMin && dlceAT < 0))            
    %Clamp the fiber length along the tendon
    lceAT  = lceATMin;
    dlceAT = 0;    
    isClamped = 1;
end

fiberState.lceAT     = lceAT;
fiberState.dlceAT    = dlceAT;
fiberState.isClamped = isClamped; 