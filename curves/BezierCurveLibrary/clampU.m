function uC = clampU(u)
%%
% This is used to ensure that the argument for the Bezier curve's Bernstein
% polynomial is limited to [0,1] and not a bit outside
%
% @param u
% @returns uC, version of u that is guaranteed to be within [0,1]
%%
uC = u;
if(u<0.0)
    uC=0;
end
if(u>1.0)
    uC=1;
end

