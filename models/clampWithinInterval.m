function [uC, isClampedToLB,isClampedToUB]= clampWithinInterval(u,umin,umax)
%%
% This function clamps u to be within the interval [0,1]
%
% @param u
% @return uC, version of u that is guaranteed to be within [0,1]
%%
uC = u;
isClampedToLB=0;
isClampedToUB=0;
if(u<=umin)
    uC=umin;
    isClampedToLB=1;
end
if(u>=umax)
    uC=umax;
    isClampedToUB=1;
end

