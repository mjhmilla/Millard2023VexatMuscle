function f = calc1DBezierCurveValue(u, pV)
%%
% This function implements De Casteljau's recursive algorithm for
% evaluating an nth order 1D Bezier curve of the form
%
%  B(u) = sum_{i=0}^n  [(n choose i)(1-u)^{n-1} u^i] p_i
%
% where 
%  u: argument of the curve
%  n: order of the curve
%  p_i: value of the ith control point
%
% For an n th order curve with (n+1) points this algorithm requires n! 
% subtractions, multiplications, and additions to terminate, and a stack 
% that is n! deep. Although this algorithm is very general faster results 
% can be obtained using optimized code if the order of the Bezier curve is 
% known ahead of time using optimized code.
%
% @params u : [0,1] the argument of the Bezier curve
% @params pV: vector of control points
% @returns f: the value of the Bezier curve evaluated at u
%
% Reference
% [1] http://en.wikipedia.org/wiki/De_Casteljau%27s_algorithm
%
%%
bV0 = zeros(size(pV,1)-1,1);


for i=1:1:length(bV0)
    bV0(i) = (pV(i+1)-pV(i))*(u) + pV(i);
end

if(length(bV0) == 1)    
    f  = bV0;
else
    f = calc1DBezierCurveValue(u,bV0);
end


