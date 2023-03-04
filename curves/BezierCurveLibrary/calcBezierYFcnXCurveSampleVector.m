function curveValues = calcBezierYFcnXCurveSampleVector(curveParams, npts, domain)
%%
% This function evaluates a Bezier spline curve (x(u), y(u)) across its
% entire domain, and and slightly beyond to get values in the extrapolated
% region. The spline is evaluated at its value, and first 3 derivatives. In
% addition if the curve has its integral function defined, then the integral
% is also evaluated.
%
% @param curveParams : a structure with the fields
%           .xpts
%           .ypts
%           .integral
%
%   xpts: n x m matrix of Bezier control points defining x(t). 
%
%              n x m matrix 
%              n: Bezier curve order + 1 
%              m: # spline sections
%                
%              Each column defines the control points used for a Bezier 
%              section. Thus a quintic Bezier spline with 3 sections will 
%              have a 6 x 3 matrix for xpts
%              
%   ypts: n x m matrix of Bezier control points defining y(t)
%
%   integral: if the integral curve has not been numerically computed, 
%             then this field will be empty. Otherwise it will have fields
%             of
%
%                      .xptsN : points the integral is evaluated at  
%                      .yptsN : numerical value of the integral
%                      .y1ptsN: '' integral's first derivative.
%                      .y2ptsN: '' integral's 2nd derivative.
%
% @param npts: the number of samples to use across the curve domain
%              
% @return curveValues, a struct with the fields
%
%	.x      : vector of x values of the Bezier spline
%	.y      : vector of y values of the Bezier spline
%	.dydx   : first derivative        dy/dx
%	.d2ydx2 : second derivative       d^2y/dx^2
%	.d3ydx3 : third derivative        d^3y/dx^3
%	.intYdx : integral of y w.r.t. x  dy/dx
%%
xmin = 0;
xmax = 0;

if(isempty(domain)==0)
   xmin = domain(1);
   xmax = domain(2);
else
    xmin  = min(min(curveParams.xpts));
    xmax  = max(max(curveParams.xpts));
    delta = 0;
    if(    sum(isnan(curveParams.dydxEnd))   == 0 ...
        && sum(isinf(curveParams.dydxEnd))   == 0 ...
        && sum(isnan(curveParams.d2ydx2End)) == 0 ...
        && sum(isinf(curveParams.d2ydx2End)) == 0)    
        delta = xmax-xmin;
    end

    xmin  = xmin-delta/10;
    xmax  = xmax+delta/10;
end

ymin  = min(min(curveParams.ypts));
ymax  = max(max(curveParams.ypts));
delta = 0;

delta = 0;
if(    sum(isnan(curveParams.dydxEnd))   == 0 ...
    && sum(isinf(curveParams.dydxEnd))   == 0 ...
    && sum(isnan(curveParams.d2ydx2End)) == 0 ...
    && sum(isinf(curveParams.d2ydx2End)) == 0)

    delta = ymax-ymin;
end

ymin  = ymin-delta/10;
ymax  = ymax+delta/10;

x = [xmin:((xmax-xmin)/npts):xmax]';  
    
y      = zeros(size(x));
dydx   = zeros(size(x));
d2ydx2 = zeros(size(x));
d3ydx3 = zeros(size(x));
intYdx = [];


if(isempty(curveParams.integral) == 0)
   intYdx = zeros(size(x)); 
end

intYdx = zeros(size(x));

for k=1:1:length(x)
   %if(isempty(curveParams.integral) == 0)
   %  intYdx(k) = calcBezierYFcnXDerivative(x(k), curveParams, -1); 
   %end 
    
   y(k)      = calcBezierYFcnXDerivative(x(k), curveParams, 0);
   dydx(k)   = calcBezierYFcnXDerivative(x(k), curveParams, 1);
   d2ydx2(k) = calcBezierYFcnXDerivative(x(k), curveParams, 2);
   d3ydx3(k) = calcBezierYFcnXDerivative(x(k), curveParams, 3);
end

curveValues.x      = x;
curveValues.y      = y;
curveValues.dydx   = dydx;
curveValues.d2ydx2 = d2ydx2;
curveValues.d3ydx3 = d3ydx3;
curveValues.intYdx = intYdx;
