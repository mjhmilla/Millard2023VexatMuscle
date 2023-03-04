function dydx = calcCentralDifferenceDataSeries(x, y)
%%
%This function will compute a numerical central difference for dy/dx. The
%endpoints of the numerical derivative are linearly extrapolated to obtain
%a reasonable level of accuracy.
%
% @param x 
%   An n x 1 domain vector
% @param y
%   An n x 1 range vector
% @returns dydx
%   An n x 1 vector of dy/dx which has been numerically calculated
%
%@author M.Millard
%@date 2013/9/27
%%
dydx = [];

if(isempty(x)==1 || isempty(y) == 1)
    disp('calcCentralDifference: X or Y is empty');
    return;
end
if(length(x) ~= length(y) || size(x,1) ~= size(y,1))
   disp('calcCentralDifference: X and Y must have the same dimensions'); 
   return
end

len = length(x);
xdiff = diff(x);    
ydiff = diff(y);
dydx = zeros(len,1);

dydx(2:len-1) = (ydiff(1:len-2) + ydiff(2:len-1))...
              ./(xdiff(1:len-2) + xdiff(2:len-1));

%Linearily extrapolate the endpoints
ddydx2 = (dydx(3)-dydx(2))/(x(3)-x(2));
dx = x(2)-x(1);    
dydx(1) = dydx(2) - ddydx2*dx;


ddydx2 = (dydx(len-1)-dydx(len-2))/(x(len-1)-x(len-2));
dx = x(len)-x(len-1);    
dydx(len) = dydx(len-1) + ddydx2*dx;
