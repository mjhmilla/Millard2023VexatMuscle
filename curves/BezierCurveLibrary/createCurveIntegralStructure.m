function integralStruct = createCurveIntegralStructure(curveParams,...
                                                       npts,... 
                                                       tol, ...
                                                       xScaling,...
                                                       flag_usingOctave ,...
                                                       flag_integrateLeftToRight)
%%
% Numerically evaluates the integral at npts over the domain of the 2D 
% Bezier curve defined by the matrix of x control points and y control.
% The first and second derivative of the curve integral are also evaluated
% at these points so that it is possible to interpolate the curve using a
% quintic hermine spline.
%
% @param x: value to 
% @param xpts: n x m matrix of Bezier control points defining x(t). 
%
%              n x m matrix 
%              n: Bezier curve order + 1 
%              m: # spline sections
%                
%              Each column defines the control points used for a Bezier 
%              section. Thus a quintic Bezier spline with 3 sections will 
%              have a 6 x 3 matrix for xpts
%              
% @param ypts: n x m matrix of Bezier control points defining y(t)
%
% @param npts: number of intermediate points between xmin and xmax to
%              evaluate the integeral.
% @param tol: relative and absolute tolerance on the integral.
%
% @return integralStruct: A structure containing the numerically calculated
%                         integral, its first and second derivative so that 
%                         the integral curve can be interpolated using a 
%                         quintic Hermite spline.
%
%           Has fields of
%                      .xptsN : points the integral is evaluated at  
%                      .yptsN : numerical value of the integral
%                      .y1ptsN: '' integral's first derivative.
%                      .y2ptsN: '' integral's 2nd derivative.
%
%%

xmin = min(min(curveParams.xpts));
xmax = max(max(curveParams.xpts));
xv = [];
if(flag_integrateLeftToRight==1)
  xv   = [xmin:((xmax-xmin)/(npts-1)):xmax];
else
  xv   = [xmax:(-(xmax-xmin)/(npts-1)):xmin];    
end    

xe = [];
ye1 = [];
ye2 = [];

if(flag_usingOctave == 0)
  fcn = @(arg,arg1)calcBezierYFcnXDerivative(arg, curveParams, 0);

  options = odeset('RelTol',tol,'AbsTol',tol);
  [xe ye] = ode45(fcn,xv,0,options);

  ye1 = zeros(size(ye));
  ye2 = zeros(size(ye));

else
  fcn = @(arg1,arg)calcBezierYFcnXDerivative(arg, curveParams, 0);

  lsode_options('relative tolerance', tol);
  lsode_options('absolute tolerance', tol);
   
  t0 = clock();          
  ye = lsode(fcn,xv,0);                                                   
  xe = xv;
  
end

for i=1:1:length(ye)
   ye1(i) =  calcBezierYFcnXDerivative(xe(i), curveParams, 0);
   ye2(i) =  calcBezierYFcnXDerivative(xe(i), curveParams, 1);
end


integralStruct.xptsN  = xe;
integralStruct.yptsN  = ye;
integralStruct.y1ptsN = ye1;
integralStruct.y2ptsN = ye2;
integralStruct.xScaling = xScaling;

