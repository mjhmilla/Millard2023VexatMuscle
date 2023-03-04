function partialDer  = ...
    calcFixedWidthPennationPartialDerivatives( alpha,...
                                       fiberLength,...
                                       optimalFiberLength,...
                                       pennationAngleAtOptimalFiberLength)
%%
%This function evaluates the partial derivative of the pennation angle 
%w.r.t. a change in fiber length.
%
% @param alpha: pennation angle of the fiber (radians)
% @param fiberLength (m)
% @param optimalFiberLength (m)
% @param pennationAngleAtOptimalFiberLength (radians)
%
% @return Dalpha_Dlce the partial derivative of the pennation angle w.r.t.
%                     the fiber length
%%
                                                   
Dalpha_Dlce = 0;
Ddalpha_Ddlce = 0;

if(pennationAngleAtOptimalFiberLength > eps^0.5)
    lce       = fiberLength;   
    lceOpt    = optimalFiberLength;         
    alphaOpt  = pennationAngleAtOptimalFiberLength;     
    cosAlpha  = cos(alpha);
    lceAT     = lce*cosAlpha;
    assert(lceAT > eps^0.5,...
           ['Impending singularity: lceAT < eps^0.5 ']);
        
    h         = lceOpt*sin(alphaOpt);

    %x         = h/lceAT;
    %dxdlce    = -h*cosAlpha/(lceAT*lceAT);

    %Incorrect: dxdlce should be
    %dxdlce    = -h*(cosAlpha - lce*sinAlpha*Dalpha_Dlce)/(lceAT*lceAT);    
    %Dalpha_Dlce_dep = (1/(1 + x*x))*dxdlce;

    Ddalpha_Ddlce = -(1/lce)*tan(alpha); 
    
    % sin(alpha) = h/lce
    %  x = h/lce
    %  dx/dlce = -h/(lce*lce)
    % alpha = asin(x)
    % dalpha/dx = 1/(1-x^2)
    % (dalpha/dx)(dx/dlce) = dalpha/dlce = (1/(1-(h/lce)^2)*(-h/(lce*lce))
    y = h/lce;
    dy = -h/(lce*lce);
    Dalpha_Dlce = dy/sqrt(1-y*y);

end

partialDer = struct(...
    'Dalpha_Dlce',  Dalpha_Dlce,...
    'Ddalpha_Ddlce',Ddalpha_Ddlce);

                                                     
