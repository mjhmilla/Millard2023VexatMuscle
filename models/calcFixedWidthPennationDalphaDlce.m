function Dalpha_Dlce  = ...
    calcFixedWidthPennationDalphaDlce( alpha,...
                                       fiberLength,...
                                       optimalFiberLength,...
                                       pennationAngleAtOptimalFiberLength,...
                                       minimumFiberLengthAlongTendon)
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

if(pennationAngleAtOptimalFiberLength > eps^0.5)
    lce       = fiberLength;   
    lceOpt    = optimalFiberLength;         
    alphaOpt  = pennationAngleAtOptimalFiberLength;     
    cosAlpha  = cos(alpha);
    lceAT     = lce*cosAlpha;
    assert(lceAT > minimumFiberLengthAlongTendon,...
           ['Impending singularity: lceAT < lceATMin ']);
        
    h         = lceOpt*sin(alphaOpt);

    %x         = h/lceAT;
    %dxdlce    = -h*cosAlpha/(lceAT*lceAT);
    %Dalpha_Dlce_dep = (1/(1 + x*x))*dxdlce;

    y = h/lce;
    dydlce = -h/(lce*lce);

    Dalpha_Dlce =  dydlce/sqrt(1-y*y) ;
    %assert(0,'Check this');

end

                                                     