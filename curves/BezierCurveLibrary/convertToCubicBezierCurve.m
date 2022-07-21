function cubicCurve = convertToCubicBezierCurve(higherOrderCurve,...
    subDivisions,cubicZeroSecondDerivative,verbose)

assert(size(higherOrderCurve.xpts,1) > 3)

rootEPS = eps^0.5;

if isempty(subDivisions)
    subDivisions=1;
end

cubicCurve = struct(...
    'xpts',zeros(4,size(higherOrderCurve.xpts,2)*subDivisions),...
    'ypts',zeros(4,size(higherOrderCurve.xpts,2)*subDivisions),...
    'xEnd',higherOrderCurve.xEnd,...
    'yEnd',higherOrderCurve.yEnd,...
    'dydxEnd',higherOrderCurve.dydxEnd,...
    'd2ydx2End',[0,0],...
    'integral',[]);


idx=1;
for i=1:1:size(higherOrderCurve.xpts,2)

    xDivision = higherOrderCurve.xpts(end,i)-higherOrderCurve.xpts(1,i);
    xSubdivision = xDivision/subDivisions;

    for j=1:1:subDivisions
        x0 = higherOrderCurve.xpts(1,i) + (j-1)*xSubdivision;
        x3 = x0+xSubdivision;

        %x0 = higherOrderCurve.xpts(1,i);
        %x1 = higherOrderCurve.xpts(end,i);
        
        y0 = calcBezierYFcnXDerivative(x0,higherOrderCurve,0);
        y3 = calcBezierYFcnXDerivative(x3,higherOrderCurve,0);
    
        dydx0 = calcBezierYFcnXDerivative(x0,higherOrderCurve,1);    
        dydx3 = calcBezierYFcnXDerivative(x3,higherOrderCurve,1);

        d2ydx20 = calcBezierYFcnXDerivative(x0,higherOrderCurve,2);    
        d2ydx23 = calcBezierYFcnXDerivative(x3,higherOrderCurve,2);

    
        p0 = [x0,y0];
        p3 = [x3,y3];
    
        %1. Calculate the location where the two lines intersect
        % (x-x0)*dydx0 + y0 = (x-x1)*dydx1 + y1
        %   x*(dydx0-dydx1) = y1-y0-x1*dydx1+x0*dydx0
        %                 x = (y1-y0-x1*dydx1+x0*dydx0)/(dydx0-dydx1);
        
        xC = 0;
        yC1 = 0;
        yC2 = 0;
        
        if(abs(dydx0-dydx3) > rootEPS)
            xC = (y3-y0-x3*dydx3+x0*dydx0)/(dydx0-dydx3);    
            yC2 = (xC-x3)*dydx3 + y3;
            yC1=  yC2;
        else
            xC = (x3+x0)/2;
            yC1 = (xC-x0)*dydx0 + y0;
            yC2 = (xC-x3)*dydx3 + y3;
        end  
        if(abs(yC1-yC2) >= rootEPS)
            here=1;
        end
        assert(abs(yC1-yC2)<rootEPS);
        yC = yC1;

        if (abs(d2ydx20) < rootEPS && abs(d2ydx23) < rootEPS) ...
                || (cubicZeroSecondDerivative==1)
            % z(u)      = [p0, p1, p2, p3]
            % z(0)      = p0; 
            % z(1)      = p3
            %
            % dzdu      = 3* diff(p)
            %           = 3 [(p1-p0), (p2-p1), (p3-p2)]
            % dzdu(0)   = 3(p1-p0)
            % dzdu(1)   = 3(p3-p2)            
            % 
            % f = dydu / dxdu
            %
            % d2zdu2    = 3*2*diff(diff(p))
            % d2zdu2    = 6[(p2-2p1+p0), (p3-2p2+p1)]
            % d2zdu2(0) = 6(p2-2p1+p0)
            % d2zdu2(1) = 6(p3-2p2+p1)
            %
            % If p2 = p1 then
            % d2zdu2(0) = 6(-p1+p0) = -2*dzdu(0)
            %               p1 = p0+(1/3)*dzdu(0)
            % d2zdu2(1) = 6(p3-p2)  =  2*dzdu(1)
            %               p2 = 6p3-2*dzdu(1) 
            %
            %    g = ( (d2ydu2) - (d2xdu2)*dydu/dxdu )(1/(dxdu^2)) 
            %      = ( (d2ydu2) - (d2xdu2)*f )(1/(dxdu^2)) 
            %
            % [1] g(0) = ( 6(y2-2y1+y0) - 6(x2-2x1+x0)*f(0) )( 1/(3*(x1-x0))^2 ) 
            %
            % If x1=x2 and y1=y2
            %
            % g(0) = ( 6(-y1+y0) - 6(-x1+x0)*f(0) )( 1/(3*(x1-x0))^2 )   
            % g(0) = ( 6(-y1+y0) + 6(x1-x0)*f(0) )( 1/(3*(x1-x0))^2 )   
            % 
            % -2dydu(0) = 6(-y1+y0)
            %  2dxdu(0) = 6(x1-x0)
            %  2dxdu(0)*f(0) = 2*dydu(0)
            %
            % g(0) = (-2dydu + 2dydu )( 1/(3*(x1-x0))^2 )
            %      = 0
            %
            % Similarily
            %
            % [2] g(1) = ( 6(y3-2y2+y1) - 6(x3-2x2+x1)*f(1) )( 1/(3*(x1-x0))^2 )
            %      = ( 6(y3-y2) - 6(x3-x2)*f(1) )( 1/(3*(x1-x0))^2 )
            %
            %  2dydu(1) = 6(y3-y2)
            %  2dxdu(1) = 6(x3-x2)
            %  2dxdu*f(1) = 2*dydu(1)
            % g(1) = ( 2*dydu - 2*dydu )( 1/(3*(x1-x0))^2 )
            % g(1) = 0
            %
            % In the more general case we can substitute
            %  x1 = d01(xC-x0) + x0
            %  y1 = d01(yC-y0) + y0
            %
            %  x2 = d32(xC-x3) + x3
            %  y2 = d32(yC-y3) + y3            
            %
            % into Eqns 1 and 2. Then we end up with a system of two 
            % quadratic equations in two unknowns. Since the solution
            % may not exist (resulting in imaginary roots), and the 
            % solution is a bit ugly in general, below we use the bisection
            % method to solve for the values of d01 and d32 to get the 
            % best valid solution possible for these two factors.
            %
            % After looking in more detail to how to analytically solve
            % two coupled quadratic equations I've come to the tentative
            % conclusion that this is, in general, very difficult. I'm 
            % going to stick with the numerical solution below. This can
            % be speeded up by using something that converges faster than
            % bisection, but for now, this does not need to be fast.

            x1=xC;
            y1=yC;
            x2=xC;
            y2=yC;
            d2ydx20=0;
            d2ydx23=0;
        else 

            %Numerically solve for the normalized distance d01 that is between
            %p0 and pC, and d32 the distance between point p3 and p2
            d01 = 0.5;
            x1 = x0 + d01*(xC-x0);
            y1 = y0 + d01*(yC-y0);
            
            d32 = 0.5;
            x2 = x3 + d32*(xC-x3);
            y2 = y3 + d32*(yC-y3);
    
            tempCubicCurve.xpts = [x0;x1;x2;x3];
            tempCubicCurve.ypts = [y0;y1;y2;y3];
            tempCubicCurve.xEnd = [x0,x3];
            tempCubicCurve.yEnd = [y0,y3];
            tempCubicCurve.dydxEnd = [dydx0,dydx3];
            tempCubicCurve.yEnd = [d2ydx20,d2ydx23];
    
            d2ydx20ErrorBest = ...
                abs(calcBezierYFcnXDerivative(x0+eps,tempCubicCurve,2) ...
                    - d2ydx20) ;
            d2ydx23ErrorBest = ...
                abs(calcBezierYFcnXDerivative(x3-eps,tempCubicCurve,2) ...
                    - d2ydx23);
    
            errBest = d2ydx20ErrorBest+d2ydx23ErrorBest;
    
            xptsBest = [x0;x1;x2;x3];
            yptsBest = [y0;y1;y2;y3];
    
            numBisections=26; 
            %This will take the final solution to a tolerance of sqrt(eps)
            %or 1.49e-8 for a double precision float.
            
    
            h=0.5;
            for indexBisection=1:1:numBisections
                d01Upd=d01;
                d32Upd=d32;
                errBestUpd=errBest;            
                for stepD01 = -1:2:1
                    d01test = d01+stepD01*h;
                    for stepD32 = -1:2:1
                        d32test = d32+stepD32*h;
                            
                        %Evaluate this set of candidate points
                        x1 = x0 + d01test*(xC-x0);
                        y1 = y0 + d01test*(yC-y0);
                        
                        x2 = x3 + d32test*(xC-x3);
                        y2 = y3 + d32test*(yC-y3);
                
                        tempCubicCurve.xpts = [x0;x1;x2;x3];
                        tempCubicCurve.ypts = [y0;y1;y2;y3];
                        tempCubicCurve.xEnd = [x0,x3];
                        tempCubicCurve.yEnd = [y0,y3];
                        tempCubicCurve.dydxEnd = [nan,nan];
                        tempCubicCurve.d2ydx2End = [nan,nan];
                
                        d2ydx20test = ...
                            abs(calcBezierYFcnXDerivative(x0+eps,tempCubicCurve,2) ...
                                -d2ydx20);
                        d2ydx23test = ...
                            abs(calcBezierYFcnXDerivative(x3-eps,tempCubicCurve,2) ...
                                -d2ydx23);
    
                        errTest = d2ydx20test+d2ydx23test;
    
                        if( errTest < errBestUpd )
                            d01Upd=d01test;
                            d32Upd=d32test;
                            errBestUpd=errTest;
                        end
    
                    end
                end
                if(errBestUpd < errBest)
                    d01=d01Upd;
                    d32=d32Upd;
                    errBest=errBestUpd;
                end
                h=h/2;
            end

            if(abs(errBest)> rootEPS && verbose==1)
                fprintf('  %1.3e\t Could not satisfy 2nd derivative\n',...
                    abs(errBest));
            end
            
            x1 = x0 + d01*(xC-x0);
            y1 = y0 + d01*(yC-y0);
            
            x2 = x3 + d32*(xC-x3);
            y2 = y3 + d32*(yC-y3);
    
            tempCubicCurve.xpts = [x0;x1;x2;x3];
            tempCubicCurve.ypts = [y0;y1;y2;y3];
            tempCubicCurve.xEnd = [x0,x3];
            tempCubicCurve.yEnd = [y0,y3];
            tempCubicCurve.dydxEnd = [nan,nan];
            tempCubicCurve.yEnd = [nan,nan];    

            d2ydx20=calcBezierYFcnXDerivative(x0+eps,tempCubicCurve,2);
            d2ydx23 = calcBezierYFcnXDerivative(x3-eps,tempCubicCurve,2);
        end        

        cubicCurve.xpts(:,idx) = [x0;x1;x2;x3];
        cubicCurve.ypts(:,idx) = [y0;y1;y2;y3];
    
        if(idx==1)
            cubicCurve.dydxEnd(1,1)=dydx0;
            cubicCurve.d2ydx2End(1,1)=d2ydx20;
        end
        if(idx==size(cubicCurve.xpts,2))
            cubicCurve.dydxEnd(1,2)=dydx3;            
            cubicCurve.d2ydx2End(1,1)=d2ydx23;            
        end
        idx=idx+1;
    end
end

here=1;