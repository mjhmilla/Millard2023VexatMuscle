function quadraticCurve = convertToQuadraticBezierCurve(higherOrderCurve,subDivisions)

assert(size(higherOrderCurve.xpts,1) > 3)

rootEPS = eps^0.5;

if isempty(subDivisions)
    subDivisions=1;
end

quadraticCurve = struct(...
    'xpts',zeros(3,size(higherOrderCurve.xpts,2)*subDivisions),...
    'ypts',zeros(3,size(higherOrderCurve.xpts,2)*subDivisions),...
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
        x1 = x0+xSubdivision;

        %x0 = higherOrderCurve.xpts(1,i);
        %x1 = higherOrderCurve.xpts(end,i);
        
        y0 = calcBezierYFcnXDerivative(x0,higherOrderCurve,0);
        y1 = calcBezierYFcnXDerivative(x1,higherOrderCurve,0);
    
        dydx0 = calcBezierYFcnXDerivative(x0,higherOrderCurve,1);    
        dydx1 = calcBezierYFcnXDerivative(x1,higherOrderCurve,1);
    
        p0 = [x0,y0];
        p2 = [x1,y1];
    
        %1. Calculate the location where the two lines intersect
        % (x-x0)*dydx0 + y0 = (x-x1)*dydx1 + y1
        %   x*(dydx0-dydx1) = y1-y0-x1*dydx1+x0*dydx0
        %                 x = (y1-y0-x1*dydx1+x0*dydx0)/(dydx0-dydx1);
        
        xC = 0;
        yC1 = 0;
        yC2 = 0;
        
        if(abs(dydx0-dydx1) > rootEPS)
            xC = (y1-y0-x1*dydx1+x0*dydx0)/(dydx0-dydx1);    
            yC2 = (xC-x1)*dydx1 + y1;
            yC1=  yC2;
        else
            xC = (x1+x0)/2;
            yC1 = (xC-x0)*dydx0 + y0;
            yC2 = (xC-x1)*dydx1 + y1;
        end  
        if(abs(yC1-yC2) >= rootEPS)
            here=1;
        end
        assert(abs(yC1-yC2)<rootEPS);
    
    
        quadraticCurve.xpts(:,idx) = [x0; xC;x1];
        quadraticCurve.ypts(:,idx) = [y0;yC1;y1];
    
        if(idx==1)
            quadraticCurve.dydxEnd(1,1)=dydx0;
        end
        if(idx==size(quadraticCurve.xpts,2))
            quadraticCurve.dydxEnd(1,2)=dydx1;
        end
        idx=idx+1;
    end
end

here=1;