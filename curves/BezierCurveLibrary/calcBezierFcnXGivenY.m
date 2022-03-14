function x = calcBezierFcnXGivenY(y, curve, x0)

iterMax=100;
iter=1;
err = 1.0;

cols = size(curve.ypts,2);
x = NaN;
if(exist('x0','var')==0)
  x = curve.xpts(1,1);
  intervalWidth = -1;
  if(y <= curve.ypts(1,1))
      x = curve.xpts(1,1);
  elseif( y >= curve.ypts(1,cols) )
      x = curve.xpts(1,cols);
  else
      for(i=1:1:(size(curve.ypts,2)-1))
          if(y > curve.ypts(1,i) && y <= curve.ypts(1,i+1))
             x = 0.5*(curve.xpts(1,i)+curve.xpts(6,i)); 
             intervalWidth = curve.xpts(6,i)-curve.xpts(1,i);
          end
      end    
  end
else
  x = x0;
  intervalWidth = 0;
  if(size(curve.xpts,2) > 1)
    for(i=1:1:(size(curve.xpts,2)-1))
        if(x > curve.xpts(1,i) && x <= curve.xpts(1,i+1))
           intervalWidth = curve.xpts(6,i)-curve.xpts(1,i);
        end
    end     
  else
    intervalWidth = curve.xpts(6,1)-curve.xpts(1,1);      
  end
end

err  = calcBezierYFcnXDerivative( x, curve, 0) ...
            -y;

tol = 1e-9;        
        
while(abs(err) > tol && iter < iterMax)
    err  = calcBezierYFcnXDerivative( x, curve, 0) ...
            -y;
    derr = calcBezierYFcnXDerivative( x, curve, 1);
    if(abs(err) > tol && abs(derr) > eps)
        delta = -err/derr;

        if(intervalWidth > 0)
           if(abs(delta) > 0.1*intervalWidth)
              delta = 0.1*intervalWidth*sign(delta); 
           end
        end
        x = x+delta;
    end
    
    iter=iter+1;
end
if(abs(err) > tol)
   here=1; 
end
assert( abs(err) <= tol, 'Error: failed to converge!');