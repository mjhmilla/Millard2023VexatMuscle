function [curveSampleECMHalf, curveSampleTitin, curveSampleTitinActive,...
          curveSampleIgp, curveSamplePevk] = ...
  sampleTitinCurves(sampleVectorHalfLength,...
                    forceLengthECMHalfCurve,...
                    forceLengthIgpCurve,...
                    forceLengthIgpInverseCurve,...
                    forceLengthPevkCurve,...
                    forceLengthPevkInverseCurve,...
                    normLengthZToT12,...
                    normLengthIgpAtOptimalFiberLength,...
                    normLengthIgdFixed)
                                
z0 = zeros(size(sampleVectorHalfLength));
curveSampleECMHalf  = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                         'd3ydx3',z0,'intYdx',[]);
curveSampleIgp  = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                         'd3ydx3',z0,'intYdx',[]);
curveSamplePevk = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                         'd3ydx3',z0,'intYdx',[]);
curveSampleTitin = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                         'd3ydx3',z0,'intYdx',[]);
curveSampleTitinActive = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                         'd3ydx3',z0,'intYdx',[]); 

                       

for i=1:1:length(sampleVectorHalfLength)
  xH = sampleVectorHalfLength(i,1);

  curveSampleECMHalf.x(i,1)=xH;
  curveSampleECMHalf.y(i,1)=...
    calcBezierYFcnXDerivative(xH, forceLengthECMHalfCurve,0);
  curveSampleECMHalf.dydx(i,1)=...
    calcBezierYFcnXDerivative(xH, forceLengthECMHalfCurve,1);
  curveSampleECMHalf.d2ydx2(i,1)=...
    calcBezierYFcnXDerivative(xH, forceLengthECMHalfCurve,2);
  curveSampleECMHalf.d3ydx3(i,1)=...
    calcBezierYFcnXDerivative(xH, forceLengthECMHalfCurve,3);

  %Solve for the length of the IgP and PEVK sections such that the 
  %two elements are in a force equilbrium (both are in series) and the 
  %total titin length is the same as xH.
  f = 0.1;
  df = 0;
  lerr = Inf;
  iter=0;
  iterMax=100;
  while(abs(lerr) > 1e-6 && iter < iterMax)
    ligpH = calcBezierYFcnXDerivative(f,forceLengthIgpInverseCurve,0);
    dligpH = calcBezierYFcnXDerivative(f,forceLengthIgpInverseCurve,1);

    lpevkH = calcBezierYFcnXDerivative(f,forceLengthPevkInverseCurve,0);
    dlpevkH = calcBezierYFcnXDerivative(f,forceLengthPevkInverseCurve,1);            

    lerr = (ligpH + lpevkH + normLengthZToT12 + normLengthIgdFixed)-xH;
    dlerr= dligpH + dlpevkH;
    df   = -lerr/dlerr;
    f    = f+df;

    iter=iter+1;
  end
  assert(abs(lerr)<=1e-6);

  curveSampleIgp.x(i,1)=ligpH;
  curveSampleIgp.y(i,1)=...
    calcBezierYFcnXDerivative(ligpH,forceLengthIgpCurve,0);
  curveSampleIgp.dydx(i,1)=...
    calcBezierYFcnXDerivative(ligpH,forceLengthIgpCurve,1);
  curveSampleIgp.d2ydx2(i,1)=...
    calcBezierYFcnXDerivative(ligpH,forceLengthIgpCurve,2);
  curveSampleIgp.d3ydx3(i,1)=...
    calcBezierYFcnXDerivative(ligpH,forceLengthIgpCurve,3);    

  curveSamplePevk.x(i,1)=lpevkH;
  curveSamplePevk.y(i,1)=...
    calcBezierYFcnXDerivative(lpevkH,forceLengthPevkCurve,0);
  curveSamplePevk.dydx(i,1)=...
    calcBezierYFcnXDerivative(lpevkH,forceLengthPevkCurve,1);
  curveSamplePevk.d2ydx2(i,1)=...
    calcBezierYFcnXDerivative(lpevkH,forceLengthPevkCurve,2);
  curveSamplePevk.d3ydx3(i,1)=...
    calcBezierYFcnXDerivative(lpevkH,forceLengthPevkCurve,3);     

  curveSampleTitin.x(i,1) = xH;
  curveSampleTitin.y(i,1) = f;

  kigp  = curveSampleIgp.dydx(i,1);
  kpevk = curveSamplePevk.dydx(i,1);
  curveSampleTitin.dydx(i,1) = ((1/kigp)+(1/kpevk))^(-1);

  %These higher derivatives are easy to derive, but I don't need them
  %for these plots.
  curveSampleTitin.d2ydx2(i,1) = NaN;
  curveSampleTitin.d3ydx3(i,1) = NaN;

  ligpHSaturated = ligpH;
  if(ligpH > normLengthIgpAtOptimalFiberLength)
    ligpHSaturated=normLengthIgpAtOptimalFiberLength;
  end

  curveSampleTitinActive.x(i,1)     = normLengthZToT12 + normLengthIgdFixed + ligpHSaturated + lpevkH;
  curveSampleTitinActive.y(i,1)     = curveSamplePevk.y(i,1);
  curveSampleTitinActive.dydx(i,1)  = curveSamplePevk.dydx(i,1);
  curveSampleTitinActive.d2ydx2(i,1)= NaN;
  curveSampleTitinActive.d3ydx3(i,1)= NaN;
end