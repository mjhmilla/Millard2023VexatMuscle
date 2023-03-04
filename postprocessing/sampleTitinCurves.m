function [curveSampleECMHalf, curveSampleTitin, curveSampleTitinActive,...
          curveSampleIgp, curveSamplePevk,curveSampleIgd,...
          curveSampleProximal, curveSampleDistal] = ...
  sampleTitinCurves(sampleVectorHalfLength,...
                    forceLengthECMHalfCurve,...
                    forceLengthProximalCurve,...
                    forceLengthProximalInverseCurve,...
                    forceLengthDistalCurve,...
                    forceLengthDistalInverseCurve,...
                    forceLengthIgpCurve,...
                    forceLengthIgpInverseCurve,...
                    forceLengthPevkCurve,...
                    forceLengthPevkInverseCurve,...
                    forceLengthIgdCurve,...
                    forceLengthIgdInverseCurve,...
                    normLengthZToT12,...
                    normLengthProximalSegmentAtOptimalFiberLength,...
                    normLengthPEVKAtOptimalFiberLength,...
                    normLengthIgdFixed, ...
                    titinModelType)
                                
z0 = zeros(size(sampleVectorHalfLength));
curveSampleECMHalf  = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                         'd3ydx3',z0,'intYdx',[]);

curveSampleProximal  = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                         'd3ydx3',z0,'intYdx',[]);
curveSampleDistal  = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                         'd3ydx3',z0,'intYdx',[]);

curveSampleIgp  = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                         'd3ydx3',z0,'intYdx',[]);
curveSamplePevk = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                         'd3ydx3',z0,'intYdx',[]);
curveSampleIgd  = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
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


  %Solve for the length of the proximal and distal sections sections such that the 
  %two elements are in a force equilbrium (all are in series) and the 
  %total titin length is the same as xH.
  f = 0.1;
  df = 0;
  lerr = Inf;
  iter=0;
  iterMax=100;
  while(abs(lerr) > 1e-6 && iter < iterMax)
    lPH   = calcBezierYFcnXDerivative(f,forceLengthProximalInverseCurve,0);
    dPH  = calcBezierYFcnXDerivative(f,forceLengthProximalInverseCurve,1);

    lDH   = calcBezierYFcnXDerivative(f,forceLengthDistalInverseCurve,0);
    dDH  = calcBezierYFcnXDerivative(f,forceLengthDistalInverseCurve,1);

    lerr = (lPH + lDH + normLengthZToT12 + normLengthIgdFixed)-xH;
    dlerr= dPH + dDH;
    df   = -lerr/dlerr;
    f    = f+df;

    iter=iter+1;
  end
  f2=f;
  assert(abs(lerr)<=1e-6);

  curveSampleProximal.x(i,1)=lPH;
  curveSampleProximal.y(i,1)=...
    calcBezierYFcnXDerivative(lPH,forceLengthProximalCurve,0);
  curveSampleProximal.dydx(i,1)=...
    calcBezierYFcnXDerivative(lPH,forceLengthProximalCurve,1);
  curveSampleProximal.d2ydx2(i,1)=...
    calcBezierYFcnXDerivative(lPH,forceLengthProximalCurve,2);
  curveSampleProximal.d3ydx3(i,1)=...
    calcBezierYFcnXDerivative(lPH,forceLengthProximalCurve,3);   

  curveSampleDistal.x(i,1)=lDH;
  curveSampleDistal.y(i,1)=...
    calcBezierYFcnXDerivative(lDH,forceLengthDistalCurve,0);
  curveSampleDistal.dydx(i,1)=...
    calcBezierYFcnXDerivative(lDH,forceLengthDistalCurve,1);
  curveSampleDistal.d2ydx2(i,1)=...
    calcBezierYFcnXDerivative(lDH,forceLengthDistalCurve,2);
  curveSampleDistal.d3ydx3(i,1)=...
    calcBezierYFcnXDerivative(lDH,forceLengthDistalCurve,3);   

  %Solve for the length of the IgP, PEVK, and IgD sections such that the 
  %three elements are in a force equilbrium (all are in series) and the 
  %total titin length is the same as xH.
  f = 0.1;
  df = 0;
  lerr = Inf;
  iter=0;
  iterMax=100;
  while(abs(lerr) > 1e-6 && iter < iterMax)
    ligpH   = calcBezierYFcnXDerivative(f,forceLengthIgpInverseCurve,0);
    dligpH  = calcBezierYFcnXDerivative(f,forceLengthIgpInverseCurve,1);

    lpevkH  = calcBezierYFcnXDerivative(f,forceLengthPevkInverseCurve,0);
    dlpevkH = calcBezierYFcnXDerivative(f,forceLengthPevkInverseCurve,1);            

    ligdH  = calcBezierYFcnXDerivative(f,forceLengthIgdInverseCurve,0);
    dligdH = calcBezierYFcnXDerivative(f,forceLengthIgdInverseCurve,1);            


    lerr = (ligpH + lpevkH + ligdH + normLengthZToT12 + normLengthIgdFixed)-xH;
    dlerr= dligpH + dlpevkH + dligdH;
    df   = -lerr/dlerr;
    f    = f+df;

    iter=iter+1;
  end
  f3 = f; 
  assert(abs(lerr)<=1e-6);

  if(abs(f2-f3)>5e-2)
      here=1;
  end


  %assert(abs(f2-f3)<1e-2);

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

  curveSampleIgd.x(i,1)=ligdH;
  curveSampleIgd.y(i,1)=...
    calcBezierYFcnXDerivative(ligdH,forceLengthIgdCurve,0);
  curveSampleIgd.dydx(i,1)=...
    calcBezierYFcnXDerivative(ligdH,forceLengthIgdCurve,1);
  curveSampleIgd.d2ydx2(i,1)=...
    calcBezierYFcnXDerivative(ligdH,forceLengthIgdCurve,2);
  curveSampleIgd.d3ydx3(i,1)=...
    calcBezierYFcnXDerivative(ligdH,forceLengthIgdCurve,3);  

  curveSampleTitin.x(i,1) = xH;
  curveSampleTitin.y(i,1) = f;

  kigp  = curveSampleIgp.dydx(i,1);
  kpevk = curveSamplePevk.dydx(i,1);
  kigd  = curveSampleIgd.dydx(i,1);
  curveSampleTitin.dydx(i,1) = ((1/kigp)+(1/kpevk)+(1/kigd))^(-1);

  %These higher derivatives are easy to derive, but I don't need them
  %for these plots.
  curveSampleTitin.d2ydx2(i,1) = NaN;
  curveSampleTitin.d3ydx3(i,1) = NaN;

  ligpHSaturated = ligpH;
  if(ligpH > normLengthProximalSegmentAtOptimalFiberLength)
    ligpHSaturated=normLengthProximalSegmentAtOptimalFiberLength;
  end

  lpevkHSaturated = lpevkH;
  if(lpevkH > normLengthPEVKAtOptimalFiberLength)
    lpevkHSaturated=normLengthPEVKAtOptimalFiberLength;
  end

  if(titinModelType==0)

      curveSampleTitinActive.x(i,1)     = normLengthZToT12 + normLengthIgdFixed + ligpHSaturated + lpevkH;
      curveSampleTitinActive.y(i,1)     = curveSamplePevk.y(i,1);
      curveSampleTitinActive.dydx(i,1)  = curveSamplePevk.dydx(i,1);
      curveSampleTitinActive.d2ydx2(i,1)= NaN;
      curveSampleTitinActive.d3ydx3(i,1)= NaN;
  end
  if(titinModelType==1)
      %With this model the prox and distal Ig segments are lumped together.
      curveSampleTitinActive.x(i,1)     = normLengthZToT12 + normLengthIgdFixed + ligpH + lpevkHSaturated;
      curveSampleTitinActive.y(i,1)     = curveSampleIgp.y(i,1);
      curveSampleTitinActive.dydx(i,1)  = curveSampleIgp.dydx(i,1);
      curveSampleTitinActive.d2ydx2(i,1)= NaN;
      curveSampleTitinActive.d3ydx3(i,1)= NaN;

  end

end