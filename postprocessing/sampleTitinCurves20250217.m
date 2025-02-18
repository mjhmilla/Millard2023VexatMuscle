%%
% SPDX-FileCopyrightText: 2023 Matthew Millard <millard.matthew@gmail.com>
%
% SPDX-License-Identifier: MIT
%
% If you use this code in your work please cite the pre-print of this paper
% or the most recent peer-reviewed version of this paper:
%
%    Matthew Millard, David W. Franklin, Walter Herzog. 
%    A three filament mechanistic model of musculotendon force and impedance. 
%    bioRxiv 2023.03.27.534347; doi: https://doi.org/10.1101/2023.03.27.534347 
%
%%

function titinCurveSample = ...
  sampleTitinCurves20250217(curves, sarcomere, npts)

xH0 = 0;
xH1 = curves.fiberForceLengthCurve.xEnd(1,2)*0.5;  
dxH = (xH1-xH0)/(npts-1);

sampleVectorHalfLength = [xH0:dxH:xH1]';


forceLengthECMHalfCurve = curves.forceLengthECMHalfCurve;

forceLengthProximalTitinCurve        = curves.forceLengthProximalTitinCurve;
forceLengthProximalTitinInverseCurve = curves.forceLengthProximalTitinInverseCurve;
forceLengthDistalTitinCurve          = curves.forceLengthDistalTitinCurve;
forceLengthDistalTitinInverseCurve   = curves.forceLengthDistalTitinInverseCurve;
forceLengthIgPTitinCurve             = curves.forceLengthIgPTitinCurve;
forceLengthIgPTitinInverseCurve      = curves.forceLengthIgPTitinInverseCurve;
forceLengthPevkTitinCurve            = curves.forceLengthPevkTitinCurve;
forceLengthPevkTitinInverseCurve     = curves.forceLengthPevkTitinInverseCurve;
forceLengthIgDTitinCurve             = curves.forceLengthIgDTitinCurve;
forceLengthIgDTitinInverseCurve      = curves.forceLengthIgDTitinInverseCurve;

normLengthZToT12    = sarcomere.ZLineToT12NormLengthAtOptimalFiberLength;
normLengthIgdFixed  = sarcomere.IGDFixedNormLengthAtOptimalFiberLength;

normLengthTitinActinBondMinimum = sarcomere.normLengthTitinActinBondMinimum;
normPevkToActinAttachmentPoint  = sarcomere.normPevkToActinAttachmentPoint;
titinModelType                  = sarcomere.titinModelType;

                                
z0 = zeros(size(sampleVectorHalfLength));
curveSampleECMHalf  = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                         'd3ydx3',z0,'intYdx',[]);

curveSampleProximal  = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                         'd3ydx3',z0,'intYdx',[]);
curveSampleDistal  = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                         'd3ydx3',z0,'intYdx',[]);

curveSampleIgP  = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                         'd3ydx3',z0,'intYdx',[]);
curveSamplePevk = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                         'd3ydx3',z0,'intYdx',[]);
curveSampleIgD  = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                         'd3ydx3',z0,'intYdx',[]);

curveSampleTitin = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                         'd3ydx3',z0,'intYdx',[]);
curveSampleTitinActive = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                         'd3ydx3',z0,'intYdx',[]); 

              
%
% Solve for the distal and proximal lengths of titin just when
% the sarcomere length reaches the minimum length at which titin
% can bond to actin
%

f = 0.1;
df = 0;
lerr = Inf;
iter=0;
iterMax=100;

xH = normLengthTitinActinBondMinimum*0.5;

while(abs(lerr) > 1e-6 && iter < iterMax)
    lPH   = calcBezierYFcnXDerivative(f,forceLengthProximalTitinInverseCurve,0);
    dPH  = calcBezierYFcnXDerivative(f,forceLengthProximalTitinInverseCurve,1);
    
    lDH   = calcBezierYFcnXDerivative(f,forceLengthDistalTitinInverseCurve,0);
    dDH  = calcBezierYFcnXDerivative(f,forceLengthDistalTitinInverseCurve,1);
    
    lerr = (lPH + lDH + normLengthZToT12 + normLengthIgdFixed)-xH;
    dlerr= dPH + dDH;
    df   = -lerr/dlerr;
    f    = f+df;
    
    iter=iter+1;
end
assert(abs(lerr)<=1e-6);

maxActiveLengthLPH = lPH;

%
%
% Numerically evaluate the 3 and 2 segment titin models
%
%

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
    lPH   = calcBezierYFcnXDerivative(f,forceLengthProximalTitinInverseCurve,0);
    dPH  = calcBezierYFcnXDerivative(f,forceLengthProximalTitinInverseCurve,1);

    lDH   = calcBezierYFcnXDerivative(f,forceLengthDistalTitinInverseCurve,0);
    dDH  = calcBezierYFcnXDerivative(f,forceLengthDistalTitinInverseCurve,1);

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
    calcBezierYFcnXDerivative(lPH,forceLengthProximalTitinCurve,0);
  curveSampleProximal.dydx(i,1)=...
    calcBezierYFcnXDerivative(lPH,forceLengthProximalTitinCurve,1);
  curveSampleProximal.d2ydx2(i,1)=...
    calcBezierYFcnXDerivative(lPH,forceLengthProximalTitinCurve,2);
  curveSampleProximal.d3ydx3(i,1)=...
    calcBezierYFcnXDerivative(lPH,forceLengthProximalTitinCurve,3);   

  curveSampleDistal.x(i,1)=lDH;
  curveSampleDistal.y(i,1)=...
    calcBezierYFcnXDerivative(lDH,forceLengthDistalTitinCurve,0);
  curveSampleDistal.dydx(i,1)=...
    calcBezierYFcnXDerivative(lDH,forceLengthDistalTitinCurve,1);
  curveSampleDistal.d2ydx2(i,1)=...
    calcBezierYFcnXDerivative(lDH,forceLengthDistalTitinCurve,2);
  curveSampleDistal.d3ydx3(i,1)=...
    calcBezierYFcnXDerivative(lDH,forceLengthDistalTitinCurve,3);   

  %Solve for the length of the IgP, PEVK, and IgD sections such that the 
  %three elements are in a force equilbrium (all are in series) and the 
  %total titin length is the same as xH.
  f = 0.1;
  df = 0;
  lerr = Inf;
  iter=0;
  iterMax=100;
  while(abs(lerr) > 1e-6 && iter < iterMax)
    ligpH   = calcBezierYFcnXDerivative(f,forceLengthIgPTitinInverseCurve,0);
    dligpH  = calcBezierYFcnXDerivative(f,forceLengthIgPTitinInverseCurve,1);

    lpevkH  = calcBezierYFcnXDerivative(f,forceLengthPevkTitinInverseCurve,0);
    dlpevkH = calcBezierYFcnXDerivative(f,forceLengthPevkTitinInverseCurve,1);            

    ligdH  = calcBezierYFcnXDerivative(f,forceLengthIgDTitinInverseCurve,0);
    dligdH = calcBezierYFcnXDerivative(f,forceLengthIgDTitinInverseCurve,1);            


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

  curveSampleIgP.x(i,1)=ligpH;
  curveSampleIgP.y(i,1)=...
    calcBezierYFcnXDerivative(ligpH,forceLengthIgPTitinCurve,0);
  curveSampleIgP.dydx(i,1)=...
    calcBezierYFcnXDerivative(ligpH,forceLengthIgPTitinCurve,1);
  curveSampleIgP.d2ydx2(i,1)=...
    calcBezierYFcnXDerivative(ligpH,forceLengthIgPTitinCurve,2);
  curveSampleIgP.d3ydx3(i,1)=...
    calcBezierYFcnXDerivative(ligpH,forceLengthIgPTitinCurve,3);    

  curveSamplePevk.x(i,1)=lpevkH;
  curveSamplePevk.y(i,1)=...
    calcBezierYFcnXDerivative(lpevkH,forceLengthPevkTitinCurve,0);
  curveSamplePevk.dydx(i,1)=...
    calcBezierYFcnXDerivative(lpevkH,forceLengthPevkTitinCurve,1);
  curveSamplePevk.d2ydx2(i,1)=...
    calcBezierYFcnXDerivative(lpevkH,forceLengthPevkTitinCurve,2);
  curveSamplePevk.d3ydx3(i,1)=...
    calcBezierYFcnXDerivative(lpevkH,forceLengthPevkTitinCurve,3);     

  curveSampleIgD.x(i,1)=ligdH;
  curveSampleIgD.y(i,1)=...
    calcBezierYFcnXDerivative(ligdH,forceLengthIgDTitinCurve,0);
  curveSampleIgD.dydx(i,1)=...
    calcBezierYFcnXDerivative(ligdH,forceLengthIgDTitinCurve,1);
  curveSampleIgD.d2ydx2(i,1)=...
    calcBezierYFcnXDerivative(ligdH,forceLengthIgDTitinCurve,2);
  curveSampleIgD.d3ydx3(i,1)=...
    calcBezierYFcnXDerivative(ligdH,forceLengthIgDTitinCurve,3);  

  curveSampleTitin.x(i,1) = xH;
  curveSampleTitin.y(i,1) = f;

  kigp  = curveSampleIgP.dydx(i,1);
  kpevk = curveSamplePevk.dydx(i,1);
  kigd  = curveSampleIgD.dydx(i,1);
  curveSampleTitin.dydx(i,1) = ((1/kigp)+(1/kpevk)+(1/kigd))^(-1);

  %These higher derivatives are easy to derive, but I don't need them
  %for these plots.
  curveSampleTitin.d2ydx2(i,1) = NaN;
  curveSampleTitin.d3ydx3(i,1) = NaN;

  


  if(titinModelType==0)
      lPHa = lPH;
      if(xH >=  (normLengthTitinActinBondMinimum*0.5))
        lPHa = maxActiveLengthLPH;
      end


      lDHa = xH-(lPHa + normLengthZToT12 + normLengthIgdFixed);

      curveSampleTitinActive.x(i,1) = xH;
      curveSampleTitinActive.y(i,1) = ...
           calcBezierYFcnXDerivative(lDHa,forceLengthDistalTitinCurve,0);
      curveSampleTitinActive.dydx(i,1) = ...
           calcBezierYFcnXDerivative(lDHa,forceLengthDistalTitinCurve,1);
      curveSampleTitinActive.dydx(i,2) = ...
           calcBezierYFcnXDerivative(lDHa,forceLengthDistalTitinCurve,2);
      curveSampleTitinActive.dydx(i,3) = ...
           calcBezierYFcnXDerivative(lDHa,forceLengthDistalTitinCurve,3);

  end
  
  assert(titinModelType==0,'Error: Have not updated this function for titinModel 1');



  titinCurveSample.curveSampleECMHalf     = curveSampleECMHalf    ; 
  titinCurveSample.curveSampleTitin       = curveSampleTitin      ;     
  titinCurveSample.curveSampleTitinActive = curveSampleTitinActive;   
  titinCurveSample.curveSampleIgP         = curveSampleIgP        ;   
  titinCurveSample.curveSamplePevk        = curveSamplePevk       ;  
  titinCurveSample.curveSampleIgD         = curveSampleIgD        ;  
  titinCurveSample.curveSampleProximal    = curveSampleProximal   ;     
  titinCurveSample.curveSampleDistal      = curveSampleDistal     ;     

%  if(titinModelType==1)
%      %With this model the prox and distal Ig segments are lumped together.
%       curveSampleTitinActive.x(i,1)     = ...
%           normLengthZToT12 ...
%           + normLengthIgdFixed ...
%           + ligpH ...
%           + lpevkHSaturated;
%       curveSampleTitinActive.y(i,1)     = curveSampleIgP.y(i,1);
%       curveSampleTitinActive.dydx(i,1)  = curveSampleIgP.dydx(i,1);
%       curveSampleTitinActive.d2ydx2(i,1)= NaN;
%       curveSampleTitinActive.d3ydx3(i,1)= NaN;
%  end

end