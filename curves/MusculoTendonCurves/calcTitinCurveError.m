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

function curveError = calcTitinCurveError( ...
                                      forceLengthPevkIgdInverseCurve, ...
                                      forceLengthIgpInverseCurve, ...
                                      fiberForceLengthData,...
                                      lT12, lmH)

curveError    = zeros(size(fiberForceLengthData,1),1);



%forceLengthIgP = zeros(size(fiberForceLengthData));
%forceLengthPevkIgD = zeros(size(fiberForceLengthData));      

f = 0;
for i=1:1:size(fiberForceLengthData,1)


  %Now solve for the force which yields this length
  iter=1;
  iterMax=100;
  err=1;
  tol = 1e-6;
  
  while(abs(err) > tol && iter < iterMax)
    lIgp     = calcBezierYFcnXDerivative(f,forceLengthIgpInverseCurve,0);
    lPevkIgd = calcBezierYFcnXDerivative(f,forceLengthPevkIgdInverseCurve,0);
    err = (lIgp+lPevkIgd+lT12+lmH)*2 - fiberForceLengthData(i,1);

    dlIgp     = calcBezierYFcnXDerivative(f,forceLengthIgpInverseCurve,1);
    dlPevkIgd = calcBezierYFcnXDerivative(f,forceLengthPevkIgdInverseCurve,1);          
    derr = (dlIgp+dlPevkIgd)*2;

    delta = -err/derr;
    f = f+delta;

    iter=iter+1;

  end

  assert(abs(err) <= tol);   
  
  %if(f > fiberForceLengthData(i,2))
  %  flag_feasible = 0;
  %end
  
  curveError(i,1) = f - fiberForceLengthData(i,2);
  
  %forceLengthIgP(i,1) = lIgp;
  %forceLengthIgP(i,2) = f;
  %forceLengthPevkIgD(i,1) = lPevkIgd;
  %forceLengthPevkIgD(i,2) = f;

     
end
        
