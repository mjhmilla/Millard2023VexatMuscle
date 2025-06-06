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

function val = interpLocallyMonotonicCurve( ...
                inputVector, ...
                outputVector, ...
                inputSample, ...
                numberOfInterpPoints)

val = NaN;              
idxA = 1;
idxB = 1;

z = (numberOfInterpPoints+1);
flag_intervalFound=0;

while flag_intervalFound == 0 && z < (length(inputVector)-numberOfInterpPoints)  
  if(    inputSample >= inputVector(z-numberOfInterpPoints,1) ...
      && inputSample <= inputVector(z+numberOfInterpPoints,1))
    idxA=z-numberOfInterpPoints;
    idxB=z+numberOfInterpPoints;    
    flag_intervalFound = 1;
  end
  z=z+1;
end

if(flag_intervalFound == 1)
  diffInput = diff(inputVector(idxA:idxB,1));
  diffOutput = diff(outputVector(idxA:idxB,1));

  idxInputNonIncreasing   = find(diffInput <= 0, 1);
  idxOutputNonIncreasing  = find(diffOutput <= 0, 1);

  if( isempty(idxInputNonIncreasing) && isempty(idxOutputNonIncreasing) )
    val = interp1( inputVector(idxA:idxB,1),...
                   outputVector(idxA:idxB,1), inputSample);
  else
    val = outputVector(idxA+numberOfInterpPoints,1);
  end
else
  
  if(inputSample >= max(inputVector))
    idxA = length(inputVector)-numberOfInterpPoints;
    idxB = length(inputVector);
    
    val = interp1(  inputVector(idxA:idxB,1),...
                   outputVector(idxA:idxB,1),...
                   inputSample, 'linear', 'extrap');
  end
  if(inputSample <= min(inputVector))
    idxA = 1;
    idxB = numberOfInterpPoints;
    
    dInputVector = abs(diff(inputVector(idxA:idxB,1)));
    idxNonZero = find(dInputVector > 0);


    if(isempty(idxNonZero)==0)
        
        val = interp1(  inputVector(idxNonZero,1),...
                       outputVector(idxNonZero,1),...
                       inputSample, 'linear', 'extrap');
    else
        val = mean(outputVector(idxA:idxB,1));
    end
  end
  
end

