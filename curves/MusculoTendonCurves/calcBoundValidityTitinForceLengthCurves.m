function boundsAreValid = ...
  calcBoundValidityTitinForceLengthCurves(arg,argLB,argUB)

boundsAreValid = 1;

for i=1:1:length(arg)
  if( ~(arg(i) >= argLB(i) && arg(i) <= argUB(i)) )
    boundsAreValid=0;
  end
end

