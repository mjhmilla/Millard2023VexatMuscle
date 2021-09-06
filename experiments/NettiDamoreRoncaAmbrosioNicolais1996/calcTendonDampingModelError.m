function errSq = calcTendonDampingModelError(arg, stiffnessData, dampingData,errSqNorm)

errSq=0;

beta0 = arg(1,1);
beta1 = arg(2,1);


for z=1:1:length(dampingData)
  beta = beta0 + beta1.*mean(stiffnessData(z).y);
  err = beta - mean(dampingData(z).y);
  errSq = errSq + err*err;
end

errSq = errSq./errSqNorm;