function  idx =  getSignalIndex(amplitudeMM, bandwidthHz, inputFunctions)

idx = 0;

for z=1:1:length(inputFunctions.amplitudeMM)
  if( abs(inputFunctions.amplitudeMM(z)-amplitudeMM<sqrt(eps)) ...
    && abs(inputFunctions.bandwidthHz(z)-bandwidthHz)<sqrt(eps) )
    idx = z;
  end          
end  