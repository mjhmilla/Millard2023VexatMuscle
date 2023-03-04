function sse = calcFrequencyDomainSquaredError(x, dataFreqRadians, dataGain,...
   dataPhaseRadians, argScaling, objScaling)

sse = 0;
assert(length(x)==2)

k    = x(1)*argScaling(1);
beta = x(2)*argScaling(2);

for i=1:1:length(dataFreqRadians)
  
  %modelResponse = k + (beta*complex(0,-1)*dataFreqRadians(i,1));      
  modelResponse = calcFrequencyModelResponse(k,beta,dataFreqRadians(i,:));
  gainError = (abs(modelResponse) - dataGain(i,1))/(dataGain(end,1));
  phaseError= (angle(modelResponse) - dataPhaseRadians(i,1))/(dataPhaseRadians(end,1));
  
  errorSq = gainError*gainError + phaseError*phaseError;
  
  sse = sse + errorSq/(length(dataFreqRadians));
end

sse = sse.*objScaling;
