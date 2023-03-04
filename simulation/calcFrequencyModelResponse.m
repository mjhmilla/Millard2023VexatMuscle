function y = calcFrequencyModelResponse(k, beta, omega)

y = k + (beta*complex(0,1)).*omega;