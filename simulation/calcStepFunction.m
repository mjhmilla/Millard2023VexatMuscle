function y = calcStepFunction(t, ton, toff,stepMagnitude)

y = 0;
if(t >= ton && t <= toff)
    y = stepMagnitude;
end