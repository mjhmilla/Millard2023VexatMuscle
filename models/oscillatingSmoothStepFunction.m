function path = oscillatingSmoothStepFunction(t, length, deltaLength, ...
                                            omega, timeStepStart, timeStepEnd,...
                                            stepHeight)


y  = length+deltaLength*sin(omega*t);
dy = (deltaLength*omega)*cos(omega*t);

if(t > timeStepStart && t < timeStepEnd)
    dt = 1/(timeStepEnd-timeStepStart);
    delta = (t-timeStepStart)*dt;
    delta2=delta*delta;
    delta3=delta2*delta;
    s = stepHeight*(3*delta2-2*delta3);

    ds = stepHeight*(6*delta-6*delta2)*dt;

    y=y+s;
    dy=dy+ds;

end

if(t >= timeStepEnd)
    y=y+stepHeight;
end

path = [y;dy];