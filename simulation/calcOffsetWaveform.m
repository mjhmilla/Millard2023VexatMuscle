function waveState = calcOffsetWaveform(t, lengthOffset, timeSeries,lengthSeries,lengthDotSeries)

waveState = zeros(2,1);

waveState(1) = interp1(timeSeries,lengthDotSeries,t);
waveState(2) = interp1(timeSeries,lengthSeries,t)+lengthOffset;
