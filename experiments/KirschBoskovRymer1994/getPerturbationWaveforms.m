function [inputFunctions] = getPerturbationWaveforms( ...
                              amplitudeMM,...
                              bandwidthHz,...
                              samplePoints,...
                              paddingPoints,...
                              sampleFrequency,...
                              workingFolder,...
                              flag_generateRandomInput,...
                              flag_processInputFunctions,...
                              flag_usingOctave)

  signalFileEnding = sprintf('_%sHz_%s',num2str(sampleFrequency),...
                                     num2str(samplePoints));
         
numberOfFunctions = length(amplitudeMM)*length(bandwidthHz);                          
                                  
totalPoints=samplePoints;

inputFunctions = struct('padding', paddingPoints,...
                        'samples', samplePoints,...
                        'totalpoints', totalPoints,...
                        'idxSignal',  [1:1:samplePoints],...
                        'time',       zeros(totalPoints,1), ...
                        'freqHz',     zeros(samplePoints,1),...
                        'freq',   zeros(samplePoints,1),...                        
                       'function',    zeros(totalPoints, numberOfFunctions),...
                       'scaling',     zeros(1, numberOfFunctions),...
                       'amplitudeMM',   zeros(1,numberOfFunctions),...
                       'bandwidthHz',   zeros(1,numberOfFunctions),...
                       'indexMaxFreq',zeros(1,numberOfFunctions),...
                       'sampleFrequency', sampleFrequency,...                       
                       'x'   ,  zeros(totalPoints, numberOfFunctions),...
                       'xdot',  zeros(totalPoints, numberOfFunctions),...
                       'y',     zeros(samplePoints, numberOfFunctions),...
                       'p',     zeros(samplePoints/2+1, numberOfFunctions),...                       
                       'xo',    zeros(samplePoints, 1),...
                       'yo',    zeros(samplePoints, 1),...
                       'po',    zeros(samplePoints/2+1, 1),...
                       'labels',    strings(numberOfFunctions,1));  
                     
if(flag_generateRandomInput==1)
  xo = randn(samplePoints,1);
  
  assert(paddingPoints*10 < samplePoints,...
    'Padding points exceed 1/10th of the sample');
  
  xo(paddingPoints:1:(samplePoints-paddingPoints)) = ...
      xo(paddingPoints:1:(samplePoints-paddingPoints)) ...
      -mean(xo(paddingPoints:1:(samplePoints-paddingPoints)));

  xo(1:1:paddingPoints) = 0;
  xo((samplePoints-paddingPoints):1:samplePoints) = 0;

  inputFunctions.xo = xo;  
  save([workingFolder,'baseFunction',signalFileEnding,'.mat'],'xo');

else
  tmp = load([workingFolder,'baseFunction',signalFileEnding,'.mat']);
  inputFunctions.xo=tmp.xo;
end

if(flag_processInputFunctions == 1)                     

  idx=1;                     
  for i=1:1:length(amplitudeMM)
    for j=1:1:length(bandwidthHz)
      if(flag_usingOctave==0)
        inputFunctions.labels(idx,:) = [num2str(amplitudeMM(i,1)),'mm ',num2str(bandwidthHz(j,1)),'Hz'];
      end
      inputFunctions.amplitudeMM(1,idx) = amplitudeMM(i);
      inputFunctions.bandwidthHz(1,idx) = bandwidthHz(j);
      idx=idx+1;
    end
  end

  inputFunctions.time = [0:1:(totalPoints-1)]' ...
        .*( (totalPoints/sampleFrequency)/(totalPoints-1));
    
  inputFunctions.yo = fft(inputFunctions.xo(:,1));
  p2 = abs(inputFunctions.yo/samplePoints);
  p1 = p2(1:samplePoints/2 + 1);
  p1(2:end-1) = 2*p1(2:end-1);
  inputFunctions.po = p1;
  inputFunctions.freqHz = ([1:(samplePoints)]' ...
                           ).*((sampleFrequency/samplePoints));
                         
  inputFunctions.freq = inputFunctions.freqHz.*(2*pi); 
                         
  idx=1;                     
  for i=1:1:length(amplitudeMM)
    for j=1:1:length(bandwidthHz)
      inputFunctions.bandwidthHz(1,idx);
      for k=1:1:length(inputFunctions.freq)
        if(inputFunctions.bandwidthHz(1,idx) >= inputFunctions.freqHz(k,1))
          inputFunctions.indexMaxFreq(1,idx) = k;
        end
      end
      idx=idx+1;
    end
  end
  

  %%
  % Create the derived input signals by filtering and scaling the input
  % wave form.
  %%

  idx=1;
  for i=1:1:length(amplitudeMM)
    for j=1:1:length(bandwidthHz)

      %Filter the base signal
      wn = bandwidthHz(j,1)/(0.5*sampleFrequency);
      [b,a] = butter(2, wn, 'low');
      inputFunctions.x(:,idx) = ...
        filter(b,a, inputFunctions.xo(:,1));
                  
      %Scale the base signal
      scaling = (amplitudeMM(i)./1000);
      xScaling = scaling/max(inputFunctions.x(:,idx));
      inputFunctions.x(:,idx) = inputFunctions.x(:,idx).*(xScaling);
      inputFunctions.scaling(1,idx) = xScaling;
      
      %Form the power spectrum in the frequency domain
      inputFunctions.y(:,idx) = fft(inputFunctions.x(inputFunctions.idxSignal,idx));
      p2 = abs(inputFunctions.y(:,idx)/samplePoints);
      p1 = p2(1:samplePoints/2 + 1);
      p1(2:end-1) = 2*p1(2:end-1);
      inputFunctions.p(:,idx) = p1;

      %Form the derivative signal: A central difference is not accurate
      %   enough - it attenuates the high frequency components.
      
      %       xEnd   = length(inputFunctions.x(:,idx));
      %       xLeft  = [1:1:xEnd-1];
      %       xRight = [2:1:xEnd];
      %       dt = 1/inputFunctions.sampleFrequency;
      %       xdotNum = zeros(size(inputFunctions.x,1),1);
      %             
      %       xdotNum(2:1:(xEnd-1),idx) = ...
      %          (diff(inputFunctions.x(xRight,idx)) ...
      %         + diff(inputFunctions.x(xLeft,idx)))./(2*dt);
      % 
      %       xdotNum(1,idx) = xdotNum(2,idx) ...
      %          -(xdotNum(3,idx)-xdotNum(2,idx));
      % 
      %       xdotNum(end,idx) = xdotNum(end-1,idx) ...
      %          +(xdotNum(end-1,idx)-xdotNum(end-2,idx));      

      s = complex(0,1).*(inputFunctions.freq(:,1));
      xdotS = ifft(inputFunctions.y(:,idx).*s,'symmetric'); 
      
      assert(abs(xdotS(1,1)) < 1e-3,...
          ['Initial xdot is ',num2str(abs(xdotS(1,1))),...
          ' but should be close to zero']);
      xdotS(1,1) = 0.;   %Muscle model currently can only be initialized
                         %with a velocity that is < sqrt(eps)
      inputFunctions.xdot(:,idx) = xdotS;
      
      idx=idx+1;
    end
  end
  save([workingFolder,'systemIdInputFunctions',signalFileEnding,'.mat'],...
          'inputFunctions');

else
  tmp=load([workingFolder,'systemIdInputFunctions',signalFileEnding,'.mat']);  
  inputFunctions = tmp.inputFunctions;
  
  %Check to make sure that the loaded functions are consistent with 
  %the above configuration

  assert( inputFunctions.sampleFrequency == sampleFrequency);
  assert( size(inputFunctions.time,1)    == totalPoints);
  
end    