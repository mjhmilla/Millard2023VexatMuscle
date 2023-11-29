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

function [success] = calcSignalGainAndPhase(...
                        normFiberLength,...
                        nominalForce,...
                        nominalForceSteps,...
                        amplitudeMM,...
                        bandwidthHz,...
                        numberOfSimulations,...
                        minFreqHz,...
                        coherenceSqThreshold,...
                        simSeriesFiles,...   
                        simSeriesNames,...
                        simSeriesColors,...
                        outputFolder,...
                        flag_plotStiffnessDamping,...
                        flag_plotDetailedSpectrumData,...
                        flag_zeroPaddingData,...
                        flag_usingOctave)
success = 0;

%All of the inputFunctions should be identical across the benchRecord
%structs. Here we load one as a reference and later check that all others
%match.
benchRecordSet = load([outputFolder,simSeriesFiles{1}]);
inputFunctions = benchRecordSet.inputFunctions;

assert(max(abs(amplitudeMM-benchRecordSet.amplitudeMM))<1e-10,...
    'Error: amplitudeMM does not match benchRecord.amplitudeMM');
assert(max(abs(bandwidthHz-benchRecordSet.bandwidthHz))<1e-10,...
    'Error: bandwidthHz does not match benchRecord.bandwidthHz');

samplePoints=inputFunctions.samples;


freqSimData = struct('force',  zeros(samplePoints, numberOfSimulations),...      
                     'gain',   zeros(samplePoints*2, numberOfSimulations),...
                     'phase',  zeros(samplePoints*2, numberOfSimulations),...
                     'gainKD',  zeros(samplePoints*2, numberOfSimulations),...
                     'phaseKD',         zeros(samplePoints*2, numberOfSimulations),...                       
                     'coherenceSq',    zeros(samplePoints, numberOfSimulations),... 
                     'coherenceSqFrequency', zeros(samplePoints, numberOfSimulations),... 
                     'freqHz',        zeros(samplePoints, numberOfSimulations),...
                     'freq',          zeros(samplePoints, numberOfSimulations),...
                     'idxFreqRange',  zeros(2, numberOfSimulations),...                
                     'stiffness',     zeros(1,numberOfSimulations),...
                     'damping',       zeros(1,numberOfSimulations),...
                     'vafTime',           zeros(1,numberOfSimulations),...
                     'vafGain',           zeros(1,numberOfSimulations),...          
                     'vafPhase',          zeros(1,numberOfSimulations),...                               
                     'forceKD',       zeros(samplePoints,numberOfSimulations),...
                     'amplitudeMM',   zeros(1,numberOfSimulations),...
                     'bandwidthHz',   zeros(1,numberOfSimulations),...
                     'normFiberLength', zeros(1,numberOfSimulations),...
                     'activation',      zeros(1,numberOfSimulations),...
                     'nominalForceDesired', zeros(1,numberOfSimulations),...
                     'nominalForce'      ,zeros(1,numberOfSimulations)); 


 

flag_useKoopmansTransformEstimator=1;

if(flag_plotDetailedSpectrumData == 1)
  fig_phaseMdl        = figure;
  fig_magMdl          = figure;
  fig_coherenceSqMdl  = figure;
  fig_timeResponseMdl = figure;
  fig_kin             = figure;      
end
if(flag_plotDetailedSpectrumData == 1 || ...
   flag_plotStiffnessDamping == 1)
  fig_stiffness = figure;
  fig_damping   = figure;     
end

freqMax = inputFunctions.sampleFrequency;




for idxModel = 1:1:length(simSeriesFiles)
  benchRecordSet = load([outputFolder,simSeriesFiles{idxModel}]);
  benchRecord = benchRecordSet.benchRecord;
  idx=1;           

  %Check for consistency across all of the benchRecords
  assert(max(abs(amplitudeMM-benchRecordSet.amplitudeMM))<1e-10,...
      'Error: amplitudeMM does not match benchRecord.amplitudeMM');
  assert(max(abs(bandwidthHz-benchRecordSet.bandwidthHz))<1e-10,...
      'Error: bandwidthHz does not match benchRecord.bandwidthHz');

  flag_Hill = 0;
  if(isempty(strfind(simSeriesFiles{idxModel},'Hill'))==0)
    flag_Hill=1;
  end

  tag = 'benchRecord';
  z = strfind(simSeriesFiles{idxModel},tag);
  outputFileName = ['freqResponse',...
                    simSeriesFiles{idxModel}(1,(z+length(tag)):end)];

  for idxNormFiberLength = 1:1:length(normFiberLength)
    for idxActivation = 1:1:nominalForceSteps     
      for i=1:1:length(amplitudeMM)        
        for j=1:1:length(bandwidthHz)

          here=0;
          if(idx==21)
            here=1;
          end
          % Save all of the configuration data
          freqSimData.amplitudeMM(1,idx)      = amplitudeMM(i);
          freqSimData.bandwidthHz(1,idx)      = bandwidthHz(j);
          freqSimData.normFiberLength(1,idx)  = normFiberLength(idxNormFiberLength);
          freqSimData.activation(1,idx)       = benchRecord.activation(1,idx);
          
          freqSimData.nominalForceDesired(1,idx) = nominalForce(1,idxActivation);
          
          %idxMidPadding = round(0.5*inputFunctions.padding);
          freqSimData.nominalForce(1,idx)     = ...
            benchRecord.tendonForce(end,idx);


          idxWave = getSignalIndex(amplitudeMM(i),bandwidthHz(j),...
                                    inputFunctions);

          xErr = max(abs(inputFunctions.x( inputFunctions.idxSignal, idxWave) ...
                 - benchRecordSet.inputFunctions.x( inputFunctions.idxSignal, idxWave)));

          %Check for consistency across all of the benchRecords          
          assert(xErr<1e-10,...
              ['Error: inputFunction structure is inconsistent between',...
               'simulation results. This likely means that some of the ',...
               'output files are from different incomplete runs of ',...
               'main_KirschBoskovRymer.m'] );

          x  = inputFunctions.x(   inputFunctions.idxSignal, idxWave);
          
          y  = benchRecord.tendonForce(inputFunctions.idxSignal, idx);
          yo = y(round(inputFunctions.padding*0.5),1);
          y  = y - yo;
          freqSimData.force(:,idx)  =  benchRecord.tendonForce(:, idx);

          %On the very last simulation in the force series Vexat is
          %poorly initialized and this means
          %that the model spends most of the padding zone going from an 
          %activation of 0 to 1. If this data is included in the subsequent
          %analysis it throws the results off. Until I've figured out what
          %the initialization problem is I am just zeroing out the padding
          %points.
          if(flag_zeroPaddingData==1)
              padding = inputFunctions.padding;
              xo = x(round(padding*0.75));
              yo = y(round(padding*0.75));
    
              x(1:1:padding) = xo;
              y(1:1:padding) = yo;
          end

          %Evaluate the cross spectral densities
          [cpsd_Gxy,cpsd_FxyHz] = cpsd(x,y,[],[],[],freqMax,'onesided');
          [cpsd_Gxx,cpsd_FxxHz] = cpsd(x,x,[],[],[],freqMax,'onesided');
          [cpsd_Gyy,cpsd_FyyHz] = cpsd(y,y,[],[],[],freqMax,'onesided');
          [cpsd_Gyx,cpsd_FyxHz] = cpsd(y,x,[],[],[],freqMax,'onesided');
          


          %Make sure all of the cpsd results have the same length
          %which is equivalent to making sure that Matlab used the 
          %same window size for each one
          assert(length(cpsd_FyxHz)==length(cpsd_FxxHz));
          assert(length(cpsd_FyxHz)==length(cpsd_FyyHz));
          assert(length(cpsd_FyxHz)==length(cpsd_FxyHz));

          maxIdx = length(cpsd_Gxy);

          freqHz = cpsd_FyxHz;
          freq   = freqHz.*(2*pi);
          freqSimData.freq(1:maxIdx,idx)   =  cpsd_FyxHz.*(pi/180);
          freqSimData.freqHz(1:maxIdx,idx) =  cpsd_FyxHz;          
          freqSimData.gain(1:maxIdx,idx)   =  abs(  cpsd_Gyx./ cpsd_Gxx);
          freqSimData.phase(1:maxIdx,idx)  =  angle(cpsd_Gyx./ cpsd_Gxx);

          freqSimData.coherenceSq(1:maxIdx,idx) = ...
              ( abs(cpsd_Gxy).*abs(cpsd_Gxy) ) ./ (cpsd_Gxx.*cpsd_Gyy) ;
          freqSimData.coherenceSqFrequency(1:maxIdx,idx) = cpsd_FxyHz;

          idxKirschMin = find(freqSimData.freqHz(1:maxIdx,idx) >= max(0,max(minFreqHz)-1), 1);
          idxKirschMax = find(freqSimData.freqHz(1:maxIdx,idx) <= bandwidthHz(j)+1, 1,'last');



          idxLb = idxKirschMin;
          while(freqSimData.coherenceSq(idxLb,idx)<coherenceSqThreshold ...
                  && idxLb < idxKirschMax)
              idxLb=idxLb+1;
          end

          idxUb = idxKirschMax;
          while(freqSimData.coherenceSq(idxUb,idx)<coherenceSqThreshold ...
                  && idxUb > idxKirschMin)
              idxUb=idxUb-1;
          end
          assert(idxUb > idxLb, 'Error: no part of the bandwidth meets the coherence threshold');
          
          idxFreqRange = [idxLb:1:idxUb]';
          idxFreqRangeFull = [1:1:idxUb]';

          freqSimData.idxFreqRange(1,idx)=idxLb;
          freqSimData.idxFreqRange(2,idx)=idxUb;          


          %Solve for the spring damping coefficients of best fit to the
          %data
          [stiffness,damping,exitFlag] = ...
              fitSpringDamperToGainAndPhaseProfiles( ...
                             freq(idxFreqRange,1),...
                             freqSimData.gain(idxFreqRange,idx),...
                             freqSimData.phase(idxFreqRange,idx),...
                             flag_usingOctave);

          freqSimData.stiffness(1,idx) = stiffness;
          freqSimData.damping(1,idx)   = damping;

          %Evaluate the time-domain response and frequency-domain
          %response of the spring-damper model of best fit.

          modelResponseTime = ...
            (inputFunctions.x(inputFunctions.idxSignal,idxWave).*freqSimData.stiffness(1,idx) ...
            +inputFunctions.xdot(inputFunctions.idxSignal,idxWave).*freqSimData.damping(1,idx));


          freqSimData.forceKD(:,idx) = ...
             (inputFunctions.x(:,idxWave)   ).*(freqSimData.stiffness(1,idx)) ...
            +(inputFunctions.xdot(:,idxWave)).*(freqSimData.damping(1,idx));

          modelResponseFreq = ...
            calcFrequencyModelResponse( freqSimData.stiffness(1,idx),...
                                        freqSimData.damping(1,idx),freq);     

          %Evaluate the VAF or Variance Accounted For in the time
          %domain. 
          
          yVar  = var(y);
          ymVar = var(y-modelResponseTime);              
          freqSimData.vafTime(1,idx)     = (yVar-ymVar)/yVar;

          modelGain  =    abs(modelResponseFreq);
          modelPhase =  angle(modelResponseFreq);

          freqSimData.gainKD(1:maxIdx,idx) = modelGain;
          freqSimData.phaseKD(1:maxIdx,idx) = modelPhase;

          gVar  = var(freqSimData.gain(idxFreqRange,idx));
          gmVar = var(freqSimData.gain(idxFreqRange,idx) ...
                     - modelGain(idxFreqRange,1));                                                
          freqSimData.vafGain(1,idx)     = (gVar-gmVar)/gVar;

          pVar  = var(freqSimData.phase(idxFreqRange,idx));
          pmVar = var(freqSimData.phase(idxFreqRange,idx) ...
                     -modelPhase(idxFreqRange,1)); 

          freqSimData.vafPhase(1,idx)     = (pVar-pmVar)/pVar;              


          flag_debugW=0;
          if(flag_debugW==1)
              
              figW = figure;
              subplot(1,3,1);

                  yyaxis left;
                  plot(inputFunctions.time, ...
                       x,...
                       'b');
                  hold on;
                  xlabel('Time (s)');
                  ylabel('Distance (mm)');
    
                  yyaxis right;
                  plot(inputFunctions.time, ...
                       y,...
                       'b');
                  plot(inputFunctions.time,...
                       modelResponseTime, 'm');
                  hold on;
                  ylabel('Force (N)');
                  box off;
    
              subplot(1,3,2);
                  plot(freqHz(idxFreqRange,1),...
                          freqSimData.gain(idxFreqRange,idx),'b');
                  hold on;
                  plot(freqHz(idxFreqRange,1),...
                       abs(modelResponseFreq(idxFreqRange,1)),'--k');
                  hold on;
                  box off
                  xlabel('Frequency (Hz)');
                  ylabel('Gain (N/m)');

              subplot(1,3,3);            
                  plot( freqHz(idxFreqRange,1),...
                        freqSimData.phase(idxFreqRange,idx).*(180/pi),'b');
                  hold on;
                  plot(freqHz(idxFreqRange,1),...
                       angle(modelResponseFreq(idxFreqRange,1)).*(180/pi),'--k');
                  hold on;
                  box off;
                  xlabel('Frequency (Hz)');
                  ylabel('Phase (degrees)');
              fprintf('%1.2f\t%1.1f\t%1.1f : VAF, k, d\n',freqSimData.vafTime(1,idx), stiffness, damping);
              here=1;
              pause(0.1);
              close(figW);
          end



         
          trialLabel = '';

          if(flag_usingOctave==0)
            trialLabel = inputFunctions.labels(idxWave,:);              
          else
            trialLabel = sprintf('%s mm %s Hz',...
              num2str(inputFunctions.amplitudeMM(1,idxWave)),...
              num2str(inputFunctions.bandwidthHz(1,idxWave)));                
          end


          fprintf('%i. %s %s\n',idx, simSeriesNames{idxModel},trialLabel);
          fprintf('  K %1.3f N/mm D %1.3f N/(mm/s) vafTime %1.1f vafGain %1.1f vafPhase %1.1f Exit %i\n',...
            freqSimData.stiffness(1,idx)/1000,...
            freqSimData.damping(1,idx)/1000,...
            freqSimData.vafTime(1,idx)*100,...
            freqSimData.vafGain(1,idx)*100,...
            freqSimData.vafPhase(1,idx)*100,...
            exitFlag);

          greyMix = 0.75;
          dataColor = simSeriesColors(idxModel,:).*(1-greyMix) ...
                             + [0.5,0.5,0.5].*greyMix;
          modelColor = simSeriesColors(idxModel,:);

          if(exitFlag < 1)
            dataColor =[1,0,0];
          end

          idxSubplot = idxWave;
          if(length(amplitudeMM)*length(bandwidthHz) < size(inputFunctions.bandwidthHz,2))
            idxSubplot = mod(idx,length(amplitudeMM)*length(bandwidthHz))+1;
          end

          if(flag_plotDetailedSpectrumData == 1)             
            figure(fig_magMdl)
            subplot(length(amplitudeMM),length(bandwidthHz),idxSubplot);    

              plot(freqHz(idxFreqRangeFull,:)./(2*pi), ...
                   freqSimData.gain(indexFreqRangeFull,idx)./1000,...
                  'Color',dataColor,'LineWidth',2);    
              hold on;
              plot(freqHz(indexFreqRangeFull,:)./(2*pi),...
                   abs(modelResponseFreq(indexFreqRangeFull,:))./1000,...
                   '-','Color',[1,1,1],'LineWidth',2);
              hold on;                
              plot(freqHz(indexFreqRangeFull,:)./(2*pi),...
                   abs(modelResponseFreq(indexFreqRangeFull,:))./1000,...
                   '--','Color',modelColor);
              hold on;

              xlabel('Frequency (Hz)');
              ylabel('Magnitude  (N/mm)');
              title(trialLabel);

            figure(fig_phaseMdl)
            subplot(length(amplitudeMM),length(bandwidthHz),idxSubplot);      

              plot(freqHz(indexFreqRangeFull,:)./(2*pi), ...
                   freqSimData.phase(indexFreqRangeFull,idx).*(180/pi),...
                  'Color',dataColor,'LineWidth',2);    
              hold on;
              plot(freqHz(indexFreqRangeFull,:)./(2*pi),...
                   angle(modelResponseFreq(indexFreqRangeFull,:)).*(180/pi),...
                   '-','Color',[1,1,1],'LineWidth',2);
              hold on;

              plot(freqHz(indexFreqRangeFull,:)./(2*pi),...
                   angle(modelResponseFreq(indexFreqRangeFull,:)).*(180/pi),...
                   '--','Color',modelColor);
              hold on;

              xlabel('Frequency (Hz)');
              ylabel('Phase');
              title(trialLabel);

            figure(fig_coherenceSqMdl)
            subplot(length(amplitudeMM),length(bandwidthHz),idxSubplot);          
              plot(freqHz(indexFreqRangeFull,:)./(2*pi), ...
                   freqSimData.coherenceSq(indexFreqRangeFull,idx),...
                  'Color',dataColor,'LineWidth',2);    
              hold on;

              ylim([min(0,min(freqSimData.coherenceSq(indexCorrFreqRange,idx))),...
                    max(1.2,max(freqSimData.coherenceSq(indexCorrFreqRange,idx)))]);
              xlabel('Frequency (Hz)');
              ylabel('CoherenceSq');
              title(trialLabel);

            figure(fig_timeResponseMdl)
            subplot(length(amplitudeMM),length(bandwidthHz),idxSubplot);      
              timeChunk = [1:1:100]';
              plot(inputFunctions.time(inputFunctions.idxSignal(1,timeChunk),1), ...
                   freqSimData.force(timeChunk,idx),...
                  'Color',dataColor,'LineWidth',2);    
              hold on;
              plot(inputFunctions.time(inputFunctions.idxSignal(1,timeChunk),1), ...
                   modelResponseTime(timeChunk,:)+yo,...
                   '-','Color',[1,1,1],'LineWidth',2);
              hold on;

              plot(inputFunctions.time(inputFunctions.idxSignal(1,timeChunk),1), ...
                   modelResponseTime(timeChunk,:)+yo.*forceNorm,...
                   '--','Color',modelColor);
              hold on;

              xlabel('Time (s)');
              ylabel('Force (N)');
              title(trialLabel);

             if(isempty(benchRecord.dstate) == 0)
               figure(fig_kin)
                subplot(length(amplitudeMM),length(bandwidthHz),idxSubplot);      
                  timeChunk = [1:1:100]';

                  plot(inputFunctions.time(inputFunctions.idxSignal(1,timeChunk),1), ...
                       benchRecord.dstate(inputFunctions.idxSignal(1,timeChunk),idx,1),...
                      'Color',dataColor,'LineWidth',2);    
                  hold on;

                  xlabel('Time (s)');

                  yLabelStr = '';
                  if(flag_Hill==1 )
                    yLabelStr = 'Velocity (m/s)';
                  else
                    yLabelStr = 'Acceleration (m/s2)';                          
                  end

                  ylabel(yLabelStr);
                  title(trialLabel);             
             end
          end
          if(flag_plotDetailedSpectrumData == 1 || ...
             flag_plotStiffnessDamping == 1)

          figure(fig_stiffness)
            subplot(length(amplitudeMM),length(bandwidthHz),idxSubplot);      

              plot(benchRecord.activation(1,idx),...
                   freqSimData.stiffness(1,idx)./1000,...
                  'x','Color',dataColor,'LineWidth',2);    
              hold on;

              xlabel('Activation');
              ylabel('Stiffness (N/mm)');
              title(trialLabel);             

          figure(fig_damping)
            subplot(length(amplitudeMM),length(bandwidthHz),idxSubplot);      

              plot(benchRecord.activation(1,idx),...
                   freqSimData.damping(1,idx)./1000,...
                  'o','Color',dataColor,'LineWidth',2);    
              hold on;

              xlabel('Activation');
              ylabel('Damping (N/(mm/s))');
              title(trialLabel);             


         end

          idx=idx+1;

        end    
      end          
    end
  end
  save([outputFolder,outputFileName],'freqSimData');
end                      
                      
success = 1;
                      