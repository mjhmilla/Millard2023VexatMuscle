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

function [success] = plotPerturbationWaveforms( inputFunctions,...
                                                amplitudeMM,...
                                                bandwidthHz,...
                                                plotLayoutSettings)

success = 0;


numberOfHorizontalPlotColumns = plotLayoutSettings.numberOfHorizontalPlotColumns;
numberOfVerticalPlotRows      = plotLayoutSettings.numberOfVerticalPlotRows;
flag_fixedPlotWidth           = plotLayoutSettings.flag_fixedPlotWidth;
plotWidth                     = plotLayoutSettings.plotWidth;
plotHeight                    = [];

flag_usingOctave              = plotLayoutSettings.flag_usingOctave;
plotConfig;

paddingPoints  = inputFunctions.padding;
sampleFrequency= inputFunctions.sampleFrequency;
samplePoints   = inputFunctions.samples;



  fig_base =figure;
  subplot(1,2,1);
    plot(inputFunctions.time(inputFunctions.idxSignal,1), inputFunctions.xo,'b');
    hold on;
    xlabel('Time (s)');
    ylabel('Magnitude');
    title('Base Function: Time Domain');

  subplot(1,2,2);
    plot(inputFunctions.freqHz(1:1:(1+samplePoints/2),1), inputFunctions.po,'b');
    hold on;
    xlabel('Frequency (Hz)');
    ylabel('Single-Sided Amplitude');
    title('Base Function: Frequency Domain');


    fig_functionsTime      = figure;
    fig_functionsFreq      = figure;
    fig_reconstructionTime = figure;

    idx=1;
    for i=1:1:length(amplitudeMM)
      for j=1:1:length(bandwidthHz)

        for z=1:1:length(inputFunctions.amplitudeMM)
          if( abs(inputFunctions.amplitudeMM(z)-amplitudeMM(i))<sqrt(eps) ...
            && abs(inputFunctions.bandwidthHz(z)-bandwidthHz(j))<sqrt(eps))
            idx = z;
          end          
        end
        
        %Scale the base signal
        xScaling = inputFunctions.scaling(1,idx);

        indexChunk =[1:1:floor(inputFunctions.samples*0.5)]';

        idxSignal = floor(0.1*inputFunctions.samples);
        
        idxX       = [1:1:idxSignal];
        idxXoTime  = [1:1:idxSignal];
        idxXo      = [1:1:idxSignal];
        
        
        %Plot the time domain signal
        figure(fig_functionsTime)
        subplot(length(amplitudeMM),length(bandwidthHz),idx);
          %xoScaling = scaling/max(inputFunctions.xo);
          plot( inputFunctions.time(idxXoTime,1),...
                inputFunctions.xo(idxXo,1).*(xScaling.*1000),'Color',[1,1,1].*0.5);
          hold on;    
          plot(inputFunctions.time(idxX,1),inputFunctions.x(idxX,idx).*1000,...
              'k','LineWidth',2);    
          box off;
          ylim([-max(amplitudeMM),max(amplitudeMM)].*(2));
          xlabel('Time (s)');
          ylabel('Length (mm)');
          title(inputFunctions.labels(idx,:));

        %Plot the frequency domain signal
        figure(fig_functionsFreq)
        subplot(length(amplitudeMM),length(bandwidthHz),idx);
          %xoScaling = scaling/max(inputFunctions.xo);
          plot( inputFunctions.freqHz(1:1:(1+samplePoints/2),1),...
                inputFunctions.po.*xScaling,'Color',[1,1,1].*0.5);
          hold on;    
          plot(inputFunctions.freqHz(1:1:(1+samplePoints/2),1),inputFunctions.p(:,idx),...
              'k','LineWidth',2);    
          box off;
          xlabel('Freq (Hz)');
          ylabel('Magnitude');
          title(inputFunctions.labels(idx,:));

        %Plot the time domain signal
        xre = ifft(inputFunctions.y(:,idx));
        
        figure(fig_reconstructionTime)
        subplot(length(amplitudeMM),length(bandwidthHz),idx);          
          plot(inputFunctions.time(inputFunctions.idxSignal,1),...
               (xre-inputFunctions.x(inputFunctions.idxSignal,idx)).*1000,...
              'r','LineWidth',0.5);              
          box off;
          ylim([-max(amplitudeMM),max(amplitudeMM)].*(2));
          xlabel('Time (s)');
          ylabel('Error (mm)');
          title(inputFunctions.labels(idx,:));

        %idx=idx+1;
          
      end
    end  

success = 1;