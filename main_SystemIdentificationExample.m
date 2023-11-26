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

clc;
close all;
clear all;

rootDir         = getRootProjectDirectory();
projectFolders  = getProjectFolders(rootDir);

flag_addNonlinearityToOutput = 0;
ampNonlinearity = 0.1;

flag_addTimeVaryingNoiseToOutput = 0;
freqTimeVaryingNoiseHz   = 0.5; %Hz
ampTimeVaryingNoiseHz    = 0.5;

flag_addRandomNoiseToOutput = 0;
randomNoiseAmplitude        = 0.1;


addpath( genpath(projectFolders.postprocessing) );
%%
%Publication configuration
%%
flag_usingOctave=0;
plotWidth = 7.5;
plotHeight= 2.0;
plotHorizMarginCm=1.5;
plotVertMarginCm = 1.5;
numberOfVerticalPlotRows = 7;
numberOfHorizontalPlotColumns=1;

pageWidth = numberOfHorizontalPlotColumns*(plotWidth+plotHorizMarginCm)...
            +plotHorizMarginCm;
pageHeight = numberOfVerticalPlotRows*(plotHeight+plotVertMarginCm) ...
            +plotVertMarginCm;
plotConfigGeneric;

%%
%Generate a bandwidth limited stochastic signal
%%
samples           = 512;
paddingWidth      = 50;
perturbationLength= 1.6;%1.6 mm in units of m
bandwidth         = 35;    
sampleFrequency   = 1/0.002;
nyquistFrequency  = sampleFrequency*0.5;

timeVec           = [0:(1/(samples-1)):1]'.*(samples/sampleFrequency);


[b,a]             =butter(2,bandwidth/nyquistFrequency);
xRandom           = rand(samples,1);
xRandom = xRandom-mean(xRandom);

%Use a rectangular window to fix the nominal length to zero
xRandom(1:paddingWidth)=0;
xRandom((samples-paddingWidth):1:samples) = 0;


%xTimeDomain is the length change: the nominal length has been removed
xTimeDomain       = filter(b,a,xRandom);
xTimeDomainScale  = perturbationLength/ (max(abs(xTimeDomain)));
xTimeDomain = xTimeDomain.*xTimeDomainScale;

%%
%Make some synthetic force data: 
%  This entire section of code is not necessary once we have experimental
%  data
%%
k=4.46;  %Stiffness: 4.46 N/mm
d=0.0089;%Damping  : feel free to adjust

%Form a nearly perfect derivative of x
% Due to numerical error, I presume, the first padding area of xDot is 
% non-zero even though this padding area in x is exactly zero. This
% bit of non-zero noise in xDot is the reason why y doesn't start at zero.
frequencyHz = [0:(1/(samples)): (1-(1/samples)) ]' .* (sampleFrequency);
frequencyRadians = frequencyHz.*(2*pi);
s    = complex(0,1).*frequencyRadians;
xFreqDomain    = fft(xTimeDomain);
xDotFreqDomain = xFreqDomain.*s;
xDotTimeDomain = ifft(xDotFreqDomain,'symmetric');



%Finally evaluate the time domain force response of the spring damper
%where the nominal force has been removed.
%To test a model (or experimental data) you would replace yTimeDomain with
%the simulated (or measured) force with the nominal value subracted off.
yTimeDomain = k.*xTimeDomain + d.*xDotTimeDomain;


if(flag_addNonlinearityToOutput==1)
    yTimeDomain = yTimeDomain + (ampNonlinearity*k).*(xTimeDomain.^2);
end

%Add some time variant noise: a cosine wave that goes from trough to peak
%over the time window
if(flag_addTimeVaryingNoiseToOutput==1)
    yTimeDomain = yTimeDomain ...
        + cos((0.5*pi/max(timeVec)).*timeVec).*ampTimeVaryingNoiseHz;
end

if(flag_addRandomNoiseToOutput==1)
    yTimeDomain = yTimeDomain + randomNoiseAmplitude.*randn(length(yTimeDomain),1);
end    


%%
%Use these two signals to evaluate the gain and phase shift between
%the length and force perturbation
%%


%Evaluate the cross-spectral density between x and y using Welch's method
%Welch's method breaks up the time domain signals into overlapping blocks
%Each block is tranformed into the frequency domain, and the final
%result is the average of the each block in the frequency domain.
%The resulting signal has a lower frequency resolution but is not
%so sensitive to noise
[cpsd_Gxy,cpsd_Fxy] = cpsd(xTimeDomain,yTimeDomain,[],[],[],sampleFrequency,'onesided');
[cpsd_Gxx,cpsd_Fxx] = cpsd(xTimeDomain,xTimeDomain,[],[],[],sampleFrequency,'onesided');
[cpsd_Gyy,cpsd_Fyy] = cpsd(yTimeDomain,yTimeDomain,[],[],[],sampleFrequency,'onesided');
[cpsd_Gyx,cpsd_Fyx] = cpsd(yTimeDomain,xTimeDomain,[],[],[],sampleFrequency,'onesided');

coherenceSq     = ( abs(cpsd_Gyx).*abs(cpsd_Gyx) ) ./ (cpsd_Gxx.*cpsd_Gyy) ;
freqHz          = cpsd_Fyx;
freqRadians     = freqHz.*(2*pi);
idxBW         = find(freqHz <= max(bandwidth));

gain  = abs(cpsd_Gyx./cpsd_Gxx);
phase = angle(cpsd_Gyx./cpsd_Gxx);


%Check this evaluation with Matlab's own internal function
[coherenceSqCheck,freqCpsdCheck] = mscohere(xTimeDomain,yTimeDomain,[],[],[],sampleFrequency);
assert( max(abs(coherenceSqCheck-coherenceSq)) < 1e-6);



figRough = figure;
subplot(3,1,1);
    yyaxis left;
    plot(timeVec, xTimeDomain,'b','DisplayName','x(t)'); 
    hold on;
    ylabel('Length (m)');

    yyaxis right
    plot(timeVec, yTimeDomain,'r','DisplayName','y(t)'); 
    hold on;
    ylabel('Force (N)');
    
    xlabel('Time (s)');
    legend;
    box off;

    title('Time domain signals');
subplot(3,1,2);
    yyaxis left;
    plot(freqHz(idxBW,1),...
         gain(idxBW,1),'-','Color',[0.75,0.75,1],'LineWidth',2,...
         'DisplayName','gain');
    
    yGainLim = ylim;
    if(min(yGainLim)>0)
        ylim([0,max(yGainLim)]);
    end

    ylabel('Gain (N/m)');


    yyaxis right;
    plot(freqHz(idxBW,1),...
         phase(idxBW,1).*(180/pi),'-',...
         'Color',[1,0.75,0.75],'LineWidth',2,...
         'DisplayName','phase');
    
    ylabel('Phase (deg)');

    xlabel('Frequency (Hz)');
    legend; 
    box off;

    title('Frequency domain response');

subplot(3,1,3);
    plot(freqHz(idxBW,1),...
         coherenceSq(idxBW,1),'k');
    hold on;

    if(flag_addTimeVaryingNoiseToOutput==1)
        [d, idx] = min(abs(freqHz-freqTimeVaryingNoiseHz));
        
        plot(freqHz(idx,1),...
             coherenceSq(idx,1),'xb');
        hold on;
        
        plot([freqHz(idx,1),freqHz(idx,1)+4],...
             [coherenceSq(idx,1),coherenceSq(idx,1)-0.05],'b');
        hold on;
    
        text(freqTimeVaryingNoiseHz+5,coherenceSq(idx,1)-0.05,'Frequency of time-varying noise',...
             'HorizontalAlignment','left');
        hold on;
    end
  
    yticks([0:0.25:1]');
    ylim([-0.01,1.01])

    box off;
    xlabel('Frequency (Hz)');
    ylabel('Coherence$$^2$$');
    title('Coherence between input and output signals');

disp('Given the gain and phase profiles you can identify');
disp('what topology fits the data usually by inspection.');
disp('here the pattern follows that of the blue line in ');
disp('Fig. 2 in the pdf (spring-damper in parallel) though');
disp('the parameters differ.');
disp('With the topology identified, it is possible to fit ');
disp('coefficients to the data in the frequency domain by ');
disp('minimizing the squared differences between the gain ');
disp('and phase response of the model and that of the ');
disp('experimental data.');


%%
%Plot the results
%%

figPub=figure;
ax=axes(figPub);
set(ax,'TitleHorizontalAlignment','Left');
subplot('Position',reshape(subPlotPanel(1,1,:),1,4));
%subplot(2,3,1); 

  yyaxis right;
  plot(timeVec, yTimeDomain,'-','Color',[217,81,23]./256,'LineWidth',1,...
       'DisplayName','Output');
  hold on;
  ylabel('Output Amplitude');  
  yticks([-1,0,1].*round(max(abs(yTimeDomain)),1));
  ylim([-1.01,1.01].*round(max(abs(yTimeDomain)),1));

  yyaxis left;
  plot(timeVec, xTimeDomain,'-','Color',[1,1,1],'LineWidth',2,...
       'DisplayName','');
  hold on;
  plot(timeVec, xTimeDomain,'-','Color',[0,112,189]./256,'LineWidth',1,...
       'DisplayName','Input');
  hold on;
  ylim([-1.01,1.01].*round(max(abs(yTimeDomain)),1));
  ylabel('Input Amplitude');

  yticks([-1,0,1].*round(max(abs(xTimeDomain)),1));
  ylim([-1.01,1.01].*round(max(abs(yTimeDomain)),1));
  
%  yticks([-1,0,1].*round(max(abs(xTimeDomain)),2));
%  ylim([-1.01,1.01].*round(max(abs(xTimeDomain)),1))
%  ylim([-1.01,1.01].*round(max(abs(yTimeDomain)),1))


  xlabel('Time (s)');  
  xlim([0,samples/sampleFrequency]);
  xticks([0:0.1:round(samples/sampleFrequency)]);

  limitsX = xlim; 
  limitsY = ylim;

  box off;

  title('A. Input and output signal time-domain data', ...
        'HorizontalAlignment', 'left', 'position', [limitsX(1), limitsY(2)]);
  
  
%subplot('Position',reshape(subPlotPanel(2,1,:),1,4));  
%subplot(2,3,2); 


fullPanel  =   reshape(subPlotPanel(2,1,:),1,4);
leftPanel  = [  fullPanel(1,1),...
                fullPanel(1,2),...
                fullPanel(1,3)*0.4,...
                fullPanel(1,4)];
rightPanel = [fullPanel(1,1)+fullPanel(1,3)*0.6,...
              fullPanel(1,2),...
              fullPanel(1,3)*0.4,...
              fullPanel(1,4)];

subplot('Position',leftPanel);  
  ts = 0.005;
  timeVecExample = [0:ts:(1-ts)]';
  n = length(timeVecExample);

  signalExample    = 0.25*sin((2*pi*30).*timeVecExample) ...
                  + 0.75*sin((2*pi*3).*timeVecExample);
  plot(timeVecExample,signalExample,'-','Color',[0,0,0],'LineWidth',1),
  hold on;
  box off;
  
  text(0.05,1.8,'$$y(t)= 0.25 \sin(2\pi \, 30t)$$','FontSize',6);
  hold on;
  text(0.205,1.5,'$$+ 0.75 \sin(2\pi \, 3t)$$','FontSize',6);
  hold on;
  
  xticks([0,0.2,0.4,0.6,0.8,1]);
  ylim([-1.3,2]);

  limitsX = xlim; 
  limitsY = ylim;
  xlabel('Time (s)');
  ylabel('Magnitude');

  title('B. Time-domain signal', ...
        'HorizontalAlignment', 'left', 'position', [limitsX(1), limitsY(2)]);  

subplot('Position',rightPanel);  
  fftOfExample = fft(signalExample);
  fs = 1/ts;
  freqVec = [0:n-1]'*(fs/n);
  idxNyquistBandwidth = find(freqVec < (0.5*fs));
  magnitudeOfExample = abs(fftOfExample(idxNyquistBandwidth,1))./ (n/2);

  %Plot the vertical lines
  for i=1:1:length(idxNyquistBandwidth)
    x0 = freqVec(idxNyquistBandwidth(i,1),1);
    x1 = freqVec(idxNyquistBandwidth(i,1),1);
    y0 = 0;
    y1 = magnitudeOfExample(i,1);
    
    plot([x0;x1],[y0;y1],'-','Color',[0.5,0.5,1],'LineWidth',1);
    hold on;
    plot([x1],[y1],'.','Color',[0.5,0.5,1],'LineWidth',1);
    hold on;
  end
  text(3+1,0.75,'(3 Hz, 0.75)','FontSize',6,'HorizontalAlignment','left',...
              'VerticalAlignment','top');
  hold on;
  text(30+1,0.25,'(30 Hz, 0.25)','FontSize',6,'HorizontalAlignment','left',...
       'VerticalAlignment','top');
  hold on;

  box off; 
  xlabel('Frequency (Hz)');
  ylabel('Magnitude');
  xticks([0:0.25:1].*(fs/2));
  xlim([0,fs*0.5].*1.01);

  yticks([0.25, 0.75]);
  ylim([0,max(magnitudeOfExample)]);
  limitsX = xlim; 
  limitsY = ylim;

  title('C. Freq.-domain signal', ...
        'HorizontalAlignment', 'left', 'position', [limitsX(1), limitsY(2)]);  


subplot('Position',reshape(subPlotPanel(3,1,:),1,4));  
  sinA = 0.8.*sin(timeVec.*((90/32)*2*pi));
  sinB = 1.6.*sin(timeVec.*((90/32)*2*pi)+pi/3 );
  plot(timeVec,sinA,'-','Color',[1,1,1].*0.5,'LineWidth',1);
  hold on;
  plot(timeVec,sinB,'-','Color',[1,1,1].*0,'LineWidth',1);
  hold on;

  text( 0.2,1.6,'Gain: 2','FontSize',6);
  hold on;
  text( 0.35,-1.6,'Phase: $\pi/3$','FontSize',6);
  hold on;
  
  box off
  yticks(round([-1.6,-0.8,0,0.8,1.6],2));

  xlim([0,samples/sampleFrequency]);
  xticks([0:0.1:round(samples/sampleFrequency)]);
  limitsX = xlim; 
  limitsY = ylim;
  title('D. Gain and phase between two sinusoids', ...
        'HorizontalAlignment', 'left', 'position', [limitsX(1), limitsY(2)]);  


subplot('Position',reshape(subPlotPanel(4,1,:),1,4));    
%subplot(2,3,4); 
  plot(freqHz(idxBW,1),...
       gain(idxBW,1),'-','Color',[.5,0.5,1],'LineWidth',1);
  hold on;
  plot(max(freqHz(idxBW,1)),...
       max(gain(idxBW,1)),'.','Color',[.5,0.5,1]);
  hold on;
  text(max(freqHz(idxBW,1)),...
       max(gain(idxBW,1)), ...
       sprintf('%1.1f',max(gain(idxBW,1))),...
       'HorizontalAlignment','right',...
       'VerticalAlignment','top');
  hold on;
  xlim([0,bandwidth]);
  xticks([0:5:bandwidth]);
  xlabel('Frequency (Hz)');
  ylabel('Gain (N/mm)');

  limitsX = xlim; 
  limitsY = ylim;
  title('E. Gain response (frequency-domain)', ...
        'HorizontalAlignment', 'left', 'position', [limitsX(1), limitsY(2)]);  
  
  box off;
  %xticklabels({'','','','','','','','','',''});
  yticks(round([0,k],1));
  ylim([0,1.01*max(gain(idxBW,1))]);

subplot('Position',reshape(subPlotPanel(5,1,:),1,4));    
%subplot(2,3,5); 
  plot(freqHz(idxBW,1),...
       phase(idxBW,1).*(180/pi),'-','Color',[.5,0.5,1],'LineWidth',1);
  xlabel('Frequency (Hz)');
  ylabel('Phase (degrees)');
  xlim([0,bandwidth]);
  xticks([0:5:bandwidth]);
  yticks(round([0,max(phase(idxBW,1).*(180/pi))],0));
  ylim([0,1.01*round(max(phase(idxBW,1).*(180/pi)),0)]);
  box off;

  limitsX = xlim; 
  limitsY = ylim;
  title('F. Phase response (frequency-domain)', ...
        'HorizontalAlignment', 'left', 'position', [limitsX(1), limitsY(2)]);  





figure(figPub);  
set(figPub,'Units','centimeters',...
'PaperUnits','centimeters',...
'PaperSize',[pageWidth pageHeight],...
'PaperPositionMode','manual',...
'PaperPosition',[0 0 pageWidth pageHeight]);     
%set(findall(figList(i).h,'-property','FontSize'),'FontSize',10);     
set(figPub,'renderer','painters');     
print('-dpdf', ...
    fullfile(projectFolders.output_plots_SystemIdentificationExample,...
               'fig_Pub_SystemIdentificationExample.pdf')); 

  