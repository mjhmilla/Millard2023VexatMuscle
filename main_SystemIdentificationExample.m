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

addpath( genpath(projectFolders.postprocessing) );
%%
%Publication configuration
%%
flag_usingOctave=0;
plotWidth = 7.5;
plotHeight= 2.0;
plotHorizMarginCm=1.5;
plotVertMarginCm = 1.5;
numberOfVerticalPlotRows = 4;
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
paddingWidth      = 100;
perturbationLength= 1.6;%1.6 mm in units of m
bandwidth         = 35;    
sampleFrequency   = 1/0.002;
nyquistFrequency  = sampleFrequency*0.5;



[b,a]             =butter(2,bandwidth/nyquistFrequency);
xRandom           = rand(samples,1);
xRandom = xRandom-mean(xRandom);
xRandom(1:paddingWidth)=0;
xRandom((samples-paddingWidth):1:samples) = 0;

%xTimeDomain is the length change: the nominal length has been removed
xTimeDomain       = filter(b,a,xRandom);
xTimeDomainScale  = perturbationLength/ (max(abs(xTimeDomain)));
xTimeDomain = xTimeDomain.*xTimeDomainScale;
timeVec           = [0:(1/(samples-1)):1].*(samples/sampleFrequency);

%%
%Make some synthetic force data: 
%  This entire section of code is not necessary once we have experimental
%  data
%%
k=4.46;  %Stiffness: 4.46 N/mm
d=0.0089;%Damping  : feel free to adjust

%Form the perfect derivative of x
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
yTimeDomain = k.*xTimeDomain - d.*xDotTimeDomain;

%%
%Use these two signals to evaluate the gain and phase shift between
%the length and force perturbation
%%

frequencyConvHz = [0:1/(samples*2):( 1-(1/(samples*2)) )]' .* (sampleFrequency);
frequencyConvRadians = frequencyConvHz.*(2*pi);

%We only analyze the frequency spectrum that is within the bandwith of the
%perturbation signal: everything else will be too weak to produce coherent
%data
idxFreqInBandwidth = find(frequencyConvHz <= bandwidth);

xyTimeDomain = conv(xTimeDomain,yTimeDomain);
xxTimeDomain = conv(xTimeDomain,xTimeDomain);
yyTimeDomain = conv(yTimeDomain,yTimeDomain);

xyFreqDomain = fft(xyTimeDomain);
xxFreqDomain = fft(xxTimeDomain);
yyFreqDomain = fft(yyTimeDomain);

xFreqDomain = fft(xTimeDomain);
yFreqDomain = fft(yTimeDomain);

gain = abs(xyFreqDomain./xxFreqDomain);
phase=-angle(xyFreqDomain./xxFreqDomain);
coherenceSq = (xyFreqDomain.*xyFreqDomain) ...
  ./ (xxFreqDomain.*yyFreqDomain);


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

fig=figure;
ax=axes(fig);
set(ax,'TitleHorizontalAlignment','Left');
subplot('Position',reshape(subPlotPanel(1,1,:),1,4));
%subplot(2,3,1); 
  plot(timeVec, xTimeDomain,'-','Color',[1,1,1].*0.5,'LineWidth',1);
  xlabel('Time (s)');
  ylabel('Length(mm)');
  box off;
  xlim([0,samples/sampleFrequency]);
  xticks([0:0.1:round(samples/sampleFrequency)]);
  xticklabels({'','','','','','','','','',''});
  yticks([-1,0,1].*round(max(abs(xTimeDomain)),2));
  ylim([-1.01,1.01].*round(max(abs(xTimeDomain)),1))
%  ylim([-1.01,1.01].*round(max(abs(yTimeDomain)),1))

  limitsX = xlim; 
  limitsY = ylim;
  title('A. Length perturbation (time domain)', ...
        'HorizontalAlignment', 'left', 'position', [limitsX(1), limitsY(2)]);
  
  
subplot('Position',reshape(subPlotPanel(2,1,:),1,4));  
%subplot(2,3,2); 
  plot(timeVec, yTimeDomain,'-','Color',[1,1,1].*0.,'LineWidth',1);
  xlabel('Time (s)');
  ylabel('Force (N)');
  title('B. Force response (time domain)');  
  box off;
  xlim([0,samples/sampleFrequency]);
  xticks([0:0.1:round(samples/sampleFrequency)]);

  yticks([-1,0,1].*round(max(abs(yTimeDomain)),1));
  ylim([-1.01,1.01].*round(max(abs(yTimeDomain)),1))
  limitsX = xlim; 
  limitsY = ylim;
  title('B. Force response (time domain)', ...
        'HorizontalAlignment', 'left', 'position', [limitsX(1), limitsY(2)]);  
% 
% subplot('Position',reshape(subPlotPanel(3,1,:),1,4));  
%   sinA = sin(timeVec.*((90/32)*2*pi));
%   sinB = 1.6.*sin(timeVec.*((90/32)*2*pi)+pi/4 );
%   plot(timeVec,sinA,'-','Color',[1,1,1].*0.5,'LineWidth',1);
%   hold on;
%   plot(timeVec,sinB,'-','Color',[1,1,1].*0,'LineWidth',1);
%   hold on;
%   box off
%   yticks(round([-1.6,-1,0,1,1.6],2));
% 
%   xlim([0,samples/sampleFrequency]);
%   xticks([0:0.1:round(samples/sampleFrequency)]);
% 
%   limitsX = xlim; 
%   limitsY = ylim;
%   title('C. Gain and phase between two sinusoids', ...
%         'HorizontalAlignment', 'left', 'position', [limitsX(1), limitsY(2)]);  
% 
% subplot('Position',reshape(subPlotPanel(4,1,:),1,4));  
% 
% 
%   limitsX = xlim; 
%   limitsY = ylim;
%   title('D. Frequency domain transformation', ...
%         'HorizontalAlignment', 'left', 'position', [limitsX(1), limitsY(2)]);  



subplot('Position',reshape(subPlotPanel(3,1,:),1,4));    
%subplot(2,3,4); 
  plot(frequencyConvHz(idxFreqInBandwidth,1),...
       gain(idxFreqInBandwidth,1),'-','Color',[0,0,0],'LineWidth',1);
  xlabel('Frequency (Hz)');
  ylabel('Gain (N/mm)');

  limitsX = xlim; 
  limitsY = ylim;
  title('E. Gain response (frequency domain)', ...
        'HorizontalAlignment', 'left', 'position', [limitsX(1), limitsY(2)]);  
  
  box off;
  xlim([0,bandwidth]);  
  xticks([0:10:bandwidth]);
  xticklabels({'','','','','','','','','',''});
  yticks(round([0,k,max(gain(idxFreqInBandwidth,1))],1));
  ylim([0,1.01*max(gain(idxFreqInBandwidth,1))]);

subplot('Position',reshape(subPlotPanel(4,1,:),1,4));    
%subplot(2,3,5); 
  plot(frequencyConvHz(idxFreqInBandwidth,1),...
       phase(idxFreqInBandwidth,1).*(180/pi),'-','Color',[0,0,0],'LineWidth',1);
  xlabel('Frequency (Hz)');
  ylabel('Phase (degrees)');
  xlim([0,bandwidth]);
  xticks([0:10:bandwidth]);
  yticks(round([0,max(phase(idxFreqInBandwidth,1).*(180/pi))],0));
  ylim([0,1.01*round(max(phase(idxFreqInBandwidth,1).*(180/pi)),0)]);
  box off;

  limitsX = xlim; 
  limitsY = ylim;
  title('F. Phase response (frequency domain)', ...
        'HorizontalAlignment', 'left', 'position', [limitsX(1), limitsY(2)]);  





  
set(fig,'Units','centimeters',...
'PaperUnits','centimeters',...
'PaperSize',[pageWidth pageHeight],...
'PaperPositionMode','manual',...
'PaperPosition',[0 0 pageWidth pageHeight]);     
%set(findall(figList(i).h,'-property','FontSize'),'FontSize',10);     
set(fig,'renderer','painters');     
print('-dpdf', ...
    fullfile(projectFolders.output_plots_SystemIdentificationExample,...
               'fig_Pub_SystemIdentificationExample.pdf')); 

  