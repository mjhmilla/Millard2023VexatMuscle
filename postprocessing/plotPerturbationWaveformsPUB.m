function [success] = plotPerturbationWaveformsPUB(inputFunctions,...
  workingFolder,...
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




%%
% Generate pub-quality plots
%%

  fig_input = figure;
  
  ampPlot   = [1.6,1.6,1.6];
  bwPlot    = [90,35,15];
      
  idxPlot   = zeros(size(ampPlot));
  colorPlot = [0,0,0; 0.25,0.25,0.25; 0.5,0.5,0.5];
  lineWidth = [1,1,1];

  %%
  %Time-domain plots
  %%
  
  idxChunk = [];
  offset = 0;
  
  currentSubPlot = subPlotRectangleSmall;
  subPlotHeight  = subPlotRectangleSmall(1,4);
  subPlotOffsetY = subPlotHeight*1.1;
  
  for i=1:1:length(ampPlot)
    idxPlot(1,i) = getSignalIndex(ampPlot(1,i),bwPlot(1,i),inputFunctions);
    idx = idxPlot(1,i);

    chunkDuration = 0.2;  
        
    timePadding    = paddingPoints/sampleFrequency;
    timeChunkStart = (paddingPoints*0.9)/sampleFrequency;
    timeChunkEnd   = (timePadding+chunkDuration);    
    idxChunk       = [round(paddingPoints*0.9):1:...
                     (paddingPoints+sampleFrequency*chunkDuration)];    
    timeTicks = [round(timeChunkStart,2), round(timeChunkEnd,2)];
    
    
    %Plot the time domain signal
    figure(fig_input);
    subplot('Position', [ currentSubPlot(1,1),...
                          currentSubPlot(1,2),...
                          currentSubPlot(1,3),...
                          currentSubPlot(1,4)]);
                        
    if(i < length(ampPlot))
      currentSubPlot(1,2) = currentSubPlot(1,2)-subPlotOffsetY;
    end
    
    plot(inputFunctions.time(idxChunk,1),...
         inputFunctions.x(idxChunk,idx).*1000,...
        'Color',colorPlot(i,:),'LineWidth',lineWidth(1,i));    
      
    hold on;

    if(i==1)
      x0    = min(inputFunctions.time(idxChunk,1));
      xSpan = max(inputFunctions.time(idxChunk,1))-min(inputFunctions.time(idxChunk,1));
      
      y0    = max(inputFunctions.amplitudeMM)*1.2;
      ySpan = y0*2;
      
      text( x0-0.1*xSpan, (y0 + ySpan*(0.5/(currentSubPlot(1,4)*pageHeight))),...
            'A. Example perturbation waveforms','FontSize',11);
      hold on;
    end
    
    box off;
    set(gca,'color','none')
    
    xlim([min(inputFunctions.time(idxChunk,1)),max(inputFunctions.time(idxChunk,1))*1.01])
    ylim([-max(inputFunctions.amplitudeMM),max(inputFunctions.amplitudeMM)].*(1.2));   
       
       
    if(i==2)
      ylabel('Length (mm)');            
    end
    
    if(flag_usingOctave == 0)       
      yticks([-max(ampPlot),0,max(ampPlot)]);
      xticks(timeTicks);
    
      hL = legend( inputFunctions.labels(idx,:) ,'Location','NorthEast');
      pos = get(hL,'Position');
      pos(1,2) = pos(1,2)+0.0275;
      set(hL,'Position',pos);
      legend boxoff;
    end    
    %title(inputFunctions.labels(idx,:));
    
    
    if(i==length(ampPlot))
      xlabel('Time (s)');      
    else
      timeTickLabels = [];
      for z=1:1:length(timeTicks)
        timeTickLabels = [num2str(timeTickLabels),''];
      end
      if(flag_usingOctave == 0)
        xticklabels(timeTickLabels);
        set(gca,'xcolor','None');              
      end
    end
    
  end
  
  %%
  %Frequency domain plots
  %%
  subPlotHeight = currentSubPlot(1,4)*3;
  subPlotOffsetY= subPlotHeight+0.05;
  currentSubPlot(1,4) = subPlotHeight;
  currentSubPlot(1,2) = currentSubPlot(1,2) - subPlotOffsetY;
  for i=1:1:length(ampPlot)
    idxPlot(1,i) = getSignalIndex(ampPlot(1,i),bwPlot(1,i),inputFunctions);
    idx = idxPlot(1,i);    
    
    figure(fig_input);
    subplot('Position', [ currentSubPlot(1,1),...
                          currentSubPlot(1,2),...
                          currentSubPlot(1,3),...
                          currentSubPlot(1,4)]);
    
    plot( inputFunctions.freqHz(1:(1+samplePoints/2),:),...
          inputFunctions.p(:,idx)./max(inputFunctions.p(:,idx)),...
          'Color',colorPlot(i,:),'LineWidth',lineWidth(1,i)*2);    
          
    hold on;
    
    if(i==length(ampPlot))
      
      xlim([(-1),(max(bwPlot)+10)]);
      ylim([-0.01, 1.01]);   
       
      if(flag_usingOctave==0)
        yticks([0,1/sqrt(2),1]);
        yticklabels({'0','$$\frac{1}{\sqrt{2}}$$','1'});
        xticks([0, sort(bwPlot)]);      
      end      
      xlabel('Frequency (Hz)')
      ylabel('Norm. Power');
      
      box off;
      set(gca,'color','none')
      
      x0    = 0;
      xSpan = 100;
      
      y0    = 1;
      ySpan = 1;
      
      text( x0-0.1*xSpan, ...
            (y0+ySpan*(0.5/(currentSubPlot(1,4)*pageHeight))),...
            'B. Example perturbation power spectrum','FontSize',11);
      hold on;      
      
    end
  end
  
  

   set(fig_input,'Units','centimeters',...
   'PaperUnits','centimeters',...
   'PaperSize',[pageWidth pageHeight],...
   'PaperPositionMode','manual',...
   'PaperPosition',[0 0 pageWidth pageHeight]);     
   %set(findall(figList(i).h,'-property','FontSize'),'FontSize',10);     
   set(fig_input,'renderer','painters');     
   print('-dpdf', [workingFolder,'fig_Pub_StochasticInput.pdf']); 

   here=1;



success = 1;