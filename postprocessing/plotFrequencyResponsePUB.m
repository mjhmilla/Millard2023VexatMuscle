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

function [success] = plotFrequencyResponsePUB(...
                      dataInputFolder,...
                      freqSeriesFiles,...
                      freqSeriesName,...
                      freqSeriesColor,...
                      inputFunctions,...                      
                      normFiberLength,...
                      nominalForce,...
                      nominalForceTargetIndex,...
                      dataKBR1994Fig3Gain,...
                      dataKBR1994Fig3Phase,...                      
                      flag_useElasticTendon,...
                      plotNameEnding,...
                      plotLayoutSettings,...
                      pubOutputFolder)
                      
success = 0;


numberOfHorizontalPlotColumns = plotLayoutSettings.numberOfHorizontalPlotColumns;
numberOfVerticalPlotRows      = plotLayoutSettings.numberOfVerticalPlotRows;
flag_fixedPlotWidth           = plotLayoutSettings.flag_fixedPlotWidth;
plotWidth                     = plotLayoutSettings.plotWidth;
plotHeight                    = [];

flag_usingOctave              = plotLayoutSettings.flag_usingOctave;
plotConfig;          
                    
samplePoints    = inputFunctions.samples;  
paddingPoints   = inputFunctions.padding;
sampleFrequency = inputFunctions.sampleFrequency;


  fig_freqResponse = figure;

  expAmpPlot   = [1.6,1.6];
  expBWPlot    = [90 ,15];

  expModelPlotColor     = [0.25,  0.25, 1;...
                            0.5,   0.5, 1];

  expPlotColor     = [0.25,  0.25, 0.25;...
                      0.5, 0.5, 0.5];
  expLineWidth = [1,1];
  expMarkType  = {'-','-'};
  expMarkFaceColor = expPlotColor;
  expLegendEntry = {'Kirsch: 1.6mm 90Hz', 'Kirsch: 1.6mm 15Hz'};

  chunkDuration = 0.2;        
  timeChunkStart = round((paddingPoints)*0.5)/sampleFrequency;  
  timeChunkEnd   = paddingPoints/sampleFrequency+chunkDuration;
  idxChunk       = [round((paddingPoints)*0.5):1: ...
                    (paddingPoints + sampleFrequency*chunkDuration)];                  
  timeTicks = [round(paddingPoints/sampleFrequency,2),...
               round(timeChunkEnd,2)];

  subPlotList = zeros(6,4);

  ySpace = 0.05;
  xSpace = 0.075;
  
  %%
  % Gain
  %%
  figure(fig_freqResponse);
    currentSubPlot  = subPlotRectangleMedium;
    currentSubPlot(1,1) = currentSubPlot(1,1)*0.5; 
    subPlotHeight   = subPlotRectangleMedium(1,4);
    subPlotWidth    = subPlotRectangleMedium(1,3);
    subPlotOffsetY  = subPlotHeight+ySpace;
    subPlotOffsetX  = subPlotWidth+xSpace;

  %%
  % Model vs spring-damper of best fit: 1.6 mm 90 Hz
  %%    
    
  subPlotList(1,:) = currentSubPlot(1,:);

  %%
  % Hill vs spring-damper of best fit: 1.6 mm 90 Hz
  %%    
  
  currentSubPlot(1,1) = currentSubPlot(1,1)+subPlotOffsetX;
  subPlotList(2,:)    = currentSubPlot(1,:);  

  %%
  % Phase
  %%
  currentSubPlot(1,1) = subPlotRectangle(1,1)*0.5;
  currentSubPlot(1,3) = subPlotRectangle(1,3);
  currentSubPlot(1,4) = subPlotRectangle(1,4);

  subPlotHeight  = subPlotRectangle(1,4);
  subPlotOffsetY = subPlotHeight + ySpace;

  %%
  % Gain
  %%  
  currentSubPlot(1,2) = currentSubPlot(1,2)-subPlotOffsetY-ySpace*0.5;
  subPlotList(3,:)    = currentSubPlot(1,:);  

  %%
  % Phase
  %%
  currentSubPlot(1,2) = currentSubPlot(1,2)-subPlotOffsetY;%-0.25*ySpace;
  subPlotList(4,:) = currentSubPlot(1,:);

  subPlotList(5,:) = subPlotList(3,:);
  subPlotList(5,1) = subPlotList(5,1) + subPlotOffsetX;
  
  subPlotList(6,:) = subPlotList(4,:);
  subPlotList(6,1) = subPlotList(6,1) + subPlotOffsetX;

  %%
  %Plot blank axis to check the layout
  %%
  for i=1:1:size(subPlotList,1)
    subplot('Position', [ subPlotList(i,1),...
                          subPlotList(i,2),...
                          subPlotList(i,3),...
                          subPlotList(i,4)]);
    box off;
    set(gca,'color','none')
    
  end
  

  
  %%
  %Plot model data
  %%
  modelAmpPlot =[1.6];
  modelBWPlot  =[ 90];
  modelForce   =[1,1,1].*nominalForce(1,nominalForceTargetIndex);
  modelNormFiberLength = [1,1,1].*normFiberLength;

  colorPlot = [0.,0.,1; 0.25,0.25,1; 0.5,0.5,1];
  lineWidth = [1,1,1].*0.5;

  %Here 'Full' means a 'Full' muscle model: Hill or otherwise  
  %      KD    means the spring of best fit.
  
  lineColorKD   = zeros(length(modelAmpPlot),3);
  lineColorFull = zeros(length(modelAmpPlot),3); 
  lineWidthFull = [ones(length(modelAmpPlot),1).*1 ...
                   ones(length(modelAmpPlot),1).*0.5];
  
  plotForceMax = 17.2;  
  vafInput     = zeros(length(freqSeriesFiles),length(modelAmpPlot));
    
  lineHandlesModelGain = [];
  lineHandlesModelPhase = [];  
  
  for z=1:1:length(freqSeriesFiles)

    modelColor = freqSeriesColor(z,:);
    tmp = load([dataInputFolder, freqSeriesFiles{1,z}]);
    freqSimData = tmp.freqSimData;

    fprintf('%s VAF: %1.1f pm %1.1f, [%1.1f, %1.1f]\n',...
             freqSeriesName{1,z}, ...
             mean(freqSimData.vafTime)*100,...
             std(freqSimData.vafTime)*100,...
             min(freqSimData.vafTime)*100,...
             max(freqSimData.vafTime)*100);
    
    flag_Hill = 0;
    if(isempty(strfind(freqSeriesFiles{1,z},'Hill'))==0)
      flag_Hill = 1;
    end
    
    modelLineWidth = expLineWidth;

    for k=1:1:length(modelAmpPlot)

      k01                = 0.75*(k-1)/(length(modelAmpPlot));
      lineColorFull(k,:) = modelColor;
      lineColorKD(k,:)   = [0.75,0.75,0.75];
      lineWidthKD        = 2;
      
      idxSim = 0;
      tol = 1e-6;
      for m=1:1:size(freqSimData.force,2)     
        if( abs(freqSimData.amplitudeMM(1,m)     - modelAmpPlot(1,k)        ) <= tol && ...
            abs(freqSimData.bandwidthHz(1,m)     - modelBWPlot(1,k)         ) <= tol && ...
            abs(freqSimData.nominalForceDesired(1,m)    - modelForce(1,k)          ) <= tol && ...
            abs(freqSimData.normFiberLength(1,m) - modelNormFiberLength(1,k)) <= tol)

          if(idxSim == 0)
            idxSim = m;
          else
            assert(0); %Error condition: there should not be 2 simulations with 
                       %the same configuration
          end
        end
      end
      
      idxWave = getSignalIndex(modelAmpPlot(1,k),...
                               modelBWPlot(1,k),...
                               inputFunctions);                            
                             
      idxForce = 1;
      idxGain  = 3;
      idxPhase = 4;
      if(flag_Hill)
        idxForce = 2;
        idxGain  = 5;
        idxPhase = 6;
      end
      
      subplot('Position', [ subPlotList(idxForce,1),...
                            subPlotList(idxForce,2),...
                            subPlotList(idxForce,3),...
                            subPlotList(idxForce,4)]);                          

      yo = freqSimData.nominalForce(1,idxSim);
     
      vafInput(z,k) = freqSimData.vafTime(1,idxSim);
      
      spring = freqSimData.stiffness(1,idxSim)./1000;
      damper = freqSimData.damping(1,idxSim)./1000;
      kLabel = sprintf('%1.1f',spring);
      dLabel = sprintf('%1.3f',damper);
      
      kdLabel = ['K: ',kLabel,'N/mm, $$\beta$$: ',dLabel,'N/(mm/s)',...
                 ' VAF ', num2str(round(freqSimData.vafTime(1,idxSim)*100)),'\%'];
      
          
      lineLabel = [ freqSeriesName{z}];

      pidKD = plot( inputFunctions.time(idxChunk,1),...
             freqSimData.forceKD(idxChunk,idxSim)+yo,...
             'Color',lineColorKD(k,:),...
             'LineWidth',lineWidthKD);
      hold on;
      
      
      pidMdl = plot( inputFunctions.time(idxChunk,1),...
            freqSimData.force(idxChunk,idxSim),...
            'Color',lineColorFull(k,:),...
            'LineWidth',lineWidthFull(k,z));
      hold on;  
      
      

                      
      box off;
      set(gca,'color','none')
      
      lh = legend([pidMdl,pidKD],lineLabel,kdLabel,'Location','NorthWest');
      %lh.NumColumns=2;
      lh.Position(1,1) = lh.Position(1,1)-0.04;
      lh.Position(1,2) = lh.Position(1,2)+0.055;
      legend boxoff;
      
      %legend boxoff;
      
      xlim([inputFunctions.time(idxChunk(1),1  )-0.01,...
            inputFunctions.time(idxChunk(end),1)+0.01]);
      ylim([0,plotForceMax]);
      xticks(timeTicks);
      
      f0 = freqSimData.force(idxChunk(1),idxSim);
      fmax=max(freqSimData.force(idxChunk,idxSim));
      
      yticks([0,round(f0,1),round(fmax,1)]);
      set(gca,'color','none')

      tmin = min(inputFunctions.time(idxChunk,1));
      tmax = max(inputFunctions.time(idxChunk,1));

      if(flag_Hill==0)
        %ta0 = text(tmin-0.15*(tmax-tmin), 1.45*plotForceMax,...
        %     'A. Example time-domain force response to the 1.6 mm 90 Hz perturbation',...
        %     'FontSize',11);

        ta0 = text(tmin-0.15*(tmax-tmin), 1.45*plotForceMax,...
             'A.','FontSize',11);        
        
      else
        tb0 = text(tmin-0.15*(tmax-tmin), 1.45*plotForceMax,...
             'B.','FontSize',11);                
      end
             
      
      ylabel('Force (N)');
      xlabel('Time (s)');


    end


    
%%
% Plot the experimental data
%%
    for k=1:1:length(expAmpPlot)

      idxSim = 0;
      tol = 1e-6;
      for m=1:1:size(freqSimData.force,2)     
        if( abs(freqSimData.amplitudeMM(1,m) - expAmpPlot(1,k)) <= tol && ...
            abs(freqSimData.bandwidthHz(1,m) - expBWPlot(1,k))  <= tol && ...
            abs(freqSimData.nominalForceDesired(1,m)- modelForce(1,k)) <= tol && ...
            abs(freqSimData.normFiberLength(1,m)-modelNormFiberLength(1,k)) <= tol)

          if(idxSim == 0)
            idxSim = m;
          else
            assert(0); %Error condition: there should not be 2 simulations with 
                       %the same configuration
          end
        end
      end      

      fprintf('%i. K: %1.3f D: %1.3f\n',idxSim,...
                freqSimData.stiffness(1,idxSim)./1000,...
                freqSimData.damping(1,idxSim)./1000);
      
      %
      idxRange = [freqSimData.idxFreqRange(1,idxSim):1: ...
                  freqSimData.idxFreqRange(2,idxSim)];

      idxCutoff = freqSimData.idxFreqRange(1,idxSim);
                
      
      subplot('Position', subPlotList(idxGain,:));

      markType = '.';
      markFaceColor = freqSeriesColor(z,:);
      markSize = 2;
      lineWidth = 0.5;
      kdMarkType = '-';
      kdLineWidth= 1;
      kdWhiteLineWidth = 2;
      kdLineColor = freqSeriesColor(z,:);
      if(k==2)        
        markType = 'o';
        markFaceColor = [1,1,1];
        markSize = 3;
        lineWidth=0.1;
        kdMarkType = '-';
        kdLineWidth= 1;
        kdWhiteLineWidth = 3;
        kdLineColor = freqSeriesColor(z,:).*0.5 + [1,1,1].*0.5;
        
      end

      trialName = sprintf('%1.1fmm %1.0fHz',...
                          freqSimData.amplitudeMM(1,idxSim),...
                          freqSimData.bandwidthHz(1,idxSim));

      plot(freqSimData.freqHz(idxRange,1), ...
           freqSimData.gain(idxRange,idxSim)./1000,...
           markType,'Color',kdLineColor, ...
           'MarkerFaceColor',markFaceColor,'LineWidth',lineWidth,...
           'MarkerSize',markSize,'DisplayName',[freqSeriesName{z},trialName]);  
      hold on;

      pidW=plot(freqSimData.freqHz(idxRange,1), ...
           freqSimData.gainKD(idxRange,idxSim)./1000,...
           '-','Color',[1,1,1], ...
           'LineWidth',kdWhiteLineWidth);  
      hold on;
      set(get(get(pidW,'Annotation'),...
                'LegendInformation'),...
                'IconDisplayStyle','off');      
              
      pidL=plot(freqSimData.freqHz(idxRange,1), ...
           freqSimData.gainKD(idxRange,idxSim)./1000,...
           kdMarkType,'Color',kdLineColor, ...
           'LineWidth',kdLineWidth,'DisplayName', ['K-$$\beta$$: ',trialName]);
      hold on;
      
      lineHandlesModelGain = [lineHandlesModelGain,...
                              pidW,...
                              pidL];


      
      fmin = 0.;
      fmax = 90.01;%max(freqSimData.freqHz(idxRange,1));
      gmin = 0.;%min(freqSimData.gain(idxRange,idxSim)./1000);
      gmax = 8.01;
      
      if(flag_Hill==0 && k == 1)
        tc = text(fmin-0.1*(fmax-fmin), 1.1*gmax,...
               'C. ',...
               'FontSize',11);   
        hold on;
      elseif(flag_Hill==1 && k == 1)
        td = text(fmin-0.1*(fmax-fmin), 1.1*gmax,...
               'D. ',...
               'FontSize',11);   
        hold on;
        
      end
      
      box off;
      set(gca,'color','none')        
      
      
      hold on;
      
      
      ylim([gmin,gmax]);
      xlim([fmin,fmax]);
      xticks([0,15,90])
      yticks([0,2,4,6,8]);
     
            
      %tc0 = text(fmin-0.15*(fmax-fmin), 1.45*(gmax),...
      %       'C.','FontSize',11);       
      
      subplot('Position', subPlotList(idxPhase,:));
      
      plot(freqSimData.freqHz(idxRange,1), ...
           freqSimData.phase(idxRange,idxSim).*(180/pi),...
           markType,'Color',kdLineColor, ...
           'MarkerFaceColor',markFaceColor,'LineWidth',0.5,...
           'MarkerSize',markSize,...
           'DisplayName',[freqSeriesName{z},trialName]);   
      hold on;
      
      pidW=plot(freqSimData.freqHz(idxRange,1), ...
           freqSimData.phaseKD(idxRange,idxSim).*(180/pi),...
           '-','Color',[1,1,1], ...
           'LineWidth',kdWhiteLineWidth);  
      hold on;
      set(get(get(pidW,'Annotation'),...
                'LegendInformation'),...
                'IconDisplayStyle','off');      
      
      pidL=plot(freqSimData.freqHz(idxRange,1), ...
           freqSimData.phaseKD(idxRange,idxSim).*(180/pi),...
           kdMarkType,'Color',kdLineColor, ...
           'LineWidth',kdLineWidth,'DisplayName',['K-$$\beta$$: ',trialName]);  
      hold on;      
      
      lineHandlesModelPhase = [lineHandlesModelPhase,...
                              pidW,...
                              pidL];
      
      
      box off;
      set(gca,'color','none')
      
      pmin = 0;
      pmax = 90.01;
      
      ylim([pmin,pmax]);
      xlim([fmin,fmax]);
      xticks([0,15,90])
      yticks([0,45,90]);

      if(flag_Hill==0 && k == 1)
        tb = text(fmin-0.1*(fmax-fmin), 1.1*pmax,...
               'E. ',...
               'FontSize',11);   
        hold on;
      elseif(flag_Hill==1 && k == 1)
        tb = text(fmin-0.1*(fmax-fmin), 1.1*pmax,...
               'F. ',...
               'FontSize',11);   
        hold on;
        
      end
      
      
    end

  end
  
  
  for z=1:1:length(freqSeriesFiles)
    
    flag_Hill = 0;
    if(isempty(strfind(freqSeriesFiles{1,z},'Hill'))==0)
      flag_Hill = 1;
    end    
    
    idxForce = 1;
    idxGain  = 3;
    idxPhase = 4;
    if(flag_Hill)
      idxForce = 2;
      idxGain  = 5;
      idxPhase = 6;
    end

  
    subplot('Position', subPlotList(idxGain,:));

    for k=1:1:2
      idxMin = find(dataKBR1994Fig3Phase(k).x >= 4, 1 );
      pid = plot( dataKBR1994Fig3Gain(k).x(idxMin:1:end),...
            dataKBR1994Fig3Gain(k).y(idxMin:1:end),...
            '-', 'Color', [1,1,1],...
            'LineWidth', expLineWidth(1,k)*2,...
            'DisplayName','');
      hold on;

      set(get(get(pid,'Annotation'),...
                  'LegendInformation'),...
                  'IconDisplayStyle','off');

      plot( dataKBR1994Fig3Gain(k).x(idxMin:1:end),...
            dataKBR1994Fig3Gain(k).y(idxMin:1:end),...
            expMarkType{1,k}, 'Color', expPlotColor(k,:),...
            'LineWidth',expLineWidth(1,k),...
            'DisplayName', expLegendEntry{k});

      hold on;
    end

    box off;
    set(gca,'color','none')  
    xlabel(dataKBR1994Fig3Gain(1).xName);
    ylabel('Gain (N/mm)');

    %if(flag_useElasticTendon==0 )
    %  lh = legend('Location','SouthEast');
    %  lh.Position(1,1) = lh.Position(1,1)+0.03;
    %  lh.Position(1,2) = lh.Position(1,2)-0.03;
    %  legend boxoff;
    %end

    subplot('Position', subPlotList(idxPhase,:));

    for k=1:1:2
      
      idxMin = find(dataKBR1994Fig3Phase(k).x >= 4, 1 );
      
      pid = plot( dataKBR1994Fig3Phase(k).x(idxMin:1:end),...
            dataKBR1994Fig3Phase(k).y(idxMin:1:end),...
            '-', 'Color', [1,1,1],...
            'LineWidth',expLineWidth(1,k)*2); 
      hold on;

      set(get(get(pid,'Annotation'),...
                  'LegendInformation'),...
                  'IconDisplayStyle','off');

      plot( dataKBR1994Fig3Phase(k).x(idxMin:1:end),...
            dataKBR1994Fig3Phase(k).y(idxMin:1:end),...
            expMarkType{1,k}, 'Color', expPlotColor(k,:),...
            'LineWidth',expLineWidth(1,k),...
            'DisplayName', expLegendEntry{k}); 
      hold on;    
    end


    xlabel(dataKBR1994Fig3Phase(1).xName);
    ylabel(dataKBR1994Fig3Phase(1).yName);
    box off;
    set(gca,'color','none')

    
    %if(flag_useElasticTendon==1 )
      if(z==1)
        lh = legend('Location','SouthWest');        
        lh.Position(1,1) = lh.Position(1,1) -0.09;    
        lh.Position(1,2) = lh.Position(1,2) -0.06;           
      else
        lh = legend('Location','NorthEast');        
        lh.Position(1,1) = lh.Position(1,1) + 0.07;    
        lh.Position(1,2) = subPlotList(idxGain,2)+0.06;                   
      end
      lh.NumColumns=3;      
      lh.NumColumnsMode='manual';    
%       llenMax = 0;      
%       for idxIcon=1:1:length(icons)
%         if( contains(class(icons(idxIcon)),'Line'))
%           if(length(icons(idxIcon).XData) == 2)
%             xd = icons(idxIcon).XData; 
%             llen = diff(xd);
%             if(llen > llenMax)
%               llenMax = llen;
%             end
%           end
%         end        
%       end
%       for idxIcon=1:1:length(icons)
%         if( contains(class(icons(idxIcon)),'Line'))
%           xd = icons(idxIcon).XData;
%           if(length(xd)==2)
%             icons(idxIcon).XData = [xd(1,1), xd(1,1)+(0.5*llenMax)];            
%           else
%             icons(idxIcon).XData = [xd(1,1)-(0.25*llenMax)];                        
%           end
%         end 
%         if( contains(class(icons(idxIcon)),'Text'))
%           icons(idxIcon).Position(1,1) = icons(idxIcon).Position(1,1)-0.5*llenMax;
%         end
%       end
      
      

      legend boxoff;
    %end
%     if(flag_useElasticTendon==1)
%       lh = legend('Location','NorthEast');  
%       lh.Position(1,1) = lh.Position(1,1)+0.035;    
%       lh.Position(1,2) = lh.Position(1,2)+0.09;
%       lh.EdgeColor = lh.Color;
% 
%     end
  end

% for i=1:1:length(lineHandlesModelPhase)
%   uistack(lineHandlesModelPhase(i),'top');
% end
% for i=1:1:length(lineHandlesModelGain)
%   uistack(lineHandlesModelGain(i),'top');
% end
%   
  
set(fig_freqResponse,'Units','centimeters',...
'PaperUnits','centimeters',...
'PaperSize',[pageWidth pageHeight],...
'PaperPositionMode','manual',...
'PaperPosition',[0 0 pageWidth pageHeight]);     
%set(findall(figList(i).h,'-property','FontSize'),'FontSize',10);     
set(fig_freqResponse,'renderer','painters');     
set(gcf,'InvertHardCopy','off')

tendonTag = '_ElasticTendon';
if(flag_useElasticTendon==0)
  tendonTag = '_RigidTendon';
end

print('-dpdf', [pubOutputFolder,'fig_Pub_ModelFrequencyResponse',...
                                tendonTag,'_',plotNameEnding,'.pdf']); 


success = 1;