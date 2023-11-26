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

function [success] = plotStiffnessDampingPUB(...
                      dataInputFolder,...
                      freqSeriesFiles,...
                      freqSeriesName,...
                      freqSeriesColor,...
                      inputFunctions,...                      
                      normFiberLength,...
                      nominalForce,...                        
                      dataKBR1994Fig12K,...
                      dataKBR1994Fig12D,... 
                      flag_useElasticTendon,...
                      plotNameEnding,...
                      plotLayoutSettings,...
                      pubOutputFolder)
                      
success = 0;

%Writes text of the VAF near the data point
flag_addVAF=0;

numberOfHorizontalPlotColumns = plotLayoutSettings.numberOfHorizontalPlotColumns;
numberOfVerticalPlotRows      = plotLayoutSettings.numberOfVerticalPlotRows;
flag_fixedPlotWidth           = plotLayoutSettings.flag_fixedPlotWidth;
plotWidth                     = plotLayoutSettings.plotWidth;
plotHeight                    = plotWidth;
flag_usingOctave              = plotLayoutSettings.flag_usingOctave;
plotConfig;          
                    
samplePoints    = inputFunctions.samples;  
paddingPoints   = inputFunctions.padding;
sampleFrequency = inputFunctions.sampleFrequency;

fig_fig12 = figure;

  targetAmplitude   = [0.8];
  targetBandwidth    = [35];

  expPlotColor      = [0.5,  0.5, 0.5;...
                       0.5,  0.5, 0.5];
  expLineWidth      = [1,1];
  expMarkType       = {'o','o'};
  expMarkSize       = [5,5];
  expMarkFaceColor  = [0.5,  0.5, 0.5;...
                       0.75, 0.75,0.75];
  expLegendEntry    = {'Kirsch: 0.8mm 35Hz Soleus',...
                       'Kirsch: 0.8mm 35Hz MG'};
  perturbationName = ' 0.8mm 35Hz';

  subPlotList = zeros(2,4);

  ySpace = 0.075;
  xSpace = 0.075;
  
  kMin = -0.01;
  kMax = 12.01;
  dMin = -0.001;
  dMax = 0.1201;
  
  fMin = -0.01;
  fMax = 12.51;
  
  %%
  % Gain
  %%
  figure(fig_fig12);
    currentSubPlot  = subPlotSquare;
    currentSubPlot(1,1) = currentSubPlot(1,1)*0.5; 
    
    subPlotHeight   = currentSubPlot(1,4);
    subPlotWidth    = currentSubPlot(1,3);
    subPlotOffsetY  = subPlotHeight + ySpace;
    subPlotOffsetX  = subPlotWidth  + xSpace;

  %%
  % Stiffness vs. force
  %%    
    
  subPlotList(1,:) = currentSubPlot(1,:);

  %%
  % Damping vs. force
  %%    
  
  currentSubPlot(1,2) = currentSubPlot(1,2)-subPlotOffsetY;
  subPlotList(2,:)    = currentSubPlot(1,:);  

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
  %Plot the experimental data
  %%
  idxK = 1;
  idxD = 2;
  
  figure(fig_fig12);
  
  for(i=1:1:length(dataKBR1994Fig12K))
    
    subplot('Position', [ subPlotList(idxK,1),...
                          subPlotList(idxK,2),...
                          subPlotList(idxK,3),...
                          subPlotList(idxK,4)]);
    plot(dataKBR1994Fig12K(i).x,...
         dataKBR1994Fig12K(i).y,...
         expMarkType{i},...
         'Color',     expPlotColor(i,:),...
         'LineWidth', expLineWidth(1,i),...
         'MarkerFaceColor',expMarkFaceColor(i,:),...
         'MarkerSize', expMarkSize(1,i),...
         'DisplayName', expLegendEntry{i});
    hold on;   
    
    if(i== length(dataKBR1994Fig12K))
      xlabel('Force (N)');
      ylabel('Stiffness (N/mm)');   
      box off;
      set(gca,'color','none');     
    end
    
                        
    subplot('Position', [ subPlotList(idxD,1),...
                          subPlotList(idxD,2),...
                          subPlotList(idxD,3),...
                          subPlotList(idxD,4)]);

    plot(dataKBR1994Fig12D(i).x,...
         dataKBR1994Fig12D(i).y,...
         expMarkType{i},...
         'Color',     expPlotColor(i,:),...
         'LineWidth', expLineWidth(1,i),...
         'MarkerFaceColor',expMarkFaceColor(i,:),...
         'MarkerSize', expMarkSize(1,i),...
         'DisplayName', expLegendEntry{i});
    hold on;   
    
    if(i== length(dataKBR1994Fig12K))
      xlabel('Force (N)');
      ylabel('Damping (N/(mm/s))'); 
      box off;
      set(gca,'color','none');
    end
                        
  end

  %%
  %Plot the model data
  %%  
  
  for z=1:1:length(freqSeriesFiles)

    modelColor = freqSeriesColor(z,:);
    tmp = load([dataInputFolder, freqSeriesFiles{1,z}]);
    freqSimData = tmp.freqSimData;

    modelLineColor = [1,1,1];
    modelMarkType = 'o';
    modelMarkSize = 7;
    modelLineWidth= 1;
    modelFaceColor = modelColor;%[1,1,1];
    
    if(z==2)
      modelMarkSize = modelMarkSize -1;
      modelMarkType = 'd';
    end
    
    
    flag_Hill = 0;
    if(isempty(strfind(freqSeriesFiles{1,z},'Hill'))==0)
      flag_Hill = 1;
    end


    for i=1:1:length(nominalForce)


  
      idxSim = getFrequencySimulationIndex(targetAmplitude,  ...
                                           targetBandwidth, ...
                                           nominalForce(1,i),...
                                           normFiberLength,...
                                           freqSimData);          
%       idxSim = 0;
%       tol = 1e-6;
%       for m=1:1:size(freqSimData.force,2)     
%         if( abs(freqSimData.amplitudeMM(1,m)     - targetAmplitude   ) <= tol && ...
%             abs(freqSimData.bandwidthHz(1,m)     - targetBandwidth   ) <= tol && ...
%             abs(freqSimData.nominalForceDesired(1,m) - nominalForce(1,i) ) <= tol && ...
%             abs(freqSimData.normFiberLength(1,m) - normFiberLength   ) <= tol)
%           if(idxSim == 0)
%             idxSim = m;
%           else
%             assert(0); %Error condition: there should not be 2 simulations with 
%                        %the same configuration
%           end
%         end
%       end
%       
      if(contains(freqSeriesFiles{1,z},'Vexat')==1 ...
              && flag_useElasticTendon == 1 ...
              && i==length(nominalForce))
        here=1;
      end

      subplot('Position', [ subPlotList(idxK,1),...
                            subPlotList(idxK,2),...
                            subPlotList(idxK,3),...
                            subPlotList(idxK,4)]);    

      
      
      pid = plot(freqSimData.nominalForce(1,idxSim),...
           freqSimData.stiffness(1,idxSim)/1000,...
           modelMarkType, ...
           'Color', modelLineColor,...
           'LineWidth',modelLineWidth,...
           'MarkerSize',modelMarkSize,...
           'MarkerFaceColor',modelFaceColor,...
           'DisplayName',freqSeriesName{z});
      hold on;
      if(i > 1)
        set(get(get(pid,'Annotation'),...
                  'LegendInformation'),...
                  'IconDisplayStyle','off'); 
      end
      
      posText = [freqSimData.nominalForce(1,idxSim),...
                freqSimData.stiffness(1,idxSim)/1000];
      vafText = sprintf('%1.0f',freqSimData.vafTime(1,idxSim)*100); 
      
      textDeltaX  = 0;
      textDeltaY  = 0;
      textAlign   = '';      
      if(flag_useElasticTendon == 1)
        textDeltaX = 0.25;
        textDeltaY = 0;
        textAlign = 'left';
        if(flag_Hill == 1)
          textDeltaX = -0.25;
          textDeltaY = 0;
          textAlign = 'right';
        end
      else
        textDeltaX = -0.5;
        textDeltaY = 0.25;
        textAlign = 'right';
        if(flag_Hill == 1)
          textDeltaX = 0.1;
          textDeltaY = -0.35;
          textAlign = 'left';
        end
        
      end
      if(flag_addVAF==1)
          text(posText(1,1)+textDeltaX,...
               posText(1,2)+textDeltaY,[vafText,'\%'],...
                'HorizontalAlignment',textAlign);
          hold on;
      end
      xlim([fMin,fMax]);
      ylim([kMin,kMax]);      

      
      subplot('Position', [ subPlotList(idxD,1),...
                            subPlotList(idxD,2),...
                            subPlotList(idxD,3),...
                            subPlotList(idxD,4)]);      

                          
      pid = plot(freqSimData.nominalForce(1,idxSim),...
           freqSimData.damping(1,idxSim)/1000,...
           modelMarkType, ...
           'Color', modelLineColor,...
           'LineWidth',modelLineWidth,...
           'MarkerSize',modelMarkSize,...
           'MarkerFaceColor',modelFaceColor,...
           'DisplayName',[freqSeriesName{z},perturbationName]);
      hold on;
      if( i > 1)
        set(get(get(pid,'Annotation'),...
                  'LegendInformation'),...
                  'IconDisplayStyle','off'); 
      end 
      xlim([fMin,fMax]);
      ylim([dMin,dMax]);
    end
  end

  subplot('Position', [ subPlotList(idxD,1),...
                      subPlotList(idxD,2),...
                      subPlotList(idxD,3),...
                      subPlotList(idxD,4)]);   
  legend('Location','Northwest');
  
  
  set(fig_fig12,'Units','centimeters',...
  'PaperUnits','centimeters',...
  'PaperSize',[pageWidth pageHeight],...
  'PaperPositionMode','manual',...
  'PaperPosition',[0 0 pageWidth pageHeight]);     
  %set(findall(figList(i).h,'-property','FontSize'),'FontSize',10);     
  set(fig_fig12,'renderer','painters');     
  set(gcf,'InvertHardCopy','off')

  tendonTag = '_ElasticTendon';
  if(flag_useElasticTendon==0)
    tendonTag = '_RigidTendon';
  end

  print('-dpdf', [pubOutputFolder,'fig_Pub_StiffnessDampingKBRFig12',...
                                  tendonTag,'_',plotNameEnding,'.pdf']);   
  

success = 1;