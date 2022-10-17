function [success] = plotStiffnessDampingPUB20191227(...
                      dataInputFolder,...
                      freqSeriesFiles_RT,...
                      freqSeriesName_RT,...
                      freqSeriesColor_RT,...
                      freqSeriesFiles_ET,...
                      freqSeriesName_ET,...
                      freqSeriesColor_ET,...       
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


assert(length(freqSeriesFiles_RT)==2);
assert(length(freqSeriesFiles_ET)==2);

assert( isempty(strfind(freqSeriesFiles_RT{1,1},'Hill'))==0 ...
     && isempty(strfind(freqSeriesFiles_ET{1,1},'Hill'))==0);

freqSeriesFiles = {freqSeriesFiles_ET{1,2},...
                   freqSeriesFiles_RT{1,2},...
                   freqSeriesFiles_ET{1,1},...
                   freqSeriesFiles_RT{1,1}};
freqSeriesName =  {freqSeriesName_ET{1,2},...
                   freqSeriesName_RT{1,2},...
                   freqSeriesName_ET{1,1},...
                   freqSeriesName_RT{1,1}};
freqSeriesColor = [freqSeriesColor_ET(2,:);...
                   freqSeriesColor_RT(2,:);...
                   freqSeriesColor_ET(1,:);...
                   freqSeriesColor_RT(1,:)];
freqSeriesSubPlotOffsetIdx = [0,2,0,2];

subPlotLabel = {'A. Models with elastic tendons','B.','C. Models with rigid tendons',... 
                'D.',};
figLabelFontSize = 8*1.2;

numberOfHorizontalPlotColumns = plotLayoutSettings.numberOfHorizontalPlotColumns;
numberOfVerticalPlotRows      = plotLayoutSettings.numberOfVerticalPlotRows;
flag_fixedPlotWidth           = plotLayoutSettings.flag_fixedPlotWidth;
plotWidth                     = plotLayoutSettings.plotWidth;
plotWidth43                   = plotLayoutSettings.plotWidth43;
plotHeight43                  = plotLayoutSettings.plotHeight43;
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

  subPlotList = zeros(4,4);

  ySpace = 0.075;
  xSpace = 0.075;
  
  kMin = -0.01;
  kMax = 15.01;
  dMin = -0.001;
  dMax = 0.1501;
  
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
  % Rigid Tendon
  %%
  % Stiffness vs. force
  subPlotList(1,:) = currentSubPlot(1,:);

  % Damping vs. force
  subPlotList(2,:)    = subPlotList(1,:);
  subPlotList(2,1)    = subPlotList(2,1)+ subPlotOffsetX;  

  %%
  % Elastic Tendon
  %%
  % Stiffness vs. force
  subPlotList(3,:) = subPlotList(1,:);
  subPlotList(3,2) = subPlotList(1,2)- subPlotOffsetY;

  % Damping vs. force
  subPlotList(4,:)    = subPlotList(3,:);
  subPlotList(4,1)    = subPlotList(4,1) + subPlotOffsetX;  
  
  
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
  idxK = 0;
  idxD = 0;
  
  figure(fig_fig12);
  
  for z=1:1:2
    idxK = 1 + 2*(z-1);
    idxD = 2 + 2*(z-1);

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
    
    %if(z==2)
    %  modelMarkSize = modelMarkSize -1;
    %  modelMarkType = 'd';
    %end
    
    
    flag_Hill = 0;
    if(isempty(strfind(freqSeriesFiles{1,z},'Hill'))==0)
      flag_Hill = 1;
    end

    idxK = 1+freqSeriesSubPlotOffsetIdx(1,z);
    idxD = 2+freqSeriesSubPlotOffsetIdx(1,z);
    
    for i=1:1:length(nominalForce)
    
      idxSim = 0;
      tol = 1e-6;
      for m=1:1:size(freqSimData.force,2)     
        if( abs(freqSimData.amplitudeMM(1,m)     - targetAmplitude   ) <= tol && ...
            abs(freqSimData.bandwidthHz(1,m)     - targetBandwidth   ) <= tol && ...
            abs(freqSimData.nominalForceDesired(1,m) - nominalForce(1,i) ) <= tol && ...
            abs(freqSimData.normFiberLength(1,m) - normFiberLength   ) <= tol)
          if(idxSim == 0)
            idxSim = m;
          else
            assert(0); %Error condition: there should not be 2 simulations with 
                       %the same configuration
          end
        end
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
      
      posTextK = [freqSimData.nominalForce(1,idxSim),...
                freqSimData.stiffness(1,idxSim)/1000];
      vafText = sprintf('%1.0f',freqSimData.vafTime(1,idxSim)*100); 
      if(i==3)
        vafText = ['VAF:',vafText,'\%'];
      else
        vafText=[vafText,'\%'];
      end

      
      textDeltaX  = 0;
      textDeltaY  = 0;
      textAlign   = '';      
      if(flag_useElasticTendon == 1)
        textDeltaX = -0.25;
        textDeltaY = 0;
        textAlign = 'right';
        if(flag_Hill == 1)
          textDeltaX = 0.25;%-0.25;
          textDeltaY = 0;
          textAlign = 'left';
        end
      else
        textDeltaX = -0.25;
        textDeltaY = 0;%0.25;
        textAlign = 'right';
        if(flag_Hill == 1)
          textDeltaX = 0.25;
          textDeltaY = 0;%-0.35;
          textAlign = 'left';
        end
        
      end
      text(posTextK(1,1)+textDeltaX,...
           posTextK(1,2)+textDeltaY,vafText,...
            'HorizontalAlignment',textAlign);
      hold on;
      if(i==3)
         xyText =[sprintf('(%1.2f,%1.2f)',...
                freqSimData.nominalForce(1,idxSim),...
                freqSimData.stiffness(1,idxSim)/1000)];
          text(posTextK(1,1)+textDeltaX,...
               posTextK(1,2)+textDeltaY - 0.05*(kMax-kMin),xyText,...
                'HorizontalAlignment',textAlign);
         
      end
      xlim([fMin,fMax]);
      ylim([kMin,kMax]);      

      if(i == 1 && z==1 || i==1 && z==2)

        tc = text(fMin-0.1*(fMax-fMin), 1.1*kMax,...
               subPlotLabel{1,idxK},...
               'FontSize',figLabelFontSize);   
        hold on;
      end              
      
      
      if(z==1 && i==1)
        x0 = subPlotList(idxK,1);
        y0 = subPlotList(idxK,2);
        dx = subPlotList(idxK,3);
        dy = subPlotList(idxK,4);

        text(fMin,kMax*1.2,'Simulation of Kirsch, Boskov, \& Rymer 1994',...
        'FontSize',8*1.2,...
        'HorizontalAlignment','left',...
        'VerticalAlignment','bottom');
        hold on;     

      end
      
      subplot('Position', [ subPlotList(idxD,1),...
                            subPlotList(idxD,2),...
                            subPlotList(idxD,3),...
                            subPlotList(idxD,4)]);      

                          
      pid = plot( freqSimData.nominalForce(1,idxSim),...
                  freqSimData.damping(1,idxSim)/1000,...
                  modelMarkType, ...
                  'Color', modelLineColor,...
                  'LineWidth',modelLineWidth,...
                  'MarkerSize',modelMarkSize,...
                  'MarkerFaceColor',modelFaceColor,...
                  'DisplayName',[freqSeriesName{z},perturbationName]);
      hold on;
      posTextD = [freqSimData.nominalForce(1,idxSim),...
                freqSimData.damping(1,idxSim)/1000];
              
              
      text(posTextD(1,1)+textDeltaX,...
           posTextD(1,2)+textDeltaY,[vafText],...
            'HorizontalAlignment',textAlign);
      hold on;
      if(i==3)
         xyText =[sprintf('(%1.2f,%1.2f)',...
                freqSimData.nominalForce(1,idxSim),...
                freqSimData.damping(1,idxSim)/1000)];
          text(posTextD(1,1)+textDeltaX,...
               posTextD(1,2)+textDeltaY - 0.05*(dMax-dMin),xyText,...
                'HorizontalAlignment',textAlign);
         
      end

      
      if( i > 1)
        set(get(get(pid,'Annotation'),...
                  'LegendInformation'),...
                  'IconDisplayStyle','off'); 
      end 
      xlim([fMin,fMax]);
      ylim([dMin,dMax]);

      if(i == 1 && z==1 || i==1 && z==2)
        tc = text(fMin-0.1*(fMax-fMin), 1.1*dMax,...
               subPlotLabel{1,idxD},...
               'FontSize',figLabelFontSize);   
        hold on;
      end              
            
    end
  end

  subplot('Position', [ subPlotList(1,1),...
                      subPlotList(1,2),...
                      subPlotList(1,3),...
                      subPlotList(1,4)]);   
  legend('Location','Northwest');
  legend boxoff;
  
  subplot('Position', [ subPlotList(3,1),...
                      subPlotList(3,2),...
                      subPlotList(3,3),...
                      subPlotList(3,4)]);   
  legend('Location','Northeast');
  legend boxoff;
    
  
  
  set(fig_fig12,'Units','centimeters',...
  'PaperUnits','centimeters',...
  'PaperSize',[pageWidth pageHeight],...
  'PaperPositionMode','manual',...
  'PaperPosition',[0 0 pageWidth pageHeight]);     
  %set(findall(figList(i).h,'-property','FontSize'),'FontSize',10);     
  set(fig_fig12,'renderer','painters');     
  set(gcf,'InvertHardCopy','off')

  tendonTag = '_RigidAndElasticTendon';

  print('-dpdf', [pubOutputFolder,'fig_Pub_StiffnessDampingKBRFig12_',...
                                  tendonTag,'_',plotNameEnding,'.pdf']);   
  

success = 1;