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

function [success] = plotStiffnessDampingVariationPUB(dataFolder,...
                                              freqSeriesFiles,...
                                              freqSeriesName,...
                                              freqSeriesColor,...
                                              inputFunctions,...                      
                                              normFiberLength,...
                                              nominalForce,...  
                                              dataKBR1994Fig9A,...
                                              dataKBR1994Fig9B,... 
                                              dataKBR1994Fig10,...
                                              flag_useElasticTendon,...
                                              plotNameEnding,...
                                              plotLayoutSettings,...
                                              plotFolder)
                                            

success = 0;

assert(length(normFiberLength)==1);

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


fig_fig9AB10 = figure;


subPlotList = zeros(2,4);

ySpace = 0.05;
xSpace = 0.075;

currentSubPlot  = subPlotSquare;
currentSubPlot(1,1) = currentSubPlot(1,1)*0.5; 

subPlotHeight   = currentSubPlot(1,4);
subPlotWidth    = currentSubPlot(1,3);
subPlotOffsetY  = subPlotHeight + ySpace;
subPlotOffsetX  = subPlotWidth  + xSpace;

figure(fig_fig9AB10);

subPlotList(1,:)    = currentSubPlot(1,:);
currentSubPlot(1,2) = currentSubPlot(1,2) - subPlotOffsetY;
subPlotList(2,:)    = currentSubPlot(1,:);
currentSubPlot(1,2) = currentSubPlot(1,2) - subPlotOffsetY;
subPlotList(3,:)    = currentSubPlot(1,:);

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

idx9A = 1;
idx9B = 2;
idx10 = 3;

mdlAmp9A = [0.4,0.8,1.6];
mdlFreq9A= [ 15, 15, 15];
mdl9APlotColorMod = [0,0.5,1];

exp9APlotColor      = [0.5,  0.5, 0.5;...
                       0.5,  0.5, 0.5;...
                       0.5,  0.5, 0.5];
exp9ALineWidth      = [1,1,1].*0.5;
exp9AMarkType       = {'o', 'o', 'o'};
exp9AMarkSize       = [ 4,   4,  4];
exp9AMarkFaceColor  = [ 0.5, 0.5, 0.5;...
                        0.75, 0.75, 0.75;...
                        1, 1, 1];
exp9ALegendEntry    = {'Kirsch: 0.4mm 15Hz',...
                       'Kirsch: 0.8mm 15Hz',...
                       'Kirsch: 1.6mm 15Hz'};

                     
                     
exp9BPlotColor      = [0.5,  0.5, 0.5;...
                       0.5,  0.5, 0.5;...
                       0.5,  0.5, 0.5;...
                       0.5,  0.5, 0.5;...
                       0.5,  0.5, 0.5;...
                       0.5,  0.5, 0.5];
                     
exp9BLineWidth      = [1,1,1, 1,1,1].*0.5;
exp9BMarkType       = {'o', 'o', 'o', 's','s','s'};
exp9BMarkSize       = [  4,   4,   4,  5,  5,  5];
exp9BMarkFaceColor  = [ 0.5, 0.5, 0.5;...
                        0.75, 0.75, 0.75;...
                        1, 1, 1;...
                        0.5, 0.5, 0.5;...
                        0.75, 0.75, 0.75;...
                        1, 1, 1];
exp9BLegendEntry    = {'Kirsch: 0.4mm 15Hz',...
                       'Kirsch: 0.4mm 35Hz',...
                       'Kirsch: 0.4mm 90Hz',...
                       'Kirsch: 1.6mm 15Hz',...
                       'Kirsch: 1.6mm 35Hz',...
                       'Kirsch: 1.6mm 90Hz'};
                     
mdlAmp9B          = [0.4,0.4,0.4,1.6,1.6,1.6];
mdlFreq9B         = [ 15, 35, 90, 15, 35, 90];                     
mdl9BPlotColorMod = [  0,0.5,  1,  0,0.5,  1];                     
                     

exp10PlotColor      = [0.5,  0.5, 0.5;...
                       0.5,  0.5, 0.5;...
                       0.5,  0.5, 0.5];
exp10LineWidth      = [1,1,1].*0.5;

exp10MarkType       = {'o', 'o', 'o'};
exp10MarkSize       = [ 4,   4,  4];
exp10MarkFaceColor  = [ 0.5, 0.5, 0.5;...
                        0.75, 0.75, 0.75;...
                        1, 1, 1];
exp10LegendEntry    = {'Kirsch: 15Hz',...
                       'Kirsch: 35Hz',...
                       'Kirsch: 90Hz'};
mdlAmp10         = [0.4, 0.4, 0.4, 0.8, 0.8, 0.8, 1.6, 1.6, 1.6];
mdlFreq10        = [ 15,  35,  90,  15,  35,  90,  15,  35,  90];      
mdl10PlotColorMod= [  0, 0.5,   1,   0, 0.5,   1,   0, 0.5,   1];
mdlData10Idx = [1,2,3,1,2,3,1,2,3];
                     
%%
%9A
%%
subplot('Position', [ subPlotList(idx9A,1),...
                      subPlotList(idx9A,2),...
                      subPlotList(idx9A,3),...
                      subPlotList(idx9A,4)]);
  
                    
for i=1:1:length(dataKBR1994Fig9A)
  plot(dataKBR1994Fig9A(i).x,...
       dataKBR1994Fig9A(i).y,...
       exp9AMarkType{i},...
       'Color',     exp9APlotColor(i,:),...
       'LineWidth', exp9ALineWidth(1,i),...
       'MarkerFaceColor',exp9AMarkFaceColor(i,:),...
       'MarkerSize', exp9AMarkSize(1,i),...
       'DisplayName', exp9ALegendEntry{i});
    hold on;         
       
end

for z=1:1:length(freqSeriesFiles)

  modelColor = freqSeriesColor(z,:);
  tmp = load([dataFolder, freqSeriesFiles{1,z}]);
  freqSimData = tmp.freqSimData;

  flag_Hill = 0;
  if(isempty(strfind(freqSeriesFiles{1,z},'Hill'))==0)
    flag_Hill = 1;
  end


  
  for i=1:1:length(mdlAmp9A)
    modelMarkType   = exp9AMarkType{i};
    modelFaceColor  = (1-mdl9APlotColorMod(1,i)).*modelColor ...
                     + mdl9APlotColorMod(1,i).*[1,1,1];
    modelMarkSize   = exp9AMarkSize(1,i);
    modelLineColor  = modelColor;
    modelLineWidth  = exp9ALineWidth(1,i);
    
    targetAmplitudeMM = mdlAmp9A(1,i);
    targetBandwidthHz = mdlFreq9A(1,i);
    targetNormFiberLength = normFiberLength;
  
    for j=1:1:length(nominalForce)
      targetNominalForceN = nominalForce(1,j);
      
      idxSim = getFrequencySimulationIndex(targetAmplitudeMM,  ...
                                           targetBandwidthHz, ...
                                           targetNominalForceN,...
                                           targetNormFiberLength,...
                                           freqSimData); 
      
      pid = plot(freqSimData.nominalForce(1,idxSim),...
           freqSimData.stiffness(1,idxSim)/1000,...
           modelMarkType, ...
           'Color', modelLineColor,...
           'LineWidth',modelLineWidth,...
           'MarkerSize',modelMarkSize,...
           'MarkerFaceColor',modelFaceColor,...
           'DisplayName',freqSeriesName{z});
      hold on;
      if(i > 1 || j > 1)
        set(get(get(pid,'Annotation'),...
                  'LegendInformation'),...
                  'IconDisplayStyle','off'); 
      end                                         
    end
  end
  
end


lh=legend('location','east');
legendPosition = get(lh,'position');
legendPosition(1,1) = legendPosition(1,1)+subPlotWidth;
set(lh,'position',legendPosition);                 

box off;
set(gca,'color','none')   


xlabel('Force (N)');
ylabel('Stiffness (N/mm)');

%%
%9B
%%
subplot('Position', [ subPlotList(idx9B,1),...
                      subPlotList(idx9B,2),...
                      subPlotList(idx9B,3),...
                      subPlotList(idx9B,4)]);
  
                    
for i=1:1:length(dataKBR1994Fig9B)
  plot(dataKBR1994Fig9B(i).x,...
       dataKBR1994Fig9B(i).y,...
       exp9BMarkType{i},...
       'Color',     exp9BPlotColor(i,:),...
       'LineWidth', exp9BLineWidth(1,i),...
       'MarkerFaceColor',exp9BMarkFaceColor(i,:),...
       'MarkerSize', exp9BMarkSize(1,i),...
       'DisplayName', exp9BLegendEntry{i});
    hold on;         
       
end

for z=1:1:length(freqSeriesFiles)

  modelColor = freqSeriesColor(z,:);
  tmp = load([dataFolder, freqSeriesFiles{1,z}]);
  freqSimData = tmp.freqSimData;

  flag_Hill = 0;
  if(isempty(strfind(freqSeriesFiles{1,z},'Hill'))==0)
    flag_Hill = 1;
  end


  
  for i=1:1:length(mdlAmp9B)
    modelMarkType   = exp9BMarkType{i};
    modelFaceColor  = (1-mdl9BPlotColorMod(1,i)).*modelColor ...
                     + mdl9BPlotColorMod(1,i).*[1,1,1];
    modelMarkSize   = exp9BMarkSize(1,i);
    modelLineColor  = modelColor;
    modelLineWidth  = exp9BLineWidth(1,i);
    
    targetAmplitudeMM = mdlAmp9B(1,i);
    targetBandwidthHz = mdlFreq9B(1,i);
    targetNormFiberLength = normFiberLength;
  
    for j=1:1:length(nominalForce)
      targetNominalForceN = nominalForce(1,j);
      
      idxSim = getFrequencySimulationIndex(targetAmplitudeMM,  ...
                                           targetBandwidthHz, ...
                                           targetNominalForceN,...
                                           targetNormFiberLength,...
                                           freqSimData); 
      
      pid = plot(freqSimData.nominalForce(1,idxSim),...
           freqSimData.stiffness(1,idxSim)/1000,...
           modelMarkType, ...
           'Color', modelLineColor,...
           'LineWidth',modelLineWidth,...
           'MarkerSize',modelMarkSize,...
           'MarkerFaceColor',modelFaceColor,...
           'DisplayName',freqSeriesName{z});
      hold on;
      if(i > 1 || j > 1)
        set(get(get(pid,'Annotation'),...
                  'LegendInformation'),...
                  'IconDisplayStyle','off'); 
      end                                         
    end
  end
  
end

lh=legend('location','east');
legendPosition = get(lh,'position');
legendPosition(1,1) = legendPosition(1,1)+subPlotWidth;
set(lh,'position',legendPosition);                 

box off;
set(gca,'color','none')   


xlabel('Force (N)');
ylabel('Stiffness (N/mm)');

%%
%9C
%%
subplot('Position', [ subPlotList(idx10,1),...
                      subPlotList(idx10,2),...
                      subPlotList(idx10,3),...
                      subPlotList(idx10,4)]);

for i=1:1:length(dataKBR1994Fig10)
  plot(dataKBR1994Fig10(i).x,...
       dataKBR1994Fig10(i).y,...
       exp10MarkType{i},...
       'Color',     exp10PlotColor(i,:),...
       'LineWidth', exp10LineWidth(1,i),...
       'MarkerFaceColor',exp10MarkFaceColor(i,:),...
       'MarkerSize', exp10MarkSize(1,i),...
       'DisplayName', exp10LegendEntry{i});
    hold on;         
       
end

               
       
%%
%Model
%%
for z=1:1:length(freqSeriesFiles)

  modelColor = freqSeriesColor(z,:);
  tmp = load([dataFolder, freqSeriesFiles{1,z}]);
  freqSimData = tmp.freqSimData;

  flag_Hill = 0;
  if(isempty(strfind(freqSeriesFiles{1,z},'Hill'))==0)
    flag_Hill = 1;
  end


  
  for i=1:1:length(mdlAmp10)
    idxData = mdlData10Idx(1,i);
    modelMarkType   = exp10MarkType{idxData};
    modelFaceColor  = (1-mdl10PlotColorMod(1,i)).*modelColor ...
                     + mdl10PlotColorMod(1,i).*[1,1,1];
    modelMarkSize   = exp10MarkSize(1,idxData);
    modelLineColor  = modelColor;
    modelLineWidth  = exp10LineWidth(1,idxData);
    
    targetAmplitudeMM = mdlAmp10(1,i);
    targetBandwidthHz = mdlFreq10(1,i);
    targetNormFiberLength = normFiberLength;
  
    for j=1:1:length(nominalForce)
      targetNominalForceN = nominalForce(1,j);
      
      idxSim = getFrequencySimulationIndex(targetAmplitudeMM,  ...
                                           targetBandwidthHz, ...
                                           targetNominalForceN,...
                                           targetNormFiberLength,...
                                           freqSimData); 
      
      pid = plot(freqSimData.nominalForce(1,idxSim),...
           freqSimData.damping(1,idxSim)/1000,...
           modelMarkType, ...
           'Color', modelLineColor,...
           'LineWidth',modelLineWidth,...
           'MarkerSize',modelMarkSize,...
           'MarkerFaceColor',modelFaceColor,...
           'DisplayName',freqSeriesName{z});
      hold on;
      if(i > 1 || j > 1)
        set(get(get(pid,'Annotation'),...
                  'LegendInformation'),...
                  'IconDisplayStyle','off'); 
      end                                         
    end
  end
  
end

                 
lh=legend('location','east');
legendPosition = get(lh,'position');
legendPosition(1,1) = legendPosition(1,1)+subPlotWidth;
set(lh,'position',legendPosition);                 
                 
box off;
set(gca,'color','none')   


xlabel('Force (N)');
ylabel('Damping (N/(mm/s))');     


set(fig_fig9AB10,'Units','centimeters',...
'PaperUnits','centimeters',...
'PaperSize',[pageWidth pageHeight],...
'PaperPositionMode','manual',...
'PaperPosition',[0 0 pageWidth pageHeight]);     
%set(findall(figList(i).h,'-property','FontSize'),'FontSize',10);     
set(fig_fig9AB10,'renderer','painters');     
set(gcf,'InvertHardCopy','off')

tendonTag = '_ElasticTendon';
if(flag_useElasticTendon==0)
  tendonTag = '_RigidTendon';
end

print('-dpdf', [plotFolder,'fig_Pub_StiffnessDampingFrequencyVariation',...
                                tendonTag,'_',plotNameEnding,'.pdf']); 

success = 1;
                                            
                                            