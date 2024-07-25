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

%%
% Setup plot parameters
%%
%totalWidth = 177.13668/10; %Frontiers journal text width.
totalWidth = 21.0;

if(exist('pageWidth','var')==0)
    pageWidth  = 21.0;
end 
if(exist('pageHeight','var')==0)
    pageHeight = 29.7;
end

if(flag_usingOctave == 0)
  set(groot, 'defaultAxesFontSize',8);
  set(groot, 'defaultTextFontSize',8);
  set(groot, 'defaultAxesLabelFontSizeMultiplier',1.2);
  set(groot, 'defaultAxesTitleFontSizeMultiplier',1.2);
  set(groot, 'defaultAxesTickLabelInterpreter','latex');
  set(groot, 'defaultLegendInterpreter','latex');
  set(groot, 'defaultTextInterpreter','latex');
  set(groot, 'defaultAxesTitleFontWeight','bold');  
  set(groot, 'defaultFigurePaperUnits','centimeters');
  set(groot, 'defaultFigurePaperSize',[pageWidth pageHeight]);
  set(groot,'defaultFigurePaperType','A4');
end

if(exist('plotHorizMarginCm','var')==0)
  plotHorizMarginCm = 7;  
end

if(exist('plotVertMarginCm','var')==0)
    plotVertMarginCm  = 1;
end

if(exist('plotHorizMarginSplitCm','var')==0)
    plotHorizMarginSplitCm = 2.0;
end
if(exist('plotVertMarginSplitCm','var')==0)
    plotVertMarginSplitCm  = 2.0;
end

%plotWidth = 7; %((pageWidth-plotHorizMarginCm)/numberOfHorizontalPlotColumns);
if(isempty(plotWidth)==1)
  plotWidth=((pageWidth-plotHorizMarginCm)/numberOfHorizontalPlotColumns);
end
if(isempty(plotHeight)==1)
  plotHeight= ((pageHeight-plotVertMarginCm)/numberOfVerticalPlotRows);
  if(plotWidth < plotHeight)
    plotHeight=plotWidth; 
  else
    plotWidth=plotHeight;
  end  
end


if(flag_fixedPlotWidth == 0)
  plotWidth = ((pageWidth-2*plotHorizMarginCm)/numberOfHorizontalPlotColumns);
end

plotWidth  = plotWidth/pageWidth;
plotHeight = plotHeight/pageHeight;
halfPlotHeight=0.5*plotHeight;
fullPlotHeight=plotHeight;

plotHeightSquare = plotWidth*(pageWidth/pageHeight);

plotHorizMargin = plotHorizMarginCm/pageWidth;
plotVertMargin  = plotVertMarginCm/pageHeight;

plotHorizMarginSplit = plotHorizMarginSplitCm/pageWidth;
plotVertMarginSplit  = plotVertMarginSplitCm/pageHeight;

oneCmVertical = 1/pageHeight;
oneCmHorizontal = 1/pageWidth;


topLeft = [0/pageWidth pageHeight/pageHeight];

subPlotSquare = zeros(1,4);
subPlotSquare(1,1) = topLeft(1)  + plotHorizMargin;
subPlotSquare(1,2) = topLeft(2)  - 3*plotVertMargin - plotHeightSquare;
subPlotSquare(1,3) = plotWidth;
subPlotSquare(1,4) = plotHeightSquare;


subPlotRectangle = zeros(1,4);
subPlotRectangle(1,1) = topLeft(1)  + plotHorizMargin;
subPlotRectangle(1,2) = topLeft(2)  - 3*plotVertMargin - plotHeight;
subPlotRectangle(1,3) = plotWidth;
subPlotRectangle(1,4) = 0.5*plotHeight;

subPlotRectangleSmall = zeros(1,4);
subPlotRectangleSmall(1,1) = topLeft(1)  + plotHorizMargin;
subPlotRectangleSmall(1,2) = topLeft(2)  - plotVertMargin/3 - plotHeight;
subPlotRectangleSmall(1,3) = plotWidth;
subPlotRectangleSmall(1,4) = 0.125*plotHeight;

subPlotRectangleMedium = zeros(1,4);
subPlotRectangleMedium(1,1) = topLeft(1)  + plotHorizMargin;
subPlotRectangleMedium(1,2) = topLeft(2)  - plotVertMargin/3 - plotHeight;
subPlotRectangleMedium(1,3) = plotWidth;
subPlotRectangleMedium(1,4) = 0.25*plotHeight;

if(exist('plotWidth43')==0)
  plotWidth43 = 5.5;
end
if(exist('plotHeight43')==0)
  plotHeight43 = 4.0;  
end

subPlotRectangle43 = zeros(1,4);
subPlotRectangle43(1,1) = topLeft(1)  + 2/pageWidth;
subPlotRectangle43(1,2) = topLeft(2)  - 2/pageHeight - plotHeight;
subPlotRectangle43(1,3) = plotWidth43/pageWidth;
subPlotRectangle43(1,4) = plotHeight43/pageHeight;



subPlotPanel=zeros(numberOfVerticalPlotRows,numberOfHorizontalPlotColumns,4);

flag_short22Plot = 0;
idx=1;
scale=1;
for(ai=1:1:numberOfVerticalPlotRows)
  for(aj=1:1:numberOfHorizontalPlotColumns)
      if(idx==4 && flag_short22Plot==1)
         scale = 0.5;
      else
          scale=1;
      end
      subPlotPanel(ai,aj,1) = topLeft(1) + plotHorizMargin...
                            + (aj-1)*(plotWidth);
      subPlotPanel(ai,aj,2) = topLeft(2) -plotHeight-plotVertMargin ...                                       
                            + (ai-1)*(-plotHeight);
      subPlotPanel(ai,aj,3) = (plotWidth-plotHorizMarginSplit);
      subPlotPanel(ai,aj,4) = (plotHeight*scale-plotHorizMarginSplit);
      idx=idx+1;
  end
end

subPlotPairPanel=zeros(numberOfVerticalPlotRows*2,numberOfHorizontalPlotColumns,4);

idx=1;
scale=1;

prevCorner = topLeft + [oneCmHorizontal*3, -oneCmVertical*2] ;

plotSmallVertMargin  = 0.5/pageHeight;
plotMedVertMargin    = 1.0/pageHeight;


for(ai=1:1:numberOfVerticalPlotRows)
  for(aj=1:1:numberOfHorizontalPlotColumns)

      rowIdx = (ai-1)*2+1;
    
      subPlotPairPanel(rowIdx,aj,1) = prevCorner(1);
      subPlotPairPanel(rowIdx,aj,2) = prevCorner(2) - 0.25*plotHeight-plotMedVertMargin;
      subPlotPairPanel(rowIdx,aj,3) = (plotWidth);
      subPlotPairPanel(rowIdx,aj,4) = (plotHeight*0.25);      

      prevCorner = [subPlotPairPanel(rowIdx,aj,1), subPlotPairPanel(rowIdx,aj,2)];
      
      rowIdx=rowIdx+1;
      
      subPlotPairPanel(rowIdx,aj,1) = prevCorner(1);
      subPlotPairPanel(rowIdx,aj,2) = prevCorner(2) - plotHeight - plotSmallVertMargin;
      subPlotPairPanel(rowIdx,aj,3) = (plotWidth);
      subPlotPairPanel(rowIdx,aj,4) = (plotHeight);

      
      prevCorner = [subPlotPairPanel(rowIdx,aj,1), subPlotPairPanel(rowIdx,aj,2)];
      prevCorner = prevCorner + [0, -plotVertMargin];
      
  end
end

subPlotHerzogLeonard2000Stability    = zeros(4,4);
prevCorner = topLeft + [oneCmHorizontal*3, -oneCmVertical*2];

rowIdx=1;
colIdx=1;
plotBigVertMargin    = 1.5/pageHeight;

subPlotHerzogLeonard2000Stability(rowIdx,colIdx,1) = prevCorner(1);
subPlotHerzogLeonard2000Stability(rowIdx,colIdx,2) = prevCorner(2) - 0.25*plotHeight-plotBigVertMargin;
subPlotHerzogLeonard2000Stability(rowIdx,colIdx,3) = (plotWidth);
subPlotHerzogLeonard2000Stability(rowIdx,colIdx,4) = (plotHeight*0.25); 

prevCorner = [subPlotHerzogLeonard2000Stability(rowIdx,colIdx,1),...
              subPlotHerzogLeonard2000Stability(rowIdx,colIdx,2)];

rowIdx=rowIdx+1;
colIdx=1;
subPlotHerzogLeonard2000Stability(rowIdx,colIdx,1) = prevCorner(1);
subPlotHerzogLeonard2000Stability(rowIdx,colIdx,2) = prevCorner(2) - plotHeight-plotBigVertMargin;
subPlotHerzogLeonard2000Stability(rowIdx,colIdx,3) = (plotWidth);
subPlotHerzogLeonard2000Stability(rowIdx,colIdx,4) = (plotHeight); 

prevCorner = [subPlotHerzogLeonard2000Stability(rowIdx,colIdx,1),...
              subPlotHerzogLeonard2000Stability(rowIdx,colIdx,2)];

rowIdx=rowIdx+1;
colIdx=1;
subPlotHerzogLeonard2000Stability(rowIdx,colIdx,1) = prevCorner(1);
subPlotHerzogLeonard2000Stability(rowIdx,colIdx,2) = prevCorner(2) - plotHeight-plotBigVertMargin;
subPlotHerzogLeonard2000Stability(rowIdx,colIdx,3) = (plotWidth);
subPlotHerzogLeonard2000Stability(rowIdx,colIdx,4) = (plotHeight); 

prevCorner = [subPlotHerzogLeonard2000Stability(rowIdx,colIdx,1),...
              subPlotHerzogLeonard2000Stability(rowIdx,colIdx,2)];
            
rowIdx=rowIdx+1;
colIdx=1;
subPlotHerzogLeonard2000Stability(rowIdx,colIdx,1) = prevCorner(1);
subPlotHerzogLeonard2000Stability(rowIdx,colIdx,2) = prevCorner(2) - plotHeight-plotBigVertMargin;
subPlotHerzogLeonard2000Stability(rowIdx,colIdx,3) = (plotWidth);
subPlotHerzogLeonard2000Stability(rowIdx,colIdx,4) = (plotHeight); 

% Row

subPlotHerzogLeonard2000StabilityRow    = zeros(4,4);
prevCorner = topLeft + [oneCmHorizontal*3, -oneCmVertical*2];
plotWidthRow = plotWidth*0.7;
plotHeightRow= plotWidthRow*0.7;

rowIdx=1;
colIdx=1;
plotBigVertMargin    = 1.25/pageHeight;
plotBigHorizMargin   = 1.25/pageWidth;

subPlotHerzogLeonard2000StabilityRow(rowIdx,colIdx,1) = prevCorner(1);
subPlotHerzogLeonard2000StabilityRow(rowIdx,colIdx,2) = prevCorner(2) ...
                                    - 0.25*plotHeightRow-plotBigVertMargin;
subPlotHerzogLeonard2000StabilityRow(rowIdx,colIdx,3) = (plotWidthRow);
subPlotHerzogLeonard2000StabilityRow(rowIdx,colIdx,4) = (plotHeightRow*0.25); 

prevCorner = [subPlotHerzogLeonard2000StabilityRow(rowIdx,colIdx,1),...
              subPlotHerzogLeonard2000StabilityRow(rowIdx,colIdx,2)];

rowIdx=rowIdx+1;
colIdx=1;
subPlotHerzogLeonard2000StabilityRow(rowIdx,colIdx,1) = prevCorner(1);
subPlotHerzogLeonard2000StabilityRow(rowIdx,colIdx,2) = prevCorner(2) ...
                                      - plotHeightRow - plotBigVertMargin;
subPlotHerzogLeonard2000StabilityRow(rowIdx,colIdx,3) = (plotWidthRow);
subPlotHerzogLeonard2000StabilityRow(rowIdx,colIdx,4) = (plotHeightRow); 

prevCorner = [subPlotHerzogLeonard2000StabilityRow(rowIdx,colIdx,1),...
              subPlotHerzogLeonard2000StabilityRow(rowIdx,colIdx,2)];

rowIdx=rowIdx+1;
colIdx=1;
subPlotHerzogLeonard2000StabilityRow(rowIdx,colIdx,1) = prevCorner(1) ...
                                     +plotWidthRow + plotBigHorizMargin;
subPlotHerzogLeonard2000StabilityRow(rowIdx,colIdx,2) = prevCorner(2);
                                        
subPlotHerzogLeonard2000StabilityRow(rowIdx,colIdx,3) = (plotWidthRow);
subPlotHerzogLeonard2000StabilityRow(rowIdx,colIdx,4) = (plotHeightRow); 

prevCorner = [subPlotHerzogLeonard2000StabilityRow(rowIdx,colIdx,1),...
              subPlotHerzogLeonard2000StabilityRow(rowIdx,colIdx,2)];
            
rowIdx=rowIdx+1;
colIdx=1;
subPlotHerzogLeonard2000StabilityRow(rowIdx,colIdx,1) = prevCorner(1) ...
                                     +plotWidthRow + plotBigHorizMargin;    
subPlotHerzogLeonard2000StabilityRow(rowIdx,colIdx,2) = prevCorner(2);
subPlotHerzogLeonard2000StabilityRow(rowIdx,colIdx,3) = (plotWidthRow);
subPlotHerzogLeonard2000StabilityRow(rowIdx,colIdx,4) = (plotHeightRow); 


%%
%Line config
%%

lineColor             = [255 0  0;...
                         0   255   0;...
                         0   0   255]./255;                     
lineWidth             = [2,1,1];
lineType              ={'-','-','-'};

xTickTime = [0,0.5,1,1.5,2,2.5,3,3.5];


