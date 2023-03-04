%%
% Setup plot parameters
%%

plotFontName = 'latex';

if(flag_usingOctave == 0)
set(groot, 'defaultAxesFontSize',8);
set(groot, 'defaultTextFontSize',8);
set(groot, 'defaultAxesLabelFontSizeMultiplier',1.2);
set(groot, 'defaultAxesTitleFontSizeMultiplier',1.2);
set(groot, 'defaultAxesTickLabelInterpreter','latex');
%set(groot, 'defaultAxesFontName',plotFontName);
%set(groot, 'defaultTextFontName',plotFontName);
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTitleFontWeight','bold');  
set(groot, 'defaultFigurePaperUnits','centimeters');
set(groot, 'defaultFigurePaperSize',[pageWidth pageHeight]);
set(groot,'defaultFigurePaperType','A4');
end



plotWidth  = plotWidth/pageWidth;
plotHeight = plotHeight/pageHeight;

plotHorizMargin = plotHorizMarginCm/pageWidth;
plotVertMargin  = plotVertMarginCm/pageHeight;

topLeft = [0/pageWidth pageHeight/pageHeight];

subPlotPanel=zeros(numberOfVerticalPlotRows,numberOfHorizontalPlotColumns,4);
subPlotPanelIndex = zeros(numberOfVerticalPlotRows,numberOfHorizontalPlotColumns);

idx=1;
scaleVerticalMargin = 0.;
for(ai=1:1:numberOfVerticalPlotRows)
  if(ai > 1)
    scaleVerticalMargin = 1;
  end
  for(aj=1:1:numberOfHorizontalPlotColumns)
      subPlotPanelIndex(ai,aj) = idx;
      scaleHorizMargin=1;
      subPlotPanel(ai,aj,1) = topLeft(1) + plotHorizMargin...
                            + (aj-1)*(plotWidth + plotHorizMargin);
      %-plotVertMargin*scaleVerticalMargin ...                             
      subPlotPanel(ai,aj,2) = topLeft(2) -plotHeight -plotVertMargin...                            
                            + (ai-1)*(-plotHeight -plotVertMargin);
      subPlotPanel(ai,aj,3) = (plotWidth);
      subPlotPanel(ai,aj,4) = (plotHeight);
      idx=idx+1;
  end
end




