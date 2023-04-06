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




