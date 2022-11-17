function [success] = plotStiffnessDampingVAFBarChart( ...
                            opus31Table_RT,...
                            hillTable_RT,...
                            opus31Table_ET,...
                            hillTable_ET,...
                            kbr1994Table,...                            
                            opus31Color,...
                            hillColor,...
                            kbr1994Color,...
                            flag_useElasticTendon,... 
                            strFittingBandwidth,...
                            plotFolder,...
                            flag_usingOctave)

success=0;

freqTable = [15,35,90];
ampTable = [0.4,0.8,1.6];
noDataCode = -42;

freq    = [ 15,  15,  15,  35,  35,  35,  90,  90,  90];
amp     = [0.4, 0.8, 1.6, 0.4, 0.8, 1.6, 0.4, 0.8, 1.6];

tableFreqCol=[  1,   1,   1,   2,   2,   2,   3,   3,   3];
tableAmpRow =[  1,   2,   3,   1,   2,   3,   1,   2,   3];

fig_barChart = figure;

numberOfHorizontalPlotColumns = 3;
numberOfVerticalPlotRows      = 3;
plotHorizMarginCm             = 2;
plotVertMarginCm              = 2;
plotHorizMarginSplitCm        = 1;
plotVertMarginSplitCm         = 2;
flag_fixedPlotWidth           = 1;
plotWidth                     = 5;
plotWidth43                   = 4;
plotHeight43                  = 3;
plotHeight                    = plotWidth;
pageWidth   = numberOfHorizontalPlotColumns*plotWidth...
    +(numberOfHorizontalPlotColumns+2)*plotHorizMarginSplitCm...
    +plotHorizMarginCm;
pageHeight  = numberOfVerticalPlotRows*plotHeight...
    +(numberOfVerticalPlotRows+2)*plotVertMarginSplitCm...
    +plotVertMarginCm;

flag_usingOctave              = flag_usingOctave;
plotConfig;


barWidth   = 0.04;
skinnySpaceWidth = 0.02;
spaceWidth = 0.05;

subPlotTypes = {'stiffness','damping','vaf'};

type_rigidTendon=2;
type_elasticTendon=1;
yLimVals = [];
for idxTendon=1:1:2

    for idxType = 1:1:length(subPlotTypes)
    
        row=nan;
        switch subPlotTypes{idxType}
            case 'stiffness'
                if(idxTendon==type_elasticTendon)
                    row=1;            
                else
                    row=1;
                end
                yLimVals = [0.0;1.55];
            case 'damping'
                if(idxTendon==type_elasticTendon)
                    row=2;            
                else
                    row=2;
                end
                yLimVals = [0.000;0.025];
                
            case 'vaf'
                if(idxTendon==type_elasticTendon)
                    row=3;            
                else
                    row=3;
                end
                yLimVals = [0;100];

            otherwise
                assert(0,'Type not in list');            
        end
    
        for idxTrial=1:1:length(freq)    
            col = nan;
            switch freq(idxTrial)
                case 15
                    col=1;
                case 35
                    col=2;
                case 90
                    col=3;
                otherwise
                    assert(0,'Frequency not in list');
            end
            subplot('Position', reshape(subPlotPanel(row,col,:),1,4) );
    
            i = tableFreqCol(1,idxTrial);
            j = tableAmpRow(1,idxTrial);
            k = 1;
            
            seriesMean = [];
            seriesMin  = [];
            seriesMax   = [];
            
            fieldName = subPlotTypes{idxType};
            if(contains(fieldName,'vaf'))
                k=2;
            end
            if(idxTendon == type_elasticTendon)
                seriesMean = [kbr1994Table.(fieldName).pMean(i,j,k),...
                               opus31Table_ET.(fieldName).pMean(i,j,k),...
                               hillTable_ET.(fieldName).pMean(i,j,k)];
                seriesMin = [kbr1994Table.(fieldName).p95CIMin(i,j,k),...
                               opus31Table_ET.(fieldName).p95CIMin(i,j,k),...
                               hillTable_ET.(fieldName).p95CIMin(i,j,k)];
                seriesMax = [kbr1994Table.(fieldName).p95CIMax(i,j,k),...
                               opus31Table_ET.(fieldName).p95CIMax(i,j,k),...
                               hillTable_ET.(fieldName).p95CIMax(i,j,k)];
            else
                seriesMean = [kbr1994Table.(fieldName).pMean(i,j,k),...
                               opus31Table_RT.(fieldName).pMean(i,j,k),...
                               hillTable_RT.(fieldName).pMean(i,j,k)];
                seriesMin = [kbr1994Table.(fieldName).p95CIMin(i,j,k),...
                               opus31Table_RT.(fieldName).p95CIMin(i,j,k),...
                               hillTable_RT.(fieldName).p95CIMin(i,j,k)];
                seriesMax = [kbr1994Table.(fieldName).p95CIMax(i,j,k),...
                               opus31Table_RT.(fieldName).p95CIMax(i,j,k),...
                               hillTable_RT.(fieldName).p95CIMax(i,j,k)];
            end

            if(contains(fieldName,'vaf'))
                seriesMean  = seriesMean.*100;
                seriesMin   = seriesMin.*100;
                seriesMax   = seriesMax.*100;                
            end

            seriesColors = [];
            edgeColors = [];

            if(idxTendon == type_elasticTendon)
                seriesColors = [kbr1994Color;...
                                opus31Color;...
                                hillColor];
                edgeColors = [kbr1994Color;...
                                opus31Color;...
                                hillColor];
%                 edgeColors = [-1,-1,-1;
%                               -1,-1,-1;
%                               -1,-1,-1];
            else
%                 seriesColors = [kbr1994Color;...
%                                 1,1,1;...
%                                 1,1,1];
                opus31ColorLight    = opus31Color.*(.5) + [1,1,1].*0.5;
                hillColorLight      = hillColor.*(.5) + [1,1,1].*0.5;
                
                seriesColors = [kbr1994Color;...
                                opus31ColorLight;...
                                hillColorLight];

                    
                edgeColors = [kbr1994Color;...
                                opus31Color;...
                                hillColor];
                
            end
    
            x0 = (tableAmpRow(1,idxTrial)-1)*(6*barWidth+2*skinnySpaceWidth+spaceWidth);
            if(idxTendon == type_rigidTendon)
                x0 = x0+barWidth;
            end
            
            y0 = 0;
    
            
            for k=1:1:length(seriesMean)
                
                x1 = x0+barWidth;
                y1 = seriesMean(1,k);
                if(y1 == noDataCode)
                    y1=0;
                end
                
                if(sum(edgeColors(k,:))<0)
                    fill([x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0],seriesColors(k,:),...
                          'EdgeColor','none');
                    hold on;
                else
                    fill([x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0],seriesColors(k,:),...
                          'EdgeColor',edgeColors(k,:));
                    hold on;
                end

                if(k==2 && idxTendon == type_elasticTendon)  
                    xM = 0.5*(x0+x1);
                    text(xM,-0.1,sprintf('%1.1fmm',amp(1,idxTrial)),...
                        'Units','normalized',...
                        'HorizontalAlignment','center');
                    hold on;
                end
                x0 = x1+barWidth+skinnySpaceWidth;
            end
            set(gca,'xticklabel',{[]});
            xlim([0,1]);
            ylim(yLimVals);
            box off;
    
            if(i==1)
                switch idxType
                    case 1
                        ylabel('Norm. Stiffness ((N/mm)/N)');
                    case 2 
                        ylabel('Norm. Damping ((N/mm)/N)');                
                    case 3
                        ylabel('VAF (\%)');                
                    otherwise
                        assert(0);
                end
            end
            if(j==1)
                title(sprintf('%iHz',freq(idxTrial)));
                hold on;
            end
    
        end
    end
end


set(fig_barChart,'Units','centimeters',...
    'PaperUnits','centimeters',...
    'PaperSize',[pageWidth pageHeight],...
    'PaperPositionMode','manual',...
    'PaperPosition',[0 0 pageWidth pageHeight]);     
set(fig_barChart,'renderer','painters');     
set(gcf,'InvertHardCopy','off')

tendonTag = '';

if(flag_useElasticTendon==1)
    tendonTag = '_ElasticTendon';
else
    tendonTag = '_RigidTendon';
end


print('-dpdf', [plotFolder,'fig_Pub_StiffnessDampingVAFBarChart.pdf']);


success=1;