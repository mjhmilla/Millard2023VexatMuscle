function [success] = plotStiffnessDampingVAFBarChart( ...
                            opus31Table_RT,...
                            hillTable_RT,...
                            opus31Table_ET,...
                            hillTable_ET,...
                            kbr1994Table,...                            
                            opus31Color,...
                            hillColor,...
                            kbr1994Color,...
                            inputFunctions,...
                            strFittingBandwidth,...
                            plotFolder,...
                            flag_usingOctave)

success=0;

freqTable = [15,35,90];
ampTable = [0.4,0.8,1.6];
noDataCode = -42;


%tableFreqCol=[  1,   1,   1,   2,   2,   2,   3,   3,   3];
%tableAmpRow =[  1,   2,   3,   1,   2,   3,   1,   2,   3];

fig_barChart = figure;

numberOfHorizontalPlotColumns = 3;
numberOfVerticalPlotRows      = 3;
plotHorizMarginCm             = 2;
plotVertMarginCm              = 2;
plotHorizMarginSplitCm        = 1;
plotVertMarginSplitCm         = 3;
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
spaceWidth = 0.1;

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
                yLimVals = [0.0;4];
            case 'damping'
                if(idxTendon==type_elasticTendon)
                    row=2;            
                else
                    row=2;
                end
                yLimVals = [0.000;0.05];
                
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
    
        for idxTrial=1:1:length(inputFunctions.bandwidthHz)    
            col = nan;
            titleLabel = '';
            switch inputFunctions.bandwidthHz(idxTrial)
                case 15
                    col=1;
                    titleLabel = 'A. ';
                case 35
                    col=2;
                    titleLabel = 'B. ';                    
                case 90
                    col=3;
                    titleLabel = 'C. ';                    
                otherwise
                    assert(0,'Frequency not in list');
            end
            subplot('Position', reshape(subPlotPanel(row,col,:),1,4) );
    
            %i = tableFreqCol(1,idxTrial);
            %j = tableAmpRow(1,idxTrial);
            if(idxTrial==4)
                here=1;
            end

            [val, idxAmp] = min(abs(ampTable-inputFunctions.amplitudeMM(1,idxTrial)));
            [val, idxFreq] = min(abs(freqTable-inputFunctions.bandwidthHz(1,idxTrial)));            
            k = 1;
            
            seriesMean = [];
            seriesMin  = [];
            seriesMax   = [];
            
            fieldName = subPlotTypes{idxType};

            if(idxTendon == type_elasticTendon)
                if(contains(fieldName,'vaf')==1)
                    seriesMean = [ mean(kbr1994Table.(fieldName).data(idxTrial).y),...
                                   mean(opus31Table_ET.(fieldName).data(idxTrial).y),...
                                   mean(hillTable_ET.(fieldName).data(idxTrial).y)];
                    seriesMin =  [ min(kbr1994Table.(fieldName).data(idxTrial).y),...
                                   min(opus31Table_ET.(fieldName).data(idxTrial).y),...
                                   min(hillTable_ET.(fieldName).data(idxTrial).y) ];
                    seriesMax =  [ max(kbr1994Table.(fieldName).data(idxTrial).y),...
                                   max(opus31Table_ET.(fieldName).data(idxTrial).y),...
                                   max(hillTable_ET.(fieldName).data(idxTrial).y )];
                else
                    seriesMean = [ mean(kbr1994Table.(fieldName).data(idxTrial).yN),...
                                   mean(opus31Table_ET.(fieldName).data(idxTrial).yN),...
                                   mean(hillTable_ET.(fieldName).data(idxTrial).yN)];
                    seriesMin =  [ min(kbr1994Table.(fieldName).data(idxTrial).yN),...
                                   min(opus31Table_ET.(fieldName).data(idxTrial).yN),...
                                   min(hillTable_ET.(fieldName).data(idxTrial).yN) ];
                    seriesMax =  [ max(kbr1994Table.(fieldName).data(idxTrial).yN),...
                                   max(opus31Table_ET.(fieldName).data(idxTrial).yN),...
                                   max(hillTable_ET.(fieldName).data(idxTrial).yN )];
                end
            else
                if(contains(fieldName,'vaf')==1)
                    seriesMean = [ mean(kbr1994Table.(fieldName).data(idxTrial).y),...
                                   mean(opus31Table_RT.(fieldName).data(idxTrial).y),...
                                   mean(hillTable_RT.(fieldName).data(idxTrial).y)];
                    seriesMin =  [ min(kbr1994Table.(fieldName).data(idxTrial).y),...
                                   min(opus31Table_RT.(fieldName).data(idxTrial).y),...
                                   min(hillTable_RT.(fieldName).data(idxTrial).y) ];
                    seriesMax =  [ max(kbr1994Table.(fieldName).data(idxTrial).y),...
                                   max(opus31Table_RT.(fieldName).data(idxTrial).y),...
                                   max(hillTable_RT.(fieldName).data(idxTrial).y )];
                else
                    seriesMean = [ mean(kbr1994Table.(fieldName).data(idxTrial).yN),...
                                   mean(opus31Table_RT.(fieldName).data(idxTrial).yN),...
                                   mean(hillTable_RT.(fieldName).data(idxTrial).yN)];
                    seriesMin =  [ min(kbr1994Table.(fieldName).data(idxTrial).yN),...
                                   min(opus31Table_RT.(fieldName).data(idxTrial).yN),...
                                   min(hillTable_RT.(fieldName).data(idxTrial).yN) ];
                    seriesMax =  [ max(kbr1994Table.(fieldName).data(idxTrial).yN),...
                                   max(opus31Table_RT.(fieldName).data(idxTrial).yN),...
                                   max(hillTable_RT.(fieldName).data(idxTrial).yN )];                    

                end
            end

            idxNoData = find(seriesMean == noDataCode);
            seriesMean(idxNoData) = nan;
            idxNoData = find(seriesMax == noDataCode);
            seriesMax(idxNoData) = nan;
            idxNoData = find(seriesMin == noDataCode);
            seriesMin(idxNoData) = nan;


            if(length(seriesMin)==2)
                here=1;
            end

            if(contains(fieldName,'vaf'))
                seriesMean  = seriesMean.*100;
                seriesMin   = seriesMin.*100;
                seriesMax   = seriesMax.*100;                
            end

            seriesColors = [];
            edgeColors = [];
            seriesLabels = [];

            if(idxTendon == type_elasticTendon)
                seriesColors = [kbr1994Color;...
                                opus31Color;...
                                hillColor];
                edgeColors = [kbr1994Color;...
                                opus31Color;...
                                hillColor];
                seriesLabels = {'Kirsch 1994',...
                                'Model',...
                                'Hill'};
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

                seriesLabels = {'Kirsch 1994',...
                                'Model (RT)',...
                                'Hill (RT)'};
                
                
            end
    
            x0 = (idxAmp-1)*(6*barWidth+2*skinnySpaceWidth+spaceWidth);
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

                    pid = fill([x0,x1,x1,x0,x0],...
                               [y0,y0,y1,y1,y0],...
                               seriesColors(k,:),...
                               'EdgeColor','none',...
                               'DisplayName',seriesLabels{k});
                    hold on;
                else
                    pid = fill([x0,x1,x1,x0,x0],...
                               [y0,y0,y1,y1,y0],...
                               seriesColors(k,:),...
                               'EdgeColor',edgeColors(k,:),...
                               'DisplayName',seriesLabels{k});
                    hold on;
                end

                if( (idxType == 3 && idxTrial == 6)==0)
                    set(get(get(pid,'Annotation'),...
                        'LegendInformation'),...
                        'IconDisplayStyle','off');
                end

                lineColors = edgeColors(k,:).*0.5 + [0,0,0].*0.5;
                yMin = seriesMin(1,k);
                yMax = seriesMax(1,k);
                xM = 0.5*(x0+x1);
                xL = xM - 0.5*(x1-x0);
                xR = xM + 0.5*(x1-x0);

                annotationColor = edgeColors(k,:).*0.5 + [1,1,1].*0.5;

                pid=plot([xM;xM],[yMin;y1],'Color',annotationColor);
                hold on;
                set(get(get(pid,'Annotation'),...
                        'LegendInformation'),...
                        'IconDisplayStyle','off');


                pid=plot([xM;xM],[y1;yMax],'Color',annotationColor);
                hold on;
                set(get(get(pid,'Annotation'),...
                        'LegendInformation'),...
                        'IconDisplayStyle','off');

                pid=plot([xL;xR],[yMin;yMin],'Color',annotationColor);
                hold on;
                set(get(get(pid,'Annotation'),...
                        'LegendInformation'),...
                        'IconDisplayStyle','off');
                
                pid=plot([xL;xR],[yMax;yMax],'Color',annotationColor);
                hold on;
                set(get(get(pid,'Annotation'),...
                        'LegendInformation'),...
                        'IconDisplayStyle','off');

                if(k==2 && idxTendon == type_elasticTendon)  
                    text(xM,-0.1,...
                        sprintf('%1.1fmm',inputFunctions.amplitudeMM(idxTrial)),...
                        'Units','normalized',...
                        'HorizontalAlignment','center');
                    hold on;
                end
                x0 = x1+barWidth+skinnySpaceWidth;
            end
            set(gca,'xticklabel',{[]});
           %xlim([0,1]);
            axis tight;
            ylim(yLimVals);
            box off;

            if(idxType == 3 && idxTrial == 9 && idxTendon == type_rigidTendon)
                lh=legend('Location','South');
                lh.Position(1,1) = lh.Position(1,1);% -0.15625;        
                lh.Position(1,2) = lh.Position(1,2) -0.2;%40625;   
                lh.NumColumns=6;
                legend boxoff;

            end
            if(idxTrial==1)
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
            if(idxAmp==2)
                title([titleLabel, sprintf('%iHz',inputFunctions.bandwidthHz(idxTrial))]);
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



print('-dpdf', [plotFolder,'fig_Pub_StiffnessDampingVAFBarChart.pdf']);


success=1;