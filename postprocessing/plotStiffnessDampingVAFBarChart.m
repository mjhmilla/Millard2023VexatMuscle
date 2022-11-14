function [success] = plotStiffnessDampingVAFBarChart( ...
                            opus31Table,...
                            hillTable,...
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
flag_fixedPlotWidth           = 1;
plotWidth                     = 5;
plotWidth43                   = 5;
plotHeight43                  = 3;
plotHeight                    = [];

flag_usingOctave              = flag_usingOctave;
plotConfig;


barWidth   = 0.1;
spaceWidth = 0.05;

subPlotTypes = {'stiffness','damping','vaf'};

for idxType = 1:1:length(subPlotTypes)

    row=nan;
    switch subPlotTypes{idxType}
        case 'stiffness'
            row=1;            
        case 'damping'
            row=2;            
        case 'vaf'
            row=3;            
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

        seriesMeans = [];
        seriesMins  = [];
        seriesMax   = [];

        switch subPlotTypes{idxType}
            case 'stiffness'
                seriesMean = [kbr1994Table.stiffness.pMean(i,j,1),...
                               opus31Table.stiffness.pMean(i,j,1),...
                               hillTable.stiffness.pMean(i,j,1)];
                seriesMin = [kbr1994Table.stiffness.p95CIMin(i,j,1),...
                               opus31Table.stiffness.p95CIMin(i,j,1),...
                               hillTable.stiffness.p95CIMin(i,j,1)];
                seriesMax = [kbr1994Table.stiffness.p95CIMax(i,j,1),...
                               opus31Table.stiffness.p95CIMax(i,j,1),...
                               hillTable.stiffness.p95CIMax(i,j,1)];                                
            case 'damping'
                seriesMean = [kbr1994Table.damping.pMean(i,j,1),...
                               opus31Table.damping.pMean(i,j,1),...
                               hillTable.damping.pMean(i,j,1)];
                seriesMin = [kbr1994Table.damping.p95CIMin(i,j,1),...
                               opus31Table.damping.p95CIMin(i,j,1),...
                               hillTable.damping.p95CIMin(i,j,1)];
                seriesMax = [kbr1994Table.damping.p95CIMax(i,j,1),...
                               opus31Table.damping.p95CIMax(i,j,1),...
                               hillTable.damping.p95CIMax(i,j,1)];                
            case 'vaf'
                 seriesMean = [kbr1994Table.vaf.pMean(i,j,2),...
                               opus31Table.vaf.pMean(i,j,2),...
                               hillTable.vaf.pMean(i,j,2)];
                seriesMin = [kbr1994Table.vaf.p95CIMin(i,j,2),...
                               opus31Table.vaf.p95CIMin(i,j,2),...
                               hillTable.vaf.p95CIMin(i,j,2)];
                seriesMax = [kbr1994Table.vaf.p95CIMax(i,j,2),...
                               opus31Table.vaf.p95CIMax(i,j,2),...
                               hillTable.vaf.p95CIMax(i,j,2)];  
                seriesMean  = seriesMean.*100;
                seriesMin   = seriesMean.*100;
                seriesMax   = seriesMean.*100;                
            otherwise
                assert(0,'Type not in list');            
        end

        seriesColors = [kbr1994Color;...
                        opus31Color;...
                        hillColor];

        x0 = (tableAmpRow(1,idxTrial)-1)*(3*barWidth+spaceWidth);
        y0 = 0;

        
        for k=1:1:length(seriesMean)
            x1 = x0+barWidth;
            y1 = seriesMean(1,k);
            if(y1 == noDataCode)
                y1=0;
            end
            fill([x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0],seriesColors(k,:),...
                  'EdgeColor','none');
            hold on;
            if(k==2)  
                xM = 0.5*(x0+x1);
                text(xM,-0.1,sprintf('%1.1fmm',amp(1,idxTrial)),...
                    'Units','normalized',...
                    'HorizontalAlignment','center');
                hold on;
            end
            x0 = x1;
        end
        set(gca,'xticklabel',{[]});
        xlim([0,1]);
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


print('-dpdf', [plotFolder,'fig_Pub_StiffnessDampingVAFBarChart',...
                tendonTag,'.pdf']);


success=1;