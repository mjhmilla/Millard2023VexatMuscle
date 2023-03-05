function figH = plotHerzogLeonardFig7Comparision(...
                  expConfigHerzogLeonard2002,dataOpus31, dataDampedEq,  ...
                  figureNumber,subFigureNumber,trialNumber,...
                  figH,subPlotPairPanel)
               

figure(figH);

lineColorExp   = [1,1,1].*0.75;
fillColorExp   = [1,1,1].*0.9;
lineColorModel = [0,0,1];
lineColorHill  = [1,0,0];

expLineWidth = 1;

idxSubplot = (trialNumber-1)*2 + 1;
subplot('Position',reshape(subPlotPairPanel(idxSubplot,1,:),1,4));  

plot( expConfigHerzogLeonard2002.dataRamp.time,...
      expConfigHerzogLeonard2002.dataRamp.length,...
      'Color',[1,1,1].*0.75,'LineWidth',2);
hold on;

box off;

subFigTitle = '';
switch subFigureNumber
  case 1
    subFigTitle = 'A';
  case 2
    subFigTitle = 'B';
  case 3
    subFigTitle = 'C';
  otherwise
    assert(0);
end  
ylabel('Length (mm)');
set(gca,'XTickLabel',[]);

title(['Herzog and Leonard 2002: ',num2str(figureNumber),...
        subFigTitle,num2str(trialNumber)]);  


idxSubplot = idxSubplot+1;
subplot('Position',reshape(subPlotPairPanel(idxSubplot,1,:),1,4));  


%n=(1-1)/2;
%lineColor = [0,0,0.1].*(1-n) + [1,1,1].*(n*0.75);
%lineColor = lineColor.*0.5 + [1,1,1].*0.5;  

plot( expConfigHerzogLeonard2002.dataRamp.time,...
      expConfigHerzogLeonard2002.dataRamp.force,...
      'Color',lineColorExp,'LineWidth',expLineWidth);
hold on;

%n=(2-1)/2;
%lineColor = [0,0,0.1].*(1-n) + [1,1,1].*(n*0.75);
%lineColor = lineColor.*0.5 + [1,1,1].*0.5;  
%plot( expConfigHerzogLeonard2002.dataStatic.time,...
%      expConfigHerzogLeonard2002.dataStatic.force,...
%      'Color',[1,1,1].*0.75,'LineWidth',2);
%hold on;

%n=(3-1)/2;
%lineColor = [0,0,0.1].*(1-n) + [1,1,1].*(n*0.75);
%lineColor = lineColor.*0.5 + [1,1,1].*0.5;  

if(trialNumber==3)
  plot( expConfigHerzogLeonard2002.dataPassive.time,...
        expConfigHerzogLeonard2002.dataPassive.force,...
        'Color',lineColorExp,'LineWidth',expLineWidth);
  hold on;
end  
lineOpus31 = [];
lineHill = [];

if(isempty(dataOpus31.benchRecord)==0)
    for i=1:1:size(dataOpus31.benchRecord.time,2)
      n=0;
      if(size(dataOpus31.benchRecord.time,2)>1)
        n = (i-1)/(size(dataOpus31.benchRecord.time,2)-1);
      end
      lineColor = [1,0,0].*(1-n) + [0.75,0.25,0.25].*(n);
      if(i==1)
        lineOpus31 = plot(dataOpus31.benchRecord.time(:,i), ...
                          dataOpus31.benchRecord.tendonForce(:,i),...
                          'Color',lineColorModel,'LineWidth',1);      
      else
        plot(dataOpus31.benchRecord.time(:,i), ...
                          dataOpus31.benchRecord.tendonForce(:,i),...
                          'Color',lineColorModel,'LineWidth',1);
      end
      hold on;
    end
end

if(isempty(dataDampedEq.benchRecord)==0)
    for i=1:1:size(dataDampedEq.benchRecord.time,2)
      n=0;
      if(size(dataDampedEq.benchRecord.time,2)>1)
        n = (i-1)/(size(dataDampedEq.benchRecord.time,2)-1);
      end
      lineColor = [0,0,1].*(1-n) + [0.25,0.25,0.75].*(n);
      if(i==1)
        lineHill = plot(dataDampedEq.benchRecord.time(:,i), ...
                        dataDampedEq.benchRecord.tendonForce(:,i),...
                       '-','Color',lineColorHill,'LineWidth',1);
      else
        plot(dataDampedEq.benchRecord.time(:,i), ...
                        dataDampedEq.benchRecord.tendonForce(:,i),...
                       '-','Color',lineColorHill,'LineWidth',1);
      end
      hold on;
    end
end
if(trialNumber == 3)
  legend([lineOpus31,lineHill],'Model', 'Hill-type');
  legend boxoff;
end

box off;

ylim([0,40]);
xlim([0,max(expConfigHerzogLeonard2002.dataRamp.time)]);
xlabel('Time (s)');
ylabel('Tendon Force (N)');
                  