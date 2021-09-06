function figH = plotHerzogLeonardDescendingLimb(figH,...
                  expConfigHerzogLeonard2002,dataOpus31, dataDampedEq,  ...
                  nominalNormalizedFiberLength,flag_useElasticTendon,...
                  figureNumber,subFigureNumber,trialNumber,...
                  subPlotPosition, outputFolder, fileName,...
                  pageWidth,pageHeight)
                
success=0;
%figH=figure;

figure(figH);

lineColorExp   = [1,1,1].*0.75;
fillColorExp   = [1,1,1].*0.9;

lineColorModel = [0,0,1];
lineColorHill  = [1,0,0];
lineWidthModel = 1;
lineWidthHill  = 1;
lineTypeModel = '-';
lineTypeHill = '-';

flag_addNormFiberLengthLabel = 0;

nameModel = 'Model: ET';
nameHill  = 'Hill: ET';

if(flag_useElasticTendon==0)
  lineColorModel = lineColorModel.*0.25 + [1,1,1].*0.75;
  lineColorHill  = lineColorHill.*0.25 + [1,1,1].*0.75;  
  lineWidthModel = 1.5;
  lineWidthHill  = 1.5;  
  nameModel = 'Model: RT';
  nameHill  = 'Hill: RT';  
  lineTypeModel = '-';
  lineTypeHill = '-';
  
end

expLineWidth = 1;



subPlotTitleFontSize = 8*1.2;
%%
% Extract the key times from the data
%%

ta0=expConfigHerzogLeonard2002.stimulationKeyTimes(1,1);
ta1=expConfigHerzogLeonard2002.stimulationKeyTimes(1,2);

tr0=expConfigHerzogLeonard2002.lengthRampKeyPoints(1,1);
tr1=expConfigHerzogLeonard2002.lengthRampKeyPoints(2,1);

t0 = expConfigHerzogLeonard2002.timeSpan(1,1);
t1 = expConfigHerzogLeonard2002.timeSpan(1,2);

tvec =round([t0,ta0,tr0,tr1,ta1,t1],1);

f0 = interp1(expConfigHerzogLeonard2002.dataRamp.time,...
             expConfigHerzogLeonard2002.dataRamp.force,...
             t0);
f1 = interp1(expConfigHerzogLeonard2002.dataRamp.time,...
             expConfigHerzogLeonard2002.dataRamp.force,...
             t1);
fp1 = interp1(expConfigHerzogLeonard2002.dataPassive.time,...
             expConfigHerzogLeonard2002.dataPassive.force,...
             t1);
           
fr0 = interp1(expConfigHerzogLeonard2002.dataRamp.time,...
             expConfigHerzogLeonard2002.dataRamp.force,...
             tr0-0.05);

fr1 = interp1(expConfigHerzogLeonard2002.dataRamp.time,...
             expConfigHerzogLeonard2002.dataRamp.force,...
             tr1);
fa1 = interp1(expConfigHerzogLeonard2002.dataRamp.time,...
             expConfigHerzogLeonard2002.dataRamp.force,...
             ta1);
           
           
           
fvec = round([f0,fp1,f1,fr0,fr1],1);
           
%%
% Length
%%

idxSubplot = 1;
subplot('Position',reshape(subPlotPosition(idxSubplot,1,:),1,4));  


if(flag_useElasticTendon==0)

  plot( expConfigHerzogLeonard2002.dataRamp.time,...
        expConfigHerzogLeonard2002.dataRamp.length,...
        'Color',lineColorExp,'LineWidth',expLineWidth);
  hold on;

  box off;

end

dyNorm = 0.5/pageHeight;


xTicksVector = [];

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
xlabel('Time (s)');

ylim([-0.5,9.5]);
xticks(tvec);
yticks([0,9]);

%set(gca,'XTickLabel',[]);
x01 = xlim;
y01 = ylim;

dyPlot = dyNorm/subPlotPosition(idxSubplot,1,4);  
dy01 = (dyPlot)*diff(y01);

if(flag_useElasticTendon==1)
  text(x01(1,1),y01(1,2),'A. Length profile',...
    'FontSize',subPlotTitleFontSize,...
    'HorizontalAlignment','left',...
    'VerticalAlignment','bottom');

  hold on;
end

%if(flag_useElasticTendon==1)
if(flag_useElasticTendon==1)
  text(mean(x01),y01(1,2)+1.75*dy01,'Simulation of Herzog \& Leonard 2002',...
    'FontSize',8*1.2,...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom');
  hold on;
end
%   text(mean(x01),y01(1,2)+dy01,'(models with elastic tendons)',...
%     'FontSize',8*1.2,...
%     'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom');
%   hold on;
%   
% else
%   text(mean(x01),y01(1,2)+1.75*dy01,'Simulation of Herzog \& Leonard 2002',...
%     'FontSize',8*1.2,...
%     'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom');
%   hold on;  
%   text(mean(x01),y01(1,2)+dy01,'(models with rigid tendons)',...
%     'FontSize',8*1.2,...
%     'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom');
%   hold on;  
% end

% text(mean(x01),y01(1,2)+dy01,'Simulation of Herzog and Leonard 2002',...
%   'FontSize',8*1.2,...
%   'HorizontalAlignment','center',...
%   'VerticalAlignment','bottom');
% hold on;

normLength31r0 = interp1(dataOpus31.benchRecord.time(:,1),...
                         dataOpus31.benchRecord.normFiberLength(:,1),tr0);

normLength31r1 = interp1(dataOpus31.benchRecord.time(:,1),...
                         dataOpus31.benchRecord.normFiberLength(:,1),tr1);

normLengthHr0 = interp1( dataDampedEq.benchRecord.time(:,1),...
                         dataDampedEq.benchRecord.normFiberLength(:,1),tr0);

normLengthHr1 = interp1(dataDampedEq.benchRecord.time(:,1),...
                         dataDampedEq.benchRecord.normFiberLength(:,1),tr1);
                       
                  
normLength0 = 0.5*(normLength31r0+normLengthHr0);
normLength1 = 0.5*(normLength31r1+normLengthHr1);

normDiffLength0 = 0.5*abs(normLength31r0-normLengthHr0);
normDiffLength1 = 0.5*abs(normLength31r1-normLengthHr1);

strNormDiffLength0 = num2str(round(normDiffLength0,4));
strNormDiffLength1 = num2str(round(normDiffLength1,4));

strNormLength0 = num2str(round(normLength0,2));
strNormLength1 = num2str(round(normLength1,2));


if(flag_addNormFiberLengthLabel==1)
  text(tr0+1,0,['Sim. $',strNormLength0,'\pm',strNormDiffLength0,'\, \ell^{M}_{\circ}$'],...
              'FontSize',6,...
              'HorizontalAlignment','left',...
              'VerticalAlignment','bottom');
  hold on;
  text(tr1+1,9,['Sim. $',strNormLength1,'\pm',strNormDiffLength1,'\, \ell^{M}_{\circ}$'],...
              'FontSize',6,...
              'HorizontalAlignment','left',...
              'VerticalAlignment','top');

  hold on;
end
%title(['Herzog and Leonard 2002: ',num2str(figureNumber),...
%        subFigTitle,num2str(trialNumber)]);  

%%
% Force
%%
      
      
idxSubplot = idxSubplot+1;
subplot('Position',reshape(subPlotPosition(idxSubplot,1,:),1,4));        

lineExp     = [];
lineOpus31  = [];
lineHill    = [];

if(flag_useElasticTendon==0)

  fill( expConfigHerzogLeonard2002.dataStatic.time,...
        expConfigHerzogLeonard2002.dataStatic.force,...
        fillColorExp,'EdgeColor','none');%);%,'DisplayName',[]);
  hold on;

  text( ta1-0.1,fr0*0.95,...
        'Fixed-length activation (9mm)',...
        'FontSize',6,...
        'HorizontalAlignment','right',...
        'VerticalAlignment','top');

  plot( expConfigHerzogLeonard2002.dataPassive.time,...
        expConfigHerzogLeonard2002.dataPassive.force,...
        '-','Color',[1,1,1],'LineWidth',expLineWidth*2);%);%,'DisplayName',[]);
  hold on;

  text(t1,0.8*max(expConfigHerzogLeonard2002.dataPassive.force),...
       'Passive lengthening (0-9mm)',...
       'FontSize',6,...
       'HorizontalAlignment','right',...
       'VerticalAlignment','top');

  plot( expConfigHerzogLeonard2002.dataPassive.time,...
        expConfigHerzogLeonard2002.dataPassive.force,...
        '--','Color',lineColorExp,'LineWidth',expLineWidth);%,...
%        'DisplayName',[]);
  hold on;

  expName ='Exp.';
  if(flag_useElasticTendon==1)
    expName = 'none';
  end

  lineExp = plot( expConfigHerzogLeonard2002.dataRamp.time,...
                  expConfigHerzogLeonard2002.dataRamp.force,...
                  'Color',lineColorExp,'LineWidth',expLineWidth,...
                  'DisplayName',expName);
  hold on;

end



idxMapH = [];
switch(size(dataDampedEq.benchRecord.time,2))
  case 1
    idxMapH=[1];
  case 2
    idxMapH=[1,2];    
  case 3
    idxMapH=[1,2];    
  otherwise
    assert(0);
end

for i=1:1:size(idxMapH,2)


  if(i==1)
    lineHill = plot(dataDampedEq.benchRecord.time(:,idxMapH(1,i)), ...
                    dataDampedEq.benchRecord.tendonForce(:,idxMapH(1,i)),...
                   lineTypeHill,'Color',lineColorHill,'LineWidth',lineWidthHill,...
                   'DisplayName',nameHill);
  else
    plot(dataDampedEq.benchRecord.time(:,idxMapH(1,i)), ...
                    dataDampedEq.benchRecord.tendonForce(:,idxMapH(1,i)),...
                   lineTypeHill,'Color',lineColorHill,'LineWidth',lineWidthHill,...
                   'DisplayName',nameHill);
  end
  hold on;
end


idxMap31 = [];
switch(size(dataOpus31.benchRecord.time,2))
  case 1
    idxMap31=[1];
  case 2
    idxMap31=[1,2];    
  case 3
    idxMap31=[1,3];    
  otherwise
    assert(0);
end

for i=1:1:size(idxMap31,2)
  if(i==1)
    lineOpus31 = plot(dataOpus31.benchRecord.time(:,idxMap31(1,i)), ...
                      dataOpus31.benchRecord.tendonForce(:,idxMap31(1,i)),...
                      lineTypeModel,'Color',lineColorModel,'LineWidth',lineWidthModel,...
                      'DisplayName',nameModel);      
  else
    plot(dataOpus31.benchRecord.time(:,idxMap31(1,i)), ...
                      dataOpus31.benchRecord.tendonForce(:,idxMap31(1,i)),...
                      lineTypeModel,'Color',lineColorModel,'LineWidth',lineWidthModel,...
                      'DisplayName',nameModel);
  end
  hold on;
end


expLabel=['Exp. Fig.',num2str(figureNumber),...
         subFigTitle,num2str(trialNumber)];
%expLabel = ['Exp. HL2002 ',num2str(figureNumber),...
%        subFigTitle,num2str(trialNumber)];



box off;

%ylim([0,40]);
if(flag_useElasticTendon==1)
  xlim([0,max(expConfigHerzogLeonard2002.dataRamp.time)]);
  xlabel('Time (s)');
  ylabel('Tendon Force (N)');

  xticks(tvec);
  yticks(fvec);
end

x01 = xlim;
y01 = ylim;


dyPlot = dyNorm/subPlotPosition(idxSubplot,1,4);  
dy01 = (dyPlot)*diff(y01);

if(flag_useElasticTendon==1)
  text(x01(1,1),y01(1,2),...
    'B. Active lengthening force profile',...
    'FontSize',subPlotTitleFontSize,...
    'HorizontalAlignment','left',...
    'VerticalAlignment','bottom');
  hold on;
end


%%
% Stiffness
%%

idxSubplot = idxSubplot+1;
subplot('Position',reshape(subPlotPosition(idxSubplot,1,:),1,4));        
    
plot(expConfigHerzogLeonard2002.timeSpan,[0,0],'k','LineWidth',0.5);%,...
     %'DisplayName',[]);
hold on;

i=1;

%for i=1:1:size(idxMap31,2)

  plot(dataOpus31.benchRecord.time(:,idxMap31(1,i)), ...
      dataOpus31.benchRecord.musculotendonStiffness(:,idxMap31(1,i))./1000,...
      lineTypeModel,'Color',lineColorModel,'LineWidth',lineWidthModel,...
      'DisplayName',nameModel);
  hold on;
%end


%for i=1:1:size(idxMapH,2)


  plot(dataDampedEq.benchRecord.time(:,idxMapH(1,i)), ...
                  dataDampedEq.benchRecord.musculotendonStiffness(:,idxMapH(1,i))./1000,...
                 lineTypeHill,'Color',lineColorHill,'LineWidth',lineWidthHill,...
                 'DisplayName',nameHill);
  hold on;
%end


box off;

k31r1 = interp1(dataOpus31.benchRecord.time(:,idxMap31(1,i)),...
              dataOpus31.benchRecord.musculotendonStiffness(:,idxMap31(1,i)),...
              tr1)/1000;
 
kHr1 = interp1(dataDampedEq.benchRecord.time(:,idxMapH(1,i)),...
               dataDampedEq.benchRecord.musculotendonStiffness(:,idxMapH(1,i)),...
              tr1)/1000;
kH1 = interp1(dataDampedEq.benchRecord.time(:,idxMapH(1,i)),...
               dataDampedEq.benchRecord.musculotendonStiffness(:,idxMapH(1,i)),...
              t1)/1000;
            
kH  = max(dataDampedEq.benchRecord.musculotendonStiffness(:,idxMapH(1,i))./1000);


% plot(tr1,k31r1,'o',...
%   'Color',lineColorModel,'MarkerFaceColor',[1,1,1],...
%   'MarkerSize',4);
% 
% hold on;
% text(tr1,k31r1,num2str(round(k31r1,2)),'FontSize',6,...
%      'VerticalAlignment','bottom',...
%      'HorizontalAlignment','left');
% hold on;

% plot(tr1,kHr1,'o','Color',lineColorHill,...
%      'MarkerFaceColor',[1,1,1],...
%      'MarkerSize',4);
% hold on;
% text(tr1,kHr1,num2str(round(kHr1,2)),'FontSize',6,...
%      'VerticalAlignment','top',...
%      'HorizontalAlignment','left');
% hold on;


kMax = max(dataOpus31.benchRecord.musculotendonStiffness(:,idxMap31(1,i))./1000);
[kMin, idxKmin] = min(dataDampedEq.benchRecord.musculotendonStiffness(:,idxMapH(1,i))./1000);

tkMin = dataDampedEq.benchRecord.time(idxKmin,idxMapH(1,i));

plot([tkMin;tkMin],[kMin+2;kMin+0.5],'k','LineWidth',0.5);%,'DisplayName',[]);
hold on;
plot([tkMin],[kMin+0.5],'v','Color',[0,0,0],'LineWidth',0.5,...
      'MarkerSize',2,'MarkerFaceColor',[0,0,0]);%,'DisplayName',[]);
hold on;
if(flag_useElasticTendon==1)
  text(tkMin, kMin+2.5, 'Negative stiffness!');
end
kvec = round([kMin,0,kH1,kMax],1);
kDelta = (kMax-kMin).*0.05;

if(flag_useElasticTendon==0)
  ylim([(kMin-kDelta),(kMax+kDelta)]);
  yticks(kvec);
else
  yt = yticks;
  %yts = sort([yt(:);round(kMax,1)]);
  %yticks(yts);
end

xlim([0,max(expConfigHerzogLeonard2002.dataRamp.time)]);
xticks(tvec);

if(flag_useElasticTendon==1)
  xlabel('Time (s)');
  ylabel('Stiffness (N/mm)');
  xticks(tvec);
end

x01 = xlim;
y01 = ylim;


dyPlot = dyNorm/subPlotPosition(idxSubplot,1,4);  
dy01 = (dyPlot)*diff(y01);

if(flag_useElasticTendon==1)
  text(x01(1,1),y01(1,2),...
    'C. Stiffness profile',...
    'FontSize',subPlotTitleFontSize,...
    'HorizontalAlignment','left',...
    'VerticalAlignment','bottom');
  hold on;
end


%%
% Damping
%%

idxSubplot = idxSubplot+1;
subplot('Position',reshape(subPlotPosition(idxSubplot,1,:),1,4));        

% if(flag_useElasticTendon==0)
%   plot(expConfigHerzogLeonard2002.timeSpan,[0,0],'k','LineWidth',0.5);%,...
%        %'DisplayName',[]);
%   hold on;
% end
%for i=1:1:size(idxMap31,2)
i=1;
  

  plot(dataDampedEq.benchRecord.time(:,idxMapH(1,i)), ...
                  dataDampedEq.benchRecord.musculotendonDamping(:,idxMapH(1,i))./1000,...
                 lineTypeHill,'Color',lineColorHill,'LineWidth',lineWidthHill,...
                 'DisplayName',nameHill);
  hold on;

  plot(dataOpus31.benchRecord.time(:,idxMap31(1,i)), ...
      dataOpus31.benchRecord.musculotendonDamping(:,idxMap31(1,i))./1000,...
      lineTypeModel,'Color',lineColorModel,'LineWidth',lineWidthModel,...
      'DisplayName',nameModel);
  hold on;
%end

%for i=1:1:size(idxMapH,2)


%end

d310 = interp1(dataOpus31.benchRecord.time(:,idxMap31(1,i)),...
              dataOpus31.benchRecord.musculotendonDamping(:,idxMap31(1,i)),...
              t0)/1000;
d31r1 = interp1(dataOpus31.benchRecord.time(:,idxMap31(1,i)),...
              dataOpus31.benchRecord.musculotendonDamping(:,idxMap31(1,i)),...
              tr1)/1000;
 
dHr1 = interp1(dataDampedEq.benchRecord.time(:,idxMapH(1,i)),...
               dataDampedEq.benchRecord.musculotendonDamping(:,idxMapH(1,i)),...
              tr1)/1000;
            
dH  = max(dataDampedEq.benchRecord.musculotendonDamping(:,idxMapH(1,i))./1000);


% plot(tr1,d31r1,'o',...
%   'Color',lineColorModel,'MarkerFaceColor',[1,1,1],...
%   'MarkerSize',4);

% hold on;
% text(tr1,d31r1,num2str(round(d31r1,2)),'FontSize',6,...
%      'VerticalAlignment','bottom',...
%      'HorizontalAlignment','left');
% hold on;
% 
% plot(tr1,dHr1,'o','Color',lineColorHill,...
%      'MarkerFaceColor',[1,1,1],...
%      'MarkerSize',4);
% hold on;
% text(tr1,dHr1,num2str(round(dHr1,2)),'FontSize',6,...
%      'VerticalAlignment','top',...
%      'HorizontalAlignment','left');
% hold on;

dvec = round(sort([0,d31r1,d310,dH]),2);
dDelta = (dH).*0.05;

box off;

if(flag_useElasticTendon==0)
  ylim([(0-dDelta),(dH+dDelta)]);
  yticks(dvec);
end

%ylim([0,40]);
xlim([0,max(expConfigHerzogLeonard2002.dataRamp.time)]);

if(flag_useElasticTendon==1)
  xlabel('Time (s)');
  ylabel('Damping (N/(mm/s))');
else
  xticks(tvec);  
end

x01 = xlim;
y01 = ylim;

dyPlot = dyNorm/subPlotPosition(idxSubplot,1,4);  
dy01 = (dyPlot)*diff(y01);

if(flag_useElasticTendon==1)
  text(x01(1,1),y01(1,2),...
     'D. Damping profile',...
     'FontSize',subPlotTitleFontSize,...
    'HorizontalAlignment','left',...
    'VerticalAlignment','bottom');
  hold on;
end

if(trialNumber == 3 && flag_useElasticTendon==1)
  %legend([lineOpus31,lineHill,lineExp],'Model', 'Hill-type',expLabel);
  lh = legend('Location','NorthEast');
  
  lh.Position(1,2) = lh.Position(1,2)*1.2;
  legend boxoff;
end


%%
%
%%
if(flag_useElasticTendon==1)
  figure(figH);
  set(figH,'Units','centimeters',...
  'PaperUnits','centimeters',...
  'PaperSize',[pageWidth pageHeight],...
  'PaperPositionMode','manual',...
  'PaperPosition',[0 0 pageWidth pageHeight]);     
  %set(findall(figList(i).h,'-property','FontSize'),'FontSize',10);     
  set(figH,'renderer','painters');     
  set(gcf,'InvertHardCopy','off')



  print('-dpdf', [outputFolder,fileName]);                     
  %benchRecord.musculotendonStiffness
  %benchRecord.musculotendonDamping
end
success=1;
