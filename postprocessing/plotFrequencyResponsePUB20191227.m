function [success] = plotFrequencyResponsePUB20191227(...
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
                      nominalForceTargetIndex,...
                      dataKBR1994Fig3Gain,...
                      dataKBR1994Fig3Phase,...                      
                      flag_useElasticTendon,...  
                      flag_Mode15Hz90HzBoth,...
                      plotNameEnding,...
                      plotLayoutSettings,...
                      pubOutputFolder)
                      
success = 0;

%flag_Mode15Hz90HzBoth = 0; %0:15, 1:90, 2:Both
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


numberOfHorizontalPlotColumns = plotLayoutSettings.numberOfHorizontalPlotColumns;
numberOfVerticalPlotRows      = plotLayoutSettings.numberOfVerticalPlotRows;
flag_fixedPlotWidth           = plotLayoutSettings.flag_fixedPlotWidth;
plotWidth                     = plotLayoutSettings.plotWidth;
plotWidth43                   = plotLayoutSettings.plotWidth43;
plotHeight43                  = plotLayoutSettings.plotHeight43;
plotHeight                    = [];

flag_usingOctave              = plotLayoutSettings.flag_usingOctave;
plotConfig;          

figLabelFontSize = 9;

samplePoints    = inputFunctions.samples;  
paddingPoints   = inputFunctions.padding;
sampleFrequency = inputFunctions.sampleFrequency;



fig_freqResponse43 = figure;



switch flag_Mode15Hz90HzBoth
  case 0
    expAmpPlot   = [1.6];
    expBWPlot    = [15];
    expLineWidth = [1];
    expMarkType  = {'-'};
    expLegendEntry = {'Kirsch: 1.6mm 15Hz'};   
    expModelPlotColor = [0.25,  0.25, 1;...
                         0.5,   0.5,  1];
    expPlotColor     = [0.25,  0.25, 0.25;...
                         0.5,   0.5, 0.5];   
    expMarkFaceColor = expPlotColor;
    xTicksVector = [0,4,15];
    
                       
  case 1
    expAmpPlot   = [1.6];
    expBWPlot    = [90];
    expLineWidth = [1];
    expMarkType  = {'-'};
    expLegendEntry = {'Kirsch: 1.6mm 90Hz'};   
    expModelPlotColor = [0.25,  0.25, 1;...
                         0.5,   0.5,  1];
    expPlotColor     = [0.25,  0.25, 0.25;...
                         0.5,   0.5, 0.5];
    expMarkFaceColor = expPlotColor;
    xTicksVector = [0,4,15,90];                       
  case 2
    expAmpPlot   = [1.6,1.6];
    expBWPlot    = [90 ,15];    
    expLineWidth = [1,1];
    expMarkType  = {'-','-'};
    expLegendEntry = {'Kirsch: 1.6mm 90Hz', 'Kirsch: 1.6mm 15Hz'};    
    expModelPlotColor = [0.25,  0.25, 1;...
                         0.5,   0.5,  1];
    expPlotColor      = [0.25,  0.25, 0.25;...
                         0.5,   0.5, 0.5];
    expMarkFaceColor = expPlotColor;
    xTicksVector = [0,4,15,90];
    
  otherwise
    assert(0)
end




chunkDuration = 0.2;        
timeChunkStart = round((paddingPoints)*0.5)/sampleFrequency;  
timeChunkEnd   = paddingPoints/sampleFrequency+chunkDuration;
idxChunk       = [round((paddingPoints)*0.5):1: ...
                  (paddingPoints + sampleFrequency*chunkDuration)];                  
timeTicks = [round(paddingPoints/sampleFrequency,2),...
             round(timeChunkEnd,2)];

subPlotList = zeros(3+3+3+3,4);

ySpace = 0.05;
xSpace = 0.04;

%%
% Plot Layout
%%
figure(fig_freqResponse43);
  currentSubPlot  = subPlotRectangle43;
  currentSubPlot(1,1) = currentSubPlot(1,1)*0.5; 
  subPlotHeight   = subPlotRectangle43(1,4);
  subPlotWidth    = subPlotRectangle43(1,3);
  subPlotOffsetY  = subPlotHeight+ySpace;
  subPlotOffsetX  = subPlotWidth+xSpace;

%%
% Model 1/4
%%    
row=0;
for i=1:1:length(freqSeriesFiles)
  col = 0;
  
  
  deltaY = subPlotOffsetY*0.5;
  if( contains(freqSeriesFiles{i},'Hill') == 1)
    deltaY = -subPlotOffsetY*0.5;
  end
    
  
  %Time
  j = (i-1)*3 + 1;
  
  subPlotList(j,:)    = currentSubPlot+[0,deltaY,0,0];
  subPlotList(j,1)    = subPlotList(j,1) + subPlotOffsetX*(col); 
  subPlotList(j,2)    = subPlotList(j,2) - subPlotOffsetY*(row);  

  %Gain  
  j = (i-1)*3 + 2;  
  col = col+1;  
  subPlotList(j,:)    = currentSubPlot+[0,deltaY,0,0];
  subPlotList(j,1)    = subPlotList(j,1) + subPlotOffsetX*(col); 
  subPlotList(j,2)    = subPlotList(j,2) - subPlotOffsetY*(row);  


  %Phase
  j = (i-1)*3 + 3;  
  col = col+1;  
  subPlotList(j,:)    = currentSubPlot+[0,deltaY,0,0];
  subPlotList(j,1)    = subPlotList(j,1) + subPlotOffsetX*(col); 
  subPlotList(j,2)    = subPlotList(j,2) - subPlotOffsetY*(row);  
  
  row = row+1;

end

subPlotLabel = {'A. Model with a damped-elastic tendon','B.','C.',... 
                'D. Model with a rigid tendon','E.','F.',...
                'A. Hill model with an elastic tendon','B.','C.',...
                'D. Hill model with a rigid tendon','E.','F.',...
                'X.','Y.','Z.'};

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
%Plot model data
%%
switch flag_Mode15Hz90HzBoth
  case 0
    modelAmpPlot =[1.6];
    modelBWPlot  =[15];    
  case 1
    modelAmpPlot =[1.6];
    modelBWPlot  =[90];        
  case 2
    modelAmpPlot =[1.6, 1.6];
    modelBWPlot  =[15,  90];    
    
  otherwise
    assert(0);
end

modelForce   =[1,1,1].*nominalForce(1,nominalForceTargetIndex);
modelNormFiberLength = [1,1,1].*normFiberLength;

colorPlot = [0.,0.,1; 0.25,0.25,1; 0.5,0.5,1];
lineWidth = [1,1,1].*0.5;

%Here 'Full' means a 'Full' muscle model: Hill or otherwise  
%      KD    means the spring of best fit.

lineColorKD   = zeros(length(modelAmpPlot),3);
lineColorFull = zeros(length(modelAmpPlot),3); 
lineWidthFull = ones(length(modelAmpPlot),3);

plotForceMax = 17.2;  
vafInput     = zeros(length(freqSeriesFiles),length(modelAmpPlot));

lineHandlesModelGain = [];
lineHandlesModelPhase = [];  

for z=1:1:length(freqSeriesFiles)

  modelColor = freqSeriesColor(z,:);
  tmp = load([dataInputFolder, freqSeriesFiles{1,z}]);
  freqSimData = tmp.freqSimData;

  fprintf('%s VAF: %1.1f pm %1.1f, [%1.1f, %1.1f]\n',...
           freqSeriesName{1,z}, ...
           mean(freqSimData.vafTime)*100,...
           std(freqSimData.vafTime)*100,...
           min(freqSimData.vafTime)*100,...
           max(freqSimData.vafTime)*100);

  flag_Hill = 0;
  if(isempty(strfind(freqSeriesFiles{1,z},'Hill'))==0)
    flag_Hill = 1;
  end

  modelLineWidth = expLineWidth;

  flagPlotEmpty  = ones(4,1);
  
  for k=1:1:length(modelAmpPlot)

    k01                = 0.75*(k-1)/(length(modelAmpPlot));
    lineColorFull(k,:) = modelColor;
    lineColorKD(k,:)   = [0.75,0.75,0.75];
    lineWidthKD        = 2;

    idxSim = 0;
    tol = 1e-6;
    for m=1:1:size(freqSimData.force,2)     
      if( abs(freqSimData.amplitudeMM(1,m)     - modelAmpPlot(1,k)        ) <= tol && ...
          abs(freqSimData.bandwidthHz(1,m)     - modelBWPlot(1,k)         ) <= tol && ...
          abs(freqSimData.nominalForceDesired(1,m)    - modelForce(1,k)          ) <= tol && ...
          abs(freqSimData.normFiberLength(1,m) - modelNormFiberLength(1,k)) <= tol)

        if(idxSim == 0)
          idxSim = m;
        else
          assert(0); %Error condition: there should not be 2 simulations with 
                     %the same configuration
        end
      end
    end

    idxWave = getSignalIndex(modelAmpPlot(1,k),...
                             modelBWPlot(1,k),...
                             inputFunctions);                            

    idxForce = 3*(z-1)+1;
    idxGain  = 3*(z-1)+2;
    idxPhase = 3*(z-1)+3;
    %if(flag_Hill)
    %  idxForce = 2;
    %  idxGain  = 5;
    %  idxPhase = 6;
    %end

    subplot('Position', [ subPlotList(idxForce,1),...
                          subPlotList(idxForce,2),...
                          subPlotList(idxForce,3),...
                          subPlotList(idxForce,4)]);    
   
                      
                        
    if(flagPlotEmpty(z,1)==1 && k==length(modelAmpPlot))

      yo = freqSimData.nominalForce(1,idxSim);

      vafInput(z,k) = freqSimData.vafTime(1,idxSim);

      spring = freqSimData.stiffness(1,idxSim)./1000;
      damper = freqSimData.damping(1,idxSim)./1000;
      kLabel = sprintf('%1.1f',spring);
      dLabel = sprintf('%1.3f',damper);

      kdLabel = ['K: ',kLabel,'N/mm, $$\beta$$: ',dLabel,'N/(mm/s)',...
                 ' VAF ', num2str(round(freqSimData.vafTime(1,idxSim)*100)),'\%'];

      kLabel = ['  K: ',kLabel,'N/mm'];
      dLabel = ['  $$\beta$$: ',dLabel,'N/(mm/s)'];
      vafLabel = [' VAF ', num2str(round(freqSimData.vafTime(1,idxSim)*100)),'\%'];

      lineLabel = [ freqSeriesName{z}];



      pidKD = plot( inputFunctions.time(idxChunk,1),...
             freqSimData.forceKD(idxChunk,idxSim)+yo,...
             'Color',lineColorKD(k,:),...
             'LineWidth',lineWidthKD);
      hold on;


      pidMdl = plot( inputFunctions.time(idxChunk,1),...
            freqSimData.force(idxChunk,idxSim),...
            'Color',lineColorFull(k,:),...
            'LineWidth',lineWidthFull(k,1));
      hold on;  

      simText = [num2str(modelAmpPlot(1,k)),'mm',' ',...
                 num2str(modelBWPlot(1,k)),'Hz'];

      tlabel = inputFunctions.time(idxChunk(1,1),1) ...
               + 0.0*(max(inputFunctions.time(idxChunk,1))...
                     -min(inputFunctions.time(idxChunk,1)));             

      text( tlabel,...
            plotForceMax*0.9,...
            simText );        
      hold on;

      tlabel = inputFunctions.time(idxChunk(1,1),1) ...
               + 0.1*(max(inputFunctions.time(idxChunk,1))...
                     -min(inputFunctions.time(idxChunk,1)));
      text( tlabel,...
            plotForceMax*0.8,...
            kLabel);        
      hold on;
      text( tlabel,...
            plotForceMax*0.7,...
            dLabel);        
      hold on;
      text( tlabel,...
            plotForceMax*0.6,...
            vafLabel);        
      hold on;

      box off;
      set(gca,'color','none')

      %lh = legend([pidMdl,pidKD],lineLabel,kdLabel,'Location','NorthWest');
      %lh.Position(1,1) = lh.Position(1,1)-0.04;
      %lh.Position(1,2) = lh.Position(1,2)+0.055;


      %legend boxoff;
      tmin = inputFunctions.time(idxChunk(1),1  )-0.01;
      tmax = inputFunctions.time(idxChunk(end),1)+0.01;
      xlim([tmin,...
            tmax]);
      ylim([0,plotForceMax]);
      xticks(timeTicks);

      f0 = freqSimData.force(idxChunk(1),idxSim);
      fmax =max(freqSimData.force(idxChunk,idxSim));
      fmin = 0;
      yticks([0,round(f0,1),round(fmax,1)]);
      set(gca,'color','none')

      tmin = min(inputFunctions.time(idxChunk,1));
      tmax = max(inputFunctions.time(idxChunk,1));

      %if(k == 1)
        idxSubPlot = idxForce;
        tc = text(tmin-0.1*(tmax-tmin), 1.1*plotForceMax,...
               subPlotLabel{1,idxSubPlot},...
               'FontSize',figLabelFontSize);   
        hold on;
      %end      
      
      ylabel('Force (N)');
      xlabel('Time (s)');
      
      if(z==1)
        x0 = subPlotList(idxForce,1);
        y0 = subPlotList(idxForce,2);
        dx = subPlotList(idxForce,3);
        dy = subPlotList(idxForce,4);

        text(tmin,plotForceMax*1.2,'Simulation of Kirsch, Boskov, \& Rymer 1994',...
        'FontSize',8*1.2,...
        'HorizontalAlignment','left',...
        'VerticalAlignment','bottom');
        hold on;     

      end

      if(z==3)
        x0 = subPlotList(idxForce,1);
        y0 = subPlotList(idxForce,2);
        dx = subPlotList(idxForce,3);
        dy = subPlotList(idxForce,4);

        text(tmin,plotForceMax*1.2,'Simulation of Kirsch, Boskov, \& Rymer 1994',...
        'FontSize',8*1.2,...
        'HorizontalAlignment','left',...
        'VerticalAlignment','bottom');
        hold on;     

      end        
    end
    

  end



%%
% Plot the experimental data
%%
  for k=1:1:length(expAmpPlot)

    idxSim = 0;
    tol = 1e-6;
    for m=1:1:size(freqSimData.force,2)     
      if( abs(freqSimData.amplitudeMM(1,m) - expAmpPlot(1,k)) <= tol && ...
          abs(freqSimData.bandwidthHz(1,m) - expBWPlot(1,k))  <= tol && ...
          abs(freqSimData.nominalForceDesired(1,m)- modelForce(1,k)) <= tol && ...
          abs(freqSimData.normFiberLength(1,m)-modelNormFiberLength(1,k)) <= tol)

        if(idxSim == 0)
          idxSim = m;
        else
          assert(0); %Error condition: there should not be 2 simulations with 
                     %the same configuration
        end
      end
    end      

    fprintf('%i. K: %1.3f D: %1.3f\n',idxSim,...
              freqSimData.stiffness(1,idxSim)./1000,...
              freqSimData.damping(1,idxSim)./1000);

    %
    idxRange = [freqSimData.idxFreqRange(1,idxSim):1: ...
                freqSimData.idxFreqRange(2,idxSim)];

    idxCutoff = freqSimData.idxFreqRange(1,idxSim);


    subplot('Position', subPlotList(idxGain,:));

    markType = '.';
    markFaceColor = freqSeriesColor(z,:);
    markSize = 2;
    lineWidth = 0.5;
    kdMarkType = '-';
    kdLineWidth= 1;
    kdWhiteLineWidth = 2;
    kdLineColor = freqSeriesColor(z,:);
    if(k==2)        
      markType = 'o';
      markFaceColor = [1,1,1];
      markSize = 3;
      lineWidth=0.5;
      kdMarkType = '-';
      kdLineWidth= 1;
      kdWhiteLineWidth = 3;
      kdLineColor = freqSeriesColor(z,:).*0.5 + [1,1,1].*0.5;

    end

    trialName = sprintf(': %1.1fmm %1.0fHz',...
                        freqSimData.amplitudeMM(1,idxSim),...
                        freqSimData.bandwidthHz(1,idxSim));

    plot(freqSimData.freqHz(idxRange,1), ...
         freqSimData.gain(idxRange,idxSim)./1000,...
         markType,'Color',kdLineColor, ...
         'MarkerFaceColor',markFaceColor,'LineWidth',lineWidth,...
         'MarkerSize',markSize,'DisplayName',[freqSeriesName{z},trialName]);  
    hold on;

    pidW=plot(freqSimData.freqHz(idxRange,1), ...
         freqSimData.gainKD(idxRange,idxSim)./1000,...
         '-','Color',[1,1,1], ...
         'LineWidth',kdWhiteLineWidth);  
    hold on;
    set(get(get(pidW,'Annotation'),...
              'LegendInformation'),...
              'IconDisplayStyle','off');      

    pidL=plot(freqSimData.freqHz(idxRange,1), ...
         freqSimData.gainKD(idxRange,idxSim)./1000,...
         kdMarkType,'Color',kdLineColor, ...
         'LineWidth',kdLineWidth,'DisplayName', ['K-$$\beta$$',trialName]);
    hold on;

    lineHandlesModelGain = [lineHandlesModelGain,...
                            pidW,...
                            pidL];



    fmin = 0.;
    fmax = max(modelBWPlot);
    gmin = 0.;%min(freqSimData.gain(idxRange,idxSim)./1000);
    gmax = 8.01;

    if(k == 1)
      idxSubPlot = idxGain;
      text(fmin-0.1*(fmax-fmin), 1.1*gmax,...
             subPlotLabel{1,idxSubPlot},...
             'FontSize',figLabelFontSize);   
      hold on;
    end

    box off;
    set(gca,'color','none')        


    hold on;


    ylim([gmin,gmax]);
    xlim([fmin,fmax+0.01]);
    xticks(xTicksVector);
    yticks([0,2,4,6,8]);


    %tc0 = text(fmin-0.15*(fmax-fmin), 1.45*(gmax),...
    %       'C.','FontSize',11);       

    subplot('Position', subPlotList(idxPhase,:));

    plot(freqSimData.freqHz(idxRange,1), ...
         freqSimData.phase(idxRange,idxSim).*(180/pi),...
         markType,'Color',kdLineColor, ...
         'MarkerFaceColor',markFaceColor,'LineWidth',0.5,...
         'MarkerSize',markSize,...
         'DisplayName',[freqSeriesName{z},trialName]);   
    hold on;

    pidW=plot(freqSimData.freqHz(idxRange,1), ...
         freqSimData.phaseKD(idxRange,idxSim).*(180/pi),...
         '-','Color',[1,1,1], ...
         'LineWidth',kdWhiteLineWidth);  
    hold on;
    set(get(get(pidW,'Annotation'),...
              'LegendInformation'),...
              'IconDisplayStyle','off');      

    pidL=plot(freqSimData.freqHz(idxRange,1), ...
         freqSimData.phaseKD(idxRange,idxSim).*(180/pi),...
         kdMarkType,'Color',kdLineColor, ...
         'LineWidth',kdLineWidth,'DisplayName',['K-$$\beta$$',trialName]);  
    hold on;      

    lineHandlesModelPhase = [lineHandlesModelPhase,...
                            pidW,...
                            pidL];


    box off;
    set(gca,'color','none')

    pmin = 0;
    pmax = 90.01;

    ylim([pmin,pmax]);
    xlim([fmin,fmax+0.01]);
    xticks(xTicksVector);
    yticks([0,45,90]);

    if(k == 1)
      idxSubPlot = idxPhase;      
      text(fmin-0.1*(fmax-fmin), 1.1*pmax,...
             subPlotLabel{1,idxSubPlot},...
             'FontSize',figLabelFontSize);   
      hold on;
    end


  end

end


for z=1:1:length(freqSeriesFiles)

  flag_Hill = 0;
  if(isempty(strfind(freqSeriesFiles{1,z},'Hill'))==0)
    flag_Hill = 1;
  end    

  idxForce = 3*(z-1)+1;
  idxGain  = 3*(z-1)+2;
  idxPhase = 3*(z-1)+3;  
  
%   idxForce = 1;
%   idxGain  = 3;
%   idxPhase = 4;
%   if(flag_Hill)
%     idxForce = 2;
%     idxGain  = 5;
%     idxPhase = 6;
%   end


  subplot('Position', subPlotList(idxGain,:));

  for k=1:1:length(expAmpPlot)

    strAmp = num2str(expAmpPlot(1,k));
    strFreq= num2str(expBWPlot(1,k));
    
    idxExp = 0;
    tol = 1e-6;
    for m=1:1:length(dataKBR1994Fig3Gain)           
      if( contains(dataKBR1994Fig3Gain(m).seriesName,strFreq) && ...
          contains(dataKBR1994Fig3Gain(m).seriesName,strAmp))        
        if(idxExp == 0)
          idxExp = m;
        else
          assert(0); %Error condition: there should not be 2 simulations with 
                     %the same configuration
        end
      end
    end      
    
    idxMin = find(dataKBR1994Fig3Gain(idxExp).x >= 4 ,1, 'first' );
    idxMax = find(dataKBR1994Fig3Gain(idxExp).x <= expBWPlot(1,k) ,1, 'last' );
    pid = plot( dataKBR1994Fig3Gain(idxExp).x(idxMin:1:idxMax),...
          dataKBR1994Fig3Gain(idxExp).y(idxMin:1:idxMax),...
          '-', 'Color', [1,1,1],...
          'LineWidth', expLineWidth(1,k)*3,...
          'DisplayName','');
    hold on;

    set(get(get(pid,'Annotation'),...
                'LegendInformation'),...
                'IconDisplayStyle','off');

    plot( dataKBR1994Fig3Gain(idxExp).x(idxMin:1:idxMax),...
          dataKBR1994Fig3Gain(idxExp).y(idxMin:1:idxMax),...
          expMarkType{1,k}, 'Color', expPlotColor(k,:),...
          'LineWidth',expLineWidth(1,k),...
          'DisplayName', expLegendEntry{k});

    hold on;
  end

  box off;
  set(gca,'color','none')  
  xlabel(dataKBR1994Fig3Gain(idxExp).xName);
  ylabel('Gain (N/mm)');

  %if(flag_useElasticTendon==0 )
    
    %lh = legend('Location','SouthEast');
    %lh.Position(1,1) = lh.Position(1,1)+0.03;
    %lh.Position(1,2) = lh.Position(1,2)-0.03;
    
  %end

  subplot('Position', subPlotList(idxPhase,:));

  for k=1:1:length(expAmpPlot)
    
    strAmp = num2str(expAmpPlot(1,k));
    strFreq= num2str(expBWPlot(1,k));
    
    idxExp = 0;
    tol = 1e-6;
    for m=1:1:length(dataKBR1994Fig3Phase)           
      if( contains(dataKBR1994Fig3Phase(m).seriesName,strFreq) && ...
          contains(dataKBR1994Fig3Phase(m).seriesName,strAmp))        
        if(idxExp == 0)
          idxExp = m;
        else
          assert(0); %Error condition: there should not be 2 simulations with 
                     %the same configuration
        end
      end
    end       
    

    idxMin = find(dataKBR1994Fig3Phase(idxExp).x >= 4,1, 'first' );
    idxMax = find(dataKBR1994Fig3Phase(idxExp).x <= expBWPlot(1,k),1, 'last' );
    
    pid = plot( dataKBR1994Fig3Phase(idxExp).x(idxMin:1:idxMax),...
          dataKBR1994Fig3Phase(idxExp).y(idxMin:1:idxMax),...
          '-', 'Color', [1,1,1],...
          'LineWidth',expLineWidth(1,k)*3); 
    hold on;

    set(get(get(pid,'Annotation'),...
                'LegendInformation'),...
                'IconDisplayStyle','off');

    plot( dataKBR1994Fig3Phase(idxExp).x(idxMin:1:idxMax),...
          dataKBR1994Fig3Phase(idxExp).y(idxMin:1:idxMax),...
          expMarkType{1,k}, 'Color', expPlotColor(k,:),...
          'LineWidth',expLineWidth(1,k),...
          'DisplayName', expLegendEntry{k}); 
    hold on;    
  end


  xlabel(dataKBR1994Fig3Phase(idxExp).xName);
  ylabel(dataKBR1994Fig3Phase(idxExp).yName);
  box off;
  set(gca,'color','none')


  %if(flag_useElasticTendon==1 )
  if(z == 2)
    lh = legend('Location','South');    
    lh.Position(1,1) = lh.Position(1,1) -0.15625;        
    lh.Position(1,2) = lh.Position(1,2) -0.140625;   
    lh.NumColumns=3;
    legend boxoff;        
  end
  if(z == 4)%length(freqSeriesFiles))
    lh = legend('Location','South');    
    lh.Position(1,1) = lh.Position(1,1) -0.15625;            
    lh.Position(1,2) = lh.Position(1,2) -0.140625;    
    lh.NumColumns=3;    
    legend boxoff;    
  end
  

  
  % end
%     if(flag_useElasticTendon==1)
%       lh = legend('Location','NorthEast');  
%       lh.Position(1,1) = lh.Position(1,1)+0.035;    
%       lh.Position(1,2) = lh.Position(1,2)+0.09;
%       lh.EdgeColor = lh.Color;
% 
%     end
end

%for i=1:1:length(lineHandlesModelPhase)
%uistack(lineHandlesModelPhase(i),'top');
%end
%for i=1:1:length(lineHandlesModelGain)
%uistack(lineHandlesModelGain(i),'top');
%end


set(fig_freqResponse43,'Units','centimeters',...
'PaperUnits','centimeters',...
'PaperSize',[pageWidth pageHeight],...
'PaperPositionMode','manual',...
'PaperPosition',[0 0 pageWidth pageHeight]);     
%set(findall(figList(i).h,'-property','FontSize'),'FontSize',10);     
set(fig_freqResponse43,'renderer','painters');     
set(gcf,'InvertHardCopy','off')

tendonTag = '_RigidAndElasticTendon';

plotBW = '';
switch flag_Mode15Hz90HzBoth
  case 0
    plotBW = '_Plot15Hz';
  case 1
    plotBW = '_Plot90Hz';    
  case 2
    plotBW = '_Plot15Hz90Hz';
  otherwise
    assert(0);
end

print('-dpdf', [pubOutputFolder,'fig_Pub_ModelFrequencyResponse',...
                tendonTag,'_',plotBW,'_',plotNameEnding,'.pdf']); 


success = 1;