function [success] = plotAccelerationEquationFactorsOpus31( ...
                          dataInputFolder,...
                            freqSeriesFiles_RT,...
                            freqSeriesName_RT,...
                            freqSeriesColor_RT,...
                            freqSeriesFiles_ET,...
                            freqSeriesName_ET,...
                            freqSeriesColor_ET,...                            
                          simSeriesFiles_RT,...
                          simSeriesName_RT,...
                          simSeriesColor_RT,...
                          simSeriesFiles_ET,...
                          simSeriesName_ET,...
                          simSeriesColor_ET,...                      
                            inputFunctions,...                      
                            normFiberLength,...
                            nominalForce,...
                            nominalForceTargetIndex,...                 
                            flag_useElasticTendon,...  
                            flag_Mode15Hz90HzBoth,...
                          plotNameEnding,...
                          plotLayoutSettings,...
                          pubOutputFolder)
%dataKBR1994Fig3Gain,...
%dataKBR1994Fig3Phase,...     
                      
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

simSeriesFiles =  {simSeriesFiles_ET{1,2},...
                   simSeriesFiles_RT{1,2},...
                   simSeriesFiles_ET{1,1},...
                   simSeriesFiles_RT{1,1}};
simSeriesName =   {simSeriesName_ET{1,2},...
                   simSeriesName_RT{1,2},...
                   simSeriesName_ET{1,1},...
                   simSeriesName_RT{1,1}};
simSeriesColor =  [simSeriesColor_ET(2,:);...
                   simSeriesColor_RT(2,:);...
                   simSeriesColor_ET(1,:);...
                   simSeriesColor_RT(1,:)];

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

fig_AccelerationTerms = figure;

%Time window to plot
chunkDuration = 0.2;        
timeChunkStart = round((paddingPoints)*0.5)/sampleFrequency;  
timeChunkEnd   = paddingPoints/sampleFrequency+chunkDuration;
idxChunk       = [round((paddingPoints)*0.5):1: ...
                  (paddingPoints + sampleFrequency*chunkDuration)];                  
timeTicks = [round(paddingPoints/sampleFrequency,2),...
             round(timeChunkEnd,2)];

%Count the number of non-Hill models
numberOfModels=0;
for i=1:1:length(freqSeriesFiles)
    if( contains(freqSeriesFiles{i},'Hill') == 0)
        numberOfModels=numberOfModels+1;
    end
end

subPlotList = zeros(numberOfModels*6,4);

ySpace = 0.05;
xSpace = 0.04;

%%
% Plot Layout
%%
figure(fig_AccelerationTerms);
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
numberOfPlotsPerModel=6;
for i=1:1:numberOfModels
  col = 0;
  
  deltaY=0;
  if(i==1)
    deltaY = subPlotOffsetY*0.5;
  end
    
  %la
  j = (i-1)*numberOfPlotsPerModel + 1;
  
  subPlotList(j,:)    = currentSubPlot+[0,deltaY,0,0];
  subPlotList(j,1)    = subPlotList(j,1) + subPlotOffsetX*(col); 
  subPlotList(j,2)    = subPlotList(j,2) - subPlotOffsetY*(row);  

  %dla 
  j = (i-1)*numberOfPlotsPerModel + 2;  
  col = col+1;  
  subPlotList(j,:)    = currentSubPlot+[0,deltaY,0,0];
  subPlotList(j,1)    = subPlotList(j,1) + subPlotOffsetX*(col); 
  subPlotList(j,2)    = subPlotList(j,2) - subPlotOffsetY*(row);  


  %ddla
  j = (i-1)*numberOfPlotsPerModel + 3;  
  col = col+1;  
  subPlotList(j,:)    = currentSubPlot+[0,deltaY,0,0];
  subPlotList(j,1)    = subPlotList(j,1) + subPlotOffsetX*(col); 
  subPlotList(j,2)    = subPlotList(j,2) - subPlotOffsetY*(row);  

  row = row+1;

  %term 1
  j = (i-1)*numberOfPlotsPerModel + 4;
  col = 0;
  subPlotList(j,:)    = currentSubPlot+[0,deltaY,0,0];
  subPlotList(j,1)    = subPlotList(j,1) + subPlotOffsetX*(col); 
  subPlotList(j,2)    = subPlotList(j,2) - subPlotOffsetY*(row);  

  %term 2 
  j = (i-1)*numberOfPlotsPerModel + 5;  
  col = col+1;  
  subPlotList(j,:)    = currentSubPlot+[0,deltaY,0,0];
  subPlotList(j,1)    = subPlotList(j,1) + subPlotOffsetX*(col); 
  subPlotList(j,2)    = subPlotList(j,2) - subPlotOffsetY*(row);  

  %term 3
  j = (i-1)*numberOfPlotsPerModel + 6;  
  col = col+1;  
  subPlotList(j,:)    = currentSubPlot+[0,deltaY,0,0];
  subPlotList(j,1)    = subPlotList(j,1) + subPlotOffsetX*(col); 
  subPlotList(j,2)    = subPlotList(j,2) - subPlotOffsetY*(row);  

  row = row+1;


end

subPlotLabel = {'1A. Model with a damped-elastic tendon','B.','C.','D.','E.','F.',... 
                '2A. Model with a rigid tendon','B.','C.','D.','E.','F.'};

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



for z=1:1:length(freqSeriesFiles)


  tmp = load([dataInputFolder, freqSeriesFiles{1,z}]);
  freqSimData = tmp.freqSimData;
  
  benchRecordSet = load([dataInputFolder, simSeriesFiles{1,z}]);

  flag_Hill = 0;
  if(isempty(strfind(freqSeriesFiles{1,z},'Hill'))==0)
    flag_Hill = 1;
  end

  
  for k=1:1:length(modelAmpPlot)

    k01        = (k-1)/(length(modelAmpPlot)-1);
    modelColor = freqSeriesColor(z,:).*(0.75-0.5*k01) ...
                +[0.75,0.75,0.75].*(0.25+0.5*k01); 
    lineWidth = 1*(1-k01) + 0.5*k01;

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
                        

    idxLength       = numberOfPlotsPerModel*(z-1)+1;
    idxVelocity     = numberOfPlotsPerModel*(z-1)+2;
    idxAcceleration = numberOfPlotsPerModel*(z-1)+3;
    idxTerm1        = numberOfPlotsPerModel*(z-1)+4;
    idxTerm2        = numberOfPlotsPerModel*(z-1)+5;
    idxTerm3        = numberOfPlotsPerModel*(z-1)+6;
    



                      
                        
    if(idxSim ~= 0 && flag_Hill==0)



      numberOfSamples     = size(benchRecordSet.benchRecord.state,1);
      numberOfSimulations = size(benchRecordSet.benchRecord.state,2);
      numberOfStates      = size(benchRecordSet.benchRecord.state,3);
      numberOfExtras      = size(benchRecordSet.benchRecord.extra,3);

      extra = reshape( benchRecordSet.benchRecord.extra(:,idxSim,:),...
                       numberOfSamples, numberOfExtras);      
      state = reshape( benchRecordSet.benchRecord.state(:,idxSim,:),...
                       numberOfSamples, numberOfStates);
      dstate = reshape( benchRecordSet.benchRecord.dstate(:,idxSim,:),...
                       numberOfSamples, numberOfStates);
      
      laN    = [];
      dlaN   = [];
      ddlaN  = [];      
      l1N    = [];
      dl1N   = [];
      hillTerm     = [];
      dampingTerm  = [];
      trackingTerm = [];
      hillTermLabel     = {};
      dampingTermLabel  = {};
      trackingTermLabel = {};
      
      
      optimalFiberLength=benchRecordSet.muscleArchitecture.optimalFiberLength;

      switch numberOfStates
          case 3
              laN   =  state(:,2)./optimalFiberLength;
              dlaN  =  state(:,1)./optimalFiberLength;
              ddlaN = dstate(:,1)./optimalFiberLength;      
              l1N   =  state(:,3)./optimalFiberLength;
              dl1N  = dstate(:,3)./optimalFiberLength;
              hillTerm     = extra(:,3);
              dampingTerm  = extra(:,4);
              trackingTerm = extra(:,5);              
              hillTermLabel     = benchRecordSet.benchRecord.extraLabels{3};
              dampingTermLabel  = benchRecordSet.benchRecord.extraLabels{4};
              trackingTermLabel = benchRecordSet.benchRecord.extraLabels{5};
              
          case 4
              laN   =  state(:,3)./optimalFiberLength;
              dlaN  =  state(:,2)./optimalFiberLength;
              ddlaN = dstate(:,2)./optimalFiberLength;      
              l1N   =  state(:,4)./optimalFiberLength;
              dl1N  = dstate(:,4)./optimalFiberLength;
              hillTerm     = extra(:,3);
              dampingTerm  = extra(:,4);
              trackingTerm = extra(:,5);              
              hillTermLabel     = benchRecordSet.benchRecord.extraLabels{3};
              dampingTermLabel  = benchRecordSet.benchRecord.extraLabels{4};
              trackingTermLabel = benchRecordSet.benchRecord.extraLabels{5};
          otherwise
              assert(0,'Model has unexpected number of states')
      end

      tmin = inputFunctions.time(idxChunk(1),1  )-0.01;
      tmax = inputFunctions.time(idxChunk(end),1)+0.01;        
      
      %%
      %Length
      %%
      subplot('Position', [ subPlotList(idxLength,1),...
                            subPlotList(idxLength,2),...
                            subPlotList(idxLength,3),...
                            subPlotList(idxLength,4)]);  

      pidLength = plot(  inputFunctions.time(idxChunk,1),...
                         laN(idxChunk,:),...
                         'Color',modelColor,...
                         'LineWidth',lineWidth);
      hold on;
      
      box off;
      set(gca,'color','none')
      
      ylabel('$$\ell^A / \ell_\circ^{\mathrm{M}}$$');
      xlabel('Time (s)');
      title(subPlotLabel{idxLength});
      xlim([tmin,tmax]);
      xticks(timeTicks);
      

      %%
      %Velocity
      %%

      subplot('Position', [ subPlotList(idxVelocity,1),...
                            subPlotList(idxVelocity,2),...
                            subPlotList(idxVelocity,3),...
                            subPlotList(idxVelocity,4)]); 

      pidVelocity = plot(inputFunctions.time(idxChunk,1),...
                         dlaN(idxChunk,:),...
                         'Color',modelColor,...
                         'LineWidth',lineWidth);
      hold on;
      
      box off;
      set(gca,'color','none')
      
      ylabel('$$\dot{\ell}^A / \ell_\circ^{\mathrm{M}}$$');
      xlabel('Time (s)');
      title(subPlotLabel{idxVelocity});

      xlim([tmin,tmax]);
      xticks(timeTicks);

      %%
      %Acceleration
      %%
      subplot('Position', [ subPlotList(idxAcceleration,1),...
                            subPlotList(idxAcceleration,2),...
                            subPlotList(idxAcceleration,3),...
                            subPlotList(idxAcceleration,4)]); 

      pidAcceleration = plot(inputFunctions.time(idxChunk,1),...
                             ddlaN(idxChunk,:),...
                             'Color',modelColor,...
                             'LineWidth',lineWidth);
      hold on;
      
      box off;
      set(gca,'color','none')
      
      ylabel('$$\ddot{\ell}^A / \ell_\circ^{\mathrm{M}}$$');
      xlabel('Time (s)');
      title(subPlotLabel{idxAcceleration});
      xlim([tmin,tmax]);
      xticks(timeTicks);

      %%
      %Acceleration term 1
      %%
      subplot('Position', [ subPlotList(idxTerm1,1),...
                            subPlotList(idxTerm1,2),...
                            subPlotList(idxTerm1,3),...
                            subPlotList(idxTerm1,4)]); 

      pidTerm1 = plot(inputFunctions.time(idxChunk,1),...
                      hillTerm(idxChunk,:),...
                      'Color',modelColor,...
                      'LineWidth',lineWidth);
      hold on;
      
      box off;
      set(gca,'color','none')
      
      ylabel(hillTermLabel);
      xlabel('Time (s)');
      title(subPlotLabel{idxTerm1});

      xlim([tmin,tmax]);
      xticks(timeTicks);

      %%
      %Acceleration term 2
      %%
      subplot('Position', [ subPlotList(idxTerm2,1),...
                            subPlotList(idxTerm2,2),...
                            subPlotList(idxTerm2,3),...
                            subPlotList(idxTerm2,4)]); 

      pidTerm2 = plot(inputFunctions.time(idxChunk,1),...
                      dampingTerm(idxChunk,:),...
                      'Color',modelColor,...
                      'LineWidth',lineWidth);
      hold on;
      
      box off;
      set(gca,'color','none')
      
      ylabel(dampingTermLabel);
      xlabel('Time (s)');
      title(subPlotLabel{idxTerm2});

      xlim([tmin,tmax]);
      xticks(timeTicks);

      %%
      %Acceleration term 3
      %%
      subplot('Position', [ subPlotList(idxTerm3,1),...
                            subPlotList(idxTerm3,2),...
                            subPlotList(idxTerm3,3),...
                            subPlotList(idxTerm3,4)]); 

      pidTerm3 = plot(inputFunctions.time(idxChunk,1),...
                      trackingTerm(idxChunk,:),...
                      'Color',modelColor,...
                      'LineWidth',lineWidth);
      hold on;
      
      box off;
      set(gca,'color','none')
      
      ylabel(trackingTermLabel);
      xlabel('Time (s)');
      title(subPlotLabel{idxTerm3});

      xlim([tmin,tmax]);
      xticks(timeTicks);

    end


    

  end


end




set(fig_AccelerationTerms,'Units','centimeters',...
'PaperUnits','centimeters',...
'PaperSize',[pageWidth pageHeight],...
'PaperPositionMode','manual',...
'PaperPosition',[0 0 pageWidth pageHeight]);     
%set(findall(figList(i).h,'-property','FontSize'),'FontSize',10);     
set(fig_AccelerationTerms,'renderer','painters');     
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

print('-dpdf', [pubOutputFolder,'fig_Pub_ModelAccelerationTerms',...
                tendonTag,'_',plotBW,'_',plotNameEnding,'.pdf']); 


success = 1;