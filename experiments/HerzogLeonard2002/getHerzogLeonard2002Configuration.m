function expConfigHerzogLeonard2002 = getHerzogLeonard2002Configuration(...
                        figureNumber, subFigureNumber, trialNumber,...
                        projectFolders)

idxFitMax = 1500; %All experimental data except series 50 stops at 1500 samples

rampSlope      = 0;
rampStartTime  = 0;
rampStretch    = 0;

actT0          = 0;
actT1          = NaN;

colLength =0;
colForce = 0;
colNoStretch=0;

colLengthPassive = 0;
colForcePassive = 0;

%lceOpt   = felineMusculotendonProperties.lceOpt;
%alphaOpt = felineMusculotendonProperties.pennationAngle;
%ltSlk    = felineMusculotendonProperties.tendonSlackLength;
%eIso     = felineMusculotendonProperties.tendonStrainAtOneNormForce;

%normTendonLength = 1;
%if(useElasticTendon == 1)
%  normTendonLength = 1+eIso;
%end

%lmOpt = lceOptExp*cos(alphaOpt) + ltSlk*normTendonLength;

    
switch(figureNumber)

  case 7
    switch(subFigureNumber)
      case 1
        expDataFile = fullfile(projectFolders.experiments_HL2002,...
                        'data','dataHerzogLeonard2002Figure6and7A.csv');     
        tmp       = csvread(expDataFile,1,0);
        idxMax = 0;
        for(z=1:1:size(tmp,1))
           if(tmp(z,1) > 0)
              idxMax = z; 
           end
        end
        expData = tmp(1:1:idxMax,:);        
        rampSlope      = (3/1000); 
        colNoStretch = 2;
        switch trialNumber
            case 1
                colLength = 7;
                colForce  = colLength-1; 
            case 2
                colLength = 9;
                colForce  = colLength-1;        
            case 3
                colLength = 13;
                colForce  = colLength-1;        
            otherwise
                assert(0,'trialNumber must be 1,2, or 3');
        end
        
        colForcePassive  = 10;
        colLengthPassive = colForcePassive+1;

      case 2
        expDataFile = fullfile(projectFolders.experiments_HL2002,...
                        'data','dataHerzogLeonard2002Figure7B.csv');       
        expData       = csvread(expDataFile,1,0);        
        rampSlope      = (9/1000);  
        idxFitMax = size(expData,1);
        colNoStretch = 2;
        switch trialNumber
            case 1
                colLength = 7;
                colForce  = colLength-1; 
            case 2
                colLength = 9;
                colForce  = colLength-1;        
            case 3
                colLength = 11;
                colForce  = colLength-1;        
            otherwise
                assert(0,'trialNumber must be 1,2, or 3');
        end
        
        colForcePassive  = 12;
        colLengthPassive = colForcePassive+1;

      case 3
        expDataFile = fullfile(projectFolders.experiments_HL2002,...
                        'data','dataHerzogLeonard2002Figure7C.csv');

        expData       = csvread(expDataFile,1,0);        
        rampSlope      = (27/1000);  
        idxFitMax = size(expData,1);
        colNoStretch = 2;
        switch trialNumber
            case 1
                colLength = 7;
                colForce  = colLength-1; 
            case 2
                colLength = 9;
                colForce  = colLength-1;        
            case 3
                colLength = 11;
                colForce  = colLength-1;        
            otherwise
                assert(0,'trialNumber must be 1,2, or 3');
        end
        
        colForcePassive  = 12;
        colLengthPassive = colForcePassive+1;
        
      otherwise 
        msg = sprintf('Fig %i (%i) : Does not exist',...
                      figureNumber,subFigureNumber);
        assert(0,msg);    
    end
  otherwise
    msg = sprintf('Fig %i (%i) : Does not exist',...
                figureNumber,subFigureNumber);
    assert(0,msg);    
end

switch(trialNumber)
    case 1
        rampStretch=3/1000;
    case 2
        rampStretch=6/1000;
    case 3
        rampStretch=9/1000;
    otherwise
        assert(0,'trialNumber must be 1, 2, or 3');
end


expLengthData = [ expData(1:idxFitMax,1),...
                  expData(1:idxFitMax,colLength)];
%expDeltaLengthData = expLengthData;                    
%expLengthData(:,2) = expLengthData(:,2)./1000 ...
%                   + lmOpt.*(ones(size(expLengthData,1),1)); 

expForceData       = [  expData(1:idxFitMax,1),...
                        expData(1:idxFitMax,colForce )];

expForceStaticData = [  expData(1:idxFitMax,1),...
                        expData(1:idxFitMax,colNoStretch )];
                    
%%
% Extract the times the activation signal turns on/off
% Extract the time the ramps starts and ends.
%%
diffTime = diff(expData(:,1));
meanSampleTime = mean(diffTime);
stdSampleTime  = std(diffTime);
assert( abs(stdSampleTime/meanSampleTime) < 1e-3, ...
        'Sample time is not regular');

freq      = 1/meanSampleTime;
filtFreq  = freq*0.125;
[b,a]     = butter(2,filtFreq/freq,'low');

expLengthDataFilt       = zeros(size(expLengthData,1),3);
expForceDataFilt        = zeros(size(expLengthData,1),3);
expForceStaticDataFilt  = zeros(size(expLengthData,1),3);


expLengthDataFilt(:,1)  = expLengthData(:,1);
expForceDataFilt(:,1)   = expForceData(:,1);
expForceStaticDataFilt(:,1)   = expForceData(:,1);

expLengthDataFilt(:,2)      = filtfilt(b,a,expLengthData(:,2));
expForceDataFilt(:,2)       = filtfilt(b,a,expForceData(:,2));
expForceStaticDataFilt(:,2) = filtfilt(b,a,expForceStaticData(:,2));

%Take a central difference to get df/dt
dt = meanSampleTime;
for k=2:1:(size(expLengthData,1)-1)
    expLengthDataFilt(k,3) = ...
        ((expLengthDataFilt(k,2)-expLengthDataFilt(k-1,2)) ...
      + (expLengthDataFilt(k+1,2)-expLengthDataFilt(k,2))) / (2*dt);

    expForceDataFilt(k,3) = ...
        ((expForceDataFilt(k,2) - expForceDataFilt(k-1,2)) ...
      + (expForceDataFilt(k+1,2)- expForceDataFilt(k,2))) / (2*dt);  
  
    expForceStaticDataFilt(k,3) =...
        ((expForceStaticDataFilt(k,2) - expForceStaticDataFilt(k-1,2)) ...
      + (expForceStaticDataFilt(k+1,2)- expForceStaticDataFilt(k,2)) )/(2*dt);  
end


%%
% Get the activation signal start and end times
%%


[dFMaxVal, idxDfMax] = max(expForceStaticDataFilt(:,3));
[dFMaxVal, idxDfMin] = min(expForceStaticDataFilt(:,3));
[maxFVal  , idxFMax] = max(expForceData(:,2));

actT0 = NaN;
actT1 = NaN;
idxActT0 = NaN;
idxActT1 = NaN;
flagSet = 0;
k       = idxDfMax;
while(flagSet == 0 && k > 2)
   if(expForceData(k,2) <= expForceData(k-1,2) && flagSet == 0)
       actT0 = expForceData(k+1,1);
       idxActT0 = k+1;
       flagSet =1;
   end
   k = k-1;
end
flagSet = 0;
k       = idxDfMin;
while(flagSet==0 && k > idxFMax)
   if(expForceData(k,2) >= expForceData(k-1,2) && flagSet==0)
       actT1 = expForceData(k+1,1);
       idxActT1 = k+1;
       flagSet =1;
   end
   k = k-1;
end

%%
% Get the ramp start time
%%
[maxDlVal, idxMaxDL] = max(expLengthDataFilt(:,3));

clusterDL = kmeans(expLengthDataFilt(:,3),2);

k         = idxMaxDL;
flagSet   = 0;
lenT0     = 0;
idxRampT0 = 0;
cluster0 = clusterDL(idxMaxDL,1);

while(flagSet==0 && k > 2)
   if(clusterDL(k,1) ~= cluster0 && flagSet == 0)
    lenT0     = expLengthDataFilt(k,1);
    idxRampT0 = k;
    flagSet   = 1;
   end
   k=k-1;
end

k         = idxMaxDL;
flagSet   = 0;
lenT1     = 0;
idxRampT1 = 0;
while(flagSet==0 && k < size(expLengthData,1))
   if(clusterDL(k,1) ~= cluster0 && flagSet == 0)
    lenT1     = expLengthDataFilt(k,1);
    idxRampT1 = k;
    flagSet   = 1;
   end
   k=k+1;
end

dl = (expLengthData(idxRampT1,2)-expLengthData(idxRampT0,2))/1000;
dt =  expLengthData(idxRampT1,1)-expLengthData(idxRampT0,1);

assert( ((dl/dt) - rampSlope)/rampSlope < 0.05);

flagPlotFilteredSignals =0;
if(flagPlotFilteredSignals==1)
   figFilt = figure;
   subplot(2,1,1);
       plot(expForceData(:,1),expForceData(:,2),...
            'Color',[0.25,0.25,0.25],'LineWidth',2);
       hold on;
       plot(expData(:,1),expData(:,colForcePassive),...
            'Color',[1,0,1],'LineWidth',2);
       hold on;       
       plot(expForceDataFilt(:,1),expForceDataFilt(:,2),'r');
       hold on;
       plot(expForceDataFilt(:,1),expForceDataFilt(:,3),'b');
       hold on;
       plot(expForceData(idxActT0,1),expForceData(idxActT0,2),'x');
       hold on;
       plot(expForceData(idxActT1,1),expForceData(idxActT1,2),'x');
       hold on;
       plot(expForceData(idxRampT0-1,1),expForceData(idxRampT0-1,2),'x');
       hold on;

       xlabel('Time (s)');
       ylabel('Force (N)');
   subplot(2,1,2);
       plot(expLengthData(:,1),expLengthData(:,2),...
            'Color',[0.25,0.25,0.25],'LineWidth',2);
       hold on;
       plot(expData(:,1),expData(:,colLengthPassive),...
            'Color',[1,0,1],'LineWidth',2);
       hold on;       
       plot(expLengthDataFilt(:,1),expLengthDataFilt(:,2),'r');
       hold on;
       plot(expLengthDataFilt(:,1),expLengthDataFilt(:,3),'b');
       hold on;      
       plot([expLengthData(idxRampT0,1);expLengthData(idxRampT0,1)],...
            [-1;10],'-k');
       hold on;
       plot([expLengthData(idxRampT1,1);expLengthData(idxRampT1,1)],...
            [-1;10],'-k');
       hold on;

       xlabel('Time (s)');
       ylabel('Length (mm)');          
end


%%
%
%%

expConfigHerzogLeonard2002 = ...
  struct('figureNumber'   ,     figureNumber,...
         'subFigureNumber',     subFigureNumber,...
         'trialNumber'    ,     trialNumber,...
         'timeSpan'       ,     [min(expLengthData(:,1)),max(expLengthData(:,1))],...
         'lengthRampKeyPoints', [ expLengthData(idxRampT0,1),  expLengthData(idxRampT0,2)./1000;...
                                  expLengthData(idxRampT1,1),  expLengthData(idxRampT1,2)./1000 ],...
         'stimulationKeyTimes', [actT0, actT1],...
         'nominalForce' ,    expForceData(idxRampT0-1,2),...
         'dataRamp', [],...
         'dataStatic',[],...
         'dataPassive',[]);
       
expConfigHerzogLeonard2002.dataRamp = ...
    struct('time',    expLengthData(:,1),...
           'length',  expLengthData(:,2),...
           'force',   expForceData(:,2));
         
expConfigHerzogLeonard2002.dataStatic = ...
    struct('time',    expForceStaticData(:,1),...
           'length',  ones(size(expForceStaticData(:,2))).*expLengthData(idxRampT1,2),...
           'force',   expForceStaticData(:,2));

expConfigHerzogLeonard2002.dataPassive = ...
    struct('time',    expData(:,1),...
           'length',  expData(:,colLengthPassive),...
           'force',   expData(:,colForcePassive));         




