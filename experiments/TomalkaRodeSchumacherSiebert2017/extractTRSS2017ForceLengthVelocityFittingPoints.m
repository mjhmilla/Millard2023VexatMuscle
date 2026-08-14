function [flN, fvN]=extractTRSS2017ForceLengthVelocityFittingPoints(...
                       expTRSS2017,...
                       expTRSS2017Sets,...
                       filterFrequencyBeforeNumericalDerivative,...
                       sampleFrequencyBeforeNumericalDerivative,...
                       flagGenerateDebuggingPlot,...
                       verbose)


[b,a] = butter(2,filterFrequencyBeforeNumericalDerivative...
               /(0.5*sampleFrequencyBeforeNumericalDerivative),'low');


for i=1:1:length(expTRSS2017Sets)

  vN = filtfilt(b,a,expTRSS2017.(expTRSS2017Sets{i}).lN);
  dt =1/expTRSS2017.sampleFrequencyHz;
  timeSeries = [0:1:(length(vN)-1)]' .* dt;

  vN = calcCentralDifferenceDataSeries(timeSeries,vN);

  idxA = expTRSS2017.keyIndices.(expTRSS2017Sets{i})(1);
  idxB = expTRSS2017.keyIndices.(expTRSS2017Sets{i})(3);
  if(i==1)
    flN.x = expTRSS2017.(expTRSS2017Sets{i}).lN(idxA); 
    flN.y = expTRSS2017.(expTRSS2017Sets{i}).fNavg(idxA);    
    fvN.x = vN(idxB);

    %This is an approximation (because the length has changed we should
    % account for the change in flN, which we ignore because the change 
    % is so small) but it works in this case 
    fvN.y = expTRSS2017.(expTRSS2017Sets{i}).fNavg(idxB) ...
           /expTRSS2017.(expTRSS2017Sets{i}).fNavg(idxA);
  else
    flN.x = [flN.x;expTRSS2017.(expTRSS2017Sets{i}).lN(idxA)]; 
    flN.y = [flN.y;expTRSS2017.(expTRSS2017Sets{i}).fNavg(idxA)];

    fvNVal= expTRSS2017.(expTRSS2017Sets{i}).fNavg(idxB) ...
           /expTRSS2017.(expTRSS2017Sets{i}).fNavg(idxA);
    fvN.x = [fvN.x;vN(idxB)];
    fvN.y = [fvN.y;fvNVal];    
  end


  if(flagGenerateDebuggingPlot==1)
    if(~exist('figDebugFv','var'))
      figDebugFv=figure;
    end
    n=(i-1)/2;
    seriesColor = [0,0,0].*n + [0,0,1].*(1-n);
    subplot(2,1,1);
      plot(expTRSS2017.(expTRSS2017Sets{i}).time,...
        expTRSS2017.(expTRSS2017Sets{i}).lN,'Color',seriesColor);
      hold on;
      plot(expTRSS2017.(expTRSS2017Sets{i}).time(idxA),...
        expTRSS2017.(expTRSS2017Sets{i}).lN(idxA),...
        'o','Color',seriesColor,...
        'MarkerFaceColor',[1,1,1]);
      hold on;
      plot(expTRSS2017.(expTRSS2017Sets{i}).time(idxB),...
        expTRSS2017.(expTRSS2017Sets{i}).lN(idxB),...
        'o','Color',seriesColor,...
        'MarkerFaceColor',seriesColor);
      hold on;
      xlabel('Time (s)');
      ylabel('Norm. Length');
    subplot(2,1,2);
      plot(expTRSS2017.(expTRSS2017Sets{i}).time,...
           vN,'-','Color',seriesColor);
      hold on;
      plot(expTRSS2017.(expTRSS2017Sets{i}).time(idxA),...
           vN(idxA),'o','Color',seriesColor,...
           'MarkerFaceColor',[1,1,1]);
      hold on;
      plot(expTRSS2017.(expTRSS2017Sets{i}).time(idxB),...
           vN(idxB),'o','Color',seriesColor,...
           'MarkerFaceColor',seriesColor);
      hold on;      
      xlabel('Time (s)');
      ylabel('Norm. Velocity');

      paddingPoints = 100;

      vNTrim = vN(paddingPoints:(end-paddingPoints));
      vNMean = mean(vNTrim);
      vNStd  = std(vNTrim);
      vNP = getPercentiles(vNTrim,[0.05,0.25,0.5,0.75,0.95]);
      if(verbose==1)
        if(i==1)
          fprintf('Variation of the ramp velocity\n');
        end
        fprintf('\t(%1.3e, %1.3e)\t(%s, %s) in %s\n',...
                vNMean, vNStd,'mean','std','lps');
        fprintf('\t(%1.3e, %1.3e, %1.3e)\t(%s, %s, %s) in %s\n',...
                vNP(2), vNP(4),(vN(4)-vN(2))*0.5,'25p','75p','0.5*(75-25)','lps');
        if(i==1)
          vNMeanGroup=vNMean;
          vNStdGroup=vNStd;
          vNGroup = vNP;
        else
          vNMeanGroup=vNMeanGroup+vNMean;
          vNStdGroup=vNStdGroup+vNStd;
          vNGroup = vNGroup + vNP;
          if(i==3)
            fprintf('\t(%1.3e, %1.3e)\t%s:(%s, %s) in %s\n',...
                    vNMeanGroup/3, vNStdGroup/3,'group','mean','std', 'lps');          
            fprintf('\t(%1.3e, %1.3e)\t%s:(%s, %s) in %s\n',...
                    (vNMeanGroup/3)/2.25, (vNStdGroup/3)/2.25,'group','norm','std', 'vmax');          
            fprintf('\t(%1.3e, %1.3e)\t%s:(%s, %s) in %s\n',...
                    vNGroup(1)/3, vNGroup(5)/3,'group','5','95','lps');          
            fprintf('\t(%1.3e, %1.3e)\t%s:(%s, %s) in %s\n',...
                    (vNGroup(1)/3)/2.25, (vNGroup(5)/3)/2.25,'group','5','95','vmax');          
            
          end
        end
      end
  end

end