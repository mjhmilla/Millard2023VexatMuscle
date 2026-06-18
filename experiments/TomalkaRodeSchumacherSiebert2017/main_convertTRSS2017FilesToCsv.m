clc;
clear all;
close all;

flag_plot=1;
flag_checkMean=0; %Will plot the means stored in the mat files in magenta
                  %If everything is correct (and it is) the computed
                  %mean (in black) will not be visible.

fileList = {'Matt_0_85_1_3_L0','085','0.85','mixedLength','shortOutput';...
            'Matt_1_1_45_L0','1','1.45','uniformLength','shortOutput';...
            'Matt_07_1_15_L0' ,'07','0.70','mixedLength','shortOutput';...
            'Matt_1_1_45_L0','1','1.45','uniformLength','longOutput'};

dirTRSS2017 = pwd;
assert(contains(dirTRSS2017,'TomalkaRodeSchumacherSiebert2017'));

vmax              = 2.25;
sampleFrequencyHz = 1000; % Check with Andre
timeRampStart     = 5;    % Check with Andre
sampleRampStart   = timeRampStart*sampleFrequencyHz;
sampleRampTotal   = 1801;

set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');



if(flag_plot==1)
  figH = figure;
end

for i=1:1:size(fileList,1)
  data = load([fileList{i,1},'.mat']);
  fieldList = fields(data);

  nMeasurments=length(fieldList)-2;

  csvData =[];  
  forceColumnNames = [];
  for j=1:1:nMeasurments
    idxStr = num2str(j);
    if(length(idxStr)<2)
      idxStr=['0',idxStr];
    end
    forceColumnNames = [forceColumnNames; {['fN',idxStr]} ];
  end
  
  csvDataHeader = [{'time'}; {'lN'}; forceColumnNames;...
                    {'fNavg'} ;{'fNstd'};{'fNmin'};{'fNmax'}];  

  lengthField = ['FL_',fileList{i,2}];

  idxA = 0;
  idxB = 0;
  switch fileList{i,5}
    case 'shortOutput'
      idxA = sampleRampStart;
      if(strcmp(fileList{i,4},'uniformLength'))
        idxB = idxA+sampleRampTotal-1;
      else
        idxB = idxA+length(data.(lengthField))-1;
      end      
    case 'longOutput'
      idxA = 1;
      idxB = idxA+length(data.(lengthField))-1;
    otherwise
      assert(0,'Error: unrecognized data length option');
  end

  csvData = zeros((idxB-idxA+1),length(csvDataHeader));


  idxData=1;
  idxTime=idxData;
  csvData(:,idxData) = [idxA:idxB]./sampleFrequencyHz;
  idxData=idxData+1;
  idxLength=idxData;
  if(strcmp(fileList{i,4},'uniformLength'))
    csvData(:,idxData) = data.(lengthField)(idxA:idxB);
  else
    csvData(:,idxData) = data.(lengthField);
  end
  n=0;
  fNavg = zeros(size(csvData,1),1);
  idx_fN01=nan;
  idx_fN24=nan;
  for j=1:1:length(forceColumnNames)
    idxData=idxData+1;
    if(j==1)
      idx_fN01=idxData;
    end
    if(j==length(forceColumnNames))
      idx_fN24=idxData;
    end
    
    n = n+1;
    fieldName = ['LR',num2str(j)];
    csvData(:,idxData) = data.(fieldName)(idxA:1:idxB);
    fNavg = fNavg + csvData(:,idxData); 

    if(flag_plot==1)
      figure(figH);
      subplot(size(fileList,1),2,(i-1)*2+1);
        yyaxis left
          if(j==1)
            lineColor = [0,0,1];
            plot(csvData(:,idxTime),csvData(:,idxLength),'-','Color',lineColor);
            hold on;
            ylabel('Length ($$\ell_o^M$$)'); 
            yaxis=ylim;
            ylim([0,3]);
          end
        yyaxis right;
          nColor    = (j-1)/(24-1);
          lineColor = [1,0.5,0.5].*nColor + [1,0.75,0.75].*(1-nColor);
          plot(csvData(:,idxTime),csvData(:,idxData),'-','Color',lineColor);
          hold on;
          text(csvData(end,idxTime),csvData(end,idxData),sprintf('%i',j),...
               'FontSize',6,'HorizontalAlignment','left');
          hold on;
          yaxis = ylim;
          ylim([0,yaxis(2)]);
          xlabel('Time (s)');
          ylabel('Force ($$f_o^M$$)');
          box off;
          title(sprintf('TRSS2017 %s %s',fileList{i,3},'$$\ell_o^M$$'));

      subplot(size(fileList,1),2,(i-1)*2+2);
          plot(csvData(:,idxLength),csvData(:,idxData),'-','Color',lineColor);
          hold on;
          text(csvData(end,idxLength),csvData(end,idxData),sprintf('%i',j),...
               'FontSize',6,'HorizontalAlignment','left');
          hold on;
          yaxis = ylim;
          ylim([0,yaxis(2)]);
          
          xlabel('Length ($$\ell_o^M$$)');
          ylabel('Force ($$f_o^M$$)');
          box off;   
          title(sprintf('TRSS2017 %s %s',fileList{i,3},'$$\ell_o^M$$'));

    end
  end

  fNavg = fNavg .* (1/n);

  idxData=idxData+1;
  idxAvg=idxData;
  csvData(:,idxData)=fNavg;

  if(flag_plot==1)
    figure(figH);
    subplot(size(fileList,1),2,(i-1)*2+1);
      lineColor = [0,0,0];
      plot(csvData(:,idxTime),csvData(:,idxAvg),'-','Color',lineColor);
      hold on;
  
    subplot(size(fileList,1),2,(i-1)*2+2);
      lineColor = [0,0,0];
      plot(csvData(:,idxLength),csvData(:,idxAvg),'-','Color',lineColor);
      hold on;
  end

  fNstd = zeros(size(fNavg));
  for j=1:1:length(forceColumnNames)
    idx = idx_fN01 + (j-1);
    fNstd = fNstd + (csvData(:,idx)-fNavg).^2;
  end
  fNstd = fNstd .* (1/n);
  fNstd = sqrt(fNstd);

  idxData=idxData+1;
  idxStd=idxData;
  csvData(:,idxData)=fNstd;

  idxData=idxData+1;
  idxMin = idxData;
  idxData=idxData+1;
  idxMax = idxData;
  
  fNmin = zeros(size(fNavg));
  fNmax = zeros(size(fNavg));
  
  for j=1:1:length(forceColumnNames)
    idx = idx_fN01 + (j-1);
    if(j==1)
      fNmin=csvData(:,idx);
      fNmax=csvData(:,idx);
    end
    fNmin = min(fNmin,csvData(:,idx));
    fNmax = max(fNmax,csvData(:,idx));    
  end

  csvData(:,idxMin)=fNmin;
  csvData(:,idxMax)=fNmax;

  if(flag_plot==1)
    figure(figH);
    subplot(size(fileList,1),2,(i-1)*2+1);
      yyaxis right;
      lineColor = [0,0,0];
      plot(csvData(:,idxTime),csvData(:,idxAvg)+csvData(:,idxStd),'--','Color',lineColor);
      hold on;
      plot(csvData(:,idxTime),csvData(:,idxAvg)-csvData(:,idxStd),'--','Color',lineColor);
      hold on;
      plot(csvData(:,idxTime),csvData(:,idxMin),'-','Color',lineColor,'LineWidth',2);
      hold on;
      plot(csvData(:,idxTime),csvData(:,idxMax),'-','Color',lineColor,'LineWidth',2);
      hold on;

      if(flag_checkMean==1)
        if(strcmp(fileList{i,4},'uniformLength'))
          plot(csvData(:,idxTime),data.(['mean_LR_',fileList{i,2}])(idxA:idxB),'-','Color',[1,0,1]);
          hold on;        
        else
          plot(csvData(:,idxTime),data.(['mean_LR_',fileList{i,2}]),'-','Color',[1,0,1]);
          hold on;
        end      
      end
      ylim([0,3]);


    subplot(size(fileList,1),2,(i-1)*2+2);
      lineColor = [0,0,0];
      plot(csvData(:,idxLength),csvData(:,idxAvg)+csvData(:,idxStd),'--','Color',lineColor);
      hold on;
      plot(csvData(:,idxLength),csvData(:,idxAvg)-csvData(:,idxStd),'--','Color',lineColor);
      hold on;
      plot(csvData(:,idxLength),csvData(:,idxMin),'-','Color',lineColor,'LineWidth',2);
      hold on;
      plot(csvData(:,idxLength),csvData(:,idxMax),'-','Color',lineColor,'LineWidth',2);
      hold on;
      
      if(flag_checkMean==1)
        if(strcmp(fileList{i,4},'uniformLength'))
          plot(csvData(:,idxLength),data.(['mean_LR_',fileList{i,2}])(idxA:idxB),'-','Color',[1,0,1]);
          hold on;        
        else
          plot(csvData(:,idxLength),data.(['mean_LR_',fileList{i,2}]),'-','Color',[1,0,1]);
          hold on;
        end      
      end
      xlim([0.7,1.45]);
      ylim([0,3]);

  end

  fileName = strrep(fileList{i,3},'.','p');
  fileName = ['data_TRSS2017_',fileName,'_',fileList{i,5},'.csv'];
  fid = fopen(fullfile(dirTRSS2017,fileName),'w');
  for j=1:1:length(csvDataHeader)
    if(j==1)
      fprintf(fid,'%s',csvDataHeader{j});
    else
      fprintf(fid,',%s',csvDataHeader{j});      
    end
  end
  fprintf(fid,'\n');

  for rows=1:1:size(csvData,1)
    for cols=1:1:size(csvData,2)
      if(cols==1)
        fprintf(fid,'%1.6e',csvData(rows,cols));
      else
        fprintf(fid,',%1.6e',csvData(rows,cols));        
      end      
    end  
    fprintf(fid,'\n');
  end

  fclose(fid);

end