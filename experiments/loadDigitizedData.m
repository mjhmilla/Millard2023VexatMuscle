function data = loadDigitizedData(fileName, xLabel, yLabel, seriesLabels, dataTitle)

data(length(seriesLabels)) = struct('xName','','yName','',...
                                    'x',zeros(1,1),'y',zeros(1,1),...
                                    'seriesName','','title',dataTitle);

fid = fopen(fileName);

for i=1:1:length(seriesLabels)
  labels = textscan(fid,'%s,%s\n');
  series = textscan(fid,'%f,%f\n');
  
  labelStr = labels{1}{1};
  idComma = strfind(labelStr,',');
  
  data(i).xName = xLabel;
  data(i).yName = yLabel;
  data(i).x = series{1};
  data(i).y = series{2};
  data(i).seriesName = seriesLabels{i};
  
end

fclose(fid);