function data = loadDigitizedData(fileName, xLabel, yLabel, seriesLabels, dataTitle)

data(length(seriesLabels)) = struct('xName','','yName','',...
                                    'x',[],'y',[],...
                                    'seriesName','','title',dataTitle);

fid = fopen(fileName);

for i=1:1:length(seriesLabels)
  data(i).xName = xLabel;
  data(i).yName = yLabel;
  data(i).seriesName = seriesLabels{i};

  line =fgetl(fid);  
  while( ischar(line) && contains(line,','))
      
    idxComma = strfind(line,',');
    assert(length(idxComma)==1,...
        'Error: loadDigitized data can only handle files with 2 data columns');

    fieldAStr = line(1,1:1:(idxComma-1));
    fieldBStr = line((idxComma+1):1:end);

    fieldA = str2double(fieldAStr);
    fieldB = str2double(fieldBStr);

    if(~isnan(fieldA) && ~isnan(fieldB))
        data(i).x = [data(i).x;fieldA];
        data(i).y = [data(i).y;fieldB];        
    end
   
    line =fgetl(fid);
    
  end


end
% MM 2023/04/05
% This works on an installation of Matlab running from Ubuntu, but
% not in windows: textscan on the windows installation only returns the
% first data point; on ubuntu it returns all of the data points. 
%
% for i=1:1:length(seriesLabels)
%   labels = textscan(fid,'%s,%s\n');
%   series = textscan(fid,'%f,%f\n');
%   
%   labelStr = labels{1}{1};
%   idComma = strfind(labelStr,',');
%   
%   data(i).xName = xLabel;
%   data(i).yName = yLabel;
%   data(i).x = series{1};
%   data(i).y = series{2};
%   data(i).seriesName = seriesLabels{i};
%   
% end

fclose(fid);