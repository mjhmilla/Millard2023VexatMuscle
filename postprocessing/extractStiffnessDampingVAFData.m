function [kbr1994Table,opus31Table,hillTable] ...
    = extractStiffnessDampingVAFData( dataFolder,...
                                      freqSeriesFiles,...                                       
                                      dataKBR1994Fig9A,...
                                      dataKBR1994Fig9B,... 
                                      dataKBR1994Fig10,...
                                      dataKBR1994Fig12K,...
                                      dataKBR1994Fig12D,...
                                      inputFunctions)

freqTable = [15,35,90];
ampTable = [0.4,0.8,1.6];
noDataCode = -42;

kbr1994Data(9) = struct('amplitude',      0,...
                        'frequency',      0,...
                        'stiffness_x',    zeros(1,1),...
                        'stiffness_y',    zeros(1,1),...                        
                        'damping_x',      zeros(1,1),...
                        'damping_y',      zeros(1,1),...                                                
                        'warning',        0,...
                        'warningMessage' , '');
                      
kbr1994Stiffness(9) = struct(  'fo', [],...
                                'g', [],...
                                'x', [],...
                                'y', [],...
                                'yN',[]);

kbr1994Damping(9)   = struct(  'fo', [],...
                                'g', [],...
                                'x', [],...
                                'y', [],...
                                'yN',[]);

                      
modelData(9) = struct(  'amplitude',      0,...
                        'frequency',      0,...
                        'stiffness_x',    [],...
                        'stiffness_y',    [],...                        
                        'damping_x',      [],...
                        'damping_y',      [],...
                        'vaf_x',          [],...
                        'vaf_y',          []);
                      
modelStiffness(9) = struct(  'fo', [],...
                              'g', [],...
                              'x',[],...
                              'y',[],...
                              'yN',[]);    
                            
modelDamping(9)   = struct(  'fo', [],...
                              'g', [],...
                              'x', [],...
                              'y', [],...
                              'yN',[]);
                            
modelVaf(9)   = struct(  'fo',[],...
                         'g', [], ...
                         'x', [],...
                         'y', [],...
                         'yN',[]);
                                        
hillTable   = struct('stiffness',[],'damping',[],'vaf',[],'vafData',[]);
opus31Table = struct('stiffness',[],'damping',[],'vaf',[],'vafData',[]);
                      
                  
for i=1:1:length(inputFunctions.amplitudeMM)
  kbr1994(i).amplitude = inputFunctions.amplitudeMM(1,i);
  kbr1994(i).frequency = inputFunctions.bandwidthHz(1,i);
  
  model(i).amplitude = inputFunctions.amplitudeMM(1,i);
  model(i).frequency = inputFunctions.bandwidthHz(1,i);  
  
  hillModel(i).amplitude = inputFunctions.amplitudeMM(1,i);
  hillModel(i).frequency = inputFunctions.bandwidthHz(1,i);    
end
  
for z=1:1:length(freqSeriesFiles)

  tmp = load([dataFolder, freqSeriesFiles{1,z}]);
  freqSimData = tmp.freqSimData;


  flag_Hill = 0;
  if(isempty(strfind(freqSeriesFiles{1,z},'Hill'))==0)
    flag_Hill = 1;
  end
  
  rawData      = modelData;
  currStiffness = modelStiffness;
  currDamping   = modelDamping;
  currVaf       = modelVaf;

 

  for i = 1:1:length(freqSimData.amplitudeMM)
    idx = getIndexIntoVectors(freqSimData.amplitudeMM(1,i),...
                              freqSimData.bandwidthHz(1,i),...
                              inputFunctions.amplitudeMM,...
                              inputFunctions.bandwidthHz);
                            
    %If data is being concatenated make sure it is to the correct trial
    if( rawData(idx).amplitude ~= 0)
      assert(rawData(idx).amplitude == freqSimData.amplitudeMM(1,i));
      assert(rawData(idx).frequency == freqSimData.bandwidthHz(1,i));
    else
      rawData(idx).amplitude = freqSimData.amplitudeMM(1,i);
      rawData(idx).frequency = freqSimData.bandwidthHz(1,i);      
    end
                            
    rawData(idx).stiffness_x = [rawData(idx).stiffness_x;...
                                 freqSimData.nominalForce(1,i)];
    rawData(idx).stiffness_y = [rawData(idx).stiffness_y;...
                                 freqSimData.stiffness(1,i)/1000];
                               
    rawData(idx).damping_x = [rawData(idx).damping_x;...
                                 freqSimData.nominalForce(1,i)];
    rawData(idx).damping_y = [rawData(idx).damping_y;...
                                 freqSimData.damping(1,i)/1000];   
    
    rawData(idx).vaf_x = [rawData(idx).vaf_x;...
                           freqSimData.nominalForce(1,i)];
    rawData(idx).vaf_y = [rawData(idx).vaf_y;...
                           freqSimData.vafTime(1,i)];
  end
  
  %Order the data and fit a line to it
  for idx=1:1:length(rawData)
    % Stiffness
    data_x = rawData(idx).stiffness_x;
    data_y = rawData(idx).stiffness_y;
    
    [k_x,map] = sort(data_x);
    k_y       = data_y(map);    
    

    
    [fo, g] = fit(k_x,k_y,'poly1');
    
    currStiffness(idx).fo = fo;
    currStiffness(idx).g = g;
    currStiffness(idx).x = k_x;
    currStiffness(idx).y = k_y;
    currStiffness(idx).yN= k_y./k_x;
    
    % Damping    
    data_x = rawData(idx).damping_x;
    data_y = rawData(idx).damping_y;
    
    [d_x,map] = sort(data_x);
    d_y       = data_y(map);    
    
       
    [fo, g] = fit(d_x,d_y,'poly1');
    currDamping(idx).fo = fo;
    currDamping(idx).g = g;
    currDamping(idx).x = d_x;
    currDamping(idx).y = d_y;
    currDamping(idx).yN = d_y./d_x;
    %VAF
    data_x = rawData(idx).vaf_x;
    data_y = rawData(idx).vaf_y;
    
    [vaf_x,map] = sort(data_x);
    vaf_y       = data_y(map);    
    
 
    
    [fo, g] = fit(vaf_x,vaf_y,'poly1');
    currVaf(idx).fo = fo;
    currVaf(idx).g = g;
    currVaf(idx).x = vaf_x;
    currVaf(idx).y = vaf_y;
    currVaf(idx).yN=noDataCode;
    
%     currVaf(idx).mean = mean(vaf_y);
%     currVaf(idx).std = std(vaf_y);
%     currVaf(idx).min = min(vaf_y);
%     currVaf(idx).max = max(vaf_y);
%     currVaf(idx).frequency = freq;
%     currVaf(idx).amplitude = amp;
  end
  
  %Re-arrange the data into tables
  stiffness = extractTablesOfFittedParameters(...
                            inputFunctions.amplitudeMM, ...
                            inputFunctions.bandwidthHz, ...
                            currStiffness, noDataCode); 
  damping = extractTablesOfFittedParameters(...
                            inputFunctions.amplitudeMM, ...
                            inputFunctions.bandwidthHz, ...
                            currDamping, noDataCode); 
  vaf = extractTablesOfFittedParameters(...
                            inputFunctions.amplitudeMM, ...
                            inputFunctions.bandwidthHz, ...
                            currVaf, noDataCode); 
  
  fName = freqSeriesFiles{1,z};
  i = strfind(fName,'freqResponse');
  fNameTable   = ['table',fName((i+12):1:end)];
  %dataFolder,
  
  vafData = currVaf;
  save([dataFolder,fNameTable],'stiffness','damping','vaf');
  
  if(flag_Hill ==1)
    hillTable.stiffness     = stiffness;
    hillTable.damping       = damping;
    hillTable.vaf           = vaf;
  else
    opus31Table.stiffness   = stiffness;
    opus31Table.damping     = damping;
    opus31Table.vaf         = vaf;  
  end
                           
end


%%
% KBR 1994 Data: From Fig. 9a
%%

amp9A = [0.4,0.8,1.6];
freq9A = [15,15,15];
assert(length(dataKBR1994Fig9A)==3);

flag_useFig9A = 1;

if(flag_useFig9A==1)
    for i=1:1:length(freq9A)
      idx = getIndexIntoVectors(amp9A(1,i),freq9A(1,i),...
                            inputFunctions.amplitudeMM, ...
                            inputFunctions.bandwidthHz);
    
      data_x = [kbr1994Stiffness(idx).x;dataKBR1994Fig9A(i).x];
      data_y = [kbr1994Stiffness(idx).y;dataKBR1994Fig9A(i).y];
      [k_x,map] = sort(data_x);
      k_y       = data_y(map);
    
      kbr1994(idx).stiffness_x = k_x;
      kbr1994(idx).stiffness_y = k_y;
    
      [fo, g] = fit(k_x,k_y,'poly1');
      kbr1994Stiffness(idx).fo = fo;
      kbr1994Stiffness(idx).g  = g;
      kbr1994Stiffness(idx).x  = k_x;
      kbr1994Stiffness(idx).y  = k_y;
      kbr1994Stiffness(idx).yN  = k_y./k_x;
    end
end
%%
% KBR 1994 Data: From Fig. 9b
%%
amp9B =  [0.4,0.4,0.4,1.6,1.6,1.6];
freq9B = [ 15, 35, 90, 15, 35, 90];
assert(length(dataKBR1994Fig9B)==6);

flag_useFig9B = 1;

if(flag_useFig9B==1)
    for i=1:1:length(freq9B)
      idx = getIndexIntoVectors(amp9B(1,i),freq9B(1,i),...
                                inputFunctions.amplitudeMM, ...
                                inputFunctions.bandwidthHz);
    
      data_x = [kbr1994Stiffness(idx).x;dataKBR1994Fig9B(i).x];
      data_y = [kbr1994Stiffness(idx).y;dataKBR1994Fig9B(i).y];
      [k_x,map] = sort(data_x);
      k_y       = data_y(map);
    
      kbr1994(idx).stiffness_x = k_x;
      kbr1994(idx).stiffness_y = k_y;
    
      [fo, g] = fit(k_x,k_y,'poly1');
      kbr1994Stiffness(idx).fo = fo;
      kbr1994Stiffness(idx).g  = g;
      kbr1994Stiffness(idx).x  = k_x;
      kbr1994Stiffness(idx).y  = k_y;
      kbr1994Stiffness(idx).yN = k_y./k_x;
    end
end

%%
% KBR 1994 Data: From Fig. 10
%%
amp10 =  [0.4,0.4,0.4, 0.8,0.8,0.8, 1.6,1.6,1.6];
freq10 = [ 15, 35, 90,  15, 35, 90,  15, 35, 90];

% This figure is a bit of a special case: data from the different
% perturbations are identically labelled. Kirsch et al. did this to make
% a point that damping varies with purturbation frequency. However this 
% means that I have to be a bit careful: for this data structure it makes
% the most sense to copy the data from Fig. 10 into the 3 different
% perturbation sections. This is fine for generating a table, but it is
% not fine for doing statistics unless this data is also specially handled
% as a 'grouped' class.

flag_useFig10 = 1;

if(flag_useFig10==1)    
    for i=1:1:length(amp10)
      j = 0;
      switch freq10(i)
        case 15
          j=1;
        case 35
          j=2;
        case 90
          j=3;
        otherwise
          assert(0,'Invalid frequency');
      end
      
      idx = getIndexIntoVectors(amp10(1,i),freq10(1,i),...
                            inputFunctions.amplitudeMM, ...
                            inputFunctions.bandwidthHz);
    
      data_x    = [kbr1994Damping(idx).x;dataKBR1994Fig10(j).x];
      data_y    = [kbr1994Damping(idx).y;dataKBR1994Fig10(j).y];
      [d_x,map] = sort(data_x);
      d_y       = data_y(map);
    
      kbr1994(idx).damping_x = d_x;
      kbr1994(idx).damping_y = d_y;
    
      [fo, g] = fit(d_x,d_y,'poly1');
      kbr1994Damping(idx).fo = fo;
      kbr1994Damping(idx).g  = g;
      kbr1994Damping(idx).x  = d_x;
      kbr1994Damping(idx).y  = d_y;
      kbr1994Damping(idx).yN = d_y./d_x;
      
      kbr1994(idx).warning = 1;
      kbr1994(idx).warningMessage = 'Grouped data from Fig. 10';
    end
end
%%
% KBR 1994 Data: From Fig. 12
%%
flag_useFig12 = 1;

if(flag_useFig12 == 1)
  idx = getIndexIntoVectors(0.8,35,...
                            inputFunctions.amplitudeMM, ...
                            inputFunctions.bandwidthHz);

  assert(length(dataKBR1994Fig12K)==2);

  data_x = [kbr1994Stiffness(idx).x;dataKBR1994Fig12K(1).x;dataKBR1994Fig12K(2).x];
  data_y = [kbr1994Stiffness(idx).y;dataKBR1994Fig12K(1).y;dataKBR1994Fig12K(2).y];
  [k_x,map] = sort(data_x);
  k_y       = data_y(map);

  kbr1994(idx).stiffness_x = k_x;
  kbr1994(idx).stiffness_y = k_y;

  [fo, g] = fit(k_x,k_y,'poly1');
  kbr1994Stiffness(idx).fo = fo;
  kbr1994Stiffness(idx).g  = g;
  kbr1994Stiffness(idx).x = k_x;
  kbr1994Stiffness(idx).y = k_y;
  kbr1994Stiffness(idx).yN= k_y./k_x;

  data_x = [kbr1994Damping(idx).x;dataKBR1994Fig12D(1).x;dataKBR1994Fig12D(2).x];
  data_y = [kbr1994Damping(idx).y;dataKBR1994Fig12D(1).y;dataKBR1994Fig12D(2).y];
  [d_x,map] = sort(data_x);
  d_y       = data_y(map);

  kbr1994(idx).damping_x = d_x;
  kbr1994(idx).damping_y = d_y;

  [fo, g] = fit(d_x,d_y,'poly1');
  kbr1994Damping(idx).fo = fo;
  kbr1994Damping(idx).g  = g;
  kbr1994Damping(idx).x  = d_x;
  kbr1994Damping(idx).y  = d_y;
  kbr1994Damping(idx).yN = d_y./d_x;
end

for idx=1:1:length(kbr1994Stiffness)

    if(isempty(kbr1994Stiffness(idx).x))
        kbr1994Stiffness(idx).x = noDataCode;
        kbr1994Stiffness(idx).y = noDataCode;
        kbr1994Stiffness(idx).yN= noDataCode;
        
    end
    
    if(isempty(kbr1994Damping(idx).x))
        kbr1994Damping(idx).x = noDataCode;
        kbr1994Damping(idx).y = noDataCode;
        kbr1994Damping(idx).yN = noDataCode;
    end
end

%%
% Arrange the fitting information into a tabular form
%%

kbr1994Table = struct('stiffness',[],'damping',[],'vaf',[],'vafData',[]);





kbr1994Table.stiffness = extractTablesOfFittedParameters(...
                          inputFunctions.amplitudeMM, ...
                          inputFunctions.bandwidthHz, ...
                          kbr1994Stiffness, noDataCode);
                        
kbr1994Table.damping = extractTablesOfFittedParameters(...
                          inputFunctions.amplitudeMM, ...
                          inputFunctions.bandwidthHz, ...
                          kbr1994Damping, noDataCode);

kbr1994Table.vaf = opus31Table.vaf;

kbr1994Table.vaf.pMean = ones(size(kbr1994Table.vaf.pMean)).*(0.5*(0.88+0.99));
kbr1994Table.vaf.pMean(:,:,1)=ones(size(kbr1994Table.vaf.pMean(:,:,1))).*noDataCode;

kbr1994Table.vaf.p95CIMin = ones(size(kbr1994Table.vaf.p95CIMin)).*(0.88);
kbr1994Table.vaf.p95CIMin(:,:,1)=ones(size(kbr1994Table.vaf.p95CIMin(:,:,1))).*noDataCode;

kbr1994Table.vaf.p95CIMax = ones(size(kbr1994Table.vaf.p95CIMax)).*noDataCode;
kbr1994Table.vaf.p95CIMax(:,:,1)=ones(size(kbr1994Table.vaf.p95CIMax(:,:,1))).*noDataCode;

kbr1994Table.vaf.rmse = ones(size(kbr1994Table.vaf.rmse)).*noDataCode;

for i=1:1:3
    for j=1:1:3
        kbr1994Table.vaf.data(i,j).x = [1,1].*noDataCode;
        kbr1994Table.vaf.data(i,j).y = [0.88,0.99];
        kbr1994Table.vaf.data(i,j).yN = [1,1].*noDataCode;
    end
end




