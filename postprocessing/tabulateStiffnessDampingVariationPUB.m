function success = tabulateStiffnessDampingVariationPUB(dataFolder,...
                                              freqSeriesFiles,...
                                              freqSeriesName,...
                                              freqSeriesColor,...
                                              inputFunctions,...                      
                                              normFiberLength,...
                                              nominalForce,...  
                                              dataKBR1994Fig9A,...
                                              dataKBR1994Fig9B,... 
                                              dataKBR1994Fig10,...
                                              dataKBR1994Fig12K,...
                                              dataKBR1994Fig12D,...
                                              flag_useElasticTendon,... 
                                              tableNameEnding,...
                                              outputFolder)
success=0;                                            
               
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
                                'g', []);                              
kbr1994Damping(9)   = struct(  'fo', [],...
                                'g', []);

                      
modelData(9) = struct(  'amplitude',      0,...
                        'frequency',      0,...
                        'stiffness_x',    [],...
                        'stiffness_y',    [],...                        
                        'damping_x',      [],...
                        'damping_y',      [],...
                        'vaf_x',          [],...
                        'vaf_y',          []);
                      
modelStiffness(9) = struct(  'fo', [],...
                              'g', []);    
                            
modelDamping(9)   = struct(  'fo', [],...
                              'g', []);
                            
modelVaf(9)   = struct(  'fo', [],...
                         'g', [], ...
                         'mean', 0, ...
                         'std', 0,...
                         'min', 0,...
                         'max', 0);
                                        
hillTable   = struct('stiffness',[],'damping',[],'vaf',[],'vafData',[]);
opus31Table = struct('stiffness',[],'damping',[],'vaf',[],'vafData',[]);
                      
freq = [ 15,  15,  15,  35,  35,  35,  90,  90,  90];
amp  = [0.4, 0.8, 1.6, 0.4, 0.8, 1.6, 0.4, 0.8, 1.6];
                      
for i=1:1:length(freq)
  kbr1994(i).amplitude = amp(i);
  kbr1994(i).frequency = freq(i);
  
  model(i).amplitude = amp(i);
  model(i).frequency = freq(i);  
  
  hillModel(i).amplitude = amp(i);
  hillModel(i).frequency = freq(i);    
end
  
for z=1:1:length(freqSeriesFiles)

  tmp = load([dataFolder, freqSeriesFiles{1,z}]);
  freqSimData = tmp.freqSimData;


  flag_Hill = 0;
  if(isempty(strfind(freqSeriesFiles{1,z},'Hill'))==0)
    flag_Hill = 1;
  end
  
  currData      = modelData;
  currStiffness = modelStiffness;
  currDamping   = modelDamping;
  currVaf       = modelVaf;

 

  for i = 1:1:length(freqSimData.amplitudeMM)
    idx = getIndexIntoVectors(freqSimData.amplitudeMM(1,i),...
                              freqSimData.bandwidthHz(1,i),amp,freq);
                            
    %If data is being concatenated make sure it is to the correct trial
    if( currData(idx).amplitude ~= 0)
      assert(currData(idx).amplitude == freqSimData.amplitudeMM(1,i));
      assert(currData(idx).frequency == freqSimData.bandwidthHz(1,i));
    else
      currData(idx).amplitude = freqSimData.amplitudeMM(1,i);
      currData(idx).frequency = freqSimData.bandwidthHz(1,i);      
    end
                            
    currData(idx).stiffness_x = [currData(idx).stiffness_x;...
                                 freqSimData.nominalForce(1,i)];
    currData(idx).stiffness_y = [currData(idx).stiffness_y;...
                                 freqSimData.stiffness(1,i)/1000];
                               
    currData(idx).damping_x = [currData(idx).damping_x;...
                                 freqSimData.nominalForce(1,i)];
    currData(idx).damping_y = [currData(idx).damping_y;...
                                 freqSimData.damping(1,i)/1000];   
    
    currData(idx).vaf_x = [currData(idx).vaf_x;...
                           freqSimData.nominalForce(1,i)];
    currData(idx).vaf_y = [currData(idx).vaf_y;...
                           freqSimData.vafTime(1,i)];
  end
  
  %Order the data and fit a line to it
  for idx=1:1:length(currData)
    % Stiffness
    data_x = currData(idx).stiffness_x;
    data_y = currData(idx).stiffness_y;
    
    [k_x,map] = sort(data_x);
    k_y       = data_y(map);    
    
    currData(idx).stiffness_x = k_x;
    currData(idx).stiffness_y = k_y;
    
    [fo, g] = fit(k_x,k_y,'poly1');
    
    currStiffness(idx).fo = fo;
    currStiffness(idx).g = g;
    
    % Damping    
    data_x = currData(idx).damping_x;
    data_y = currData(idx).damping_y;
    
    [d_x,map] = sort(data_x);
    d_y       = data_y(map);    
    
    currData(idx).damping_x = d_x;
    currData(idx).damping_y = d_y;    
    
    [fo, g] = fit(d_x,d_y,'poly1');
    currDamping(idx).fo = fo;
    currDamping(idx).g = g;
    
    %VAF
    data_x = currData(idx).vaf_x;
    data_y = currData(idx).vaf_y;
    
    [vaf_x,map] = sort(data_x);
    vaf_y       = data_y(map);    
    
    currData(idx).vaf_x = vaf_x;
    currData(idx).vaf_y = vaf_y;    
    
    [fo, g] = fit(vaf_x,vaf_y,'poly1');
    currVaf(idx).fo = fo;
    currVaf(idx).g = g;
    currVaf(idx).mean = mean(vaf_y);
    currVaf(idx).std = std(vaf_y);
    currVaf(idx).min = min(vaf_y);
    currVaf(idx).max = max(vaf_y);
  end
  
  %Re-arrange the data into tables
  stiffness = extractTablesOfFittedParameters(...
                            ampTable, freqTable, amp, freq, ...
                            currStiffness, noDataCode); 
  damping = extractTablesOfFittedParameters(...
                            ampTable, freqTable, amp, freq, ...
                            currDamping, noDataCode); 
  vaf = extractTablesOfFittedParameters(...
                            ampTable, freqTable, amp, freq, ...
                            currVaf, noDataCode); 
  
  fName = freqSeriesFiles{1,z};
  i = strfind(fName,'freqResponse');
  fNameTable   = ['table',fName((i+12):1:end)];
  %dataFolder,
  
  vafData = currVaf;
  save([dataFolder,fNameTable],'stiffness','damping','vaf','vafData');
  
  if(flag_Hill ==1)
    hillTable.stiffness = stiffness;
    hillTable.damping   = damping;
    hillTable.vaf       = vaf;
    hillTable.vafData   = vafData;
  else
    opus31Table.stiffness = stiffness;
    opus31Table.damping   = damping;
    opus31Table.vaf       = vaf;  
    opus31Table.vafData   = vafData;
  end
                           
end


%%
% KBR 1994 Data: From Fig. 9a
%%

amp9A = [0.4,0.8,1.6];
freq9A = [15,15,15];
assert(length(dataKBR1994Fig9A)==3);

for i=1:1:length(freq9A)
  idx = getIndexIntoVectors(amp9A(1,i),freq9A(1,i),amp,freq);

  data_x = [dataKBR1994Fig9A(i).x];
  data_y = [dataKBR1994Fig9A(i).y];
  [k_x,map] = sort(data_x);
  k_y       = data_y(map);

  kbr1994(idx).stiffness_x = k_x;
  kbr1994(idx).stiffness_y = k_y;

  [fo, g] = fit(k_x,k_y,'poly1');
  kbr1994Stiffness(idx).fo = fo;
  kbr1994Stiffness(idx).g  = g;
end

%%
% KBR 1994 Data: From Fig. 9b
%%
amp9B =  [0.4,0.4,0.4,1.6,1.6,1.6];
freq9B = [ 15, 35, 90, 15, 35, 90];
assert(length(dataKBR1994Fig9B)==6);

for i=1:1:length(freq9B)
  idx = getIndexIntoVectors(amp9B(1,i),freq9B(1,i),amp,freq);

  data_x = [dataKBR1994Fig9B(i).x];
  data_y = [dataKBR1994Fig9B(i).y];
  [k_x,map] = sort(data_x);
  k_y       = data_y(map);

  kbr1994(idx).stiffness_x = k_x;
  kbr1994(idx).stiffness_y = k_y;

  [fo, g] = fit(k_x,k_y,'poly1');
  kbr1994Stiffness(idx).fo = fo;
  kbr1994Stiffness(idx).g  = g;
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
  
  idx = getIndexIntoVectors(amp10(1,i),freq10(1,i),amp,freq);

  data_x    = [dataKBR1994Fig10(j).x];
  data_y    = [dataKBR1994Fig10(j).y];
  [d_x,map] = sort(data_x);
  d_y       = data_y(map);

  kbr1994(idx).damping_x = d_x;
  kbr1994(idx).damping_y = d_y;

  [fo, g] = fit(d_x,d_y,'poly1');
  kbr1994Damping(idx).fo = fo;
  kbr1994Damping(idx).g  = g;
  kbr1994(idx).warning = 1;
  kbr1994(idx).warningMessage = 'Grouped data from Fig. 10';
end

%%
% KBR 1994 Data: From Fig. 12
%%
flag_addInFig12 = 1;

if(flag_addInFig12 == 1)
  idx = getIndexIntoVectors(0.8,35,amp,freq);

  assert(length(dataKBR1994Fig12K)==2);

  data_x = [dataKBR1994Fig12K(1).x;dataKBR1994Fig12K(2).x];
  data_y = [dataKBR1994Fig12K(1).y;dataKBR1994Fig12K(2).y];
  [k_x,map] = sort(data_x);
  k_y       = data_y(map);

  kbr1994(idx).stiffness_x = k_x;
  kbr1994(idx).stiffness_y = k_y;

  [fo, g] = fit(k_x,k_y,'poly1');
  kbr1994Stiffness(idx).fo = fo;
  kbr1994Stiffness(idx).g  = g;

  data_x = [dataKBR1994Fig12D(1).x;dataKBR1994Fig12D(2).x];
  data_y = [dataKBR1994Fig12D(1).y;dataKBR1994Fig12D(2).y];
  [d_x,map] = sort(data_x);
  d_y       = data_y(map);

  kbr1994(idx).damping_x = d_x;
  kbr1994(idx).damping_y = d_y;

  [fo, g] = fit(d_x,d_y,'poly1');
  kbr1994Damping(idx).fo = fo;
  kbr1994Damping(idx).g  = g;
end


%%
% Arrange the fitting information into a tabular form
%%


kbr1994TableStiffness = extractTablesOfFittedParameters(...
                          ampTable, freqTable, amp, freq, ...
                          kbr1994Stiffness, noDataCode);
                        
kbr1994TableDamping = extractTablesOfFittedParameters(...
                          ampTable, freqTable, amp, freq, ...
                          kbr1994Damping, noDataCode);
                        

%%
% Transform the table of fitted variables into a tex file
%
%     KBR 1994    Model      Hill
% K/F 15 35 90    15 35 90   15 35 90
% 0.4 
% 0.8
% 1.6
%
%     KBR 1994    Model      Hill
% D/F 15 35 90    15 35 90   15 35 90
% 0.4 
% 0.8
% 1.6
%                 Model      Hill
% VAF 15 35 90    15 35 90   15 35 90
% 0.4 
% 0.8 
% 1.6
%%

tendonTag = '_ElasticTendon';
if(flag_useElasticTendon==0)
  tendonTag = '_RigidTendon';
end

fid = fopen([outputFolder,'tableStiffnessDampingVaf',tendonTag,'_',tableNameEnding,'.tex'],'w');

switch flag_useElasticTendon
    case 0
        strCaption = 'Mean normalized stiffness coefficients (A.), damping coefficients (B.) and VAF (C.) for models with rigid tendons. All additional details are identical to those of Table \label{tbl:KBR1994Sim_ET} except the tendon of the model is rigid.';
        strLabel = '\label{tbl:KBR1994Sim_RT}';
    case 1
        strCaption = 'Mean normalized stiffness coefficients (A.), damping coefficients (B.) and VAF (C.) for models with elastic tendons. Here the proposed model has been fitted to Figure 12 of Kirsch et al. \cite{Kirsch1994MuscleImpedance}. The impedance experiments at each combination of perturbation amplitude and frequncy have been evaluated at 3 different nominal forces: 2.5N, 5N, and 11.5N. The normalized results presented in the table are the mean values of the 2.5N, 5.0N, and 11.5N simulations. Finally, note that the VAF is evaluated between the model and the spring-damper of best fit to the response of the model, rather than to the response of biological muscle (which was not published by Kirsch et al. \cite{Kirsch1994MuscleImpedance}).';
        strLabel = '\label{tbl:KBR1994Sim_ET}';
    otherwise
        assert(0,'flag_useElasticTendon must be 0 or 1');
end




fprintf(fid,'\\begin{table}[!h]\n');
fprintf(fid,'\\caption{%s %s}\n',...
    strCaption, ...
    strLabel);
fprintf(fid,'\\begin{center}\n');
fprintf(fid,'\\begin{tabular}{r | r r r || r r r || r r r}\n');
%            A    15   35   90   15   35   90   15   35   90  

fprintf(fid,['%s & \\multicolumn{3}{c}{Kirsch et al.} & ',...
            '\\multicolumn{3}{c}{Model} & \\multicolumn{3}{c}{Hill} \\\\ \n'],'');
fprintf(fid,'%s & %dHz & %dHz & %dHz & %dHz & %dHz & %dHz & %dHz & %dHz & %dHz \\\\ \n',...
            'A. \hfill Norm. Stiffness ($\frac{K}{F})$', 15, 35, 90, 15, 35, 90, 15, 35, 90 );

firstLineExtra = '\hline ';          
for i=1:1:length(ampTable)
  
  kbrK    = cell(1,3);
  opus31K = cell(1,3);
  hillK   = cell(1,3);

  for j=1:1:3
    strVal = '';
    if(kbr1994TableStiffness.pMean(i,j,1) == noDataCode)
      strVal = '--';      
    else
      strVal = sprintf('%1.2f',kbr1994TableStiffness.pMean(i,j,1));
    end
    kbrK(1,j) = {strVal};
    
    strVal = '';
    if(opus31Table.stiffness.pMean(i,j,1) == noDataCode)
      strVal = '--';      
    else
      strVal = sprintf('%1.2f',opus31Table.stiffness.pMean(i,j,1));
    end
    opus31K(1,j) = {strVal};
    
    strVal = '';
    if(hillTable.stiffness.pMean(i,j,1) == noDataCode)
      strVal = '--';      
    else
      strVal = sprintf('%1.2f',hillTable.stiffness.pMean(i,j,1));
    end
    hillK(1,j) = {strVal};    
  end
  
  


  fprintf(fid,['%s %1.1f mm & %s & %s & %s ',...
                        '& %s & %s & %s ',...
                        '& %s & %s & %s \\\\ \n'],...
              firstLineExtra, ampTable(i), ...
              kbrK{1,1},    kbrK{1,2},    kbrK{1,3},...
              opus31K{1,1}, opus31K{1,2}, opus31K{1,3}, ...
              hillK{1,1},   hillK{1,2},   hillK{1,3} );
  firstLineExtra = '';    
  
end

%fprintf(fid,'\\end{tabular}\n');
%fprintf(fid,'\\end{center}\n');
%fprintf(fid,'\\end{table}\n');

%fprintf(fid,'\\begin{table}[!h]\n');
%fprintf(fid,'\\caption{%s %s %s}\n',...
%    'Normalized damping coefficients of', ...
%    strCaption,...
%    ['\label{',strLabel,'D}']);
%fprintf(fid,'\\begin{center}\n');
%fprintf(fid,'\\begin{tabular}{r | r r r || r r r || r r r}\n');
%            A    15   35   90   15   35   90   15   35   90  

fprintf(fid,'\\multicolumn{10}{c}{} \\\\ \n');
%fprintf(fid,['%s & \\multicolumn{3}{c}{Kirsch et al.} & ',...
%            '\\multicolumn{3}{c}{Model} & \\multicolumn{3}{c}{Hill} \\\\ \n'],'');
fprintf(fid,'%s & %dHz & %dHz & %dHz & %dHz & %dHz & %dHz & %dHz & %dHz & %dHz \\\\ \n',...
            'B. \hfill Norm. damping $\frac{\beta}{F}$', 15, 35, 90, 15, 35, 90, 15, 35, 90 );


firstLineExtra = '\hline ';          
for i=1:1:length(ampTable)
  
  kbrD    = cell(1,3);
  opus31D = cell(1,3);
  hillD   = cell(1,3);

  for j=1:1:3
    strVal = '';
    if(kbr1994TableDamping.pMean(i,j,1) == noDataCode)
      strVal = '--';      
    else
      strVal = sprintf('%1.4f',kbr1994TableDamping.pMean(i,j,1));
    end
    kbrD(1,j) = {strVal};
    
    strVal = '';
    if(opus31Table.damping.pMean(i,j,1) == noDataCode)
      strVal = '--';      
    else
      strVal = sprintf('%1.4f',opus31Table.damping.pMean(i,j,1));
    end
    opus31D(1,j) = {strVal};
    
    strVal = '';
    if(hillTable.damping.pMean(i,j,1) == noDataCode)
      strVal = '--';      
    else
      strVal = sprintf('%1.4f',hillTable.damping.pMean(i,j,1));
    end
    hillD(1,j) = {strVal};    
  end


  fprintf(fid,['%s %1.1fmm & %s & %s & %s ',...
                        '& %s & %s & %s ',...
                        '& %s & %s & %s \\\\ \n'],...
              firstLineExtra, ampTable(i), ...
              kbrD{1,1},    kbrD{1,2},    kbrD{1,3},...
              opus31D{1,1}, opus31D{1,2}, opus31D{1,3}, ...
              hillD{1,1},   hillD{1,2},   hillD{1,3} );
  firstLineExtra = '';    
  
end

%fprintf(fid,'\\end{tabular}\n');
%fprintf(fid,'\\end{center}\n');
%fprintf(fid,'\\end{table}\n');


%fprintf(fid,'\\begin{table}[!h]\n');
%fprintf(fid,'\\caption{%s %s %s}\n','VAF of', ...
%    strCaption,...
%    ['\label{',strLabel,'VAF}']);
%fprintf(fid,'\\begin{center}\n');
%fprintf(fid,'\\begin{tabular}{r | r r r || r r r || r r r}\n');
%            A    15   35   90   15   35   90   15   35   90  

%fprintf(fid,['%s & \\multicolumn{3}{c}{Kirsch et al.} & ',...
%            '\\multicolumn{3}{c}{Model} & \\multicolumn{3}{c}{Hill} \\\\ \n'],'');

fprintf(fid,'\\multicolumn{10}{c}{} \\\\ \n');

fprintf(fid,'%s & %dHz & %dHz & %dHz & %dHz & %dHz & %dHz & %dHz & %dHz & %dHz \\\\ \n',...
            'C. VAF (\%) \hfill', 15, 35, 90, 15, 35, 90, 15, 35, 90 );
firstLineExtra = '\hline ';          
for i=1:1:length(ampTable)


  kbrVaf    = cell(1,3);
  opus31Vaf = cell(1,3);
  hillVaf   = cell(1,3);

  for j=1:1:3
    
    idx = getIndexIntoVectors(ampTable(1,i),freqTable(1,j),amp,freq);   
    
%     strVal = '';
%     if(kbr1994TableDamping.pMean(i,j,1) == noDataCode)
%       strVal = '--';      
%     else
%       strVal = sprintf('%1.3f',kbr1994TableDamping.pMean(i,j,1));
%     end
    kbrVaf(1,j) = {''};
    
    
    
    strVal = '';
    if(opus31Table.vafData(idx).mean == noDataCode)
      strVal = '--';      
    else
      strVal = sprintf('%1.1f',opus31Table.vafData(idx).mean*100);
    end
    opus31Vaf(1,j) = {strVal};
    
    strVal = '';
    if(hillTable.vafData(idx).mean == noDataCode)
      strVal = '--';      
    else
      strVal = sprintf('%1.1f',hillTable.vafData(idx).mean*100);
    end
    hillVaf(1,j) = {strVal};    
  end


  fprintf(fid,['%s %1.1fmm & %s & %s & %s ',...
                        '& %s & %s & %s ',...
                        '& %s & %s & %s \\\\ \n'],...
              firstLineExtra, ampTable(i), ...
              kbrVaf{1,1},    kbrVaf{1,2},    kbrVaf{1,3},...
              opus31Vaf{1,1}, opus31Vaf{1,2}, opus31Vaf{1,3}, ...
              hillVaf{1,1},   hillVaf{1,2},   hillVaf{1,3} );
  firstLineExtra = '';     
end
          

fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\end{center}\n');
fprintf(fid,'\\end{table}\n');

fclose(fid);
success=1;

% %%
% % Transform the table of fitted variables into a tex file
% %
% %     KBR 1994    Model      Hill
% % K/F 15 35 90    15 35 90   15 35 90
% % 0.4 
% % 0.8
% % 1.6
% %
% %     KBR 1994    Model      Hill
% % D/F 15 35 90    15 35 90   15 35 90
% % 0.4 
% % 0.8
% % 1.6
% %                 Model      Hill
% % VAF 15 35 90    15 35 90   15 35 90
% % 0.4 
% % 0.8 
% % 1.6
% %%
% 
% tendonTag = '_ElasticTendon';
% if(flag_useElasticTendon==0)
%   tendonTag = '_RigidTendon';
% end
% 
% fid = fopen([outputFolder,'tableStiffnessDampingVaf',tendonTag,'_',tableNameEnding,'.tex'],'w');
% 
% fprintf(fid,'\\begin{table}\n');
% fprintf(fid,'\\begin{tabular}{r | r r r || r r r || r r r}\n');
% %            A    15   35   90   15   35   90   15   35   90  
% 
% fprintf(fid,['%s & \\multicol{3}{c}{Kirsch et al.} & ',...
%             '\\multicol{3}{c}{Model} & \\multicol{3}{c}{Hill} \\\\ \n'],'');
% fprintf(fid,'%s & %d & %d & %d & %d & %d & %d & %d & %d & %d \\\\ \n',...
%             '$\frac{K}{F}$', 15, 35, 90, 15, 35, 90, 15, 35, 90 );
% 
% firstLineExtra = '\hline ';          
% for i=1:1:length(ampTable)
%   
%   kbrK    = cell(1,3);
%   opus31K = cell(1,3);
%   hillK   = cell(1,3);
% 
%   for j=1:1:3
%     strVal = '';
%     if(kbr1994TableStiffness.pMean(i,j,1) == noDataCode)
%       strVal = '--';      
%     else
%       strVal = sprintf('%1.2f',kbr1994TableStiffness.pMean(i,j,1));
%     end
%     kbrK(1,j) = {strVal};
%     
%     strVal = '';
%     if(opus31Table.stiffness.pMean(i,j,1) == noDataCode)
%       strVal = '--';      
%     else
%       strVal = sprintf('%1.2f',opus31Table.stiffness.pMean(i,j,1));
%     end
%     opus31K(1,j) = {strVal};
%     
%     strVal = '';
%     if(hillTable.stiffness.pMean(i,j,1) == noDataCode)
%       strVal = '--';      
%     else
%       strVal = sprintf('%1.2f',hillTable.stiffness.pMean(i,j,1));
%     end
%     hillK(1,j) = {strVal};    
%   end
%   
%   
% 
% 
%   fprintf(fid,['%s %1.1f & %s & %s & %s ',...
%                         '& %s & %s & %s ',...
%                         '& %s & %s & %s \\\\ \n'],...
%               firstLineExtra, ampTable(i), ...
%               kbrK{1,1},    kbrK{1,2},    kbrK{1,3},...
%               opus31K{1,1}, opus31K{1,2}, opus31K{1,3}, ...
%               hillK{1,1},   hillK{1,2},   hillK{1,3} );
%   firstLineExtra = '';    
%   
% end
% 
% fprintf(fid,'\\end{tabular}\n');
% fprintf(fid,'\\end{table}\n');
% 
% fprintf(fid,'\\begin{table}\n');
% fprintf(fid,'\\begin{tabular}{r | r r r || r r r || r r r}\n');
% %            A    15   35   90   15   35   90   15   35   90  
% 
% fprintf(fid,['%s & \\multicol{3}{c}{Kirsch et al.} & ',...
%             '\\multicol{3}{c}{Model} & \\multicol{3}{c}{Hill} \\\\ \n'],'');
% fprintf(fid,'%s & %d & %d & %d & %d & %d & %d & %d & %d & %d \\\\ \n',...
%             '$\frac{\beta}{F}$', 15, 35, 90, 15, 35, 90, 15, 35, 90 );
% 
% 
% firstLineExtra = '\hline ';          
% for i=1:1:length(ampTable)
%   
%   kbrD    = cell(1,3);
%   opus31D = cell(1,3);
%   hillD   = cell(1,3);
% 
%   for j=1:1:3
%     strVal = '';
%     if(kbr1994TableDamping.pMean(i,j,1) == noDataCode)
%       strVal = '--';      
%     else
%       strVal = sprintf('%1.4f',kbr1994TableDamping.pMean(i,j,1));
%     end
%     kbrD(1,j) = {strVal};
%     
%     strVal = '';
%     if(opus31Table.damping.pMean(i,j,1) == noDataCode)
%       strVal = '--';      
%     else
%       strVal = sprintf('%1.4f',opus31Table.damping.pMean(i,j,1));
%     end
%     opus31D(1,j) = {strVal};
%     
%     strVal = '';
%     if(hillTable.damping.pMean(i,j,1) == noDataCode)
%       strVal = '--';      
%     else
%       strVal = sprintf('%1.4f',hillTable.damping.pMean(i,j,1));
%     end
%     hillD(1,j) = {strVal};    
%   end
% 
% 
%   fprintf(fid,['%s %1.1f & %s & %s & %s ',...
%                         '& %s & %s & %s ',...
%                         '& %s & %s & %s \\\\ \n'],...
%               firstLineExtra, ampTable(i), ...
%               kbrD{1,1},    kbrD{1,2},    kbrD{1,3},...
%               opus31D{1,1}, opus31D{1,2}, opus31D{1,3}, ...
%               hillD{1,1},   hillD{1,2},   hillD{1,3} );
%   firstLineExtra = '';    
%   
% end
% 
% fprintf(fid,'\\end{tabular}\n');
% fprintf(fid,'\\end{table}\n');
% 
% 
% fprintf(fid,'\\begin{table}\n');
% fprintf(fid,'\\begin{tabular}{r | r r r || r r r || r r r}\n');
% %            A    15   35   90   15   35   90   15   35   90  
% 
% fprintf(fid,['%s & \\multicol{3}{c}{Kirsch et al.} & ',...
%             '\\multicol{3}{c}{Model} & \\multicol{3}{c}{Hill} \\\\ \n'],'');
% fprintf(fid,'%s & %d & %d & %d & %d & %d & %d & %d & %d & %d \\\\ \n',...
%             'avg. VAF', 15, 35, 90, 15, 35, 90, 15, 35, 90 );
% firstLineExtra = '\hline ';          
% for i=1:1:length(ampTable)
% 
% 
%   kbrVaf    = cell(1,3);
%   opus31Vaf = cell(1,3);
%   hillVaf   = cell(1,3);
% 
%   for j=1:1:3
%     
%     idx = getIndexIntoVectors(ampTable(1,i),freqTable(1,j),amp,freq);   
%     
% %     strVal = '';
% %     if(kbr1994TableDamping.pMean(i,j,1) == noDataCode)
% %       strVal = '--';      
% %     else
% %       strVal = sprintf('%1.3f',kbr1994TableDamping.pMean(i,j,1));
% %     end
%     kbrVaf(1,j) = {''};
%     
%     
%     
%     strVal = '';
%     if(opus31Table.vafData(idx).mean == noDataCode)
%       strVal = '--';      
%     else
%       strVal = sprintf('%1.1f',opus31Table.vafData(idx).mean*100);
%     end
%     opus31Vaf(1,j) = {strVal};
%     
%     strVal = '';
%     if(hillTable.vafData(idx).mean == noDataCode)
%       strVal = '--';      
%     else
%       strVal = sprintf('%1.1f',hillTable.vafData(idx).mean*100);
%     end
%     hillVaf(1,j) = {strVal};    
%   end
% 
% 
%   fprintf(fid,['%s %1.1f & %s & %s & %s ',...
%                         '& %s & %s & %s ',...
%                         '& %s & %s & %s \\\\ \n'],...
%               firstLineExtra, ampTable(i), ...
%               kbrVaf{1,1},    kbrVaf{1,2},    kbrVaf{1,3},...
%               opus31Vaf{1,1}, opus31Vaf{1,2}, opus31Vaf{1,3}, ...
%               hillVaf{1,1},   hillVaf{1,2},   hillVaf{1,3} );
%   firstLineExtra = '';     
% end
%           
% 
% fprintf(fid,'\\end{tabular}\n');
% fprintf(fid,'\\end{table}\n');
% 
% fclose(fid);
% success=1;


