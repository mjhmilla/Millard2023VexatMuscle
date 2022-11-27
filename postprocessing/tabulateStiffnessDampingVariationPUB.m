function success = tabulateStiffnessDampingVariationPUB( ...
                                              opus31Table,...
                                              hillTable,...
                                              kbr1994Table,...
                                              flag_useElasticTendon,... 
                                              tableNameEnding,...
                                              outputFolder)

%                                               dataFolder,...
%                                               freqSeriesFiles,...
%                                               freqSeriesName,...
%                                               freqSeriesColor,...
%                                               inputFunctions,...                      
%                                               normFiberLength,...
%                                               nominalForce,...  
%                                               dataKBR1994Fig9A,...
%                                               dataKBR1994Fig9B,... 
%                                               dataKBR1994Fig10,...
%                                               dataKBR1994Fig12K,...
%                                               dataKBR1994Fig12D,...
success=0;                                            
               
freqTable = [15,35,90];
ampTable = [0.4,0.8,1.6];
noDataCode = -42;

freq = [ 15,  15,  15,  35,  35,  35,  90,  90,  90];
amp  = [0.4, 0.8, 1.6, 0.4, 0.8, 1.6, 0.4, 0.8, 1.6];

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
    if(kbr1994Table.stiffness.pMean(i,j,1) == noDataCode)
      strVal = '--';      
    else
      strVal = sprintf('%1.2f',kbr1994Table.stiffness.pMean(i,j,1));
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



fprintf(fid,'\\multicolumn{10}{c}{} \\\\ \n');


fprintf(fid,'%s & %dHz & %dHz & %dHz & %dHz & %dHz & %dHz & %dHz & %dHz & %dHz \\\\ \n',...
            'B. \hfill Norm. damping $\frac{\beta}{F}$', 15, 35, 90, 15, 35, 90, 15, 35, 90 );


firstLineExtra = '\hline ';          
for i=1:1:length(ampTable)
  
  kbrD    = cell(1,3);
  opus31D = cell(1,3);
  hillD   = cell(1,3);

  for j=1:1:3
    strVal = '';
    if(kbr1994Table.damping.pMean(i,j,1) == noDataCode)
      strVal = '--';      
    else
      strVal = sprintf('%1.4f',kbr1994Table.damping.pMean(i,j,1));
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
    strVal = sprintf('%1.1f',mean(opus31Table.vaf.data(idx).y*100));

    opus31Vaf(1,j) = {strVal};
    
    strVal = '';    
    strVal = sprintf('%1.1f',mean(hillTable.vaf.data(idx).y*100));

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



