%%
% SPDX-FileCopyrightText: 2023 Matthew Millard <millard.matthew@gmail.com>
%
% SPDX-License-Identifier: MIT
%
% If you use this code in your work please cite the pre-print of this paper
% or the most recent peer-reviewed version of this paper:
%
%    Matthew Millard, David W. Franklin, Walter Herzog. 
%    A three filament mechanistic model of musculotendon force and impedance. 
%    bioRxiv 2023.03.27.534347; doi: https://doi.org/10.1101/2023.03.27.534347 
%
%%

function success = tabulateStiffnessDampingVariationPUB(opus31Table,...
                                          hillTable,...
                                          kbr1994Table,...
                                          inputFunctions,...
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
noDataCode = nan;



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
kbrK    = cell(3,3);
opus31K = cell(3,3);
hillK   = cell(3,3);
  
for idxTrial=1:1:length(kbr1994Table.stiffness.data)
  

  [val, idxAmp] = min(abs(ampTable-inputFunctions.amplitudeMM(1,idxTrial)));
  [val, idxFreq] = min(abs(freqTable-inputFunctions.bandwidthHz(1,idxTrial)));
    
  val       = mean(kbr1994Table.stiffness.data(idxTrial).yN);
  strVal    = sprintf('%1.2f',val);
  kbrK(idxAmp,idxFreq) = {strVal};

  val           = mean(opus31Table.stiffness.data(idxTrial).yN);
  strVal        = sprintf('%1.2f',val);
  opus31K(idxAmp,idxFreq)  = {strVal};

  val           = mean(hillTable.stiffness.data(idxTrial).yN);
  strVal        = sprintf('%1.2f',val);
  hillK(idxAmp,idxFreq)  = {strVal};
  
  


end


for i=1:1:3
      fprintf(fid,['%s %1.1f mm & %s & %s & %s ',...
                            '& %s & %s & %s ',...
                            '& %s & %s & %s \\\\ \n'],...
                  firstLineExtra, ampTable(1,i), ...
                     kbrK{i,1},    kbrK{i,2},    kbrK{i,3},...
                  opus31K{i,1}, opus31K{i,2}, opus31K{i,3}, ...
                    hillK{i,1},   hillK{i,2},   hillK{i,3} );
      firstLineExtra = '';    
end 



fprintf(fid,'\\multicolumn{10}{c}{} \\\\ \n');


fprintf(fid,'%s & %dHz & %dHz & %dHz & %dHz & %dHz & %dHz & %dHz & %dHz & %dHz \\\\ \n',...
            'B. \hfill Norm. damping $\frac{\beta}{F}$', 15, 35, 90, 15, 35, 90, 15, 35, 90 );


firstLineExtra = '\hline ';          

kbrD    = cell(3,3);
opus31D = cell(3,3);
hillD   = cell(3,3);
  
for idxTrial=1:1:length(kbr1994Table.damping.data)
  
  [val, idxAmp] = min(abs(ampTable-inputFunctions.amplitudeMM(1,idxTrial)));
  [val, idxFreq] = min(abs(freqTable-inputFunctions.bandwidthHz(1,idxTrial)));



  val       = mean(kbr1994Table.damping.data(idxTrial).yN);
  strVal    = sprintf('%1.4f',val);
  kbrD(idxAmp,idxFreq) = {strVal};

  val           = mean(opus31Table.damping.data(idxTrial).yN);
  strVal        = sprintf('%1.4f',val);
  opus31D(idxAmp,idxFreq)  = {strVal};

  val           = mean(hillTable.damping.data(idxTrial).yN);
  strVal        = sprintf('%1.4f',val);
  hillD(idxAmp,idxFreq)  = {strVal};
  
end

for i=1:1:3
      fprintf(fid,['%s %1.1f mm & %s & %s & %s ',...
                            '& %s & %s & %s ',...
                            '& %s & %s & %s \\\\ \n'],...
                  firstLineExtra, ampTable(1,i), ...
                     kbrD{i,1},    kbrD{i,2},    kbrD{i,3},...
                  opus31D{i,1}, opus31D{i,2}, opus31D{i,3}, ...
                    hillD{i,1},   hillD{i,2},   hillD{i,3} );
      firstLineExtra = '';    
end 

fprintf(fid,'\\multicolumn{10}{c}{} \\\\ \n');

fprintf(fid,'%s & %dHz & %dHz & %dHz & %dHz & %dHz & %dHz & %dHz & %dHz & %dHz \\\\ \n',...
            'C. VAF (\%) \hfill', 15, 35, 90, 15, 35, 90, 15, 35, 90 );
firstLineExtra = '\hline ';  


kbrVAF    = cell(3,3);
opus31VAF = cell(3,3);
hillVAF   = cell(3,3);
  
for idxTrial=1:1:length(kbr1994Table.vaf.data)
  


  [val, idxAmp] = min(abs(ampTable-inputFunctions.amplitudeMM(1,idxTrial)));
  [val, idxFreq] = min(abs(freqTable-inputFunctions.bandwidthHz(1,idxTrial)));



  val       = mean(kbr1994Table.vaf.data(idxTrial).y)*100;
  strVal    = sprintf('%1.0f',val);
  kbrVAF(idxAmp,idxFreq) = {strVal};

  val           = mean(opus31Table.vaf.data(idxTrial).y)*100;
  strVal        = sprintf('%1.0f',val);
  opus31VAF(idxAmp,idxFreq)  = {strVal};

  val           = mean(hillTable.vaf.data(idxTrial).y)*100;
  strVal        = sprintf('%1.0f',val);
  hillVAF(idxAmp,idxFreq)  = {strVal};
  
  


end 

for i=1:1:3
      fprintf(fid,['%s %1.1f mm & %s & %s & %s ',...
                            '& %s & %s & %s ',...
                            '& %s & %s & %s \\\\ \n'],...
                  firstLineExtra, ampTable(1,i), ...
                     kbrVAF{i,1},    kbrVAF{i,2},    kbrVAF{i,3},...
                  opus31VAF{i,1}, opus31VAF{i,2}, opus31VAF{i,3}, ...
                    hillVAF{i,1},   hillVAF{i,2},   hillVAF{i,3} );
      firstLineExtra = '';    
end 

fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\end{center}\n');
fprintf(fid,'\\end{table}\n');

fclose(fid);
success=1;



