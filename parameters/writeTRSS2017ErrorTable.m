function success = writeTRSS2017ErrorTable(fullFilePath, fitInfo, simConfig, projectFolders)


fidTable = fopen(fullFilePath,'w');

fidTemplate = fopen(fullfile(projectFolders.output_tables_TRSS2017,...
                            'TRSS2017TableTemplate.tex'));

if(simConfig.fitToIndividualTrials==1)
    tagList = {'<flParam>',  '<flError1>',  '<flError2>',  '<flError3>',...
               '<fvParam>',  '<fvError1>',  '<fvError2>',  '<fvError3>',...
               '<QParam1>'   ,'<QParam2>','<QParam3>',...
               '<QNrmse1>','<QNrmse2>','<QNrmse3>'};

    tagValues = [fitInfo.fl.arg, fitInfo.fl.yErr(1), fitInfo.fl.yErr(2), fitInfo.fl.yErr(3),...
                 fitInfo.fv.arg, fitInfo.fv.yErr(1), fitInfo.fv.yErr(2), fitInfo.fv.yErr(3),...
                 fitInfo.QToF.arg(1),fitInfo.QToF.arg(2),fitInfo.QToF.arg(3),...
                 fitInfo.QToF.nrmse(1),fitInfo.QToF.nrmse(2),fitInfo.QToF.nrmse(3)];
else
    tagList = {'<flParam>',  '<flError1>',  '<flError2>',  '<flError3>',...
               '<fvParam>',  '<fvError1>',  '<fvError2>',  '<fvError3>',...
               '<QParamAll>' ,'<NrmseQall1>','<NrmseQall2>','<NrmseQall3>'};

    tagValues = [fitInfo.fl.arg, fitInfo.fl.yErr(1), fitInfo.fl.yErr(2), fitInfo.fl.yErr(3),...
                 fitInfo.fv.arg, fitInfo.fv.yErr(1), fitInfo.fv.yErr(2), fitInfo.fv.yErr(3),...
                 fitInfo.QToF.arg(1),...
                 fitInfo.QToF.nrmse(1),fitInfo.QToF.nrmse(2),fitInfo.QToF.nrmse(3)];
end


assert(length(tagValues)==length(tagList));

lineTxt = fgetl(fidTemplate);
while(lineTxt ~= -1)
    
    %Replace any tags with the corresponding number
    for i=1:1:length(tagList)
        idxFound = strfind(lineTxt,tagList{i});
        assert(isempty(idxFound) || length(idxFound)==1);
        if(length(idxFound)==1)
          lineTxt = strrep(lineTxt,tagList{i},sprintf('%1.3f',tagValues(i)));
        end
    end
    fprintf(fidTable,'%s\n',lineTxt);

    lineTxt = fgetl(fidTemplate);
end

fclose(fidTable);
fclose(fidTemplate);

success = 1;

% 
% fitInfoFields=fields(fitInfo);
% fprintf(fidTable,'%s\n','\begin{tabular}{l l l l}');
% fprintf(fidTable,'%s\n','Parameter & Trial 1 & Trial 2 & Trial 3\\');
% 
% for idxF=1:1:length(fitInfoFields)
%     if(~isnan(fitInfo.(fitInfoFields{idxF}).rmse))
%         if(simConfig.fitToIndividualTrials)
%             fprintf(fidTable,...
%                     '%s Param',...
%                     fitInfoFields{idxF});
%             if(size(fitInfo.(fitInfoFields{idxF}).arg,2)>1)
%                 for i=1:1:size(fitInfo.(fitInfoFields{idxF}).arg,2)
%                     fprintf(fidTable,...
%                         ' & %1.3f',...
%                         fitInfo.(fitInfoFields{idxF}).arg(1,i));
%                 end
%             else
%                 fprintf(fidTable,...
%                     ' & \\multicol{3}{c}{%1.3f}\\\\ \n',...
%                     fitInfo.(fitInfoFields{idxF}).arg(1,1));
%             end
%             fprintf(fidTable,'\n');
%             
%             fprintf(fidTable,...
%                     '%s RMSE',...
%                     fitInfoFields{idxF});
%                 
%             for i=1:1:size(fitInfo.(fitInfoFields{idxF}).yErr,2)
%                 fprintf(fidTable,...
%                     ' & %1.3f',...
%                     sqrt(mean(fitInfo.(fitInfoFields{idxF}).yErr(:,i).^2)) );
%             end
%             fprintf(fidTable,'\n');
%         else
%             fprintf(fidTable,...
%                     '%s Param',...
%                     fitInfoFields{idxF});                    
%             fprintf(fidTable,...
%                 ' & \\multicol{3}{c}{%1.3f}\\\\ \n',...
%                 fitInfo.(fitInfoFields{idxF}).arg(1,1));
%             
%             fprintf(fidTable,...
%                     '%s RMSE',...
%                     fitInfoFields{idxF});                    
%             for i=1:1:size(fitInfo.(fitInfoFields{idxF}).yErr,2)
%                 fprintf(fidTable,...
%                     ' & %1.3f',...
%                     sqrt(mean(fitInfo.(fitInfoFields{idxF}).yErr(:,i).^2)) );
%             end
%             fprintf(fidTable,'%s\n','\\');
% 
%         end
%     end
% end
% fprintf(fidTable,'%s\n','\end{tabular}');
%
%fclose(fidTable);