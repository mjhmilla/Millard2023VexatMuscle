function success = writeTRSS2017ErrorTable(fullFilePath, fitInfo, ...
                      simConfig,fittingConfig, projectFolders)





flag_fittingMethodFound=0;
if(fittingConfig.fitQToF==1)
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
                 '<QParamAll>' ,'<QNrmseAll1>','<QNrmseAll2>','<QNrmseAll3>'};
  
      tagValues = [fitInfo.fl.arg, fitInfo.fl.yErr(1), fitInfo.fl.yErr(2), fitInfo.fl.yErr(3),...
                   fitInfo.fv.arg, fitInfo.fv.yErr(1), fitInfo.fv.yErr(2), fitInfo.fv.yErr(3),...
                   fitInfo.QToF.arg(1),...
                   fitInfo.QToF.nrmse(1),fitInfo.QToF.nrmse(2),fitInfo.QToF.nrmse(3)];
  end
  flag_fittingMethodFound=1;
  fidTemplate = fopen(fullfile(projectFolders.output_tables_TRSS2017,...
                            'TRSS2017TableTemplateA.tex'));
end

if(fittingConfig.fitFpeQToF==1)
  if(simConfig.fitToIndividualTrials==1)
      tagList = {'<flParam>',  '<flError1>',  '<flError2>',  '<flError3>',...
                 '<fvParam>',  '<fvError1>',  '<fvError2>',  '<fvError3>',...
                 '<QParam1>'   , '<QParam2>','<QParam3>',...
                 '<kToeParam1>','<normLengthToeParam1>','<curvinessParam1>',...
                 '<nrmse1>','<nrmse2>','<nrmse3>'};
  
      tagValues = [fitInfo.fl.arg, fitInfo.fl.yErr(1), fitInfo.fl.yErr(2), fitInfo.fl.yErr(3),...
                   fitInfo.fv.arg, fitInfo.fv.yErr(1), fitInfo.fv.yErr(2), fitInfo.fv.yErr(3),...
                   fitInfo.FpeQToF.arg,...
                   fitInfo.FpeQToF.nrmse(1),fitInfo.FpeQToF.nrmse(2),fitInfo.FpeQToF.nrmse(3)];
  else

      tagList = {'<flParam>',  '<flError1>',  '<flError2>',  '<flError3>',...
                 '<fvParam>',  '<fvError1>',  '<fvError2>',  '<fvError3>',...
                 '<QParamAll>'   ,'<kToeParamAll>','<normLengthToeParamAll>','<curvinessParamAll>',...            
                 '<nrmseAll1>','<nmseAll2>','<nrmseAll3>'};
  
      tagValues = [fitInfo.fl.arg, fitInfo.fl.yErr(1), fitInfo.fl.yErr(2), fitInfo.fl.yErr(3),...
                   fitInfo.fv.arg, fitInfo.fv.yErr(1), fitInfo.fv.yErr(2), fitInfo.fv.yErr(3),...
                   fitInfo.FpeQToF.arg,...
                   fitInfo.FpeQToF.nrmse(1),fitInfo.FpeQToF.nrmse(2),fitInfo.FpeQToF.nrmse(3)];      
  end
  flag_fittingMethodFound=1;
  fidTemplate = fopen(fullfile(projectFolders.output_tables_TRSS2017,...
                            'TRSS2017TableTemplateB.tex'));
end

if(flag_fittingMethodFound==1)
  fidTable = fopen(fullFilePath,'w');


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
end
success = 1;

