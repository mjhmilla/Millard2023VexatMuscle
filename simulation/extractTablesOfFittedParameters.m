function fitTable = extractTablesOfFittedParameters(ampTable,freqTable,...
                                                    ampData,freqData,...
                                                    fittedData,...
                                                    noDataCode)

fitDim = 0;

for z=1:1:length(fittedData)
  

  if(isempty(fittedData(z).fo) == 0)
      if(length(fittedData(z).fo) > 1)
        if(length(coeffvalues(fittedData(z).fo)) > fitDim)
    
          if(fitDim == 0)
            fitDim = length(coeffvalues(fittedData(z).fo));
          else
            assert(0,'Error: all fields fo must use the same underlying parametric model');
          end
    
        end   
      else
        fitDim=1;
      end
  end
end

fitTable = struct(...
            'pMean'   ,zeros(length(ampTable),length(freqTable),fitDim),...
            'p95CIMin',zeros(length(ampTable),length(freqTable),fitDim),...
            'p95CIMax',zeros(length(ampTable),length(freqTable),fitDim),...
            'pNames', {''},...
            'rmse'    ,zeros(  length(ampTable),length(freqTable)));

pNames = {};          
flag_firstData = 0;          

for i = 1:1:length(ampTable)
  for j=1:1:length(freqTable)
    idx = getIndexIntoVectors(ampTable(1,i),freqTable(1,j),ampData,freqData);
    

    if(isempty(fittedData(idx).fo) == 0) 

        if(isa(fittedData(idx).fo,'double'))

            fitTable.pMean(i,j,1)    = fittedData(idx).fo;
            fitTable.p95CIMin(i,j,1) = fittedData(idx).fo;
            fitTable.p95CIMax(i,j,1) = fittedData(idx).fo;        
           fitTable.rmse(i,j) = 0;
        else
          if(flag_firstData == 0)
            pNames = coeffnames(fittedData(idx).fo);
            fitTable.pNames = pNames;
            flag_firstData = 1;
          end
          tmpNames = coeffnames(fittedData(idx).fo);
          assert(length(tmpNames) == length(pNames));
          for z=1:1:length(tmpNames)
            assert(strcmp(tmpNames{z},pNames{z})==1);
          end
          
          coeff   = coeffvalues(fittedData(idx).fo);
          coeffCI = confint(fittedData(idx).fo);
    
          for k=1:1:fitDim
            fitTable.pMean(i,j,k) = coeff(1,k);
            fitTable.p95CIMin(i,j,k) = coeffCI(1,k);
            fitTable.p95CIMax(i,j,k) = coeffCI(2,k);        
          end
          fitTable.rmse(i,j) = fittedData(idx).g.rmse;
        end
    else
      for k=1:1:fitDim
        fitTable.pMean(i,j,k)     = noDataCode;
        fitTable.p95CIMin(i,j,k)  = noDataCode;
        fitTable.p95CIMax(i,j,k)  = noDataCode;          
      end
      fitTable.rmse(i,j)      = noDataCode;            
    end
    
  end
end
          
          
                                                   