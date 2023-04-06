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

function fitTable = extractTablesOfFittedParameters(fittedData,...
                                                    noDataCode)

fitDim = 0;

for z=1:1:length(fittedData)
  
  if(isempty(fittedData(z).fo) == 0)
    if(length(coeffvalues(fittedData(z).fo)) > fitDim)

      if(fitDim == 0)
        fitDim = length(coeffvalues(fittedData(z).fo));
      else
        assert(0,'Error: all fields fo must use the same underlying parametric model');
      end

    end   
  end
end

fitTable = struct(...
            'pMean'   ,zeros(length(fittedData),fitDim),...
            'p95CIMin',zeros(length(fittedData),fitDim),...
            'p95CIMax',zeros(length(fittedData),fitDim),...
            'pNames', {''},...
            'rmse'   , zeros(  length(fittedData),1),...
            'data'   , []);

data(length(fittedData)) = struct('x',[],'y',[],'yN',[]);

pNames = {};          
flag_firstData = 0;          

for idx = 1:1:length(fittedData)
    
    data(idx).x=fittedData(idx).x;
    data(idx).y=fittedData(idx).y;
    data(idx).yN=fittedData(idx).yN;

    if(isempty(fittedData(idx).fo) == 0)    
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
        fitTable.pMean(idx,k) = coeff(1,k);
        fitTable.p95CIMin(idx,k) = coeffCI(k,1);
        fitTable.p95CIMax(idx,k) = coeffCI(k,2);        
      end
      fitTable.rmse(idx,1) = fittedData(idx).g.rmse;
    else
      for k=1:1:fitDim
        fitTable.pMean(idx,k)     = noDataCode;
        fitTable.p95CIMin(idx,k)  = noDataCode;
        fitTable.p95CIMax(idx,k)  = noDataCode;          
      end
      fitTable.rmse(idx,1)      = noDataCode;            
    end
    
end
        
fitTable.data = data;
          
                                                   