
if(sarcomereProperties.fitCrossBridgeStiffnessDampingToKirch199490Hz==1)
  strFittingBandwidth = '90Hz';
else
  strFittingBandwidth = '15Hz';
end

nominalNormalizedFiberLengthStr = '';
if(exist('nominalNormalizedFiberLength','var'))
  nominalNormalizedFiberLengthStr = sprintf('%1.2f',nominalNormalizedFiberLength);
  k = strfind(nominalNormalizedFiberLengthStr,'.');
  if(isempty(k)==0)
    nominalNormalizedFiberLengthStr = [nominalNormalizedFiberLengthStr(1,1:(k-1)),'p',nominalNormalizedFiberLengthStr(1,(k+1):1:end)];
  end
end



titinActiveDampingStr='';
if(exist('normActiveTitinToActinDamping','var'))
  titinActiveDampingStr = sprintf('%1.2f',normActiveTitinToActinDamping);
  k = strfind(titinActiveDampingStr,'.');
  if(isempty(k)==0)
    titinActiveDampingStr = [titinActiveDampingStr(1,1:(k-1)),'p',titinActiveDampingStr(1,(k+1):1:end)];
  end
end

titinPassiveDampingStr = '';
if(exist('normPassiveTitinToActinDamping','var'))
  titinPassiveDampingStr = sprintf('%1.2f',normPassiveTitinToActinDamping);
  k = strfind(titinPassiveDampingStr,'.');
  if(isempty(k)==0)
    titinPassiveDampingStr = [titinPassiveDampingStr(1,1:(k-1)),'p',titinPassiveDampingStr(1,(k+1):1:end)];
  end
end

kScaleStr = sprintf('%1.2f',sarcomerePropertiesOpus31.normCrossBridgeStiffness);
k = strfind(kScaleStr,'.');
if(isempty(k)==0)
  kScaleStr = [kScaleStr(1,1:(k-1)),'p',kScaleStr(1,(k+1):1:end)];
end

dScaleStr = sprintf('%1.2f',sarcomerePropertiesOpus31.normCrossBridgeDamping);
d = strfind(dScaleStr,'.');
if(isempty(d)==0)
  dScaleStr = [dScaleStr(1,1:(d-1)),'p',dScaleStr(1,(d+1):1:end)];
end

%Rigid Tendon
kScaleStr_RT = '';
dScaleStr_RT = '';
if(isempty(sarcomerePropertiesOpus31_RT)==0)
  kScaleStr_RT = sprintf('%1.2f',sarcomerePropertiesOpus31_RT.normCrossBridgeStiffness);
  k = strfind(kScaleStr_RT,'.');
  if(isempty(k)==0)
    kScaleStr_RT = [kScaleStr_RT(1,1:(k-1)),'p',kScaleStr_RT(1,(k+1):1:end)];
  end

  dScaleStr_RT = sprintf('%1.2f',sarcomerePropertiesOpus31_RT.normCrossBridgeDamping);
  d = strfind(dScaleStr_RT,'.');
  if(isempty(d)==0)
    dScaleStr_RT = [dScaleStr_RT(1,1:(d-1)),'p',dScaleStr_RT(1,(d+1):1:end)];
  end
end

%Elastic tendon
kScaleStr_ET = '';
dScaleStr_ET = '';
if(isempty(sarcomerePropertiesOpus31_ET)==0)
  kScaleStr_ET = sprintf('%1.2f',sarcomerePropertiesOpus31_ET.normCrossBridgeStiffness);
  k = strfind(kScaleStr_ET,'.');
  if(isempty(k)==0)
    kScaleStr_ET = [kScaleStr_ET(1,1:(k-1)),'p',kScaleStr_ET(1,(k+1):1:end)];
  end

  dScaleStr_ET = sprintf('%1.2f',sarcomerePropertiesOpus31_ET.normCrossBridgeDamping);
  d = strfind(dScaleStr_ET,'.');
  if(isempty(d)==0)
    dScaleStr_ET = [dScaleStr_ET(1,1:(d-1)),'p',dScaleStr_ET(1,(d+1):1:end)];
  end
end

tScaleStr = '';
if(exist('scaleSlidingTimeConstant','var'))
  tScaleStr = sprintf('%1.2f',scaleSlidingTimeConstant);
  t = strfind(tScaleStr,'.');
  if(isempty(t)==0)
    tScaleStr = [tScaleStr(1,1:(t-1)),'p',tScaleStr(1,(t+1):1:end)];
  end
end

kTLinearStr = '';
if(exist('normTendonDampingLinear','var'))
  kTLinearStr = sprintf('%1.2f',normTendonDampingLinear);
  t = strfind(kTLinearStr,'.');
  if(isempty(t)==0)
    kTLinearStr = [kTLinearStr(1,1:(t-1)),'p',kTLinearStr(1,(t+1):1:end)];
  end
end


kTConstStr = '';
if(exist('normTendonDampingConstant','var'))
  kTConstStr = sprintf('%1.2f',normTendonDampingConstant);
  t = strfind(kTConstStr,'.');
  if(isempty(t)==0)
    kTConstStr = [kTConstStr(1,1:(t-1)),'p',kTConstStr(1,(t+1):1:end)];
  end
end

outputFileEndingPropOpus31_ET = '';

if(isempty(sarcomerePropertiesOpus31_ET)==0)
  outputFileEndingPropOpus31_ET = sprintf('_K%sD%sTau%s_KTC%s_KTL%s',...
    kScaleStr_ET, dScaleStr_ET, tScaleStr, kTConstStr, kTLinearStr);
end

outputFileEndingPropOpus31_RT = '';

if(isempty(sarcomerePropertiesOpus31_RT)==0)
  outputFileEndingPropOpus31_RT = sprintf('_K%sD%sTau%s',...
    kScaleStr_RT, dScaleStr_RT, tScaleStr);%'_Std_5A1L3B3D_Descending';
end

if(flag_useElasticTendon==1)
  outputFileEndingPropOpus31 = outputFileEndingPropOpus31_ET;
else
  outputFileEndingPropOpus31 = outputFileEndingPropOpus31_RT;
end


outputFileEndingPropHill = sprintf('_D%i',flag_useFiberDamping);

outputFileEndingPropHill_RT = ...
  sprintf('_D%i', 0);
  
outputFileEndingPropHill_ET = ...
  sprintf('_D%i', 1);






