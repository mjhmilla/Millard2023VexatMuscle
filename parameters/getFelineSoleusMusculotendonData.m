function [activeForceLengthData,...
          passiveForceLengthData] = getFelineSoleusMusculotendonData(...
                                      musculotendonPropertiesExp,...                              
                                      musculotendonProperties,...
                                      normPlateauShift,...
                                      useElasticTendon,...
                                      flag_useOctave)

fitData =csvread('experiments/HerzogLeonard2002/data/dataHerzogLeonard2002Figure6and7A.csv',1,0);

%Note: 'F43', 'L43', etc are the column headings in the file.
idxFitF43=2;   % Length 0mm, Activation 0-1-0 
idxFitL43=3;
idxFitF44=4;   % Length 9mm, Activation 0-1-0
idxFitL44=5;
idxFitF45=6;   % Length 6-9mm, Activation 0-1-0
idxFitL45=7;
idxFitF47=8;   % Length 3-9mm, Activation 0-1-0
idxFitL47=9;
idxFitF49=10;  % Length 0-9mm, Activation 0 (Passive)
idxFitL49=11;
idxFitF50=12;  % Length 0-9mm, Activation 0-1-0
idxFitL50=13;

idxFitMax = 1500; %All columns except 50 series stop at 1500

%
% Extract experimental data on the passive force length curves 
%
idx49LStretch       = find(fitData(:,idxFitL49) > 1.0);
idx49LStart         = idx49LStretch(1);     
[maxF49 idxMaxF49]  = max(fitData(:,idxFitF49));


%Use the parameters that we have directly from Herzog & Leonard,
%fill the remaining needed parameters from other sources

params = {'tendonStrainAtOneNormForce',...
          'optimalFiberLength',...
          'tendonSlackLength',...
          'pennationAngle',...
          'fiso'};
vals = zeros(length(params),1);
for i=1:1:length(params)
  %if(isfield(musculotendonProperties,params{i}))
  assert(isfield(musculotendonProperties,params{i})==1);
  vals(i,1) = musculotendonProperties.(params{i});
    %else
    %  vals(i,1) = musculotendonProperties.(params{i});
  %end
end

eIso      = vals(1,1);
lceOpt    = vals(2,1);
ltSlk     = vals(3,1);
alphaOpt  = vals(4,1);
fiso      = vals(5,1);


normTendonLength = 1+eIso;
if(useElasticTendon==0)
  eIso = 0;
  normTendonLength = 1;
end

tendonForceLengthCurve = [];

if(useElasticTendon==1)
  kIso            = 1.375/eIso;
  fToe            = 2.0/3.0;
  curviness       = 0.5;
  computeIntegral = 0;
  minimumSlope = sqrt(eps);
  flag_enableNumericallyNonZeroGradients = 0;
  smallNumericallyNonZeroNumber = sqrt(sqrt(eps));
  
  tendonForceLengthCurve = ...
    createTendonForceLengthCurve2019( ...
        eIso, kIso, fToe, curviness, computeIntegral, ...
        flag_enableNumericallyNonZeroGradients,...
        smallNumericallyNonZeroNumber,...
        'temporary',flag_useOctave);

end        


%In the experimental data the total MTU length is recorded. Thus the
%lmtOptZero = lmtOptFiso, since a tendon model is not used to account
%for the difference in fiber length due to tendon stretching.

lceOptLeft = lceOpt*(1-normPlateauShift);
fibKin = calcFixedWidthPennatedFiberKinematicsAlongTendon(...
                                        lceOptLeft,...
                                        0,...
                                        lceOpt,...
                                        alphaOpt);

alphaOptLeft = fibKin.pennationAngle;                                      
lmtOptFiso  = lceOptLeft*cos(alphaOptLeft) + ltSlk*(1+eIso);
lmtOptZero  = lmtOptFiso;%lceOpt*cos(alphaOpt) + ltSlk;


idx   = 1;
idxDelta = floor((idxMaxF49-idx49LStart)/10); %So that we get ~10 pts

passiveForceLengthData = zeros(floor((idxMaxF49-idx49LStart)/idxDelta),2);

for(i=idx49LStart:idxDelta:idxMaxF49)
  %Calculate the length of the tendon
  ftN = fitData(i,idxFitF49)/fiso;
  
  if(useElasticTendon==1)
    ltN = calcBezierFcnXGivenY(ftN,tendonForceLengthCurve);
  else
    ltN=1;
  end
  lt  = ltN*ltSlk;
 
  %Calculate the length of the fiber
  lmAT    = (lmtOptZero + fitData(i,idxFitL49)/1000) - lt;
  fibKin  = calcFixedWidthPennatedFiberKinematics(lmAT,0,lceOpt,alphaOpt);
  
  lceN = fibKin.fiberLength/lceOpt;
  fceN = ( fitData(i,idxFitF49)/cos(fibKin.pennationAngle) )/fiso;
  passiveForceLengthData(idx,1) = lceN;
  passiveForceLengthData(idx,2) = fceN;    
  idx = idx+1;
end    
      
%
% Extract experimental data on the active force length curves:
%    get the point just prior to the muscle being stretched 
%
colFalF =[idxFitF43; idxFitF44; idxFitF45; idxFitF47];
colFalL =[idxFitL43; idxFitL44; idxFitL45; idxFitL47];
rowFal  =zeros(length(colFalF),1);
[tmp rowFal(1)] = max(fitData(:,idxFitF43));
[tmp rowFal(2)] = max(fitData(:,idxFitF44));
rowFal(3) = 422; % Point in time just prior to the lengthening ramp
rowFal(4) = 326; % " ... "

falLpDelta    = zeros(length(colFalL),1);
falFt         = zeros(length(colFalL),1);
falFpe        = zeros(length(colFalL),1);
falFpErr      = zeros(length(colFalL),1);

for(k=1:1:length(colFalF))
  
    falLpDelta(k) = fitData(rowFal(k), colFalL(k));
    falFt(k)      = fitData(rowFal(k), colFalF(k));
    
    %Scan through the passive-force-length trial and get the passive
    %force at the current path length
    lerr = Inf;
    for(z=1:1:idxMaxF49)
        err = abs(falLpDelta(k)-fitData(z,idxFitL49));
        if(err < lerr)
           lerr = err;
           falFpe(k) = fitData(z,idxFitF49);
           falFpErr(k) = err;
        end        
    end
end

activeForceLengthData  = zeros(4, 2);

for(k=1:1:length(colFalF))  
  %Calculate the length of the tendon
  ftN = falFt(k,1)/fiso;
  if(useElasticTendon == 1)
    ltN = calcBezierFcnXGivenY(ftN,tendonForceLengthCurve);
  else
    ltN = 1;
  end
  lt  = ltN*ltSlk;  
  
  %The changes in tendon length are being ignored here
  lmAT    = (lmtOptFiso + falLpDelta(k,1)./1000) - lt;
  fibKin  = calcFixedWidthPennatedFiberKinematics(lmAT,0,lceOpt,alphaOpt);  
  lceN = fibKin.fiberLength/lceOpt;  
  fceN = ( (falFt(k,1)-falFpe(k,1)) / cos(fibKin.pennationAngle) ) / fiso;
  
  activeForceLengthData(k,1) = lceN;    
  activeForceLengthData(k,2) = fceN;
end

