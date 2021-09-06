function [forceLengthIgpCurve, forceLengthIgpInverseCurve,...
          forceLengthPevkCurve, forceLengthPevkInverseCurve] ...
          = createTitinCurves2021( fiberForceLengthCurve,...                                   
                                   forceLengthCurveSettings,...
                                   forceLengthECMHalfCurve,...
                                   sarcomereProperties,...
                                   muscleName,...
                                   flag_computeCurveIntegrals,...
                                   flag_useElasticIgD,...
                                   flag_useHumanIgDGeometry,...
                                   flag_useOctave)
                                 

%%
%Evaluate the length and stiffness of fpe when it reaches a normalized
%force of 1.
%%


lambdaECM = sarcomereProperties.extraCellularMatrixPassiveForceFraction;


%%
% Get the titin normalized segment lengths
%%
ZLineToT12NormLengthAtOptimalFiberLength = ...
  sarcomereProperties.ZLineToT12NormLengthAtOptimalFiberLength;
IGDTotalNormLengthAtOptimalFiberLength = ...
  sarcomereProperties.IGDTotalNormLengthAtOptimalFiberLength;

IGDFixedNormLengthAtOptimalFiberLength=0;

if(flag_useHumanIgDGeometry==1)
  IGDFixedNormLengthAtOptimalFiberLength=...
    sarcomereProperties.IGDFixedNormLengthAtOptimalFiberLengthHuman;
else
  IGDFixedNormLengthAtOptimalFiberLength=...
    sarcomereProperties.IGDFixedNormLengthAtOptimalFiberLength;  
end

lceZero = fiberForceLengthCurve.xEnd(1,1);

k = size(fiberForceLengthCurve.ypts,2);
x0 = 0;
if(   fiberForceLengthCurve.ypts(end,k) < 1 )
  x0 =   fiberForceLengthCurve.xEnd(1,2)...
       + diff(fiberForceLengthCurve.xEnd)*0.1;
else
  x0 = mean(fiberForceLengthCurve.xpts(:,k));
end



lceOne = ...
  calcBezierFcnXGivenY(1, fiberForceLengthCurve, x0);

lceZeroHalf  = 0.5*lceZero;
lceOneHalf   = 0.5*lceOne;

kpeOne = calcBezierYFcnXDerivative(lceOne,...
                      fiberForceLengthCurve,...
                      1);

kpeOneHalf = 2*kpeOne;
kecmOneHalf= calcBezierYFcnXDerivative(lceOneHalf,...
                      forceLengthECMHalfCurve,...
                      1);

kTitinOneHalf = kpeOneHalf - kecmOneHalf;

igpStretchRate     = sarcomereProperties.IGPNormStretchRate;

pevkStretchRate = 0;
if(flag_useElasticIgD==1)
  pevkStretchRate = sarcomereProperties.PEVKNormStretchRate ...
                      +sarcomereProperties.IGDFreeNormStretchRate;
  
else
  pevkStretchRate   = sarcomereProperties.PEVKNormStretchRate;                      
end

                                    
%%
% ka = kigp
% kb = kpevk-igd
% kc = ktitin
%
% 1/ka + 1/kb = 1/kc [1]
%
% We also have the stretch rates extracted from Trombitas 1998 (see
% parameters/felineSoleus/getMammalianSkeletalMuscleNormalizedSarcomereProperties.m
% where the variable 'normStretchRateIgP is defined (near line 218-220)
%
% sa = rate the igp region stretches as the sarcomere stretches
% sb = rate the pevk + igd (free) region stretches as the sarcomere 
%      stretches
%
% Or more concisely
%
%  sa = delta la /delta ls     [2]
%  sb = delta lb /delta ls     [3]
%
% where ls is the sarcomere length. Since these elements are in series
%
% ka = delta f / delta la      [4]
% kb = delta f / delta lb      [5]
%
% and thus
%
% ka/kb = delta lb / delta la  [6]
%
% which is equivalent to
%
% ka/kb = sb /sa               [7]
%
% And so
%
% ka = (sb/sa)*kb              [8]
%    = A*kb
% Coming back to [1]
%
% (1/A*kb) + (1/kb)       = 1/kc
% (  A+1   ) / ( A * kb) = 1/kc
% ( A * kb ) / ( A + 1 )   = kc
% 
% kb = kc * (A+1) / (A)
%%
A  =  pevkStretchRate/igpStretchRate;

kpevk = kTitinOneHalf*(A+1)/A;
kigp     = A*kpevk; 


%Geometrically scale the other stiffnesses of the force-length curve
kpevkLow = kpevk*(forceLengthCurveSettings.kLow ...
                       /forceLengthCurveSettings.kToe);
kigpLow     = kigp*(forceLengthCurveSettings.kLow ...
                    /forceLengthCurveSettings.kToe);

kpevkZero = kpevk*(forceLengthCurveSettings.kZero ...
                       /forceLengthCurveSettings.kToe);
kigpZero     = kigp*(forceLengthCurveSettings.kZero ...
                    /forceLengthCurveSettings.kToe);                    



%Evalute the normalized length of each element at lceOneHalf
%sarcomereProperties
ltitinOneHalf = 0;

if(flag_useElasticIgD==1)

    ltitinOneHalf = lceOneHalf...
      -(IGDFixedNormLengthAtOptimalFiberLength ...
      + ZLineToT12NormLengthAtOptimalFiberLength);  

else  
  ltitinOneHalf = lceOneHalf...
    -(IGDTotalNormLengthAtOptimalFiberLength ...
    + ZLineToT12NormLengthAtOptimalFiberLength);    
end

ligpOneHalf     = ltitinOneHalf /(1 + (kigp/kpevk)); 
lpevkOneHalf = ltitinOneHalf /(1 + (kpevk/kigp));


%Evaluate the proportion of the force-length curve that is of low 
%strain
strainWidth = forceLengthCurveSettings.normLengthToe-forceLengthCurveSettings.normLengthZero;  
lowStrainWidth =  (strainWidth - (1/forceLengthCurveSettings.kToe))...
                  /(1/forceLengthCurveSettings.kToe);
eZeroTest = forceLengthCurveSettings.normLengthToe ...
           -(1/forceLengthCurveSettings.kToe)...
           -(1/forceLengthCurveSettings.kToe)*lowStrainWidth;    

%Set the low strain values for each of the igp and pevk-igd to preserve
%this same proportion

ltitinZeroHalf = 0;

if(flag_useElasticIgD==1)

    ltitinZeroHalf = lceZeroHalf...
      -(IGDFixedNormLengthAtOptimalFiberLength ...
       + ZLineToT12NormLengthAtOptimalFiberLength);  

else

    ltitinZeroHalf = lceZeroHalf...
      -(IGDTotalNormLengthAtOptimalFiberLength ...
       + ZLineToT12NormLengthAtOptimalFiberLength);    

end

ligpZeroHalf     = ligpOneHalf      - ((1-lambdaECM)/kigp) ...
                 - ((1-lambdaECM)/kigp)*lowStrainWidth; 
lpevkZeroHalf = lpevkOneHalf  - ((1-lambdaECM)/kpevk) ...
                  - ((1-lambdaECM)/kpevk)*lowStrainWidth;  



%%
% Make the Igp curve
%%
forceLengthIgpCurve  = ...
  createFiberForceLengthCurve2021(ligpZeroHalf,...
                              ligpOneHalf,...
                              (1-lambdaECM),...
                              kigpZero,...
                              kigpLow,...
                              kigp,...
                              forceLengthCurveSettings.curviness,...
                              flag_computeCurveIntegrals,...
                              muscleName,...
                              flag_useOctave);  


forceLengthIgpCurve.name = sprintf('%s.%s',...
  muscleName,'forceLengthIgpCurve');
fprintf('    forceLengthIgpCurve created\n');

forceLengthIgpInverseCurve = ...
      createInverseCurve(forceLengthIgpCurve); 
fprintf('    forceLengthIgpInverseCurve created\n');

%%
% Make the PEVK-Igd curve
%%
forceLengthPevkCurve  = ...
  createFiberForceLengthCurve2021(lpevkZeroHalf,...
                              lpevkOneHalf,...
                              (1-lambdaECM),...
                              kpevkZero,...
                              kpevkLow,...
                              kpevk,...
                              forceLengthCurveSettings.curviness,...
                              flag_computeCurveIntegrals,...
                              muscleName,...
                              flag_useOctave);  

forceLengthPevkCurve.name = sprintf('%s.%s',...
  muscleName,'forceLengthPevkCurve');
fprintf('    forceLengthPevkCurve created\n');

forceLengthPevkInverseCurve = ...
      createInverseCurve(forceLengthPevkCurve); 
fprintf('    forceLengthPevkInverseCurve created\n');

here=1;                                 
                                 