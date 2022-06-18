function [forceLengthProximalTitinCurve, forceLengthProximalTitinInverseCurve,...
            forceLengthDistalTitinCurve, forceLengthDistalTitinInverseCurve] ...
          = createTitinCurves2022( fiberForceLengthCurve,...                                   
                                   forceLengthCurveSettings,...
                                   forceLengthECMHalfCurve,...
                                   sarcomereProperties,...
                                   muscleName,...
                                   flag_createTwoSidedCurves,...
                                   flag_computeCurveIntegrals,...
                                   flag_useElasticIgD,...
                                   flag_activeTitinModel,...
                                   flag_useHumanIgDGeometry,...
                                   flag_useOctave)
                                 

%%
%Evaluate the length and stiffness of fpe when it reaches a normalized
%force of 1.
%%


lambdaECM = sarcomereProperties.extraCellularMatrixPassiveForceFraction;
normPevkToActinAttachmentPoint = sarcomereProperties.normPevkToActinAttachmentPoint;

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
distalStretchRate=0;   % From the pevk-actin attachment point to myosin
proximalStretchRate=0; % From the z-line to the pevk-actin attachment point

proximalContourLength=0;
distalContourLength = 0;

if(flag_useElasticIgD==1)

  %Sticky spring: lump IgD and part of PEVK distal to attachment point 
  if(flag_activeTitinModel == 0)
    %Lump IgD with PEVK
    pevkStretchRate = sarcomereProperties.PEVKNormStretchRate ...
                        +sarcomereProperties.IGDFreeNormStretchRate;

    distalStretchRate = ...
        (1-normPevkToActinAttachmentPoint)*sarcomereProperties.PEVKNormStretchRate ...
                        + sarcomereProperties.IGDFreeNormStretchRate;  
    
    proximalStretchRate = sarcomereProperties.IGPNormStretchRate ...
      + normPevkToActinAttachmentPoint*sarcomereProperties.PEVKNormStretchRate;               

    proximalContourLength = sarcomereProperties.IGPContourLengthNorm ...
      + normPevkToActinAttachmentPoint*sarcomereProperties.PEVKContourLengthNorm;
    distalContourLength   = ...
        (1-normPevkToActinAttachmentPoint)*sarcomereProperties.PEVKContourLengthNorm ...
        +sarcomereProperties.IGDFreeContourLengthNorm;

  %Stiff spring: lump IgP and IgD together
  elseif(flag_activeTitinModel==1)

    %Lump IgD with IgP 
    proximalStretchRate = sarcomereProperties.IGPNormStretchRate ...
                        + sarcomereProperties.IGDFreeNormStretchRate;               

    distalStretchRate = sarcomereProperties.PEVKNormStretchRate;  

    proximalContourLength = sarcomereProperties.IGPContourLengthNorm;
    distalContourLength   = sarcomereProperties.PEVKContourLengthNorm ...
                          + sarcomereProperties.IGDFreeContourLengthNorm;    
  else
    assert(0,'flag_IgDLumpingMethod must be 0 or 1');  
  end

else

  if(flag_activeTitinModel == 0)
    pevkStretchRate   = sarcomereProperties.PEVKNormStretchRate;   
    distalStretchRate = normPevkToActinAttachmentPoint*sarcomereProperties.PEVKNormStretchRate;  

    proximalStretchRate = sarcomereProperties.IGPNormStretchRate ...
      + normPevkToActinAttachmentPoint*sarcomereProperties.PEVKNormStretchRate;    

    proximalContourLength = sarcomereProperties.IGPContourLengthNorm ...
      + normPevkToActinAttachmentPoint*sarcomereProperties.PEVKContourLengthNorm;
    distalContourLength   = ...
        (1-normPevkToActinAttachmentPoint)*sarcomereProperties.PEVKContourLengthNorm;    
  else
    distalStretchRate = sarcomereProperties.PEVKNormStretchRate;  
    proximalStretchRate = sarcomereProperties.IGPNormStretchRate;   

    proximalContourLength = sarcomereProperties.IGPContourLengthNorm;
    distalContourLength   = sarcomereProperties.PEVKContourLengthNorm ...
                           +sarcomereProperties.IGDFreeContourLengthNorm;    
  end
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
%
% Similarly, we can use the same approach to evaluate ka, kb, and kc if instead the 
% titin-actin attachment point is not at the N2A epitope but instead is at some point
% in the pevk section
%
%%
A  =  pevkStretchRate/igpStretchRate;
kpevk = kTitinOneHalf*(A+1)/A;
kigp     = A*kpevk; 

A = distalStretchRate/proximalStretchRate;
kD = kTitinOneHalf*(A+1)/A;
kP = A*kD;



%Geometrically scale the other stiffnesses of the force-length curve
kDLow     = kD*(forceLengthCurveSettings.kLow ...
               /forceLengthCurveSettings.kToe);

kPLow     = kP*(forceLengthCurveSettings.kLow ...
               /forceLengthCurveSettings.kToe);

kDZero    = kD*(forceLengthCurveSettings.kZero ...
               /forceLengthCurveSettings.kToe);

kPZero    = kP*(forceLengthCurveSettings.kZero ...
               /forceLengthCurveSettings.kToe); 

%Leaving this calculation in place just as a sanity check
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


lPOneHalf  = ltitinOneHalf /(1 + (kP/kD)); 
lDOneHalf  = ltitinOneHalf /(1 + (kD/kP));

%As before
ligpOneHalf  = ltitinOneHalf /(1 + (kigp/kpevk)); 
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
lpevkZeroHalf    = lpevkOneHalf  - ((1-lambdaECM)/kpevk) ...
                  - ((1-lambdaECM)/kpevk)*lowStrainWidth;  

lPZeroHalf      = lPOneHalf      - ((1-lambdaECM)/kP) ...
                 - ((1-lambdaECM)/kP)*lowStrainWidth; 
lDZeroHalf      = lDOneHalf  - ((1-lambdaECM)/kD) ...
                 - ((1-lambdaECM)/kD)*lowStrainWidth;  
%%
% Evaluate the zero, one, kzero, klow, and k values for two segments:
%   1. From the Z-line, through the IgP segment, to the desired PEVK segment fraction
%   2. From the desired PEVK segment fraction to the distal Ig/myosin boundary.
%%


fNfailure = 5.14; %Average failure force from Leonard, Joumaa, Herzog 2010


%%
% Make the P curve
%%
if flag_createTwoSidedCurves == 0
  forceLengthProximalTitinCurve  = ...
    createWLCForceLengthCurve2022(lPZeroHalf,...
                                lPOneHalf,...
                                (1-lambdaECM),...
                                proximalContourLength,...
                                fNfailure,...
                                kPZero,...
                                kPLow,...
                                kP,...
                                forceLengthCurveSettings.curviness,...
                                flag_computeCurveIntegrals,...
                                muscleName,...
                                flag_useOctave); 
else
  forceLengthProximalTitinCurve  = ...
    createTwoSidedWLCForceLengthCurve2022(lPZeroHalf,...
                                lPOneHalf,...
                                (1-lambdaECM),...
                                proximalContourLength,...
                                fNfailure,...
                                kPZero,...
                                kPLow,...
                                kP,...
                                forceLengthCurveSettings.curviness,...
                                flag_computeCurveIntegrals,...
                                muscleName,...
                                flag_useOctave); 
end


forceLengthProximalTitinCurve.name = sprintf('%s.%s',...
  muscleName,'forceLengthProximalTitinCurve');
%fprintf('    forceLengthProximalTitinCurve created\n');

forceLengthProximalTitinInverseCurve = ...
      createInverseCurve(forceLengthProximalTitinCurve); 
%fprintf('    forceLengthProximalTitinInverseCurve created\n');


%%
% Make the D curve
%%
if flag_createTwoSidedCurves == 0
  forceLengthDistalTitinCurve  = ...
    createWLCForceLengthCurve2022(lDZeroHalf,...
                                lDOneHalf,...
                                (1-lambdaECM),...
                                distalContourLength,...
                                fNfailure,...
                                kDZero,...
                                kDLow,...
                                kD,...
                                forceLengthCurveSettings.curviness,...
                                flag_computeCurveIntegrals,...
                                muscleName,...
                                flag_useOctave);  
else 
  forceLengthDistalTitinCurve  = ...
    createTwoSidedWLCForceLengthCurve2022(lDZeroHalf,...
                                lDOneHalf,...
                                (1-lambdaECM),...
                                distalContourLength,...
                                fNfailure,...
                                kDZero,...
                                kDLow,...
                                kD,...
                                forceLengthCurveSettings.curviness,...
                                flag_computeCurveIntegrals,...
                                muscleName,...
                                flag_useOctave);    
end




forceLengthDistalTitinCurve.name = sprintf('%s.%s',...
  muscleName,'forceLengthDistalTitinCurve');
%fprintf('    forceLengthDistalTitinCurve created\n');

forceLengthDistalTitinInverseCurve = ...
      createInverseCurve(forceLengthDistalTitinCurve); 



%fprintf('    forceLengthDistalTitinInverseCurve created\n');

%%
% For reference: IgP curve
%%
% forceLengthIgpCurve  = ...
%   createFiberForceLengthCurve2021(ligpZeroHalf,...
%                               ligpOneHalf,...
%                               (1-lambdaECM),...
%                               kigpZero,...
%                               kigpLow,...
%                               kigp,...
%                               forceLengthCurveSettings.curviness,...
%                               flag_computeCurveIntegrals,...
%                               muscleName,...
%                               flag_useOctave);  


flag_debugWLCExtension=0;
if(flag_debugWLCExtension == 1)

    %% Proximal curve
    % Solve for proximal curve slack length
    x    = lPOneHalf;
    y    = calcBezierYFcnXDerivative(x,forceLengthProximalTitinCurve,0);
    dydx = calcBezierYFcnXDerivative(x,forceLengthProximalTitinCurve,1);

    lPo = x-(y/dydx);

    % Solve for WLC coefficient s.t. the two curves are equal at fiso
    cPa = 1;
    fPaWLC=calcWormLikeChainModelDer(lPOneHalf-lPo,...
                                      proximalContourLength,...
                                      cPa,cPa*cPa,[0,0]);
    fPa = calcBezierYFcnXDerivative(lPOneHalf,...
                    forceLengthProximalTitinCurve,0);
    
    cPa = fPa/fPaWLC;
    % Solve for WLC length at 5.5*fiso, the length

    lPfailure = calcWormLikeChainModelInvDer(fNfailure,proximalContourLength,... 
                                          cPa,cPa*cPa,[0,0]);
    lPfailure = lPfailure + lPo;
    fNTestP=calcWormLikeChainModelDer(lPfailure-lPo,...
        proximalContourLength,cPa,cPa*cPa,[0,0]);

    %zPfailure = lPfailure/proximalContourLength;
    %zPfailure = zPfailure+zPo;
    %lPfailure = zPfailure*proximalContourLength;

    %% Distal curve    
    % Solve for distal curve slack length
    fDa = calcBezierYFcnXDerivative(lDOneHalf,...
                    forceLengthDistalTitinCurve,0);    
    x    = lDOneHalf;
    y    = calcBezierYFcnXDerivative(x,forceLengthDistalTitinCurve,0);
    dydx = calcBezierYFcnXDerivative(x,forceLengthDistalTitinCurve,1);

    lDo = x-(y/dydx);
    
    cDa = 1;
    fDaWLC=calcWormLikeChainModelDer(lDOneHalf-lDo,...
                                      distalContourLength,...
                                      cDa,cDa*cDa,[0,0]);
    cDa = fDa/fDaWLC;

    lDfailure = calcWormLikeChainModelInvDer(fNfailure,distalContourLength,... 
                                              cDa,cDa*cDa,[0,0]);
    lDfailure = lDfailure+lDo;
    fNTestD=calcWormLikeChainModelDer(lDfailure-lDo,...
        distalContourLength,cDa,cDa*cDa,[0,0]);    
    %zDfailure = lDfailure/distalContourLength;    
    %zDfailure = zDfailure+zDo;
    %lDfailure = zDfailure*distalContourLength;

    figWLCDebug = figure;
    n=100;
    n01 = [0:(1/(n-1)):1]';    
    lenP = n01.*(lPfailure);
    lenD = n01.*(lDfailure);

    fP = zeros(size(lenP));
    fD = zeros(size(lenD));
    fPwlc = zeros(size(lenP));
    fDwlc = zeros(size(lenD));

    for i=1:1:n
        fP(i,1) = calcBezierYFcnXDerivative(lenP(i,1),...
                    forceLengthProximalTitinCurve,0);
        fD(i,1) = calcBezierYFcnXDerivative(lenD(i,1),...
                    forceLengthDistalTitinCurve,0);

        if(lenP(i,1)>lPo)
            fPwlc(i,1)= calcWormLikeChainModelDer(lenP(i,1)-lPo,...
                                          proximalContourLength,...
                                          cPa,cPa*cPa,[0,0]);
        else
            fPwlc(i,1)=0;
        end

        if(lenD(i,1)>lDo)
            fDwlc(i,1)= calcWormLikeChainModelDer(lenD(i,1)-lDo,...
                                          distalContourLength,...
                                          cDa,cDa*cDa,[0,0]);
        else
            fDwlc(i,1)=0;
        end
        
    end

    zDivLP = lenP./proximalContourLength;
    zDivLD = lenD./distalContourLength;

    subplot(1,2,1);
        plot(zDivLP,fP,'b');
        hold on;
        plot(zDivLP,fPwlc,'r');
        xlabel('Norm. Length ($$\ell/L^{1}$$)')
        ylabel('Norm. Force');

    subplot(1,2,2);
        plot(zDivLD,fD,'b');
        hold on;
        plot(zDivLD,fDwlc,'r');
        xlabel('Norm. Length ($$\ell/L^{2}$$)')
        ylabel('Norm. Force');
    here=1;
end    
%%
% For reference: PEVK-Igd curve
%%
% forceLengthPevkCurve  = ...
%   createFiberForceLengthCurve2021(lpevkZeroHalf,...
%                               lpevkOneHalf,...
%                               (1-lambdaECM),...
%                               kpevkZero,...
%                               kpevkLow,...
%                               kpevk,...
%                               forceLengthCurveSettings.curviness,...
%                               flag_computeCurveIntegrals,...
%                               muscleName,...
%                               flag_useOctave);  


here=1;                                 
                                 