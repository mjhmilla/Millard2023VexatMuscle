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

function [forceLengthProximalTitinCurve, forceLengthProximalTitinInverseCurve,...
            forceLengthDistalTitinCurve, forceLengthDistalTitinInverseCurve,...
            forceLengthIgPTitinCurve, forceLengthIgPInverseTitinCurve,...
            forceLengthPevkTitinCurve, forceLengthPevkInverseTitinCurve,...
            forceLengthIgDTitinCurve, forceLengthIgDInverseTitinCurve] ...
          = createTitinCurves2025( fiberForceLengthCurve,...                                   
                                   forceLengthCurveSettings,...
                                   forceLengthECMHalfCurve,...
                                   sarcomereProperties,...
                                   muscleName,...
                                   flag_useWLCTitinModel,...                                   
                                   flag_createTwoSidedCurves,...
                                   flag_computeCurveIntegrals,...
                                   flag_useElasticIgD,...
                                   flag_activeTitinModel,...
                                   projectFolders,...
                                   flag_useOctave)
                                 


fpeNRef = fiberForceLengthCurve.yEnd(1,2);

lambdaECM = sarcomereProperties.extraCellularMatrixPassiveForceFraction;

normPevkToActinAttachmentPoint = ...
    sarcomereProperties.normPevkToActinAttachmentPoint;

ZLineToT12NormLengthAtOptimalFiberLength = ...
  sarcomereProperties.ZLineToT12NormLengthAtOptimalFiberLength;

IGDFixedNormLengthAtOptimalFiberLength=...
    sarcomereProperties.IGDFixedNormLengthAtOptimalFiberLength;  


%%
%
% Make the reference human model with the same modified titin curves
% as the updated model. Check to see if the model can still reproduce
% Trombitas
%
%%


tmp = load(fullfile(projectFolders.output_structs_FittedModels,...
                                'defaultHumanSoleus.mat'));
defaultHumanSoleus = tmp.defaultHumanSoleus;

defaultHumanSoleus.sarcomere.normPevkToActinAttachmentPoint = ...
    normPevkToActinAttachmentPoint;

lsOptHuman = defaultHumanSoleus.sarcomere.optimalSarcomereLength;

[forceLengthProximalTitinHumanCurve, forceLengthProximalTitinInverseHumanCurve,...
 forceLengthDistalTitinHumanCurve, forceLengthDistalTitinInverseHumanCurve,...
 forceLengthIgPTitinHumanCurve, forceLengthIgPInverseTitinHumanCurve,...
 forceLengthPevkTitinHumanCurve, forceLengthPevkInverseTitinHumanCurve,...
 forceLengthIgDTitinHumanCurve, forceLengthIgDInverseTitinHumanCurve] ...
          = createCandidateTitinCurves2025(...
               fiberForceLengthCurve,...                                   
               forceLengthCurveSettings,...
               forceLengthECMHalfCurve,...
               defaultHumanSoleus.sarcomere,...
               muscleName,...
               flag_useWLCTitinModel,...                                   
               flag_createTwoSidedCurves,...
               flag_computeCurveIntegrals,...
               flag_useElasticIgD,...
               flag_activeTitinModel,...                                   
               flag_useOctave);

[forceLengthProximalTitinCurve, forceLengthProximalTitinInverseCurve,...
 forceLengthDistalTitinCurve, forceLengthDistalTitinInverseCurve,...
 forceLengthIgPTitinCurve, forceLengthIgPInverseTitinCurve,...
 forceLengthPevkTitinCurve, forceLengthPevkInverseTitinCurve,...
 forceLengthIgDTitinCurve, forceLengthIgDInverseTitinCurve] ...
          = createCandidateTitinCurves2025(...
               fiberForceLengthCurve,...                                   
               forceLengthCurveSettings,...
               forceLengthECMHalfCurve,...
               sarcomereProperties,...
               muscleName,...
               flag_useWLCTitinModel,...                                   
               flag_createTwoSidedCurves,...
               flag_computeCurveIntegrals,...
               flag_useElasticIgD,...
               flag_activeTitinModel,...                                   
               flag_useOctave);

%%
% Data
%%
fileTrombitas1998Figure5 = fullfile(projectFolders.experiments_TGFG1998, 'Trombitas1998_Figure5.csv');

dataTrombitas1998Figure5 = loadDigitizedData(fileTrombitas1998Figure5,...
                    'Sarcomere Length','PEVK Width (um)',...
                    {'a','b'},'Trombitas 1998 Figure 5');  


%%
% 1. Sample the IgP, PEVK, and IdD curves 
% 2. Form the force-length relations for the proximal and distal curves
% 3. Fit Bezier splines to each one.   
%%

flag_debug=0;
if(flag_debug==1)

    fTiNSeries = [0.01:0.01:1]' .* fpeNRef;
    n = length(fTiNSeries);
    
    %Human (for Trombitas)
    igpHSoln.x    = zeros(n,1);
    igpHSoln.y    = zeros(n,1);
    igpHSoln.dydx = zeros(n,1);
    
    %Human (for Trombitas)
    pevkHSoln.x    = zeros(n,1);
    pevkHSoln.y    = zeros(n,1);
    pevkHSoln.dydx = zeros(n,1);
    
    %Human (for Trombitas)
    igdHSoln.x    = zeros(n,1);
    igdHSoln.y    = zeros(n,1);
    igdHSoln.dydx = zeros(n,1);
    
    %Human (for Trombitas)
    lengthZ2PevkH.x = zeros(n,1);
    lengthZ2PevkH.y = zeros(n,1);
    
    %Human (for Trombitas)
    lengthZ2IgdH.x = zeros(n,1);
    lengthZ2IgdH.y = zeros(n,1);
    
    %The model
    igpSoln.x    = zeros(n,1);
    igpSoln.y    = zeros(n,1);
    igpSoln.dydx = zeros(n,1);
    
    igdSoln.x    = zeros(n,1);
    igdSoln.y    = zeros(n,1);
    igdSoln.dydx = zeros(n,1);
    
    pevkSoln.x    = zeros(n,1);
    pevkSoln.y    = zeros(n,1);
    pevkSoln.dydx = zeros(n,1);
    
    pTiSoln.x    = zeros(n,1);
    pTiSoln.y    = zeros(n,1);
    pTiSoln.dydx = zeros(n,1);
    
    dTiSoln.x    = zeros(n,1);
    dTiSoln.y    = zeros(n,1);
    dTiSoln.dydx = zeros(n,1);
    
    pTiMdlSoln.x    = zeros(n,1);
    pTiMdlSoln.y    = zeros(n,1);
    pTiMdlSoln.dydx = zeros(n,1);
    
    dTiMdlSoln.x    = zeros(n,1);
    dTiMdlSoln.y    = zeros(n,1);
    dTiMdlSoln.dydx = zeros(n,1);
    
    titinSoln.x    = zeros(n,1);
    titinSoln.y    = zeros(n,1);
    titinSoln.dydx = zeros(n,1);
    
    fpeNSoln.x    = zeros(n,1);
    fpeNSoln.y    = zeros(n,1);
    fpeNSoln.dydx = zeros(n,1);
    
    lengthZ2Pevk.x = zeros(n,1);
    lengthZ2Pevk.y = zeros(n,1);
    
    lengthZ2Igd.x = zeros(n,1);
    lengthZ2Igd.y = zeros(n,1);
    
    
    for i=1:1:n
        fTiN = fTiNSeries(i,1);
    
        %Human (for Trombitas)
        igpHSoln.y(i,1)=fTiN;    
        igpHSoln.x(i,1)= calcBezierYFcnXDerivative(fTiN,...
                            forceLengthIgPInverseTitinHumanCurve,0);
        igpHSoln.dydx(i,1) = calcBezierYFcnXDerivative(igpHSoln.x(i,1),...
                            forceLengthIgPTitinHumanCurve,1);
        
        %Human (for Trombitas)
        pevkHSoln.y(i,1)=fTiN;    
        pevkHSoln.x(i,1)= calcBezierYFcnXDerivative(fTiN,...
                            forceLengthPevkInverseTitinHumanCurve,0);
        pevkHSoln.dydx(i,1) = calcBezierYFcnXDerivative(pevkHSoln.x(i,1),...
                            forceLengthPevkTitinHumanCurve,1);
    
        %Human (for Trombitas)
        igdHSoln.y(i,1)=fTiN;    
        igdHSoln.x(i,1)= calcBezierYFcnXDerivative(fTiN,...
                            forceLengthIgDInverseTitinHumanCurve,0);
        igdHSoln.dydx(i,1) = calcBezierYFcnXDerivative(igdHSoln.x(i,1),...
                            forceLengthIgDTitinHumanCurve,1);    
        %Human
        lceN_human = 2*(igpHSoln.x(i,1)+pevkHSoln.x(i,1)+igdHSoln.x(i,1) ... 
                      + defaultHumanSoleus.sarcomere.ZLineToT12NormLengthAtOptimalFiberLength ...
                      + defaultHumanSoleus.sarcomere.IGDFixedNormLengthAtOptimalFiberLength);
    
        lengthZ2PevkH.x(i,1) = lceN_human;
        lengthZ2PevkH.y(i,1) = igpHSoln.x(i,1)...
                             +defaultHumanSoleus.sarcomere.ZLineToT12NormLengthAtOptimalFiberLength;
        
        lengthZ2IgdH.x(i,1) = lceN_human;
        lengthZ2IgdH.y(i,1) = lengthZ2PevkH.y(i,1) + pevkHSoln.x(i,1);
    
        lceN = 0;
        if(i>1)
            lceN=fpeNSoln.x(i-1,1);
        else
            lceN = fiberForceLengthCurve.yEnd(1,2);
        end
    
    
        fpeNSoln.y(i,1)=fTiN;
        fpeNSoln.x(i,1)=calcBezierFcnXGivenY(fTiN,fiberForceLengthCurve,lceN);
        fpeNSoln.dydx(i,1) = calcBezierYFcnXDerivative(fpeNSoln.x(i,1),...
                                               fiberForceLengthCurve,1);
    
        igpSoln.y(i,1) = fTiN;
        igpSoln.x(i,1) = calcBezierYFcnXDerivative(fTiN,...
                            forceLengthIgPInverseTitinCurve,0);
        igpSoln.dydx(i,1) = calcBezierYFcnXDerivative(igpSoln.x(i,1),...
                            forceLengthIgPTitinCurve,1);
    
        igdSoln.y(i,1) = fTiN;
        igdSoln.x(i,1) = calcBezierYFcnXDerivative(fTiN,...
                            forceLengthIgDInverseTitinCurve,0);
        igdSoln.dydx(i,1) = calcBezierYFcnXDerivative(igdSoln.x(i,1),...
                            forceLengthIgDTitinCurve,1);
    
        pevkSoln.y(i,1) = fTiN;
        pevkSoln.x(i,1) = calcBezierYFcnXDerivative(fTiN,...
                             forceLengthPevkInverseTitinCurve,0);
        pevkSoln.dydx(i,1) = calcBezierYFcnXDerivative(pevkSoln.x(i,1),...
                            forceLengthPevkTitinCurve,1);
    
    
        lceN = 2*(igpSoln.x(i,1)+pevkSoln.x(i,1)+igdSoln.x(i,1) ... 
                      + ZLineToT12NormLengthAtOptimalFiberLength ...
                      + IGDFixedNormLengthAtOptimalFiberLength);
    
        lengthZ2Pevk.x(i,1) = lceN;
        lengthZ2Pevk.y(i,1) = igpSoln.x(i,1)...
                             +ZLineToT12NormLengthAtOptimalFiberLength;
        
        lengthZ2Igd.x(i,1) = lceN;
        lengthZ2Igd.y(i,1) = lengthZ2Pevk.y(i,1) + pevkSoln.x(i,1);
    
        titinSoln.x(i,1) = igpSoln.x(i,1)+pevkSoln.x(i,1)+igdSoln.x(i,1);
        titinSoln.y(i,1) = fTiN;
        titinSoln.dydx(i,1) = 1/((1/igpSoln.dydx(i,1))... 
                                +(1/pevkSoln.dydx(i,1))...
                                +(1/igdSoln.dydx(i,1)));
    
        pTiSoln.x(i,1) = igpSoln.x(i,1) + pevkSoln.x(i,1).*normPevkToActinAttachmentPoint;
        pTiSoln.y(i,1) = fTiN;
    
        kigp = igpSoln.dydx(i,1);
        kigd = igdSoln.dydx(i,1);
        kpevk= pevkSoln.dydx(i,1);
    
        cigp = 1/kigp;
        cpevk= 1/kpevk;
        cigd = 1/kigd;
    
        cp = cigp + cpevk*normPevkToActinAttachmentPoint;
        cd = cigd + cpevk*(1-normPevkToActinAttachmentPoint);
    
        kp = 1/cp;
        kd = 1/cd;
    
        pTiSoln.x(i,1) = igpSoln.x(i,1) + pevkSoln.x(i,1).*(normPevkToActinAttachmentPoint);
        pTiSoln.y(i,1) = fTiN;
        pTiSoln.dydx(i,1)=kp;
    
        dTiSoln.x(i,1) = igdSoln.x(i,1) + pevkSoln.x(i,1).*(1-normPevkToActinAttachmentPoint);
        dTiSoln.y(i,1) = fTiN;
        dTiSoln.dydx(i,1)=kd;
        
        pTiMdlSoln.x(i,1) = calcBezierYFcnXDerivative(fTiN,...
                            forceLengthProximalTitinInverseCurve,0);
        pTiMdlSoln.y(i,1) = fTiN;
        pTiMdlSoln.dydx(i,1)= calcBezierYFcnXDerivative(pTiMdlSoln.x(i,1),...
                               forceLengthProximalTitinCurve,1);
    
        dTiMdlSoln.x(i,1) = calcBezierYFcnXDerivative(fTiN,...
                            forceLengthDistalTitinInverseCurve,0);
        dTiMdlSoln.y(i,1) = fTiN;
        dTiMdlSoln.dydx(i,1)= calcBezierYFcnXDerivative(pTiMdlSoln.x(i,1),...
                               forceLengthDistalTitinCurve,1);
    
    end
    
    %
    %
    % Fit the proximal and distal titin segments
    %
    %






    subplot(3,2,1);
    plot(igpSoln.x,igpSoln.y,'-b','DisplayName','IgP');
    hold on;
    plot(pevkSoln.x,pevkSoln.y,'-m','DisplayName','Pevk');
    hold on;
    plot(igdSoln.x,igdSoln.y,'-r','DisplayName','IgD');
    hold on;
    legend;
    box off;
    xlabel('Norm. Length');
    ylabel('Norm. Force')
    title('Force-length relation of titin segments');

    lceNLength = 2*(titinSoln.x ... 
                  + ZLineToT12NormLengthAtOptimalFiberLength ...
                  + IGDFixedNormLengthAtOptimalFiberLength);


    subplot(3,2,2);
    plot(lceNLength,titinSoln.y,'-g','DisplayName','Titin (passive)');
    hold on;
    plot(fpeNSoln.x,fpeNSoln.y, '-k','DisplayName','fpeN');
    hold on;
    legend;
    box off;
    xlabel('Norm. Length');
    ylabel('Norm. Force');
    title('Force-length relation: Titin vs. parallel element');


    subplot(3,2,3);
    plot(igpSoln.x,igpSoln.dydx,'-b','DisplayName','IgP');
    hold on;
    plot(pevkSoln.x,pevkSoln.dydx,'-m','DisplayName','Pevk');
    hold on;
    plot(igdSoln.x,igdSoln.dydx,'-r','DisplayName','IgD');
    hold on;
    legend;
    box off;
    xlabel('Norm. Length');
    ylabel('Norm. Stiffness')
    title('Stiffness-length relation of titin segments');

    subplot(3,2,4);
    plot(lceNLength,titinSoln.dydx.*0.5,'-g','DisplayName','Titin (passive)');
    hold on;
    plot(fpeNSoln.x,fpeNSoln.dydx, '-k','DisplayName','fpeN');
    hold on;
    legend;
    box off;    
    xlabel('Norm. Length');
    ylabel('Norm. Stiffness')
    title('Stiffness-length relation: Titin vs. parallel element'); 

    subplot(3,2,5);
    plot(lengthZ2PevkH.x.*lsOptHuman, ...
         lengthZ2PevkH.y.*lsOptHuman,...
        '-r','DisplayName','Z-PEVK (human)');
    hold on;
    plot(lengthZ2IgdH.x.*lsOptHuman, ...
         lengthZ2IgdH.y.*lsOptHuman,...
        '-r','DisplayName','Z-Igd (human)');
    hold on;

    plot( dataTrombitas1998Figure5(2).x,...
            dataTrombitas1998Figure5(2).y,...
            'o','MarkerSize',3, ...
            'MarkerFaceColor',[1,1,1].*0.9,...
            'Color',[1,1,1].*0.75,...
            'DisplayName','Exp: Z-line (ZL) to N-end 9D10');
    hold on;
    
    plot( dataTrombitas1998Figure5(1).x,...
            dataTrombitas1998Figure5(1).y,...
            'o','MarkerSize',3, ...
            'MarkerFaceColor',[1,1,1].*0.5,...
            'Color',[1,1,1].*0.5,...
            'DisplayName','Exp: ZL to C-end 9D10');
    hold on;

    legend;
    box off;
    xlabel('Sarcomere Length $$\mu m$$');
    ylabel('Length $$\mu m$$');
    title('Distance: Z-IgP/PEVK and Z-PEVK/IgD');

    subplot(3,2,6);
    plot(pTiSoln.x,pTiSoln.y,'-','Color',[0.75,0.75,1],'LineWidth',2,...
        'DisplayName','Prox. Segment Data');
    hold on;
    plot(dTiSoln.x,dTiSoln.y,'-','Color',[1,0.75,1],'LineWidth',2,...
        'DisplayName','Dist. Segment Data');
    hold on;

    plot(pTiSoln.x,pTiSoln.y,'-b','DisplayName','Prox. Segment Mdl');
    hold on;
    plot(dTiSoln.x,dTiSoln.y,'-m','DisplayName','Dist. Segment Mdl');
    hold on;


    xlabel('Norm. Length $$\ell/\ell_o^M$$');
    ylabel('Norm. Force');
    title('Two-segment titin model');
    legend;
    box off;
    here=1;

   
end


                                
                                 