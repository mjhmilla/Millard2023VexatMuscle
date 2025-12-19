%%
% SPDX-FileCopyrightText: 2023 Matthew Millard <millard.matthew@gmail.com>
%
% SPDX-License-Identifier: MIT
%
%
%%

function [ratFibrilModelsFitted,...
          benchRecordFitted,...
          fitInfo] =...
             fitRatFibrilTRSS2017( ...
                    ratFibrilModelsDefault,...
                    expTRSS2017,...
                    simConfig,...
                    fittingConfig,...
                    plotConfig,...
                    projectFolders)


ratFibrilModelsFitted=ratFibrilModelsDefault;
benchRecordFitted=[];


fitInfo = struct('fl',[],'fv',[],'timeConstant',[],'Kx',[],...
                    'QToF',[],'QToK',[],'f1HNPreload',[],'l1HNOffset',[]);

fitInfoFields= fields(fitInfo);

for i=1:1:length(fitInfoFields)
    fitInfo.(fitInfoFields{i}) = ...
        struct('rmse',nan,'x',[],'yErr',[],'y',[],'yFit',[],...
               'arg',nan,'argDelta',nan);
end

fidFitting = fopen(fullfile(projectFolders.output_structs_TRSS2017,...
                ['fittingLog_',fittingConfig.trialStr,'.txt']),'w');


%%
% Fitting the active-force-length relation
%%

figDebugFitting = figure;

if(fittingConfig.fitFl==1)

    %%
    % fal fitting
    %%


    dfalN = 0.2;

    lceOptMdl   = ratFibrilModelsFitted(1).musculotendon.optimalFiberLength;
    vmax        = ratFibrilModelsFitted(1).musculotendon.maximumNormalizedFiberVelocity;
    lceOptData  = min(expTRSS2017.activeLengtheningData(3).x);

    expfl.lce = [];
    expfl.fN   = []; 
    
    mdlfl.lce = [];
    mdlfl.fN   = [];

    assert(ratFibrilModelsFitted(1).curves.useCalibratedCurves==1,...
           'Error: the calibrated curves should be used');

    curveFalN   = ratFibrilModelsFitted(1).curves.activeForceLengthCurve; 
    curveFvN    = ratFibrilModelsFitted(1).curves.fiberForceVelocityCurve; 

    for idxTrial = simConfig.trials
        lce  = expTRSS2017.activeLengtheningData(idxTrial).x(1,1);
        fN    = expTRSS2017.activeLengtheningData(idxTrial).y(1,1);
        
        expfl.lce  = [expfl.lce;lce];
        expfl.fN    = [expfl.fN;fN];   
        
        fNMdl = calcBezierYFcnXDerivative(lce./lceOptMdl,curveFalN,0);
        mdlfl.lce  = [mdlfl.lce;lce];
        mdlfl.fN    = [mdlfl.fN;fNMdl];                
    end


    argBest  = 0.2;
    flag_compensateForCrossbridgeStiffness = 0;
    [errFlN,flNFit,curveL] = calcErrorTRSS2017ForceLengthRelationAscendingLimb(argBest,...
                        expfl.lce,expfl.fN,...
                        ratFibrilModelsFitted(1).sarcomere,...
                        ratFibrilModelsFitted(1).musculotendon,...
                        flag_compensateForCrossbridgeStiffness);
    errBest = sqrt(mean(errFlN.^2));
    errValBest=errFlN;
    fprintf('%1.2e\tfitting: fal rmse (start)\n',errBest);
    fprintf('%e\tfal-asc offset (start)\n',argBest);

    fprintf(fidFitting,'%1.2e\tfitting: fal rmse (start)\n',errBest);
    fprintf(fidFitting,'%e\tfal-asc offset (start)\n',argBest);    

    argDelta = argBest/2;

    for i=1:1:fittingConfig.numberOfBisections

        argL = argBest-argDelta;
        [errL,flNFitL,curveL] = ...
            calcErrorTRSS2017ForceLengthRelationAscendingLimb(argL,...
                            expfl.lce,expfl.fN,...
                            ratFibrilModelsFitted(1).sarcomere,...
                            ratFibrilModelsFitted(1).musculotendon,...
                            flag_compensateForCrossbridgeStiffness);
        errLMag = sqrt(mean(errL.^2));

        argR = argBest+argDelta;
        [errR,flNFitR,curveR] = ...
            calcErrorTRSS2017ForceLengthRelationAscendingLimb(argR,...
                            expfl.lce,expfl.fN,...
                            ratFibrilModelsFitted(1).sarcomere,...
                            ratFibrilModelsFitted(1).musculotendon,...
                            flag_compensateForCrossbridgeStiffness);
        errRMag = sqrt(mean(errR.^2));        

        if(errLMag < errBest && errLMag <= errRMag )
            argBest=argL;
            errBest=errLMag;
            fitInfo.fl.rmse = errBest;
            fitInfo.fl.x = expfl.lce';
            fitInfo.fl.y = expfl.fN';
            fitInfo.fl.yFit = flNFitL';            
            fitInfo.fl.yErr = errL';
            fitInfo.fl.arg = argBest;
            fitInfo.fl.argDelta = argDelta;
        elseif(errRMag < errBest && errRMag < errLMag)
            argBest=argR;
            errBest=errRMag;
            fitInfo.fl.rmse = errBest;
            fitInfo.fl.x = expfl.lce';
            fitInfo.fl.y = expfl.fN';
            fitInfo.fl.yFit = flNFitR';                        
            fitInfo.fl.yErr = errR';
            fitInfo.fl.arg = argBest;
            fitInfo.fl.argDelta = argDelta;
        end

        argDelta=argDelta*0.5;

    end



    fprintf('%1.2e\tfitting: fal rmse (end)\n',errBest);
    fprintf('%e\tfal-asc offset (end)\n\n',argBest);

    fprintf(fidFitting,'%1.2e\tfitting: fal rmse (end)\n\n',errBest);
    fprintf(fidFitting,'%e\tfal-asc offset (end)\n\n',argBest);
       
    %
    % Update all of the models
    %
    for i=simConfig.trials

        flag_compensateForCrossbridgeStiffness=0;
        [optError, flNV,falCurve] = ...
            calcErrorTRSS2017ForceLengthRelationAscendingLimb(...
                argBest, expfl.lce,expfl.fN,...
                ratFibrilModelsFitted(1).sarcomere,...
                ratFibrilModelsFitted(1).musculotendon,...
                flag_compensateForCrossbridgeStiffness);

        ratFibrilModelsFitted(i).curves.activeForceLengthCurve=falCurve;

        flag_compensateForCrossbridgeStiffness=1;
        [optErrorCal, flNV, falCurveCal] = ...
            calcErrorTRSS2017ForceLengthRelationAscendingLimb(...
                argBest, expfl.lce,expfl.fN,...
                ratFibrilModelsFitted(1).sarcomere,...
                ratFibrilModelsFitted(1).musculotendon,...
                flag_compensateForCrossbridgeStiffness);

        ratFibrilModelsFitted(i).curves.activeForceLengthCalibratedCurve=falCurveCal;
        
        curveBest= falCurve;
    end


    if(simConfig.flag_debugFitting==1)

        subplot('Position',reshape(plotConfig.subPlotPanel(1,1,:),1,4));

        sampleFlN=calcBezierYFcnXCurveSampleVector(curveBest,100,curveBest.xEnd);
        plot(sampleFlN.x.*lceOptMdl,sampleFlN.y,'-',...
             'Color',[1,1,1].*0.75,'DisplayName','$$f^L(\ell^M)$$');
        hold on;


        for idxTrial = simConfig.trials
            plot(expfl.lce(idxTrial),...
                 expfl.fN(idxTrial),...
                 'x','Color',plotConfig.lineColors.exp(idxTrial,:),...
                 'MarkerFaceColor',plotConfig.lineColors.exp(idxTrial,:),...
                 'DisplayName',...
                 expTRSS2017.activeLengtheningData(idxTrial).seriesName);
            hold on;
        end
        

        xlabel(expTRSS2017.activeLengtheningData(3).xName);
        ylabel(expTRSS2017.activeLengtheningData(3).yName);
        title('Force-Length Relation');  
        box off;
        hold on;

    end    

end

if(fittingConfig.fitFv==1)
    %%
    % fv fitting
    %%

    lceOptMdl   = ratFibrilModelsFitted(1).musculotendon.optimalFiberLength;
    vmax        = ratFibrilModelsFitted(1).musculotendon.maximumNormalizedFiberVelocity;
    lceOptData  = min(expTRSS2017.activeLengtheningData(3).x);

    mdlfv.lceN = [];
    mdlfv.vceN = [];  
    mdlfv.fN   = [];
    mdlfv.fvN   = [];

    expfv.lceN = [];
    expfv.vceN = [];  
    expfv.fN   = [];
    expfv.fvN   = [];

    for idxTrial = simConfig.trials
        idxKey = expTRSS2017.activeLengtheningData(idxTrial).keyIndices(...
                    fittingConfig.idxFvKey); 
        lce  = expTRSS2017.activeLengtheningData(idxTrial).x(idxKey,1);
        lceN = lce/lceOptMdl;
        fN0  = expTRSS2017.activeLengtheningData(idxTrial).y(1,1);
        fN1  = expTRSS2017.activeLengtheningData(idxTrial).y(idxKey,1);

        fN1 = (fN1-fN0)*fittingConfig.fv.scaleEnhancement + fN0;

        fvN = fN1/fN0;
        
        expfv.lceN  = [expfv.lceN;lceN];
        expfv.vceN  = [expfv.vceN;0.11];
        expfv.fN  = [expfv.fN;fN1];
        expfv.fvN  = [expfv.fvN;fvN];   
        
        fvNMdl = calcBezierYFcnXDerivative(0.11,curveFvN,0);
        falNMdl = calcBezierYFcnXDerivative(lce/lceOptMdl,curveFalN,0);

        mdlfv.lceN = [mdlfv.lceN;lceN];
        mdlfv.vceN  = [mdlfv.vceN;0.11];
        mdlfv.fN    = [mdlfv.fN;fvNMdl*falNMdl];
        mdlfv.fvN    = [mdlfv.fvN;fvNMdl];                
    end

    arg = 0;
    argBest=arg;
    delta = 0.4*(...
        ratFibrilModelsFitted(1).curves.fiberForceVelocityCalibratedCurve.ypts(2,3)-1);
    

    [fvNErrorV,fvNfitV, fvCurveBest] = ...
        calcErrorTRSS2017ForceVelocityRelation(...
            arg, expfv, ...
            ratFibrilModelsFitted(1).curves.activeForceLengthCurve,...
            ratFibrilModelsFitted(1).musculotendon);
    errBest = sqrt(mean(fvNErrorV.^2));

    fprintf('%1.2e\tfitting: fv rmse (start)\n',errBest);
    fprintf('%e\tfv ecc offset (start)\n',argBest);

    fprintf(fidFitting,'%1.2e\tfitting: fv rmse (start)\n',errBest);
    fprintf(fidFitting,'%e\tfv ecc offset (start)\n',argBest);


    for i=1:1:fittingConfig.numberOfBisections

        arg = argBest-delta;

        [fvNErrorV,fvNfitV, fvCurve] = ...
        calcErrorTRSS2017ForceVelocityRelation(...
            arg, expfv, ...
            ratFibrilModelsFitted(1).curves.activeForceLengthCurve,...
            ratFibrilModelsFitted(1).musculotendon);
        fvNRmse = sqrt(mean(fvNErrorV.^2));

        if(fvNRmse<errBest)
            argBest=arg;
            errBest=fvNRmse;
            %errValBest=fvNErrorV;
            fvCurveBest=fvCurve;
            fitInfo.fv.rmse = errBest;
            fitInfo.fv.x    = (expfv.lceN*lceOptMdl)';
            fitInfo.fv.y    = expfv.fvN';
            fitInfo.fv.yFit = fvNfitV';
            fitInfo.fv.yErr = errValBest';
            fitInfo.fv.arg  = argBest;
            fitInfo.fv.argDelta = delta;
            
        else
            arg = argBest+delta;
            [fvNErrorV,fvNfitV, fvCurve] = ...
            calcErrorTRSS2017ForceVelocityRelation(...
                arg, expfv, ...
                ratFibrilModelsFitted(1).curves.activeForceLengthCurve,...
                ratFibrilModelsFitted(1).musculotendon);       
            fvNRmse = sqrt(mean(fvNErrorV.^2));
            
            if(fvNRmse<errBest)
                argBest=arg;
                errBest=fvNRmse;
                fvCurveBest=fvCurve;                
                fitInfo.fv.rmse = errBest;
                fitInfo.fv.x    = (expfv.lceN*lceOptMdl)';
                fitInfo.fv.y    = expfv.fvN';
                fitInfo.fv.yFit = fvNfitV';
                fitInfo.fv.yErr = errValBest';
                fitInfo.fv.arg  = argBest;
                fitInfo.fv.argDelta = delta;
            end
        end

        delta=delta*0.5;
    end

    fprintf('%1.2e\tfitting: fv rmse (end)\n',errBest);
    fprintf('%e\tfv ecc offset (end)\n',argBest);

    fprintf(fidFitting,'%1.2e\tfitting: fv rmse (end)\n',errBest);
    fprintf(fidFitting,'%e\tfv ecc offset (end)\n',argBest);
    

    %
    % Update all of the models
    %

    for i=simConfig.trials

        ratFibrilModelsFitted(i).curves.fiberForceVelocityCurve=fvCurveBest;

        ratFibrilModelsFitted(i).curves.fiberForceVelocityInverseCurve = ...
            createInverseCurve(fvCurveBest);
        
        fvCalCurve=fvCurveBest;

        fvCalCurve.xpts = fvCalCurve.xpts ...
          .*ratFibrilModelsFitted(i).sarcomere.forceVelocityCalibrationFactor;
        
        fvCalCurve.xEnd = fvCalCurve.xEnd ...
          .*ratFibrilModelsFitted(i).sarcomere.forceVelocityCalibrationFactor;
        
        fvCalCurve.dydxEnd = fvCalCurve.dydxEnd ...
          ./ratFibrilModelsFitted(i).sarcomere.forceVelocityCalibrationFactor;        

        ratFibrilModelsFitted(i).curves.fiberForceVelocityCalibratedCurve=...
            fvCalCurve;

        ratFibrilModelsFitted(i).curves.fiberForceVelocityCalibratedInverseCurve=...
            createInverseCurve(fvCalCurve);

    end   
    here=1;
end



%%
% Fitting the properties of the model for the very first 100 ms of data
%%
fittingFraction = 0.03; %Just to the first enhancement and recovery
npts = round(200*fittingFraction);


if(fittingConfig.fitTimeConstant==1)

    lceOptMdl   = ratFibrilModelsFitted(1).musculotendon.optimalFiberLength;
    lceOptData  = min(expTRSS2017.activeLengtheningData(3).x);

    if(simConfig.flag_debugFitting==1)

        for idxTrial=simConfig.trials
            if(simConfig.flag_debugFitting==1)
                figure(figDebugFitting);
                subplot('Position',reshape(plotConfig.subPlotPanel(2,1,:),1,4));
                for idxTrial = simConfig.trials
                    plot(expTRSS2017.activeLengtheningData(idxTrial).x,...
                         expTRSS2017.activeLengtheningData(idxTrial).y,...
                         '-','Color',plotConfig.lineColors.exp(idxTrial,:),...
                         'DisplayName',...
                         expTRSS2017.activeLengtheningData(idxTrial).seriesName);
                    hold on;
                end
                xlabel(expTRSS2017.activeLengtheningData(3).xName);
                ylabel(expTRSS2017.activeLengtheningData(3).yName);
                title(expTRSS2017.activeLengtheningData(3).title);  
                box off;
                hold on;
            end
        end
    end

    %
    % Fit the lengthening time constant
    %
    
    bestValue = 30;
    deltaValue=bestValue*0.5;

    optParams.name = 'responseTimeScaling';
    optParams.value=bestValue;
    simConfig.flag_debugFitting=0;

    [bestError,bestErrorValues,figDebugFitting,...
        ratFibrilModelsUpd,benchRecord] ...
        = calcErrorTRSS2017RampFraction(optParams,...
                   fittingFraction, npts, ...
                   ratFibrilModelsFitted, expTRSS2017,simConfig,...
                   figDebugFitting,plotConfig.subPlotPanel,...
                   plotConfig.lineColors.simXE);

    flag_debugTimeConstantFitting=1;
    if(simConfig.flag_debugFitting==1 || flag_debugTimeConstantFitting)
        figure(figDebugFitting);
        subplot('Position',reshape(plotConfig.subPlotPanel(2,1,:),1,4));

        for i=1:1:3
            plot(benchRecord.normFiberLength(:,i).*lceOptMdl,...
                 benchRecord.normFiberForce(:,i),'-b');
            hold on;
        end
    end
    fprintf('%1.2e\tfitting: sliding time-constant (start)\n',bestError);
    fprintf('%e\t sliding-time constant scaling (end)\n\n',bestValue);

    fprintf(fidFitting,'%1.2e\tfitting: sliding time-constant (start)\n',bestError);
    fprintf(fidFitting,'%e\t sliding-time constant scaling (start)\n\n',bestValue);

    benchRecordBest=[];

    for i=1:1:fittingConfig.numberOfBisections
        fprintf('%i/%i\n',i,fittingConfig.numberOfBisections);

        optParams.value=bestValue-deltaValue;

        [errorVal,errorValues,figDebugFitting,...
            ratFibrilModelsUpd,benchRecord] ...
            = calcErrorTRSS2017RampFraction(optParams,...
                   fittingFraction, npts, ...
                   ratFibrilModelsFitted, expTRSS2017,simConfig,...
                   figDebugFitting,plotConfig.subPlotPanel,...
                   plotConfig.lineColors.simXE);

        if(errorVal < bestError)
           bestError=errorVal;  
           bestErrorValues = errorValues;
           bestValue = optParams.value;
           ratFibrilModelsFitted=ratFibrilModelsUpd;
           benchRecordBest=benchRecord;
        else
            optParams.value=bestValue+deltaValue;
            [errorVal,errorValues,figDebugFitting,...
                ratFibrilModelsUpd,benchRecord] ...
                = calcErrorTRSS2017RampFraction(optParams,...
                       fittingFraction, npts, ...
                       ratFibrilModelsFitted, expTRSS2017,simConfig,...
                       figDebugFitting,plotConfig.subPlotPanel,...
                       plotConfig.lineColors.simXE);
            if(errorVal < bestError)
                bestError=errorVal;        
                bestErrorValues = errorValues;
                bestValue = optParams.value;               
                ratFibrilModelsFitted=ratFibrilModelsUpd;
                benchRecordBest=benchRecord;               
            end
        end
        deltaValue=deltaValue*0.5;
    end

    if(simConfig.flag_debugFitting==1 || flag_debugTimeConstantFitting)
        figure(figDebugFitting);
        subplot('Position',reshape(plotConfig.subPlotPanel(2,1,:),1,4));

        for i=1:1:3
            plot(benchRecord.normFiberLength(:,i).*lceOptMdl,...
                 benchRecord.normFiberForce(:,i),'-b');
            hold on;
        end
    end    
    fprintf('%1.2e\tfitting: sliding time-constant (end)\n',bestError);
    fprintf('%e\t sliding-time constant scaling (end)\n\n',bestValue);

    fitInfo.timeConstant.rmse = bestError;
    fitInfo.timeConstant.x  = bestErrorValues.x;
    fitInfo.timeConstant.y  = bestErrorValues.y;
    fitInfo.timeConstant.yFit  = bestErrorValues.yFit;
    fitInfo.timeConstant.yErr  = bestErrorValues.yErr;
    fitInfo.timeConstant.arg  = bestValue;
    fitInfo.timeConstant.argDelta = deltaValue*2;
    
end


if(fittingConfig.fitKx==1)
    %
    % Scale the stiffness and cross-bridge damping
    %
    bestValue = 1;
    deltaValue=bestValue*0.5;

    optParams.name = 'xeStiffnessDampingScaling';
    optParams.value=bestValue;
    simConfig.flag_debugFitting=0;

    [bestError,bestErrorValues,figDebugFitting,...
        ratFibrilModelsUpd,benchRecord]...
        = calcErrorTRSS2017RampFraction(optParams,...
                   fittingFraction, npts, ...
                   ratFibrilModelsFitted, expTRSS2017,simConfig,...
                   figDebugFitting,plotConfig.subPlotPanel,...
                   plotConfig.lineColors.simXE);

    fprintf('%1.2e\tfitting: impedance scaling (start)\n',bestError);
    fprintf('%e\t impedance scaling (start)\n\n',bestValue);
    fprintf(fidFitting,'%1.2e\tfitting: impedance scaling (start)\n',bestError);
    fprintf(fidFitting,'%e\t impedance scaling (start)\n\n',bestValue);
    
    for i=1:1:fittingConfig.numberOfBisections
        fprintf('%i/%i\n',i,fittingConfig.numberOfBisections);

        optParams.value=bestValue-deltaValue;

        [errorVal,errorValues,figDebugFitting,...
            ratFibrilModelsUpd,benchRecord] ...
            = calcErrorTRSS2017RampFraction(optParams,...
                   fittingFraction, npts, ...
                   ratFibrilModelsFitted, expTRSS2017,simConfig,...
                   figDebugFitting,plotConfig.subPlotPanel,...
                   plotConfig.lineColors.simXE);

        if(errorVal < bestError)
           bestError=errorVal;           
           bestValue = optParams.value; 
           bestErrorValues = errorValues;
           ratFibrilModelsFitted=ratFibrilModelsUpd;
        else
            optParams.value=bestValue+deltaValue;
            [errorVal,errorValues,figDebugFitting,...
                ratFibrilModelsUpd,benchRecord] ...
                = calcErrorTRSS2017RampFraction(optParams,...
                       fittingFraction, npts, ...
                       ratFibrilModelsFitted, expTRSS2017,simConfig,...
                       figDebugFitting,plotConfig.subPlotPanel,...
                       plotConfig.lineColors.simXE);
            if(errorVal < bestError)
                bestError=errorVal;      
                bestValue = optParams.value;
                bestErrorValues = errorValues;
                ratFibrilModelsFitted=ratFibrilModelsUpd;
            end
        end
        deltaValue=deltaValue*0.5;
    end

    fitInfo.Kx.rmse = bestError;
    fitInfo.Kx.x    = bestErrorValues.x;
    fitInfo.Kx.y    = bestErrorValues.y;
    fitInfo.Kx.yFit = bestErrorValues.yFit;
    fitInfo.Kx.yErr = bestErrorValues.yErr;
    fitInfo.Kx.arg  = bestValue;
    fitInfo.Kx.argDelta = deltaValue*2;
    
    fprintf('%1.2e\tfitting: impedance scaling (end)\n',bestError);
    fprintf('%e\t impedance scaling (end)\n\n',bestValue);
    fprintf(fidFitting,'%1.2e\tfitting: impedance scaling (end)\n',bestError);    
    fprintf(fidFitting,'%e\t impedance scaling (end)\n\n',bestValue);
    

end

%%
% Fit the parameters of titin to get the closest match possible during
% the trial. Here we will fit two parameters:
%
% Q : the point within the PEVK segment that attaches to actin. A value of
%     0 corresponds to the N2A epitope (the most proximal end) while a 
%     value 1 corresponds to the most distal end of the PEVK segment
%
% f1HNPreload: the preload force that the proximal segment develops
%
% Unlike the previous fitting routines, here we need to:
% 1. Fit to one of the trials, and then apply the parameters to all
% 2. Individually fit each trial
%%

optParams.exp(3) = struct('name','','value',0,'x',[],'y',[],'dydx',[],...
                          'xLine',[],'yLine',[]);

for idxTrial=simConfig.trials
    idx1 = length(expTRSS2017.activeLengtheningData(idxTrial).x);

    x0 = expTRSS2017.activeLengtheningData(idxTrial).x(1,1);
    x1 = expTRSS2017.activeLengtheningData(idxTrial).x(end);
    xStart= x0 + 0.125*(x1-x0);
    expXTmp = [0:0.1:1]'.*(x1-xStart) + xStart;

    [expLceU,iq] = unique(expTRSS2017.activeLengtheningData(idxTrial).x);
            
    expYTmp = interp1(expTRSS2017.activeLengtheningData(idxTrial).x(iq),...
                   expTRSS2017.activeLengtheningData(idxTrial).y(iq),...
                   expXTmp);

    A = [expXTmp, ones(size(expXTmp))];
    b = expYTmp;
    p = (A'*A)\(A'*b);
    expSlope = p(1,1);

    xLine = expXTmp;
    yLine = [xLine, ones(size(xLine))]*p;

    lopt=ratFibrilModelsFitted(idxTrial).musculotendon.optimalFiberLength;
    fiso=ratFibrilModelsFitted(idxTrial).musculotendon.fiso;

    optParams.exp(idxTrial).x = expXTmp;
    optParams.exp(idxTrial).y = expYTmp;
    optParams.exp(idxTrial).dydx = expSlope;
    optParams.exp(idxTrial).xLine=xLine;
    optParams.exp(idxTrial).yLine=yLine;

end


if(fittingConfig.fitQToF ==1 || fittingConfig.fitQToK == 1)
    assert(~(fittingConfig.fitQToF && fittingConfig.fitQToK),...
      'Error: fitting Q to force and also Q to stiffness does not make sense');
    
    %
    % Solve for Q to best fit the average slope of the force development
    %
    flagDebug=1;
    simConfigTmp = simConfig;
    figDebugFittingQ=figure;
    
    loops = length(fittingConfig.titin.trials);
    if(fittingConfig.titin.individuallyFit==0)
        loops=1;
    end
    
    QInit = 0.5;
    QDeltaInit = 0.5*QInit;
    
    

    for idxLoop = 1:1:loops
        idxTrial = nan;
        if(fittingConfig.titin.individuallyFit==1)
            idxTrial = simConfig.trials(1,idxLoop);
            simConfigTmp.trials = idxTrial;
        end

        if(fittingConfig.fitQToF)
            optParams.name = 'QToF';
        end
        if(fittingConfig.fitQToK)
            optParams.name = 'QToK';
        end

        optParams.value=  QInit;
        
        simConfigTmp.flag_debugFitting=0;
        fittingFraction=1;
        npts=100;
        [optError,optErrorValues,figDebugFittingQ,...
            ratFibrilModelsFittedUpd,benchRecord] =...
            calcErrorTRSS2017RampFraction(optParams,...
                       fittingFraction, npts, ...
                       ratFibrilModelsFitted, expTRSS2017,simConfigTmp,...
                       figDebugFittingQ,plotConfig.subPlotPanel,...
                       plotConfig.lineColors.simTitinK);

        optErrorBest=optError;
        optErrorValuesBest=optErrorValues;
        QBest=QInit;
        QDelta=QDeltaInit;
        dirMap = [1,-1];



        for i=1:1:fittingConfig.numberOfBisections
            fprintf('%i/%i\n',i,fittingConfig.numberOfBisections);
        
            for j=1:1:length(dirMap)
                Q = QBest + QDelta*dirMap(1,j);
                optParams.value=Q;
                [optError,optErrorValues,figDebugFittingQ,...
                    ratFibrilModelsFittedUpd,benchRecord] =...
                    calcErrorTRSS2017RampFraction(optParams,...
                       fittingFraction, npts, ...
                       ratFibrilModelsFitted, expTRSS2017,simConfigTmp,...
                       figDebugFittingQ,plotConfig.subPlotPanel,...
                       plotConfig.lineColors.simTitinK);
                if(optError<optErrorBest)
                    optErrorBest=optError;
                    optErrorValuesBest=optErrorValues;
                    QBest=Q;
                    break;
                end
            end
            
            QDelta=QDelta*0.5;

        end

        if(fittingConfig.fitQToF)
            if(fittingConfig.titin.individuallyFit==1)
                fprintf('%1.2e\tfitting trial %i: Q rmse force profile (start)\n',...
                        optErrorBest,idxTrial);
                fprintf('%e\tQ (start)\n',QBest);    
        
                fprintf(fidFitting,...
                    '%1.2e\tfitting trial %i: Q rmse force profile (start)\n',...
                    optErrorBest,idxTrial);
                fprintf(fidFitting,...
                    '%e\tQ (start)\n',QBest);    
            else
                fprintf('%1.2e\tfitting all: Q rmse force profile (start)\n',...
                        optErrorBest);
                fprintf('%e\tQ (start)\n',QBest);    
        
                fprintf(fidFitting,...
                    '%1.2e\tfitting all: Q rmse force profile (start)\n',...
                    optErrorBest);
                fprintf(fidFitting,...
                    '%e\tQ (start)\n',QBest);    
            end            
            if(fittingConfig.titin.individuallyFit==1)
                if(isnan(fitInfo.QToF.rmse))
                    fitInfo.QToF.rmse = optErrorBest;
                    fitInfo.QToF.x  = optErrorValuesBest.x(:,idxTrial);
                    fitInfo.QToF.y     = optErrorValuesBest.y(:,idxTrial);
                    fitInfo.QToF.yFit  = optErrorValuesBest.yFit(:,idxTrial);                    
                    fitInfo.QToF.yErr  = optErrorValuesBest.yErr(:,idxTrial);
                    fitInfo.QToF.arg  = QBest;            
                    fitInfo.QToF.argDelta = QDelta*2;

                else
                    fitInfo.QToF.rmse = [fitInfo.QToF.rmse, optErrorBest];
                    fitInfo.QToF.x  = [fitInfo.QToF.x,  optErrorValuesBest.x(:,idxTrial)];
                    fitInfo.QToF.y      = [fitInfo.QToF.y,      optErrorValuesBest.y(:,idxTrial)];
                    fitInfo.QToF.yFit   = [fitInfo.QToF.yFit,   optErrorValuesBest.yFit(:,idxTrial)];
                    fitInfo.QToF.yErr   = [fitInfo.QToF.yErr,   optErrorValuesBest.yErr(:,idxTrial)];
                    fitInfo.QToF.arg  = [fitInfo.QToF.arg,  QBest];            
                    fitInfo.QToF.argDelta = QDelta*2;
                end
            else
                fitInfo.QToF.rmse = optErrorBest;
                fitInfo.QToF.x  = optErrorValuesBest.x;
                fitInfo.QToF.y     = optErrorValuesBest.y;
                fitInfo.QToF.yFit  = optErrorValuesBest.yFit;
                fitInfo.QToF.yErr  = optErrorValuesBest.yErr;                
                fitInfo.QToF.arg  = QBest;           
                fitInfo.QToF.argDelta = QDelta*2;
            end
            fprintf('%1.2e\tfitting: Q rmse force profile (end)\n',optErrorBest);
            fprintf('%e\tQ (end)\n',QBest);    
            fprintf(fidFitting,'%1.2e\tfitting: Q rmse force profile (end)\n',optErrorBest);
            fprintf(fidFitting,'%e\tQ (end)\n',QBest);              
        end

        if(fittingConfig.fitQToK)
            if(fittingConfig.titin.individuallyFit==1)
                fprintf('%1.2e\tfitting trial %i: Q rmse terminal slope (start)\n',...
                        optErrorBest,idxTrial);
                fprintf('%e\tQ (start)\n',QBest);    
        
                fprintf(fidFitting,...
                    '%1.2e\tfitting trial %i: Q rmse terminal slope (start)\n',...
                    optErrorBest,idxTrial);
                fprintf(fidFitting,...
                    '%e\tQ (start)\n',QBest);    
            else
                fprintf('%1.2e\tfitting all: Q rmse terminal slope (start)\n',...
                        optErrorBest);
                fprintf('%e\tQ (start)\n',QBest);    
        
                fprintf(fidFitting,...
                    '%1.2e\tfitting all: Q rmse terminal slope (start)\n',...
                    optErrorBest);
                fprintf(fidFitting,...
                    '%e\tQ (start)\n',QBest);    
            end      
                  
            if(fittingConfig.titin.individuallyFit==1)
                if(isnan(fitInfo.QToK.rmse))
                    fitInfo.QToK.rmse = optErrorBest;
                    fitInfo.QToK.x  = optErrorValuesBest.x(:,idxTrial);
                    fitInfo.QToK.y      = optErrorValuesBest.y(:,idxTrial);
                    fitInfo.QToK.yFit   = optErrorValuesBest.yFit(:,idxTrial);
                    fitInfo.QToK.yErr   = optErrorValuesBest.yErr(:,idxTrial);
                    fitInfo.QToK.arg  = QBest;   
                    fitInfo.QToK.argDelta = QDelta*2;
                else
                    fitInfo.QToK.rmse = [fitInfo.QToK.rmse, optErrorBest];
                    fitInfo.QToK.x  = [fitInfo.QToK.x,  optErrorValuesBest.x(:,idxTrial)];
                    fitInfo.QToK.y    = [fitInfo.QToK.y,   optErrorValuesBest.y(:,idxTrial)];
                    fitInfo.QToK.yFit = [fitInfo.QToK.yFit,optErrorValuesBest.yFit(:,idxTrial)];
                    fitInfo.QToK.yErr = [fitInfo.QToK.yErr,optErrorValuesBest.yErr(:,idxTrial)];
                    fitInfo.QToK.arg  = [fitInfo.QToK.arg,  QBest];            
                    fitInfo.QToK.argDelta = QDelta*2;
                end
            else
                fitInfo.QToK.rmse = optErrorBest;
                fitInfo.QToK.x  = optErrorValuesBest.x;
                fitInfo.QToK.y      = optErrorValuesBest.y;
                fitInfo.QToK.yFit   = optErrorValuesBest.yFit;
                fitInfo.QToK.yErr   = optErrorValuesBest.yErr;
                fitInfo.QToK.arg  = QBest;   
                fitInfo.QToK.argDelta = QDelta*2;
            end
            fprintf('%1.2e\tfitting: Q rmse terminal slope (end)\n',optErrorBest);
            fprintf('%e\tQ (end)\n',QBest);    
            fprintf(fidFitting,'%1.2e\tfitting: Q rmse terminal slope (end)\n',optErrorBest);
            fprintf(fidFitting,'%e\tQ (end)\n',QBest);              
        end

  

        %
        % Update the parameter struct
        % 
        optParams.value=QBest;
        [optError,optErrorValues, figDebugFittingQ,...
            ratFibrilModelsFittedUpd,benchRecord] =...
            calcErrorTRSS2017RampFraction(optParams,...
                       fittingFraction, npts, ...
                       ratFibrilModelsFitted, expTRSS2017,simConfigTmp,...
                       figDebugFittingQ,plotConfig.subPlotPanel,...
                       plotConfig.lineColors.simTitinK);

        if(fittingConfig.titin.individuallyFit==1)
            ratFibrilModelsFitted(idxTrial)=ratFibrilModelsFittedUpd(idxTrial);
        else
            ratFibrilModelsFitted=ratFibrilModelsFittedUpd;
        end

        if(flagDebug==1)
            if(exist('figExpSlope')==0)
                figExpSlope = figure;
            end
            if(fittingConfig.titin.individuallyFit==1)
                plot(expTRSS2017.activeLengtheningData(idxTrial).x,...
                     expTRSS2017.activeLengtheningData(idxTrial).y,'-k');
                hold on;
                plot(optParams.exp(idxTrial).x,optParams.exp(idxTrial).y,'xk');
                hold on;
                plot(optParams.exp(idxTrial).xLine,optParams.exp(idxTrial).yLine,'-b');
                hold on;
                plot(benchRecord.normFiberLength(:,idxTrial).*lopt,...
                     benchRecord.normFiberForce(:,idxTrial),'-r' );
            
            else
                for idxTrial=simConfig.trials
                    plot(expTRSS2017.activeLengtheningData(idxTrial).x,...
                         expTRSS2017.activeLengtheningData(idxTrial).y,'-k');
                    hold on;
                    plot(optParams.exp(idxTrial).x,optParams.exp(idxTrial).y,'xk');
                    hold on;
                    plot(optParams.exp(idxTrial).xLine,optParams.exp(idxTrial).yLine,'-b');
                    hold on;
                    plot(benchRecord.normFiberLength(:,idxTrial).*lopt,...
                         benchRecord.normFiberForce(:,idxTrial),'-r' );
                end

            end
            xlabel('X');
            ylabel('Y');
            title(sprintf('Trial %i',idxTrial));
        end
        
    end

    
    %
    % If we are fitting just one of the trials, then update the others
    % to have the same fitted Q value
    %
    if(length(fittingConfig.titin.trials)==1 ...
            && fittingConfig.titin.applyToAllTrials==1)
        for i=1:1:length(ratFibrilModelsFitted)
            if(i ~= fittingConfig.titin.trials(1,1))
                ratFibrilModelsFitted(i)=ratFibrilModelsFitted(fittingConfig.titin.trials);
            end
        end
    end
   
end
%
% Solve for f1HNPreload: the preload that the proximal segment of titin
%                        develops to best match the observed
%                        active-lengthening force profile.
%
if(fittingConfig.fitf1HNPreload == 1)
    
    flagDebug=1;
    simConfigTmp=simConfig;
    figDebugFittingF1HNPreload=figure;

    loops = length(fittingConfig.titin.trials);
    if(fittingConfig.titin.individuallyFit==0)
        loops=1;
    end


    for idxLoop=1:1:loops

        idxTrial=nan;
        if(fittingConfig.titin.individuallyFit==1)
            idxTrial = simConfig.trials(1,idxLoop);
            simConfigTmp.trials = idxTrial;
        end


        f1HNPreload=0;
        f1HNPreloadDelta=0.25;

        optParams.name  = 'f1HNPreload';
        optParams.value = f1HNPreload;

        simConfigTmp.flag_debugFitting=0;
        fittingFraction=1;
        npts=100;
        [optError,optErrorValues,figDebugFittingF1HNPreload,...
            ratFibrilModelsFittedUpd,benchRecord] =...
            calcErrorTRSS2017RampFraction(optParams,...
                fittingFraction,npts,...
                ratFibrilModelsFitted, expTRSS2017,simConfigTmp,...
                figDebugFittingF1HNPreload,plotConfig.subPlotPanel,...
                plotConfig.lineColors.simTitinK);

        optErrorBest=optError;
        optErrorValuesBest=optErrorValues;
        f1HNPreloadBest=f1HNPreload;
        dirMap = [1,-1];

        if(fittingConfig.titin.individuallyFit==1)
            fprintf('%1.2e\tfitting trial %i: f1HNPreload rmse (start)\n',...
                    optErrorBest,idxTrial);
            fprintf('%e\tf1HNPreload (start)\n',f1HNPreloadBest);    
            fprintf(fidFitting,...
                    '%1.2e\tfitting trial %i: f1HNPreload rmse (start)\n',...
                    optErrorBest,idxTrial);
            fprintf(fidFitting,...
                    '%e\tf1HNPreload (start)\n',f1HNPreloadBest);    
        else
            fprintf('%1.2e\tfitting all: f1HNPreload rmse (start)\n',...
                    optErrorBest);
            fprintf('%e\tf1HNPreload (start)\n',f1HNPreloadBest);    
            fprintf(fidFitting,...
                    '%1.2e\tfitting all: f1HNPreload rmse (start)\n',...
                    optErrorBest);
            fprintf(fidFitting,...
                    '%e\tf1HNPreload (start)\n',f1HNPreloadBest);    

        end
        for i=1:1:fittingConfig.numberOfBisections
            fprintf('%i/%i\n',i,fittingConfig.numberOfBisections);
        
            for j=1:1:length(dirMap)
                f1HNPreload = f1HNPreloadBest + f1HNPreloadDelta*dirMap(1,j);
                optParams.value=f1HNPreload;
                [optError,optErrorValues,figDebugFittingQ,...
                    ratFibrilModelsFittedUpd,benchRecord] =...
                    calcErrorTRSS2017RampFraction(optParams,...
                       fittingFraction, npts, ...
                       ratFibrilModelsFitted, expTRSS2017,simConfigTmp,...
                       figDebugFittingQ,plotConfig.subPlotPanel,...
                       plotConfig.lineColors.simTitinK);
                if(optError<optErrorBest)
                    optErrorBest=optError;
                    f1HNPreloadBest=f1HNPreload;
                    optErrorValuesBest=optErrorValues;                    
                    break;
                end
            end
            
            f1HNPreloadDelta=f1HNPreloadDelta*0.5;

        end

        if(fittingConfig.titin.individuallyFit==1)
            if(isnan(fitInfo.f1HNPreload.rmse))
                fitInfo.f1HNPreload.rmse = optErrorBest;
                fitInfo.f1HNPreload.x  = optErrorValuesBest.x(:,idxTrial);
                fitInfo.f1HNPreload.y     = optErrorValuesBest.y(:,idxTrial);
                fitInfo.f1HNPreload.yFit  = optErrorValuesBest.yFit(:,idxTrial);                    
                fitInfo.f1HNPreload.yErr  = optErrorValuesBest.yErr(:,idxTrial);
                fitInfo.f1HNPreload.arg  = QBest;            
                fitInfo.f1HNPreload.argDelta = QDelta*2;

            else
                fitInfo.f1HNPreload.rmse = [fitInfo.f1HNPreload.rmse, optErrorBest];
                fitInfo.f1HNPreload.x  = [fitInfo.f1HNPreload.x,  optErrorValuesBest.x(:,idxTrial)];
                fitInfo.f1HNPreload.y      = [fitInfo.f1HNPreload.y,      optErrorValuesBest.y(:,idxTrial)];
                fitInfo.f1HNPreload.yFit   = [fitInfo.f1HNPreload.yFit,   optErrorValuesBest.yFit(:,idxTrial)];
                fitInfo.f1HNPreload.yErr   = [fitInfo.f1HNPreload.yErr,   optErrorValuesBest.yErr(:,idxTrial)];
                fitInfo.f1HNPreload.arg  = [fitInfo.f1HNPreload.arg,  QBest];            
                fitInfo.f1HNPreload.argDelta = QDelta*2;
            end
        else
            fitInfo.f1HNPreload.rmse = optErrorBest;
            fitInfo.f1HNPreload.x  = optErrorValuesBest.x;
            fitInfo.f1HNPreload.y     = optErrorValuesBest.y;
            fitInfo.f1HNPreload.yFit  = optErrorValuesBest.yFit;
            fitInfo.f1HNPreload.yErr  = optErrorValuesBest.yErr;                
            fitInfo.f1HNPreload.arg  = QBest;           
            fitInfo.f1HNPreload.argDelta = QDelta*2;
        end

        
        fprintf('%1.2e\tfitting trial %i: f1HNPreload rmse (end)\n',...
            optErrorBest,idxTrial);
        fprintf('%e\tf1HNPreload (end)\n',f1HNPreloadBest);    

        fprintf(fidFitting,...
            '%1.2e\tfitting trial %i: f1HNPreload rmse (end)\n',...
            optErrorBest,idxTrial);
        fprintf(fidFitting,...
            '%e\tf1HNPreload (end)\n',f1HNPreloadBest);    

        %
        % Update the parameter struct
        % 
        optParams.value=f1HNPreloadBest;
        [optError,optErrorValues,figDebugFittingF1HNPreload,...
            ratFibrilModelsFittedUpd,benchRecord] =...
            calcErrorTRSS2017RampFraction(optParams,...
               fittingFraction, npts, ...
               ratFibrilModelsFitted, expTRSS2017,simConfigTmp,...
               figDebugFittingF1HNPreload,plotConfig.subPlotPanel,...
               plotConfig.lineColors.simTitinK);


        if(fittingConfig.titin.individuallyFit==1)
            ratFibrilModelsFitted(idxTrial)=ratFibrilModelsFittedUpd(idxTrial);
        else
            ratFibrilModelsFitted=ratFibrilModelsFittedUpd;
        end



        if(flagDebug==1)
            if(exist('figExpSlope2')==0)
                figExpSlope2 = figure;
            end
            if(fittingConfig.titin.individuallyFit==1)
                plot(expTRSS2017.activeLengtheningData(idxTrial).x,...
                     expTRSS2017.activeLengtheningData(idxTrial).y,'-k');
                hold on;
                plot(optParams.exp(idxTrial).x,optParams.exp(idxTrial).y,'xk');
                hold on;
                plot(optParams.exp(idxTrial).xLine,optParams.exp(idxTrial).yLine,'-b');
                hold on;
                plot(benchRecord.normFiberLength(:,idxTrial).*lopt,...
                     benchRecord.normFiberForce(:,idxTrial),'-r' );
            
            else
                for idxTrial=simConfig.trials
                    plot(expTRSS2017.activeLengtheningData(idxTrial).x,...
                         expTRSS2017.activeLengtheningData(idxTrial).y,'-k');
                    hold on;
                    plot(optParams.exp(idxTrial).x,optParams.exp(idxTrial).y,'xk');
                    hold on;
                    plot(optParams.exp(idxTrial).xLine,optParams.exp(idxTrial).yLine,'-b');
                    hold on;
                    plot(benchRecord.normFiberLength(:,idxTrial).*lopt,...
                         benchRecord.normFiberForce(:,idxTrial),'-r' );
                end

            end
            xlabel('X');
            ylabel('Y');
            title(sprintf('Trial %i',idxTrial));
        end

    end

    %
    % If we are fitting just one of the trials, then update the others
    % to have the same fitted Q value
    %
    if(length(fittingConfig.titin.trials)==1 ...
            && fittingConfig.titin.applyToAllTrials==1)
        for i=1:1:length(ratFibrilModelsFitted)
            if(i ~= fittingConfig.titin.trials(1,1))
                ratFibrilModelsFitted(i)=ratFibrilModelsFitted(fittingConfig.titin.trials);
            end
        end
    end
    
end

if(fittingConfig.fitl1HNOffset==1)
    flagDebug=1;
    simConfigTmp=simConfig;
    figDebugFittingL1HNOffset=figure;

    loops = length(fittingConfig.titin.trials);
    if(fittingConfig.titin.individuallyFit==0)
        loops=1;
    end


    for idxLoop=1:1:loops

        idxTrial=nan;
        if(fittingConfig.titin.individuallyFit==1)
            idxTrial = simConfig.trials(1,idxLoop);
            simConfigTmp.trials = idxTrial;
        end


        l1HNOffset=0;
        l1HNOffsetDelta=0.25;

        optParams.name  = 'l1HNOffset';
        optParams.value = l1HNOffset;

        simConfigTmp.flag_debugFitting=0;
        fittingFraction=1;
        npts=100;
        [optError,optErrorValues,figDebugFittingL1HNOffset,...
            ratFibrilModelsFittedUpd,benchRecord] =...
            calcErrorTRSS2017RampFraction(optParams,...
                fittingFraction,npts,...
                ratFibrilModelsFitted, expTRSS2017,simConfigTmp,...
                figDebugFittingL1HNOffset,plotConfig.subPlotPanel,...
                plotConfig.lineColors.simTitinK);

        optErrorBest        = optError;
        optErrorValuesBest  = optErrorValues;
        l1HNOffsetBest     = l1HNOffset;
        dirMap = [1,-1];


        if(fittingConfig.titin.individuallyFit==1)
            fprintf('%1.2e\tfitting trial %i: l1HNOffset rmse (start)\n',...
                    optErrorBest,idxTrial);
            fprintf('%e\tl1HNOffset (start)\n',l1HNOffsetBest);    
            fprintf(fidFitting,...
                    '%1.2e\tfitting trial %i: l1HNOffset rmse (start)\n',...
                    optErrorBest,idxTrial);
            fprintf(fidFitting,...
                    '%e\tl1HNOffset (start)\n',l1HNOffsetBest);    
        else
            fprintf('%1.2e\tfitting all: l1HNOffset rmse (start)\n',...
                    optErrorBest);
            fprintf('%e\tl1HNOffset (start)\n',l1HNOffsetBest);    
            fprintf(fidFitting,...
                    '%1.2e\tfitting all: l1HNOffset rmse (start)\n',...
                    optErrorBest);
            fprintf(fidFitting,...
                    '%e\tl1HNOffset (start)\n',l1HNOffsetBest);    

        end
        for i=1:1:fittingConfig.numberOfBisections
            fprintf('%i/%i\n',i,fittingConfig.numberOfBisections);
        
            for j=1:1:length(dirMap)
                l1HNOffset = l1HNOffsetBest + l1HNOffsetDelta*dirMap(1,j);
                optParams.value=l1HNOffset;
                [optError,optErrorValues,figDebugFittingL1HNOffset,...
                    ratFibrilModelsFittedUpd,benchRecord] =...
                    calcErrorTRSS2017RampFraction(optParams,...
                       fittingFraction, npts, ...
                       ratFibrilModelsFitted, expTRSS2017,simConfigTmp,...
                       figDebugFittingL1HNOffset,plotConfig.subPlotPanel,...
                       plotConfig.lineColors.simTitinK);
                if(optError<optErrorBest)
                    optErrorBest=optError;
                    l1HNOffsetBest=l1HNOffset;
                    optErrorValuesBest=optErrorValues;                    
                    break;
                end
            end
            
            l1HNOffsetDelta=l1HNOffsetDelta*0.5;

        end

        if(fittingConfig.titin.individuallyFit==1)
            if(isnan(fitInfo.l1HNOffset.rmse))
                fitInfo.l1HNOffset.rmse = optErrorBest;
                fitInfo.l1HNOffset.x  = optErrorValuesBest.x(:,idxTrial);
                fitInfo.l1HNOffset.y     = optErrorValuesBest.y(:,idxTrial);
                fitInfo.l1HNOffset.yFit  = optErrorValuesBest.yFit(:,idxTrial);                    
                fitInfo.l1HNOffset.yErr  = optErrorValuesBest.yErr(:,idxTrial);
                fitInfo.l1HNOffset.arg  = l1HNOffsetBest;            
                fitInfo.l1HNOffset.argDelta = l1HNOffsetDelta*2;

            else
                fitInfo.l1HNOffset.rmse = [fitInfo.l1HNOffset.rmse, optErrorBest];
                fitInfo.l1HNOffset.x  = [fitInfo.l1HNOffset.x,  optErrorValuesBest.x(:,idxTrial)];
                fitInfo.l1HNOffset.y      = [fitInfo.l1HNOffset.y,      optErrorValuesBest.y(:,idxTrial)];
                fitInfo.l1HNOffset.yFit   = [fitInfo.l1HNOffset.yFit,   optErrorValuesBest.yFit(:,idxTrial)];
                fitInfo.l1HNOffset.yErr   = [fitInfo.l1HNOffset.yErr,   optErrorValuesBest.yErr(:,idxTrial)];
                fitInfo.l1HNOffset.arg  = [fitInfo.l1HNOffset.arg,  l1HNOffsetBest];            
                fitInfo.l1HNOffset.argDelta = l1HNOffsetDelta*2;
            end
        else
            fitInfo.l1HNOffset.rmse = optErrorBest;
            fitInfo.l1HNOffset.x  = optErrorValuesBest.x;
            fitInfo.l1HNOffset.y     = optErrorValuesBest.y;
            fitInfo.l1HNOffset.yFit  = optErrorValuesBest.yFit;
            fitInfo.l1HNOffset.yErr  = optErrorValuesBest.yErr;                
            fitInfo.l1HNOffset.arg  = l1HNOffsetBest;           
            fitInfo.l1HNOffset.argDelta = l1HNOffsetDelta*2;
        end
        
        fprintf('%1.2e\tfitting trial %i: l1HNOffset rmse (end)\n',...
            optErrorBest,idxTrial);
        fprintf('%e\tl1HNOffset (end)\n',l1HNOffsetBest);    

        fprintf(fidFitting,...
            '%1.2e\tfitting trial %i: l1HNOffset rmse (end)\n',...
            optErrorBest,idxTrial);
        fprintf(fidFitting,...
            '%e\tl1HNOffset (end)\n',l1HNOffsetBest);    

        %
        % Update the parameter struct
        % 
        optParams.value=l1HNOffsetBest;
        [optError,optErrorValues,figDebugFittingL1HNOffset,...
            ratFibrilModelsFittedUpd,benchRecord] =...
            calcErrorTRSS2017RampFraction(optParams,...
               fittingFraction, npts, ...
               ratFibrilModelsFitted, expTRSS2017,simConfigTmp,...
               figDebugFittingL1HNOffset,plotConfig.subPlotPanel,...
               plotConfig.lineColors.simTitinK);


        if(fittingConfig.titin.individuallyFit==1)
            ratFibrilModelsFitted(idxTrial)=ratFibrilModelsFittedUpd(idxTrial);
        else
            ratFibrilModelsFitted=ratFibrilModelsFittedUpd;
        end

        if(flagDebug==1)
            if(exist('figExpSlope2')==0)
                figExpSlope2 = figure;
            end
            if(fittingConfig.titin.individuallyFit==1)
                plot(expTRSS2017.activeLengtheningData(idxTrial).x,...
                     expTRSS2017.activeLengtheningData(idxTrial).y,'-k');
                hold on;
                plot(optParams.exp(idxTrial).x,optParams.exp(idxTrial).y,'xk');
                hold on;
                plot(optParams.exp(idxTrial).xLine,optParams.exp(idxTrial).yLine,'-b');
                hold on;
                plot(benchRecord.normFiberLength(:,idxTrial).*lopt,...
                     benchRecord.normFiberForce(:,idxTrial),'-r' );
            
            else
                for idxTrial=simConfig.trials
                    plot(expTRSS2017.activeLengtheningData(idxTrial).x,...
                         expTRSS2017.activeLengtheningData(idxTrial).y,'-k');
                    hold on;
                    plot(optParams.exp(idxTrial).x,optParams.exp(idxTrial).y,'xk');
                    hold on;
                    plot(optParams.exp(idxTrial).xLine,optParams.exp(idxTrial).yLine,'-b');
                    hold on;
                    plot(benchRecord.normFiberLength(:,idxTrial).*lopt,...
                         benchRecord.normFiberForce(:,idxTrial),'-r' );
                end

            end
            xlabel('X');
            ylabel('Y');
            title(sprintf('Trial %i',idxTrial));
        end        

    end

    if(length(fittingConfig.titin.trials)==1 ...
            && fittingConfig.titin.applyToAllTrials==1)
        for i=1:1:length(ratFibrilModelsFitted)
            if(i ~= fittingConfig.titin.trials(1,1))
                ratFibrilModelsFitted(i)=ratFibrilModelsFitted(fittingConfig.titin.trials);
            end
        end
    end    
end

fclose(fidFitting);

%
% Simulate all of the trials with the updated parameters and save the
% results to file
%
optParams.name  = 'simulate';
optParams.value = nan;
optParams.expX  = [];
optParams.expY  = [];
optParams.expDYDX = nan;

simConfigTmp=simConfig;
simConfigTmp.trials =simConfig.trials;
simConfigTmp.flag_debugFitting=0;
fittingFraction=1;
npts=500;

figDebugFitting=figure;

[optError,optErrorValues,figDebugFitting,...
    ratFibrilModelsFittedUpd,benchRecordFitted] =...
    calcErrorTRSS2017RampFraction(optParams,fittingFraction,npts,...
               ratFibrilModelsFitted, expTRSS2017,simConfigTmp,...
               figDebugFitting,plotConfig.subPlotPanel,...
               plotConfig.lineColors.simTitinK);


