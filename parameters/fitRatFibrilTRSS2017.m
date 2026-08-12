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
                    expTRSS2017Raw,...
                    simConfig,...
                    fittingConfig,...
                    plotConfig,...
                    projectFolders)


ratFibrilModelsFitted=...
  ratFibrilModelsDefault;


benchRecordFitted=[];
benchRecord=[];


fitInfo = struct('fl',[],'fv',[],'timeConstant',[],'Kx',[],...
                    'QToF',[],'QToK',[],'f1HNPreload',[],'l1HNOffset',[],...
                    'FpeQToF',[]);

fitInfoFields= fields(fitInfo);

for i=1:1:length(fitInfoFields)
    fitInfo.(fitInfoFields{i}) = ...
        struct('rmse',nan,'nrmse',nan,...
               'x',[],'yErr',[],'y',[],'yFit',[],'yStd',[],...
               'arg',nan,'argDelta',nan);
end

fidFitting = fopen(fullfile(projectFolders.output_structs_TRSS2017,...
                ['fittingLog_',fittingConfig.trialStr,'.txt']),'w');


%%
% Fitting the active-force-length relation
%%
if(simConfig.flag_debugFitting==1)
  figDebugFitting = figure;
end
if(fittingConfig.fitFl==1)

    %%
    % fal fitting
    %%

    for i=1:1:length(ratFibrilModelsFitted)
      lceOptMdlErr = ...
        ratFibrilModelsFitted(1).musculotendon.optimalFiberLength ...
       -ratFibrilModelsFitted(i).musculotendon.optimalFiberLength;
      vmaxMdlErr = ...
        ratFibrilModelsFitted(1).musculotendon.maximumNormalizedFiberVelocity ...
       -ratFibrilModelsFitted(i).musculotendon.maximumNormalizedFiberVelocity;
      
      assert(abs(lceOptMdlErr)<1e-6 && abs(vmaxMdlErr) < 1e-6,...
             'Error: all models should be the same');
    end

    dfalN = 0.2;

    lceOptMdl   = ratFibrilModelsFitted(1).musculotendon.optimalFiberLength;
    vmax        = ratFibrilModelsFitted(1).musculotendon.maximumNormalizedFiberVelocity;

    expfl.lce = [];
    expfl.fN   = []; 
    
    mdlfl.lce = [];
    mdlfl.fN   = [];

    assert(ratFibrilModelsFitted(1).curves.useCalibratedCurves==1,...
           'Error: the calibrated curves should be used');

    curveFalN   = ratFibrilModelsFitted(1).curves.activeForceLengthCurve; 
    curveFvN    = ratFibrilModelsFitted(1).curves.fiberForceVelocityCurve; 

    for idxTrial = simConfig.trials
        %lce  = expTRSS2017.activeLengtheningData(idxTrial).x(1,1);
        %fN   = expTRSS2017.activeLengtheningData(idxTrial).y(1,1);
        
        seriesName = expTRSS2017Raw.trials{idxTrial};        
        lceN = expTRSS2017Raw.(seriesName).lN(1,1);  
        lce  = lceN*lceOptMdl;
        fceN = expTRSS2017Raw.(seriesName).fNavg(1,1);

        expfl.lce  = [expfl.lce;lce];
        expfl.fN    = [expfl.fN;fceN];   
        
        fNMdl = calcBezierYFcnXDerivative(lce./lceOptMdl,curveFalN,0);
        mdlfl.lce  = [mdlfl.lce;lce];
        mdlfl.fN    = [mdlfl.fN;fNMdl];                
    end


    argBest  = 0.2;
    flag_compensateForCrossbridgeStiffness = 0;
    [errFlN,flNFit,curveL] = ...
      calcErrorTRSS2017ForceLengthRelationAscendingLimb(argBest,...
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
            fitInfo.fl.nrmse= nan;
            fitInfo.fl.yStd = nan;
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
            fitInfo.fl.nrmse= nan;
            fitInfo.fl.yStd = nan;                       
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
                 num2str(idxTrial));
            hold on;
        end
        

        xlabel('Length ($$\mu m$$)');
        ylabel('Norm. Force ($$f / f_o^M$$)');
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
    %lceOptData  = min(expTRSS2017.activeLengtheningData(3).x);

    mdlfv.lceN = [];
    mdlfv.vceN = [];  
    mdlfv.fN   = [];
    mdlfv.fvN   = [];

    expfv.lceN = [];
    expfv.vceN = [];  
    expfv.fN   = [];
    expfv.fvN   = [];

    for idxTrial = simConfig.trials
%         idxKey_ = expTRSS2017.activeLengtheningData(idxTrial).keyIndices(...
%                     fittingConfig.idxFvKey); 
%         lce_  = expTRSS2017.activeLengtheningData(idxTrial).x(idxKey_,1);
%         lceN_ = lce_/lceOptMdl;
%         fN0_  = expTRSS2017.activeLengtheningData(idxTrial).y(1,1);
%         fN1_  = expTRSS2017.activeLengtheningData(idxTrial).y(idxKey_,1);

        seriesName = expTRSS2017Raw.trials{idxTrial};            
        idxKey = expTRSS2017Raw.keyIndices.(seriesName)(fittingConfig.idxFvKey); 
        lceN  = expTRSS2017Raw.(seriesName).lN(idxKey,1);
        lce   = lceN*lceOptMdl;
        fN0   = expTRSS2017Raw.(seriesName).fNavg(1,1);
        fN1   = expTRSS2017Raw.(seriesName).fNavg(idxKey,1);



        fN1 = (fN1-fN0)*fittingConfig.fv.scaleEnhancement + fN0;

        fvN = fN1/fN0;

        lpfFreqHz = simConfig.numericalDiffLowPassFilterFrequencyHz;
        [b,a]     = butter(2,lpfFreqHz/(0.5*expTRSS2017Raw.sampleFrequencyHz),'low');
        vceNnum  = calcCentralDifferenceDataSeries(...
                      expTRSS2017Raw.(seriesName).time,...
                      filtfilt(b,a,expTRSS2017Raw.(seriesName).lN));
        
        vceNavg  = mean(vceNnum);
        vceNN    = vceNavg/vmax;
        
        expfv.lceN  = [expfv.lceN;lceN];
        expfv.vceN  = [expfv.vceN;vceNN];
        expfv.fN    = [expfv.fN;fN1];
        expfv.fvN   = [expfv.fvN;fvN];   
        
        fvNMdl = calcBezierYFcnXDerivative(0.11,curveFvN,0);
        falNMdl = calcBezierYFcnXDerivative(lceN,curveFalN,0);

        mdlfv.lceN  = [mdlfv.lceN;lceN];
        mdlfv.vceN  = [mdlfv.vceN;vceNN];
        mdlfv.fN    = [mdlfv.fN;fvNMdl*falNMdl];
        mdlfv.fvN   = [mdlfv.fvN;fvNMdl];                
    end

    arg       = 0;
    argBest   = arg;
    delta     = 0.4*(...
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
            fitInfo.fv.nrmse= nan;
            fitInfo.fv.yStd = nan;                        
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
                fitInfo.fv.nrmse= nan;
                fitInfo.fv.yStd = nan;            
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
    %lceOptData  = min(expTRSS2017.activeLengtheningData(3).x);

    if(simConfig.flag_debugFitting==1)

        for idxTrial=simConfig.trials
            if(simConfig.flag_debugFitting==1)
                figure(figDebugFitting);
                subplot('Position',reshape(plotConfig.subPlotPanel(2,1,:),1,4));
                for idxTrial = simConfig.trials
                  seriesName = expTRSS2017Raw.trials{idxTrial};

                  plot(expTRSS2017Raw.(seriesName).lN.*lceOptMdl,...
                       expTRSS2017Raw.(seriesName).fNavg,...
                       '-','Color',plotConfig.lineColors.exp(idxTrial,:),...
                       'DisplayName',num2str(idxTrial));
                    hold on;

%                     plot(expTRSS2017.activeLengtheningData(idxTrial).x,...
%                          expTRSS2017.activeLengtheningData(idxTrial).y,...
%                          '-','Color',plotConfig.lineColors.exp(idxTrial,:),...
%                          'DisplayName',...
%                          expTRSS2017.activeLengtheningData(idxTrial).seriesName);
%                     hold on;
                end
                xlabel('Length ($$\mu m$$)');
                ylabel('Norm. Force ($$f / f_o^M$$)');
                title('Rat fibril (EDL) active lengthening');  
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
                   ratFibrilModelsFitted, ...
                   expTRSS2017Raw,...
                   simConfig,...
                   figDebugFitting,...
                   plotConfig.subPlotPanel,...
                   plotConfig.lineColors.simXE,...
                   benchRecord,...
                   fittingConfig.fittingEvaluationMethod,...
                   projectFolders);

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
                   ratFibrilModelsFitted, ...
                   expTRSS2017Raw,...
                   simConfig,...
                   figDebugFitting,...
                   plotConfig.subPlotPanel,...
                   plotConfig.lineColors.simXE,...
                   benchRecord,...
                   fittingConfig.fittingEvaluationMethod,...
                   projectFolders);

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
                       fittingFraction, ...
                       npts, ...
                       ratFibrilModelsFitted,...
                       expTRSS2017Raw,...
                       simConfig,...
                       figDebugFitting,...
                       plotConfig.subPlotPanel,...
                       plotConfig.lineColors.simXE,...
                       benchRecord,...
                       fittingConfig.fittingEvaluationMethod,...
                       projectFolders);

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

    %
    % Update the bench record if an approximation method was being used
    %
    if(contains(fittingConfig.fittingEvaluationMethod,'approximation'))
        optParams.value=bestValue;
        [errorVal,errorValues,figDebugFitting,...
            ratFibrilModelsFitted,benchRecordFitted] ...
            = calcErrorTRSS2017RampFraction(optParams,...
                   fittingFraction, ...
                   npts, ...
                   ratFibrilModelsFitted,...
                   expTRSS2017Raw,...
                   simConfig,...
                   figDebugFitting,...
                   plotConfig.subPlotPanel,...
                   plotConfig.lineColors.simXE,...
                   benchRecordFitted,...
                   'simulateFullModel',...
                   projectFolders);  
          bestError=errorVal;        
          bestErrorValues = errorValues;                        
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

    fitInfo.timeConstant.rmse     = bestError;
    fitInfo.timeConstant.nrmse= nan;
    fitInfo.timeConstant.yStd = nan;            
    fitInfo.timeConstant.x        = bestErrorValues.x;
    fitInfo.timeConstant.y        = bestErrorValues.y;
    fitInfo.timeConstant.yFit     = bestErrorValues.yFit;
    fitInfo.timeConstant.yErr     = bestErrorValues.yErr;
    fitInfo.timeConstant.arg      = bestValue;
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
                   fittingFraction,...
                   npts, ...
                   ratFibrilModelsFitted, ...
                   expTRSS2017Raw,...
                   simConfig,...
                   figDebugFitting,...
                   plotConfig.subPlotPanel,...
                   plotConfig.lineColors.simXE,...
                   benchRecord,...
                   fittingConfig.fittingEvaluationMethod,...
                   projectFolders);

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
                   fittingFraction,...
                   npts, ...
                   ratFibrilModelsFitted, ...
                   expTRSS2017Raw,...
                   simConfig,...
                   figDebugFitting,...
                   plotConfig.subPlotPanel,...
                   plotConfig.lineColors.simXE,...
                   benchRecord,...
                   fittingConfig.fittingEvaluationMethod,...
                   projectFolders);

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
                       fittingFraction, ...
                       npts, ...
                       ratFibrilModelsFitted, ...
                       expTRSS2017Raw,...
                       simConfig,...
                       figDebugFitting,...
                       plotConfig.subPlotPanel,...
                       plotConfig.lineColors.simXE,...
                       benchRecord,...
                       fittingConfig.fittingEvaluationMethod,...
                       projectFolders);

            if(errorVal < bestError)
                bestError=errorVal;      
                bestValue = optParams.value;
                bestErrorValues = errorValues;
                ratFibrilModelsFitted=ratFibrilModelsUpd;
            end
        end
        deltaValue=deltaValue*0.5;
    end

    %
    % Update the bench record if an approximation method was being used
    %
    if(contains(fittingConfig.fittingEvaluationMethod,'approximation'))
        optParams.value=bestValue;
        [errorVal,errorValues,figDebugFitting,...
            ratFibrilModelsUpd,benchRecordFitted] ...
            = calcErrorTRSS2017RampFraction(optParams,...
                   fittingFraction, ...
                   npts, ...
                   ratFibrilModelsFitted,...
                   expTRSS2017Raw,...
                   simConfig,...
                   figDebugFitting,...
                   plotConfig.subPlotPanel,...
                   plotConfig.lineColors.simXE,...
                   benchRecordFitted,...
                   'simulateFullModel',...
                   projectFolders);  
          bestError=errorVal;        
          bestErrorValues = errorValues;
          bestValue = optParams.value;               
          ratFibrilModelsFitted=ratFibrilModelsUpd;           
    end    

    fitInfo.Kx.rmse     = bestError;
    fitInfo.Kx.nrmse= nan;
    fitInfo.Kx.yStd = nan;              
    fitInfo.Kx.x        = bestErrorValues.x;
    fitInfo.Kx.y        = bestErrorValues.y;
    fitInfo.Kx.yFit     = bestErrorValues.yFit;
    fitInfo.Kx.yErr     = bestErrorValues.yErr;
    fitInfo.Kx.arg      = bestValue;
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
    %idx1 = expTRSS2017Raw.(seriesName).lN.*lceOptMdl;
    %idx1 = length(expTRSS2017.activeLengtheningData(idxTrial).x);

%     x0_ = expTRSS2017.activeLengtheningData(idxTrial).x(1,1);
%     x1_ = expTRSS2017.activeLengtheningData(idxTrial).x(end);
%     xStart_= x0_ + 0.125*(x1_-x0_);
%     expXTmp_ = [0:0.1:1]'.*(x1_-xStart_) + xStart_;
% 
%     [expLceU_,iq_] = unique(expTRSS2017.activeLengtheningData(idxTrial).x);
%             
%     expYTmp_ = interp1(expTRSS2017.activeLengtheningData(idxTrial).x(iq_),...
%                    expTRSS2017.activeLengtheningData(idxTrial).y(iq_),...
%                    expXTmp_);

    seriesName = expTRSS2017Raw.trials{idxTrial};

    x0 = expTRSS2017Raw.(seriesName).lN(1,1)*lceOptMdl;
    x1 = expTRSS2017Raw.(seriesName).lN(end)*lceOptMdl;
    xStart= x0 + 0.125*(x1-x0);
    expXTmp = [0:0.1:1]'.*(x1-xStart) + xStart;

    [expLceU,iq] = unique(expTRSS2017Raw.(seriesName).lN.*lceOptMdl);
            
    expYTmp = interp1(expTRSS2017Raw.(seriesName).lN(iq).*lceOptMdl,...
                   expTRSS2017Raw.(seriesName).fNavg(iq),...
                   expXTmp);

    expYStdTmp = interp1(expTRSS2017Raw.(seriesName).lN(iq).*lceOptMdl,...
                   expTRSS2017Raw.(seriesName).fNstd(iq),...
                   expXTmp);

    
    A = [expXTmp, ones(size(expXTmp))];
    b = expYTmp;
    p = (A'*A)\(A'*b);
    expSlope = p(1,1);
    
    xLine = expXTmp;
    yLine = [xLine, ones(size(xLine))]*p;


    
    %
    % -std
    %

    lopt=ratFibrilModelsFitted(idxTrial).musculotendon.optimalFiberLength;
    fiso=ratFibrilModelsFitted(idxTrial).musculotendon.fiso;

    optParams.exp(idxTrial).x = expXTmp;
    optParams.exp(idxTrial).y = expYTmp;
    optParams.exp(idxTrial).yStd = expYStdTmp;    
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
    flagDebug=0;
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
        [optErrorStart,optErrorValues,figDebugFittingQ,...
            ratFibrilModelsFittedUpd,benchRecord] =...
            calcErrorTRSS2017RampFraction(optParams,...
                       fittingFraction, ...
                       npts, ...
                       ratFibrilModelsFitted, ...
                       expTRSS2017Raw,...
                       simConfigTmp,...
                       figDebugFittingQ,...
                       plotConfig.subPlotPanel,...
                       plotConfig.lineColors.simTitinK,...
                       benchRecord,...
                       fittingConfig.fittingEvaluationMethod,...
                       projectFolders);

        optErrorBest=optErrorStart;
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
                       fittingFraction, ...
                       npts, ...
                       ratFibrilModelsFitted, ...
                       expTRSS2017Raw,...
                       simConfigTmp,...
                       figDebugFittingQ,...
                       plotConfig.subPlotPanel,...
                       plotConfig.lineColors.simTitinK,...
                       benchRecord,...
                       fittingConfig.fittingEvaluationMethod,...
                       projectFolders);

                if(optError<optErrorBest)
                    optErrorBest=optError;
                    optErrorValuesBest=optErrorValues;
                    QBest=Q;
                    break;
                end
            end
            
            QDelta=QDelta*0.5;

        end


        fitType='';
        if(fittingConfig.fitQToF)
            fitType = 'QToF';
        end
        if(fittingConfig.fitQToK)
            fitType = 'QToK';
        end

        %
        % Update the parameter struct
        % 
        optParams.value=QBest;
        [optErrorBest,optErrorValuesBest, figDebugFittingQ,...
            ratFibrilModelsFittedUpd,benchRecordFitted] =...
            calcErrorTRSS2017RampFraction(optParams,...
                       fittingFraction, ...
                       npts, ...
                       ratFibrilModelsFitted, ...
                       expTRSS2017Raw,...
                       simConfigTmp,...
                       figDebugFittingQ,...
                       plotConfig.subPlotPanel,...
                       plotConfig.lineColors.simTitinK,...
                       benchRecordFitted,...
                       'simulateFullModel',...
                       projectFolders);

        if(fittingConfig.titin.individuallyFit==1)
            fprintf('fitting trial %i: Q nrmse force profile\n',idxTrial);
            fprintf('\t%e\tQ (start)\n',QInit); 
            fprintf('\t%e\tNRMSE (start)\n',optErrorStart); 
            fprintf('\t%e\tQ (end)\n',QBest); 
            fprintf('\t%e\tNRMSE (end)\n',optErrorBest);
        else
            fprintf('fitting all: Q rmse force profile\n');
            fprintf('\t%e\tQ (start)\n',QInit); 
            fprintf('\t%e\tNRMSE (start)\n',optErrorStart); 
            fprintf('\t%e\tQ (end)\n',QBest); 
            fprintf('\t%e\tNRMSE (end)\n',optErrorBest); 
        end            
        if(fittingConfig.titin.individuallyFit==1)
            if(isnan(fitInfo.(fitType).rmse))
                fitInfo.(fitType).x     = optErrorValuesBest.x(:,idxTrial);
                fitInfo.(fitType).y     = optErrorValuesBest.y(:,idxTrial);
                fitInfo.(fitType).yFit  = optErrorValuesBest.yFit(:,idxTrial);                    
                fitInfo.(fitType).yErr  = optErrorValuesBest.yErr(:,idxTrial);
                fitInfo.(fitType).rmse  = optErrorValuesBest.rmse(:,idxTrial);
                fitInfo.(fitType).yStd  = optErrorValuesBest.yStd(:,idxTrial);
                fitInfo.(fitType).yNErr = optErrorValuesBest.yNErr(:,idxTrial);
                fitInfo.(fitType).nrmse = optErrorValuesBest.nrmse(:,idxTrial);

                fitInfo.(fitType).arg  = QBest;            
                fitInfo.(fitType).argDelta = QDelta*2;

            else
                fitInfo.(fitType).x     = [fitInfo.(fitType).x,     optErrorValuesBest.x(:,idxTrial)];
                fitInfo.(fitType).y     = [fitInfo.(fitType).y,     optErrorValuesBest.y(:,idxTrial)];
                fitInfo.(fitType).yFit  = [fitInfo.(fitType).yFit,  optErrorValuesBest.yFit(:,idxTrial)];
                fitInfo.(fitType).yErr  = [fitInfo.(fitType).yErr,  optErrorValuesBest.yErr(:,idxTrial)];
                fitInfo.(fitType).rmse  = [fitInfo.(fitType).rmse,  optErrorValuesBest.rmse(:,idxTrial)];
                fitInfo.(fitType).yStd  = [fitInfo.(fitType).yStd,  optErrorValuesBest.yStd(:,idxTrial)];
                fitInfo.(fitType).yNErr = [fitInfo.(fitType).yNErr, optErrorValuesBest.yNErr(:,idxTrial)];
                fitInfo.(fitType).nrmse = [fitInfo.(fitType).nrmse, optErrorValuesBest.nrmse(:,idxTrial)];
                
                fitInfo.(fitType).arg  = [fitInfo.(fitType).arg,  QBest];            
                fitInfo.(fitType).argDelta = QDelta*2;
            end
        else
            fitInfo.(fitType).x     = optErrorValuesBest.x;
            fitInfo.(fitType).y     = optErrorValuesBest.y;
            fitInfo.(fitType).yFit  = optErrorValuesBest.yFit;
            fitInfo.(fitType).yErr  = optErrorValuesBest.yErr;
            fitInfo.(fitType).rmse  = optErrorValuesBest.rmse;
            fitInfo.(fitType).yStd  = optErrorValuesBest.yStd;
            fitInfo.(fitType).yNErr = optErrorValuesBest.yNErr;
            fitInfo.(fitType).nrmse = optErrorValuesBest.nrmse;

            fitInfo.(fitType).arg  = QBest;           
            fitInfo.(fitType).argDelta = QDelta*2;
        end

                      
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
                %plot(expTRSS2017.activeLengtheningData(idxTrial).x,...
                %     expTRSS2017.activeLengtheningData(idxTrial).y,'-k');
                %hold on;
                seriesName = expTRSS2017Raw.trials{idxTrial};
                plot(expTRSS2017Raw.(seriesName).lN.*lceOptMdl,...
                     expTRSS2017Raw.(seriesName).fNavg,'-k');
                hold on;
                
                plot(optParams.exp(idxTrial).x,optParams.exp(idxTrial).y,'xk');
                hold on;
                plot(optParams.exp(idxTrial).xLine,optParams.exp(idxTrial).yLine,'-b');
                hold on;
                plot(benchRecord.normFiberLength(:,idxTrial).*lopt,...
                     benchRecord.normFiberForce(:,idxTrial),'-r' );
            
            else
                for idxTrial=simConfig.trials
                    %plot(expTRSS2017.activeLengtheningData(idxTrial).x,...
                    %     expTRSS2017.activeLengtheningData(idxTrial).y,'-k');
                    %hold on;
                    seriesName = expTRSS2017Raw.trials{idxTrial};
                    plot(expTRSS2017Raw.(seriesName).lN.*lceOptMdl,...
                         expTRSS2017Raw.(seriesName).fNavg,'-k');
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
end

%%
%
% Multi-variate fitting
%
%%
if(fittingConfig.fitFpeQToF ==1)
    %
    % Solve for Q to best fit the average slope of the force development
    %
    flagDebug           = 0;
    simConfigTmp        = simConfig;
    figDebugFittingFpeQ = figure;
    
    loops = length(fittingConfig.titin.trials);
    if(fittingConfig.titin.individuallyFit==0)
        loops=1;
    end
    
    QInit = 0.5;
    QDeltaInit = 0.5*QInit;

    kToeMin = 1.40;    
    kToeMax = ratFibrilModelsFitted(idxTrial).curves.forceLengthCurveSettings.kToe;
    kToeInit = 0.5*(kToeMin+kToeMax);
    kToeDeltaInit = 0.25*(kToeMax-kToeMin);

    normLengthToeMax = ...
      ratFibrilModelsFitted(idxTrial).curves.forceLengthCurveSettings.normLengthToe;
    normLengthToeMin = normLengthToeMax-0.2;
    normLengthToeInit = 0.5*(normLengthToeMax+normLengthToeMin);
    normLengthToeDeltaInit = 0.25*(normLengthToeMax-normLengthToeMin);

    curvinessMax =1;
    curvinessMin =0;
    curvinessInit = 0.5;
    curvinessDeltaInit=0.25*(curvinessMax-curvinessMin);

    paramInit = [QInit;kToeInit;normLengthToeInit;curvinessInit];
    paramInitDelta = [QDeltaInit;...
                      kToeDeltaInit;...
                      normLengthToeDeltaInit;...
                      curvinessDeltaInit];

    numParams=length(paramInit);

    paramPattern = zeros(numParams^2-1,numParams);

    for i=1:1:size(paramPattern,1)
      rowVal = dec2bin(i);
      k = length(rowVal);
      for j=size(paramPattern,2):-1:(size(paramPattern,2)-length(rowVal)+1)
        paramPattern(i,j) = str2double(rowVal(k));
        k=k-1;
      end
    end


    for idxLoop = 1:1:loops

        idxTrial = nan;
        if(fittingConfig.titin.individuallyFit==1)
            idxTrial = simConfig.trials(1,idxLoop);
            simConfigTmp.trials = idxTrial;
        end
        optParams.name = 'FpeQToF';        

        optParams.value=  paramInit;
        
        simConfigTmp.flag_debugFitting=0;
        fittingFraction=1;
        npts=100;
        [optErrorStart,optErrorValues,figDebugFittingFpeQ,...
            ratFibrilModelsFittedUpd,benchRecord] =...
            calcErrorTRSS2017RampFraction(optParams,...
                       fittingFraction, ...
                       npts, ...
                       ratFibrilModelsFitted, ...
                       expTRSS2017Raw,...
                       simConfigTmp,...
                       figDebugFittingFpeQ,...
                       plotConfig.subPlotPanel,...
                       plotConfig.lineColors.simTitinK,...
                       benchRecord,...
                       fittingConfig.fittingEvaluationMethod,...
                       projectFolders);

        optErrorBest=optErrorStart;
        optErrorValuesBest=optErrorValues;
        paramBest=paramInit;
        paramDelta=paramInitDelta;
        dirMap = [1,-1];



        for i=1:1:fittingConfig.numberOfBisections
          fprintf('%i/%i\n',i,fittingConfig.numberOfBisections);
    
          for j=1:1:size(paramPattern,1)

            for k=1:1:length(dirMap)

                paramTest = paramBest ...
                  + (paramPattern(j,:)').*(paramDelta).*dirMap(1,k);
                %disp(paramTest);
                optParams.value = paramTest;
                [optError,optErrorValues,figDebugFittingFpeQ,...
                    ratFibrilModelsFittedUpd,benchRecord] =...
                    calcErrorTRSS2017RampFraction(optParams,...
                       fittingFraction, ...
                       npts, ...
                       ratFibrilModelsFitted, ...
                       expTRSS2017Raw,...
                       simConfigTmp,...
                       figDebugFittingFpeQ,...
                       plotConfig.subPlotPanel,...
                       plotConfig.lineColors.simTitinK,...
                       benchRecord,...
                       fittingConfig.fittingEvaluationMethod,...
                       projectFolders);

                if(optError<optErrorBest)
                    optErrorBest=optError;
                    optErrorValuesBest=optErrorValues;
                    paramIterationBest=optParams.value;
                end
            end

          end
          paramBest=paramIterationBest;
          paramDelta=paramDelta*0.5;

        end


        fitType='FpeQToF';

        %
        % Update the parameter struct
        % 
        optParams.value=paramBest;
        [optErrorBest,optErrorValuesBest, figDebugFittingFpeQ,...
            ratFibrilModelsFittedUpd,benchRecordFitted] =...
            calcErrorTRSS2017RampFraction(optParams,...
                       fittingFraction, ...
                       npts, ...
                       ratFibrilModelsFitted, ...
                       expTRSS2017Raw,...
                       simConfigTmp,...
                       figDebugFittingFpeQ,...
                       plotConfig.subPlotPanel,...
                       plotConfig.lineColors.simTitinK,...
                       benchRecordFitted,...
                       'simulateFullModel',...
                       projectFolders);

        if(fittingConfig.titin.individuallyFit==1)
            fprintf('fitting trial %i: Fpe+Q nrmse force profile\n',idxTrial);
            for k=1:1:length(paramInit)
              fprintf('\t%e',paramInit(k));
            end
            fprintf('\tparam (start)\n');
            fprintf('\t%e\tNRMSE (start)\n',optErrorStart); 
            for k=1:1:length(paramBest)
              fprintf('\t%e',paramBest(k));
            end
            fprintf('\tparam (end)\n');            
            fprintf('\t%e\tNRMSE (end)\n',optErrorBest);
        else
            fprintf('fitting all: Fpe+Q nrmse force profile\n');
            for k=1:1:length(paramInit)
              fprintf('\t%e',paramInit(k));
            end
            fprintf('\tparam (start)\n');
            fprintf('\t%e\tNRMSE (start)\n',optErrorStart); 
            for k=1:1:length(paramBest)
              fprintf('\t%e',paramBest(k));
            end
            fprintf('\tparam (end)\n');            
            fprintf('\t%e\tNRMSE (end)\n',optErrorBest); 
        end            
        if(fittingConfig.titin.individuallyFit==1)
            if(isnan(fitInfo.(fitType).rmse))
                fitInfo.(fitType).x     = optErrorValuesBest.x(:,idxTrial);
                fitInfo.(fitType).y     = optErrorValuesBest.y(:,idxTrial);
                fitInfo.(fitType).yFit  = optErrorValuesBest.yFit(:,idxTrial);                    
                fitInfo.(fitType).yErr  = optErrorValuesBest.yErr(:,idxTrial);
                fitInfo.(fitType).rmse  = optErrorValuesBest.rmse(:,idxTrial);
                fitInfo.(fitType).yStd  = optErrorValuesBest.yStd(:,idxTrial);
                fitInfo.(fitType).yNErr = optErrorValuesBest.yNErr(:,idxTrial);
                fitInfo.(fitType).nrmse = optErrorValuesBest.nrmse(:,idxTrial);

                fitInfo.(fitType).arg  = paramBest';            
                fitInfo.(fitType).argDelta = paramDelta*2;

            else
                fitInfo.(fitType).x     = [fitInfo.(fitType).x,     optErrorValuesBest.x(:,idxTrial)];
                fitInfo.(fitType).y     = [fitInfo.(fitType).y,     optErrorValuesBest.y(:,idxTrial)];
                fitInfo.(fitType).yFit  = [fitInfo.(fitType).yFit,  optErrorValuesBest.yFit(:,idxTrial)];
                fitInfo.(fitType).yErr  = [fitInfo.(fitType).yErr,  optErrorValuesBest.yErr(:,idxTrial)];
                fitInfo.(fitType).rmse  = [fitInfo.(fitType).rmse,  optErrorValuesBest.rmse(:,idxTrial)];
                fitInfo.(fitType).yStd  = [fitInfo.(fitType).yStd,  optErrorValuesBest.yStd(:,idxTrial)];
                fitInfo.(fitType).yNErr = [fitInfo.(fitType).yNErr, optErrorValuesBest.yNErr(:,idxTrial)];
                fitInfo.(fitType).nrmse = [fitInfo.(fitType).nrmse, optErrorValuesBest.nrmse(:,idxTrial)];
                
                fitInfo.(fitType).arg  = [fitInfo.(fitType).arg,  paramBest'];            
                fitInfo.(fitType).argDelta = paramDelta*2;
            end
        else
            fitInfo.(fitType).x     = optErrorValuesBest.x;
            fitInfo.(fitType).y     = optErrorValuesBest.y;
            fitInfo.(fitType).yFit  = optErrorValuesBest.yFit;
            fitInfo.(fitType).yErr  = optErrorValuesBest.yErr;
            fitInfo.(fitType).rmse  = optErrorValuesBest.rmse;
            fitInfo.(fitType).yStd  = optErrorValuesBest.yStd;
            fitInfo.(fitType).yNErr = optErrorValuesBest.yNErr;
            fitInfo.(fitType).nrmse = optErrorValuesBest.nrmse;

            fitInfo.(fitType).arg  = paramBest';           
            fitInfo.(fitType).argDelta = paramDelta*2;
        end

                      
        if(fittingConfig.titin.individuallyFit==1)
            ratFibrilModelsFitted(idxTrial)=ratFibrilModelsFittedUpd(idxTrial);
        else
            ratFibrilModelsFitted=ratFibrilModelsFittedUpd;
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
                ratFibrilModelsFitted(i)...
                  =ratFibrilModelsFitted(fittingConfig.titin.trials);
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
                fittingFraction,...
                npts,...
                ratFibrilModelsFitted, ...
                expTRSS2017Raw,...
                simConfigTmp,...
                figDebugFittingF1HNPreload,...
                plotConfig.subPlotPanel,...
                plotConfig.lineColors.simTitinK,...
                fittingConfig.fittingEvaluationMethod,...
                benchRecord,...
                projectFolders);

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
                       fittingFraction, ...
                       npts, ...
                       ratFibrilModelsFitted, ...
                       expTRSS2017Raw,...
                       simConfigTmp,...
                       figDebugFittingQ,...
                       plotConfig.subPlotPanel,...
                       plotConfig.lineColors.simTitinK,...
                       benchRecord,...
                       fittingConfig.fittingEvaluationMethod,...
                       projectFolders);

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
            ratFibrilModelsFittedUpd,benchRecordFitted] =...
            calcErrorTRSS2017RampFraction(optParams,...
               fittingFraction, ...
               npts, ...
               ratFibrilModelsFitted, ...
               expTRSS2017Raw,...
               simConfigTmp,...
               figDebugFittingF1HNPreload,...
               plotConfig.subPlotPanel,...
               plotConfig.lineColors.simTitinK,...
               benchRecordFitted,...
               'simulateFullModel',...
               projectFolders);
       

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
                %plot(expTRSS2017.activeLengtheningData(idxTrial).x,...
                %     expTRSS2017.activeLengtheningData(idxTrial).y,'-k');
                %hold on;              
                seriesName = expTRSS2017Raw.trials{idxTrial};
                plot(expTRSS2017Raw.(seriesName).lN.*lceOptMdl,...
                     expTRSS2017Raw.(seriesName).fNavg,'-k');
                hold on;                              
                plot(optParams.exp(idxTrial).x,optParams.exp(idxTrial).y,'xk');
                hold on;
                plot(optParams.exp(idxTrial).xLine,optParams.exp(idxTrial).yLine,'-b');
                hold on;
                plot(benchRecord.normFiberLength(:,idxTrial).*lopt,...
                     benchRecord.normFiberForce(:,idxTrial),'-r' );
            
            else
                for idxTrial=simConfig.trials
                    %plot(expTRSS2017.activeLengtheningData(idxTrial).x,...
                    %     expTRSS2017.activeLengtheningData(idxTrial).y,'-k');
                    %hold on;
                    seriesName = expTRSS2017Raw.trials{idxTrial};
                    plot(expTRSS2017Raw.(seriesName).lN.*lceOptMdl,...
                         expTRSS2017Raw.(seriesName).fNavg,'-k');
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
                fittingFraction,...
                npts,...
                ratFibrilModelsFitted, ...
                expTRSS2017Raw,...
                simConfigTmp,...
                figDebugFittingL1HNOffset,...
                plotConfig.subPlotPanel,...
                plotConfig.lineColors.simTitinK,...
                benchRecord,...
                fittingConfig.fittingEvaluationMethod,...
                projectFolders);

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
                       fittingFraction, ...
                       npts, ...
                       ratFibrilModelsFitted, ...
                       expTRSS2017Raw,...
                       simConfigTmp,...
                       figDebugFittingL1HNOffset,...
                       plotConfig.subPlotPanel,...
                       plotConfig.lineColors.simTitinK,...
                       benchRecord,...
                       fittingConfig.fittingEvaluationMethod,...
                       projectFolders);

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
            ratFibrilModelsFittedUpd,benchRecordFitted] =...
            calcErrorTRSS2017RampFraction(optParams,...
               fittingFraction, ...
               npts, ...
               ratFibrilModelsFitted, ...
               expTRSS2017Raw,...
               simConfigTmp,...
               figDebugFittingL1HNOffset,...
               plotConfig.subPlotPanel,...
               plotConfig.lineColors.simTitinK,...
               benchRecordFitted,...
               'simulateFullModel',...
               projectFolders);


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
                %plot(expTRSS2017.activeLengtheningData(idxTrial).x,...
                %     expTRSS2017.activeLengtheningData(idxTrial).y,'-k');
                %hold on;
                seriesName = expTRSS2017Raw.trials{idxTrial};
                plot(expTRSS2017Raw.(seriesName).lN.*lceOptMdl,...
                     expTRSS2017Raw.(seriesName).fNavg,'-k');
                hold on;
                plot(optParams.exp(idxTrial).x,optParams.exp(idxTrial).y,'xk');
                hold on;
                plot(optParams.exp(idxTrial).xLine,optParams.exp(idxTrial).yLine,'-b');
                hold on;
                plot(benchRecord.normFiberLength(:,idxTrial).*lopt,...
                     benchRecord.normFiberForce(:,idxTrial),'-r' );
            
            else
                for idxTrial=simConfig.trials
                    %plot(expTRSS2017.activeLengtheningData(idxTrial).x,...
                    %     expTRSS2017.activeLengtheningData(idxTrial).y,'-k');
                    %hold on;
                    seriesName = expTRSS2017Raw.trials{idxTrial};
                    plot(expTRSS2017Raw.(seriesName).lN.*lceOptMdl,...
                         expTRSS2017Raw.(seriesName).fNavg,'-k');
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




