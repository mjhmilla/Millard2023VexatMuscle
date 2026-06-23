function success = runParameterStudyRatFibrilTRSS2017_01(...
                             ratFibrilModelsFitted,...
                             expTRSS2017Raw,...
                             simConfig,...
                             fittingConfig,...
                             plotConfig)
success=1;
%%
%Plot the exp data
%%
figSandbox=figure;

lceOptMdl = ratFibrilModelsFitted(1).musculotendon.optimalFiberLength;
for idx=simConfig.trials     
    idxTrial = simConfig.trials(idx);
    hdlVis='off';
    if(idx==1)
       hdlVis = 'on';
    end
    figure(figSandbox);
    subplot(1,1,1);  

    seriesName = expTRSS2017Raw.trials{idxTrial};
    plot(expTRSS2017Raw.(seriesName).lN,...
         expTRSS2017Raw.(seriesName).fNavg,...
         '-','Color',plotConfig.lineColors.exp(idx,:),...
         'LineWidth',1,...
         'DisplayName','TRSS2017',...
         'HandleVisibility',hdlVis);
%     plot(expTRSS2017.activeLengtheningData(idx).x./lceOptMdl,...
%          expTRSS2017.activeLengtheningData(idx).y,...
%          '-','Color',plotConfig.lineColors.exp(idx,:),...
%          'LineWidth',1,...
%          'DisplayName','TRSS2017',...
%          'HandleVisibility',hdlVis);
    hold on; 
end
box off;
xlim([0.6,1.4]);
ylim([0,2.5]);
xlabel('Norm. Length');
ylabel('Norm. Force');
%%
% Configure the fitted models
%%
fittingConfig.fitTimeConstant=1;



%%
% Simulate a variety of parameter values
%%
idxTrial=2;
scalingCoeffs = [0.9,1,1.1];
fittingFraction=1;
npts=500;
numberOfSimulations=1;
for i=1:1:length(scalingCoeffs)

    %%
    % Set the parameters
    %%
    responseTimeScaling=scalingCoeffs(1,i);
    
    ratFibrilModelsFitted(idxTrial).sarcomere.slidingTimeConstantBlendingParameter = 0.01;
    
    if(isfield(ratFibrilModelsFitted(idxTrial).sarcomere,'slidingTimeConstantLengthening'))
        ratFibrilModelsFitted(idxTrial).sarcomere.slidingTimeConstantLengthening= ...
            ratFibrilModelsFitted(idxTrial).sarcomere.slidingTimeConstantLengthening...
            *responseTimeScaling;

    else
        ratFibrilModelsFitted(idxTrial).sarcomere.slidingTimeConstantLengthening= ...
            ratFibrilModelsFitted(idxTrial).sarcomere.slidingTimeConstant...
            *responseTimeScaling;

    end

    
    ratFibrilModelsFitted(idxTrial).sarcomere.slidingTimeConstantShortening= ...
        ratFibrilModelsFitted(idxTrial).sarcomere.slidingTimeConstant;
    
    ratFibrilModelsFitted(idxTrial).sarcomere.useVariableSlidingTimeConstant = 1;

    %%
    % Run the simulation
    %%
    lceOptMdl   = ratFibrilModelsFitted(idxTrial).musculotendon.optimalFiberLength;
    vmax        = ratFibrilModelsFitted(idxTrial).musculotendon.maximumNormalizedFiberVelocity;

    seriesName = expTRSS2017Raw.trials{idxTrial};

    rampLengthStart  = expTRSS2017Raw.(seriesName).lN(1,1).*lceOptMdl;
    rampLengthEnd    = expTRSS2017Raw.(seriesName).lN(end,1).*lceOptMdl; 
    
%     rampLengthStart  = expTRSS2017.activeLengtheningData(idxTrial).x(1,1);
%     rampLengthEnd    = expTRSS2017.activeLengtheningData(idxTrial).x(end,1); 
    
    timeStart       = expTRSS2017Raw.(seriesName).time(1,1);
    timeEnd         = expTRSS2017Raw.(seriesName).time(end,1)-timeStart;    
    timeRampStart   = expTRSS2017Raw.(seriesName).time(1,1)-timeStart+0.001;
    timeRampEnd     = expTRSS2017Raw.(seriesName).time(end,1)-timeStart;
    timeStart = 0;

    fsample = 1/(expTRSS2017Raw.(seriesName).time(2,1)...
                -expTRSS2017Raw.(seriesName).time(1,1));
    lpfFreqHz=simConfig.numericalDiffLowPassFilterFrequencyHz;
    [b,a]=butter(2,lpfFreqHz/(0.5*fsample),'low');

    rampVelocityFilt = calcCentralDifferenceDataSeries(...
                        expTRSS2017Raw.(seriesName).time,...
                        filtfilt(b,a,expTRSS2017Raw.(seriesName).lN).*lceOptMdl);
    rampVelocity    = mean(rampVelocityFilt);


    timeSpan = [timeStart,timeEnd];

    timeStimulation = timeStart;
    activation      = 1;

    excitationFcn = @(argT)calcStepFunction(argT,...
                      timeStimulation-1,...
                      inf,...
                      activation);
    
    pathLengthFcn = @(argT)calcRampStateSharp(argT,...
                            timeRampStart,timeRampEnd,...
                            rampLengthStart,rampVelocity);   
    
    activationFcn = ...
        @(argU,argA)calcFirstOrderActivationDerivative(argU,argA, ...
            ratFibrilModelsFitted(idxTrial).sarcomere.activationTimeConstant,...
            ratFibrilModelsFitted(idxTrial).sarcomere.deactivationTimeConstant,0);        
  
    flag_debug=1;
    if(flag_debug==1)
      timeSeries = [0:0.01:1]'.*(timeEnd-timeStart)+timeStart;
      exSeries = zeros(length(timeSeries),1);
      pathSeries=zeros(length(timeSeries),2);
      for z=1:1:length(timeSeries)
        exSeries(z) = excitationFcn(timeSeries(z));
        pathSeries(z,:)=pathLengthFcn(timeSeries(z))';
      end
      figTest=figure;
      plot(timeSeries,exSeries,'-b','DisplayName','Excitation');
      hold on;
      plot(timeSeries,pathSeries(:,1),'-','Color',[1,1,1].*0.25,'DisplayName','Path 1');
      hold on;
      plot(timeSeries,pathSeries(:,2),'-','Color',[1,1,1].*0.75,'DisplayName','Path 2');
      hold on;      
      xlabel('Time (s)');
      ylabel('Value');
      legend;
      here=1;
    end

    %
    % Bench config
    %
    benchConfig.npts                  = npts;
    benchConfig.relTol                = 1e-6;
    benchConfig.absTol                = 1e-6;
    benchConfig.minActivation         = 0;
    benchConfig.color0                = [0,0,1].*0.5;
    benchConfig.color1                = [0,0,1];
    
    nStates     = 3;
    labelStates = {'$$\dot{\ell}_{a}$$', '$$\ell_{a}$$', '$$\ell_1$$'};       
    
    benchConfig.numberOfMuscleStates  = nStates;
    benchConfig.stateLabels           = labelStates;
    benchConfig.name                  = '';
    benchConfig.initialState          = [];
    benchConfig.initialActivation     = excitationFcn(0);
    benchConfig.pathFcn               = [];
    benchConfig.excitationFcn         = [];
    benchConfig.activationFcn         = activationFcn; 
    benchConfig.tspan                 = timeSpan;  
    
    benchConfig.useFiberDamping  = 1;
    benchConfig.useElasticTendon = 0;
    benchConfig.damping          = 0.1;
    benchConfig.iterMax          = 100;
    benchConfig.tol              = 1e-6;
    
    loopTolerance = min(benchConfig.relTol,benchConfig.absTol)/100;
    
    %
    % Initialization 
    %
    
    modelConfig = struct( ...
      'iterMax'                 , 100             , ...
      'tol'                     , loopTolerance   , ... 
      'tolInit'                 , sqrt(eps)       , ...
      'minActivation'           , 0.0             , ...
      'useElasticTendon'        , 0 , ...
      'initializeState'         , 0                     );  

    modelConfig.initializeState =1;
    activationState0    = [0;excitationFcn(0)];
    pathState0          = pathLengthFcn(0);
    muscleState0        = zeros(nStates,1);
    mtInfo = calcMillard2023VexatMuscleInfo(...
                activationState0,...
                pathState0,...
                muscleState0,...
                ratFibrilModelsFitted(idxTrial).musculotendon,...
                ratFibrilModelsFitted(idxTrial).sarcomere,...
                ratFibrilModelsFitted(idxTrial).curves,...
                modelConfig);

    muscleState0                = mtInfo.state.value;
    modelConfig.initializeState = 0;           
    
    benchConfig.numberOfMuscleStates = length(muscleState0);
    benchConfig.initialState         = muscleState0;
    
    benchConfig.minimumActivation    = 0;
    benchConfig.name                 = 'Vexat';
    benchConfig.eventFcn             = [];            
    
    calcMillard2023VexatMuscleInfoFcn = ...
         @(activationState1,pathState2,mclState3) ...
         calcMillard2023VexatMuscleInfo(    ...
            activationState1,...
            pathState2,...
            mclState3,...
            ratFibrilModelsFitted(idxTrial).musculotendon,...
            ratFibrilModelsFitted(idxTrial).sarcomere,...
            ratFibrilModelsFitted(idxTrial).curves,...
            modelConfig);
    
    %
    % Run the simulation
    %
    idx                     = idxTrial;
    flag_appendEnergetics   = 0;
    flag_useOctave          = 0;
    
    benchConfig.pathFcn               = pathLengthFcn;
    benchConfig.excitationFcn         = excitationFcn; 
    benchRecord=[];
    benchRecord = runPrescribedLengthActivationSimulation2025(...
                               calcMillard2023VexatMuscleInfoFcn,...
                               [],...
                               benchConfig,...
                               benchRecord,...
                               idx, ...
                               numberOfSimulations,...
                               flag_appendEnergetics,...
                               flag_useOctave);
    
    figure(figSandbox);
    subplot(1,1,1);     
    plot(benchRecord.normFiberLength(:,idx),...
         benchRecord.normFiberForce(:,idx),...
         '-','Color',plotConfig.lineColors.simF(:,idxTrial),...
         'LineWidth',1,...
         'DisplayName',['Sim ',num2str(i)],...
         'HandleVisibility',hdlVis);
    hold on; 
    
end


