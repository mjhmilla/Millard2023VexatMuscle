clc;
close all;
clear all;

flag_useElasticTendon=1;

%%
% Simulation parameters
%%
npts    = 1000;
absTol  = 1e-8;
relTol  = absTol;
maxStep = 0.001;
initTol = sqrt(eps);
loopTol = min(absTol,relTol)/100.;
iterMax = 100;

nCycle      = 1;
omega       = 2*pi; %0 is allowed for a static path
pathStrain  = 0.05;

transientWindowMax = 0.025;


pubOutputFolder = 'output/plots/InitializationBenchmark/';
%%
%Path setup
%%
parametersDirectoryTree       = genpath('parameters');
curvesDirectoryTree           = genpath('curves');
experimentsDirectoryTree      = genpath('experiments');
simulationDirectoryTree       = genpath('simulation');
modelDirectoryTree            = genpath('models');
postprocessingDirectoryTree   = genpath('postprocessing');

addpath(parametersDirectoryTree       );
addpath(curvesDirectoryTree           );
addpath(experimentsDirectoryTree      );
addpath(simulationDirectoryTree       );
addpath(modelDirectoryTree            );
addpath(postprocessingDirectoryTree   );

%%
%Plotting setup
%%
fig_StateRecord = figure;

plotLayoutSettings = struct('numberOfHorizontalPlotColumns',  4,...
                            'numberOfVerticalPlotRows',       2,...
                            'flag_fixedPlotWidth',            1,...
                            'plotWidth',                      7,...
                            'plotHeight',                     7,...
                            'flag_usingOctave',               0);

numberOfHorizontalPlotColumns = plotLayoutSettings.numberOfHorizontalPlotColumns;
numberOfVerticalPlotRows      = plotLayoutSettings.numberOfVerticalPlotRows;
flag_fixedPlotWidth           = plotLayoutSettings.flag_fixedPlotWidth;
plotWidth                     = plotLayoutSettings.plotWidth;
plotHeight                    = plotLayoutSettings.plotHeight;
flag_usingOctave              = plotLayoutSettings.flag_usingOctave;

plotHorizMarginCm = 2;
plotVertMarginCm  = 2;

pageWidth   = (plotWidth+plotHorizMarginCm)*numberOfHorizontalPlotColumns...
                +plotHorizMarginCm;

pageHeight  = (plotHeight+plotVertMarginCm)*numberOfVerticalPlotRows...
                +plotVertMarginCm;

plotConfigGeneric;

%%
%Load the model
%%

load('output/structs/defaultFelineSoleus.mat')
musculotendonProperties   = defaultFelineSoleus.musculotendon;
sarcomereProperties       = defaultFelineSoleus.sarcomere;
normMuscleCurves          = defaultFelineSoleus.curves;

%%
%Configure the model, save the state
%%
tendonTag='';
if(flag_useElasticTendon==1)
    tendonTag = '_ElasticTendon';
else
    tendonTag = '_RigidTendon';
end


nStates = 0;
labelStates = {''};
labelStateDerivatives = {''};

if(flag_useElasticTendon==1)
  nStates = 4;
  labelStates           = {'$$\ell_{CE}$$','$$\dot{\ell}_{a}$$','$$\ell_{a}$$','$$\ell_1$$'};
  labelStateDerivatives = {'$$\dot{\ell}_{CE}$$','$$\ddot{\ell}_{a}$$','$$\dot{\ell}_{a}$$','$$\dot{\ell}_1$$'};
  
else
  nStates = 3;
  labelStates           = {'$$\dot{\ell}_{a}$$','$$\ell_{a}$$','$$\ell_1$$'};
  labelStateDerivatives = {'$$\ddot{\ell}_{a}$$','$$\dot{\ell}_{a}$$','$$\dot{\ell}_1$$'};
end


excitationFcn = @(argT)calcStepFunction(argT,...
                              0.0,...
                              Inf,...
                              1.0);

activationFcn = @(argU,argA)calcFirstOrderActivationDerivative(...
                  argU,argA, sarcomereProperties.activationTimeConstant,...
                             sarcomereProperties.deactivationTimeConstant,0);

lceOpt  = musculotendonProperties.optimalFiberLength;
alphaOpt= musculotendonProperties.pennationAngle;
ltSlk   = musculotendonProperties.tendonSlackLength;

lp0         = lceOpt*cos(alphaOpt)+ltSlk;
amp         = pathStrain*lp0;
tStart      = 0;

if(omega == 0)
    tEnd        = nCycle*(2*pi);
else
    tEnd        = nCycle*(2*pi)/omega;
end

transientWindow = min((tEnd-tStart)/20, transientWindowMax);


pathFcn = @(argT)calcSinusoidState(argT, 0, tEnd, ...
                lp0, omega, amp);


modelConfig = struct( ...
  'iterMax'                 , iterMax   , ...
  'tol'                     , loopTol   , ... 
  'tolInit'                 , initTol   , ...
  'minActivation'           , 0.0             , ...
  'useElasticTendon'        , flag_useElasticTendon , ...
  'initializeState'         , 1                     );    


%%
%Initialize the model, save the state
%%


modelConfig.initializeState = 1;
activationState0            = [0;excitationFcn(0)];
pathState0                  = pathFcn(0);
muscleState0                = zeros(nStates,1);

mtInfo = calcMillard2019MuscleInfoOpus31(activationState0,...
                                        pathState0,...
                                        muscleState0,...
                                        musculotendonProperties,...
                                        sarcomereProperties,...
                                        normMuscleCurves,...
                                        modelConfig);

%%
%Simulate the model
%%
    modelConfig.initializeState = 0;   
    state0                = mtInfo.state.value;


    tV  = ([0:(1/(npts-1)):1]').*(tEnd-tStart);

    options = odeset('RelTol',relTol,...
                     'AbsTol',absTol,...
                     'NormControl','on',...
                     'MaxStep',maxStep,...
                     'Stats','off');

    calcMillard2019MuscleInfoOpus31Fcn = ...
       @(activationState1,pathState2,mclState3) ...
       calcMillard2019MuscleInfoOpus31(   activationState1,...
                                          pathState2,...
                                          mclState3,...
                                          musculotendonProperties,...
                                          sarcomereProperties,...
                                          normMuscleCurves,...
                                          modelConfig);

    flag_appendEnergetics=0;

    dfcn = @(argt,argState)...
        calcPrescribedMusculotendonStateDerivativeWrapper(...
                          argt,...
                          argState,...                                                      
                          pathFcn,...
                          excitationFcn,...
                          activationFcn,...
                          calcMillard2019MuscleInfoOpus31Fcn,...
                          flag_appendEnergetics);

    t0 = tic;
    [xe, ye] = ode15s(dfcn,tV,[activationState0(2,1);state0],options);                                                   
    cpuTime = toc(t0);  

    activation = ye(:,1);
    state  = ye(:,2:1:(nStates+1));
    dstate = zeros(length(tV),nStates);
    forces = zeros(length(tV),1);
    forcesLabels = {'$$f^{\mathrm{T}}$$'};

    for indexTime=1:1:length(tV)
        dactivation = activationFcn(excitationFcn(tV(indexTime,1)),...
                                    activation(indexTime,1));

        activationState = [dactivation,activation(indexTime,1)];
        pathState       = pathFcn(tV(indexTime,1));
    
        dlp = pathState(1);
        lp  = pathState(2);
        
        mtInfo = calcMillard2019MuscleInfoOpus31Fcn(...
                                   activationState,...                                   
                                   pathState,...
                                   state(indexTime,:));

        dstate(indexTime,:) = mtInfo.state.derivative;
        forces(indexTime,1) = mtInfo.muscleDynamicsInfo.tendonForce;
    end

    figure(fig_StateRecord);
    for indexState=1:1:nStates
        row = floor((indexState-1)/size(subPlotPanel,2)) + 1;
        col = indexState-(row-1)*size(subPlotPanel,2);
        
        subplot('Position',reshape(subPlotPanel(row,col,:),1,4));
        plot(tV,state(:,indexState),'Color',[0,0,1]);
        hold on;
        plot(tV,dstate(:,indexState),'Color',[1,0,0]);
        hold on;

        text(tV(end),state(end,indexState),labelStates{indexState},...
             'HorizontalAlignment','left');
        hold on;
        text(tV(end),dstate(end,indexState),...
                ['d/dt ',labelStates{indexState}],...
                'HorizontalAlignment','left');


        plot([0;1].*tV(end),...
             [state(1,indexState);state(end,indexState)],...
             'Color',[0,0,0]);
        hold on;

        xTxt  = mean([0;1].*tV(end));
        yTxt  = state(1,indexState); 
        yErr  = (state(end,indexState)-state(1,indexState))...
                /( mean(state(:,indexState)) );
        text(xTxt,yTxt,...
                ['$$y^{C}_{rel}$$',sprintf('%1.2e',yErr)],...
                'HorizontalAlignment','center',...
                'VerticalAlignment','bottom');        
        hold on;

        idxStartWindow  = find(tV < transientWindow);
        idxEndWindow    = find(tV > (max(tV)-transientWindow) );

        xTxt = 0;
        [maxAbsVal, idxStartTransient]= ...
            max(abs(state(idxStartWindow,indexState)));
        [maxAbsVal, idxEndTransient]= ...
            max(abs(state(idxEndWindow,indexState)));

        x0=tV(1,1);
        x1=tV(max(idxStartWindow),1);
        y0=min(state(idxStartWindow,indexState));
        y1=max(state(idxStartWindow,indexState));
        yMean=mean(state(:,indexState));
        fill([x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0],...
            [1,1,1].*0.5,'FaceAlpha',0.5);
        hold on;

        text(x1,y0,sprintf('$$y^{T}_{rel}$$: %1.2e',(y1-y0)/yMean),...
             'HorizontalAlignment','left',...
             'VerticalAlignment','top');
        hold on;

        x0=tV(min(idxEndWindow),1);
        x1=tV(end,1);
        y0=min(state(idxEndWindow,indexState));
        y1=max(state(idxEndWindow,indexState));
        yMean=mean(state(:,indexState));        
        fill([x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0],...
            [1,1,1].*0.5,'FaceAlpha',0.5);
        hold on;

        text(x1,y0,sprintf('$$y^{T}_{rel}$$: %1.2e',(y1-y0)/yMean),...
             'HorizontalAlignment','right',...
             'VerticalAlignment','top');
        hold on;


        xlabel('Time (s)');
        title(['State: ',labelStates{indexState},' and ', labelStateDerivatives{indexState}]);
        box off;
    end

    for indexForce= 1:1:size(forces,2)

        indexPlot   = nStates+indexForce;
        row         = floor((indexPlot-1)/size(subPlotPanel,2)) + 1;
        col         = indexPlot-(row-1)*size(subPlotPanel,2);

        subplot('Position',reshape(subPlotPanel(row,col,:),1,4));
        plot(tV,forces(:,indexForce),'b');
        hold on;

        plot([tV(1,1);tV(end,1)],... 
             [forces(1,indexForce);forces(end,indexForce)],...
             'Color',[0,0,0]);
        hold on;

        xTxt  = mean([0;1].*tV(end));
        yTxt  = forces(1,indexForce); 
        yErr  = (forces(end,indexForce)-forces(1,indexForce))...
               /( mean(forces(:,indexForce)));
        text(xTxt,yTxt,...
                ['$$y^{C}_{rel}$$',sprintf('%1.2e',yErr)],...
                'HorizontalAlignment','center',...
                'VerticalAlignment','bottom');        
        hold on;

        idxStartWindow  = find(tV < transientWindow);
        idxEndWindow    = find(tV > (max(tV)-transientWindow) );

        x0=tV(1,1);
        x1=tV(max(idxStartWindow),1);
        y0=min(forces(idxStartWindow,indexForce));
        y1=max(forces(idxStartWindow,indexForce));
        yMean = mean(forces(:,indexForce));
        fill([x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0],...
            [1,1,1].*0.5,'FaceAlpha',0.5);
        hold on;

        text(x1,y0,sprintf('dy: %1.2e',(y1-y0)/yMean),...
             'HorizontalAlignment','left',...
             'VerticalAlignment','top');
        hold on;

        x0=tV(min(idxEndWindow),1);
        x1=tV(end,1);
        y0=min(forces(idxEndWindow,indexForce));
        y1=max(forces(idxEndWindow,indexForce));
        fill([x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0],...
            [1,1,1].*0.5,'FaceAlpha',0.5);
        hold on;

        text(x1,y0,sprintf('dy: %1.2e',(y1-y0)/yMean),...
             'HorizontalAlignment','left',...
             'VerticalAlignment','top');
        hold on;

        xlabel('Time (s)');
        title(['Forces: ',forcesLabels{indexForce}]);
        box off;

    end

    
        
    set(fig_StateRecord,'Units','centimeters',...
    'PaperUnits','centimeters',...
    'PaperSize',[pageWidth pageHeight],...
    'PaperPositionMode','manual',...
    'PaperPosition',[0 0 pageWidth pageHeight]);     
    %set(findall(figList(i).h,'-property','FontSize'),'FontSize',10);     
    set(fig_StateRecord,'renderer','painters');     
    set(gcf,'InvertHardCopy','off')

    print('-dpdf', [pubOutputFolder,'fig_StateRecord',tendonTag,'_.pdf']);     


