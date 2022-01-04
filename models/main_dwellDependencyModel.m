clc;
close all;
clear all;





npts = 1000;
timeStart   = 0;
timeEnd     = 1;
time        = [timeStart:((timeEnd-timeStart)/(npts-1)):timeEnd]';


%%
%Simulate the dwell-time-dependency dynamics
%%
calcLengthFcnA = @(argT)oscillatingSmoothStepFunction(argT,1,0.01,2*pi*20, 0.5,0.55, 0.2);
calcLengthFcnB = @(argT)oscillatingSmoothStepFunction(argT,1,0.0,2*pi*2 , 0.5,0.55, 0.2);

pathA = zeros(length(time),2);
pathB = zeros(length(time),2);

for i=1:1:length(time)
    pathA(i,:) = calcLengthFcnA(time(i,1));
    pathB(i,:) = calcLengthFcnB(time(i,1));
end

velocityStick = 0.01;
shortRange    = 0.01;

zeta         = 1;
omega        = 50*2*pi;
timeConstant = 0.001;


dfcnA = @(argT,argY)stateDerivativeDwellDependency(argT,argY,...
                     calcLengthFcnA,omega,zeta,timeConstant);

dfcnB = @(argT,argY)stateDerivativeDwellDependency(argT,argY,...
                     calcLengthFcnB,omega,zeta,timeConstant);

options = odeset('MaxStep',1e-3);
[tA, yA]=ode15s(dfcnA,time,[0;0],options);
[tB, yB]=ode15s(dfcnB,time,[0;0],options);

%%
% Plotting
%%
addpath('../postprocessing');
flag_usingOctave              = 0;
numberOfHorizontalPlotColumns = 2;
numberOfVerticalPlotRows      = 4;
plotWidth                     = 7;
plotHeight                    = 7.0;
plotHorizMarginCm             = 2.0;
plotVertMarginCm              = 2.0;
pageHeight                    = 29.7;
pageWidth                     = 21.0;
plotConfigGeneric;

fig=figure;
subplot('Position',reshape(subPlotPanel(1,1,:),1,4));

    plot(time, pathA(:,1),'b','DisplayName','$$\ell^{A}(t)$$');
    hold on;
    plot(time, pathB(:,1),'r','DisplayName','$$\ell^{B}(t)$$');
    hold on;
    xlabel('Time (s)');
    ylabel('Norm. Length $$(\ell / \ell_\circ)$$')
    title('Path Length');
    box off;
    legend('Location','NorthWest');

subplot('Position',reshape(subPlotPanel(1,2,:),1,4));

    plot(time, pathA(:,2),'b','DisplayName','$$\dot{\ell}^{A}(t)$$');
    hold on;
    plot(time, pathB(:,2),'r','DisplayName','$$\dot{\ell}^{B}(t)$$');
    hold on;
    xlabel('Time (s)');
    ylabel('Norm. Velocity $$(\dot{\ell} / \ell_\circ)$$')
    title('Path Velocity');
    box off;

subplot('Position',reshape(subPlotPanel(2,1,:),1,4));
    plot(tA, yA(:,1), 'b','DisplayName', '$$s^{A}$$');
    hold on;
    plot(tB, yB(:,1), 'r','DisplayName', '$$s^{B}$$');
    hold on;
    box off;
    xlabel('Time (s)');
    ylabel('Norm. Length $$(\ell / \ell_\circ)$$')
    title('Short Range Length');
    box off;

subplot('Position',reshape(subPlotPanel(2,2,:),1,4));
    plot(tA, exp(-(yA(:,1)./shortRange).^2), 'b','DisplayName', '$$\lambda^{A}$$');
    hold on;
    plot(tB, exp(-(yB(:,1)./shortRange).^2), 'r','DisplayName', '$$\lambda^{B}$$');
    hold on;
    box off;
    xlabel('Time (s)');
    ylabel('Value')
    title('Short-Range Switching Value');
    box off;

rmpath('../postprocessing');