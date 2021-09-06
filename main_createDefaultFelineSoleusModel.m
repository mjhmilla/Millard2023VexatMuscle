flag_outerLoopMode = 0;

if(flag_outerLoopMode == 0)

  clc;
  close all;
  clear all;
  fitCrossBridgeStiffnessDampingToKirch199490Hz=1;
  flag_useFixedLambdaECM    = 0;
  flag_makeAndSavePubPlots  = 1;
end

pubOutputFolder = 'output/plots/MuscleCurves/';
postprocessingDirectoryTree      = genpath('postprocessing');
addpath(postprocessingDirectoryTree   );



flag_useOctave = 0;
flag_enableNumericallyNonZeroGradients    = 1;

% 1. Active force length curve vs. data
% Solution: There were some initial descrepencies between the experimental force
%length data and a theoretical curve. These errors almost completely go
%away if it is assumed that the experimental recordings are of total 
%path length, rather than fiber length. In this case, when the elastiticy
%of the tendon is taken into account the theoretical active-force-length 
%curve and the transformed data nicely align.

%Failed attempt:
%This creates a cat soleus with an optimal fiber length of 58 mm: this
%is simply way too big to be realistic (given the data I'm working from)
flag_solveForOptimalFiberLengthOfBestFit  = 0; 

%Failed attempt:
shiftLengthActiveForceLengthCurveDescendingCurve = 0.;%...
%  (1/3)*( (1.154-1.087) + (1.23-1.162) + (1.077-1.039) );

smallNumericallyNonZeroNumber           = sqrt(sqrt(eps));

%%
% Add the directories needed to run this script
%%
parametersDirectoryTreeMTParams     = genpath('parameters');
parametersDirectoryTreeExperiments  = genpath('experiments');
parametersDirectoryTreeModels       = genpath('models');
parametersDirectoryTreeCurves       = genpath('curves');

addpath(parametersDirectoryTreeMTParams);
addpath(parametersDirectoryTreeExperiments);
addpath(parametersDirectoryTreeModels);
addpath(parametersDirectoryTreeCurves);



scaleOptimalFiberLength      = 1.0; 

scaleMaximumIsometricTension = 1;

if(exist('fitCrossBridgeStiffnessDampingToKirch199490Hz','var')==0)
  fitCrossBridgeStiffnessDampingToKirch199490Hz = 1;
end

[felineSoleusMusculotendonProperties, ...
 felineSoleusSarcomereProperties,...
 felineSoleusActiveForceLengthData,...
 felineSoleusPassiveForceLengthData] = createFelineSoleus(...                                          
                                          scaleOptimalFiberLength,...
                                          scaleMaximumIsometricTension,...
                                          fitCrossBridgeStiffnessDampingToKirch199490Hz,...
                                          flag_useOctave);

createMusculoTendonFcn = ...
  @(argScaleFiberLength,argScaleFiso)createFelineSoleus(...
                                        argScaleFiberLength,...
                                        argScaleFiso,...
                                        flag_useOctave); 
                                        
[felineSoleusNormMuscleCurves,...
 felineSoleusMusculotendonPropertiesUpd,...
 felineSoleusSarcomerePropertiesUpd,...
 activeForceLengthCurveAnnotationPoints,...
 felineSoleusActiveForceLengthDataUpd,...
 felineSoleusPassiveForceLengthDataUpd]= ...
    createFittedMuscleCurves( ...
      felineSoleusMusculotendonProperties,...
      felineSoleusSarcomereProperties,...
      felineSoleusActiveForceLengthData,...
      felineSoleusPassiveForceLengthData,...
      shiftLengthActiveForceLengthCurveDescendingCurve,...
      flag_useFixedLambdaECM,...
      flag_enableNumericallyNonZeroGradients,...
      smallNumericallyNonZeroNumber,...
      flag_solveForOptimalFiberLengthOfBestFit,...
      createMusculoTendonFcn,...
      flag_useOctave);

figH = plotStructOfBezierSplines( felineSoleusNormMuscleCurves,...
                                  'Inverse');                          

%%
% Note the average offset between the active-force-length curve and
% the transformed data
%%

xExp = felineSoleusActiveForceLengthDataUpd(2:end,1);
yExp = felineSoleusActiveForceLengthDataUpd(2:end,2);
xCurve = zeros(size(xExp));

for i=1:1:length(xExp)
xCurve(i,1) = calcBezierFcnXGivenY(yExp(i,1), ...
  felineSoleusNormMuscleCurves.activeForceLengthCurve,... 
  xExp(i,1));
end                                    
dx = mean(xCurve-xExp);
%felineSoleusActiveForceLengthDataUpd(:,1)=...
%  felineSoleusActiveForceLengthDataUpd(:,1)+dx;

disp('Normalized length offset');
fprintf('%1.6f lce/lopt \n',felineSoleusActiveForceLengthDataUpd(1,1));
disp('Average error on the descending limb');
fprintf('%1.6f lce/lopt \n',dx);
fprintf('%1.6f mm \n',dx*(felineSoleusMusculotendonPropertiesUpd.optimalFiberLength*1000));

lceNStart = felineSoleusActiveForceLengthDataUpd(1,1);
save('output/structs/normalizedFiberLengthStartHerzogLeonard2002.mat',...
     'lceNStart');

% dl = felineSoleusActiveForceLengthDataUpd(2:end,1)-1;
% A  = [dl ones(size(dl))];
% b  = felineSoleusActiveForceLengthDataUpd(2:end,2);
% 
% x     = (A'*A)\(A'*b);
% y0    = x(2,1);
% dydx0 = x(1,1);
% 
% felineSoleusActiveForceLengthLineBestFit = zeros(length(dl),1);
% 
% felineSoleusActiveForceLengthLineBestFit(:,1) = ...
%   felineSoleusActiveForceLengthDataUpd(2:end,1);
% dl = felineSoleusActiveForceLengthLineBestFit(:,1)-1;
% 
% felineSoleusActiveForceLengthLineBestFit(:,2) = dydx0.*dl + y0;
% 
% disp('Active force length line of best fit y=(dydx)*(x-1) + y0');
% fprintf('dydx: %1.3f\n',dydx0);
% fprintf('  y0: %1.3f\n',y0);



%%
% Plot the derivative of the tendon force length curve on top of the
% stiffness curve
%% 
figure(figH.tendonStiffnessCurve);

    curveSample = calcBezierYFcnXCurveSampleVector(...
                    felineSoleusNormMuscleCurves.('tendonForceLengthCurve'), 200,[]);

    xmin = min(curveSample.x);
    xmax = max(curveSample.x);
    ymin = min(curveSample.y);
    ymax = max(curveSample.y);

subplot(2,2,1);
  plot(curveSample.x, curveSample.dydx,...
    '--','Color',[1,1,1].*0.5,'LineWidth',2);
  hold on;

%%
%Plot experimental data over top of the curves where it is available.
%%
figure(figH.activeForceLengthCurve);
  subplot(2,2,1);  
  plot(  felineSoleusActiveForceLengthDataUpd(:,1),...
       felineSoleusActiveForceLengthDataUpd(:,2),'xb');
  hold on;          

figure(figH.fiberForceLengthCurve);
  subplot(2,2,1);
  plot(   felineSoleusPassiveForceLengthDataUpd(:,1),...
          felineSoleusPassiveForceLengthDataUpd(:,2),'xb');
  hold on;          
  
%%
% Plot the passive force length curve of the 2 segment titin model and
% compare it to the passive force length curve
%%

lceN0 = calcBezierFcnXGivenY(0, ...
          felineSoleusNormMuscleCurves.fiberForceLengthCurve);
lceN1 = calcBezierFcnXGivenY(1, ...
          felineSoleusNormMuscleCurves.fiberForceLengthCurve);
        
npts = 100;
lceNSeries = [lceN0:((lceN1-lceN0)/(npts-1)):lceN1]';

fNSeries = zeros(npts,4); % fpe, fecm, fIgp, fPevkIgd,
lNSeries = zeros(npts,4); % lpe, lecm, lIgp, lPevkIgd,

normLengthIgdFixed = ...
  felineSoleusSarcomerePropertiesUpd.IGDFixedNormLengthAtOptimalFiberLength;

normLengthT12ToZ = ...
  felineSoleusSarcomerePropertiesUpd.ZLineToT12NormLengthAtOptimalFiberLength;


for i=1:1:npts
  lceN  = lceNSeries(i,1);
  fpeN  = calcBezierYFcnXDerivative(lceN,...
            felineSoleusNormMuscleCurves.fiberForceLengthCurve,0);
  fecmN = calcBezierYFcnXDerivative(lceN*0.5,...
            felineSoleusNormMuscleCurves.forceLengthECMHalfCurve,0);
          
  lIgpPevkN = lceN*0.5 - normLengthIgdFixed - normLengthT12ToZ; 
          
  [lIgpN, lPevkIgdN, fTiN] = calcSeriesSpringStretch(lIgpPevkN,...
            felineSoleusNormMuscleCurves.forceLengthIgpCurve,...
            felineSoleusNormMuscleCurves.forceLengthIgpInverseCurve, ...
            felineSoleusNormMuscleCurves.forceLengthPevkIgdCurve,...
            felineSoleusNormMuscleCurves.forceLengthPevkIgdInverseCurve);          

  fNSeries(i,:) = [fpeN,fecmN,fTiN, fTiN];
  lNSeries(i,:) = [lceN, lceN,lIgpN*2, lPevkIgdN*2];
          
end


fig_forceLength = figure;
  plot(lNSeries(:,1),fNSeries(:,1),'k');
  hold on;
  plot(lNSeries(:,1), fNSeries(:,2)+fNSeries(:,3),'b');
  hold on;
  plot(lNSeries(:,1), fNSeries(:,2),'r');
  hold on;
  plot(felineSoleusPassiveForceLengthDataUpd(:,1),...
       felineSoleusPassiveForceLengthDataUpd(:,2),'xb');
  legend('fpe','fecm+fti','fecm','data');
  xlabel('Norm. Length')
  ylabel('Norm. Force');
  


defaultFelineSoleus = struct('musculotendon',...
                            felineSoleusMusculotendonPropertiesUpd,...
                            'sarcomere',...
                            felineSoleusSarcomerePropertiesUpd,...
                            'falData',...
                            felineSoleusActiveForceLengthDataUpd,...
                            'fpeData',...
                            felineSoleusPassiveForceLengthDataUpd,...
                            'curves',...
                            felineSoleusNormMuscleCurves);
                      
save('output/structs/defaultFelineSoleus.mat',...
     'defaultFelineSoleus');                      


%%
% Generate publication quality plots
%%
if(flag_makeAndSavePubPlots==1)
  plotMuscleCurves( felineSoleusNormMuscleCurves,...
                      activeForceLengthCurveAnnotationPoints,...
                      felineSoleusMusculotendonPropertiesUpd,...
                      felineSoleusSarcomerePropertiesUpd,...
                      felineSoleusActiveForceLengthDataUpd,...
                      felineSoleusPassiveForceLengthDataUpd,...
                      pubOutputFolder);
end
   
%%
% Remove the directories ...
%%
rmpath(parametersDirectoryTreeMTParams);
rmpath(parametersDirectoryTreeExperiments);
rmpath(parametersDirectoryTreeModels);
rmpath(parametersDirectoryTreeCurves);
