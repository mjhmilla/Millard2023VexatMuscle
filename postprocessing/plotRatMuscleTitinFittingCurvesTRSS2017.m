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

function [success] = plotRatMuscleTitinFittingCurvesTRSS2017( ...
                                          muscleModel,...
                                          humanSoleusMuscle,...                      
                                          activeForceLengthCurveAnnotationPoints,...
                                          fittingData,...
                                          activeForceLengthData,...
                                          passiveForceLengthData,...   
                                          normFiberLengthAtOneNormPassiveForce,...
                                          trialId,...
                                          previousFilePathAndName,...
                                          filePathAndName,...
                                          projectFolders)
%%
% 
%%
   assert(humanSoleusMuscle.sarcomere.titinModelType ...
          == muscleModel.sarcomere.titinModelType,...
          ['Error: the muscleModel and humanSoleusModel need to have' ...
          'the same titinModelType']);

scaleLength = muscleModel.musculotendon.optimalFiberLength;   
scaleVelocity = muscleModel.musculotendon.maximumNormalizedFiberVelocity;
%%
% Series colors
%%
cs = getPaulTolColourSchemes('vibrant');

n=0.75;

lineColors.orig(1,:)=[0,0,0];
lineColors.orig(2,:)=[44,133,229]./255;
lineColors.orig(3,:)=[63,181,175]./255;

for i=1:1:3
    n0 = 0.25*(i-1);
    n1 = n0;%n0+0.25;
    lineColors.exp(i,:)         = [0,0,0].*(1-n1)+[1,1,1].*n1;  
    lineColors.simXE(i,:)       = cs.grey;
    lineColors.calcTitinF(i,:)  = lineColors.exp(i,:);%cs.blue.*(1-n0)+[1,1,1].*n0;
    lineColors.calcTitinK(i,:)  = lineColors.exp(i,:);%lineColors.calcTitinF(i,:);
    lineColors.simF(i,:)        = cs.blue.*(1-n0)+[1,1,1].*n0;%cs.red.*(1-n0)+[1,1,1].*n0;
    lineColors.simTitinF(i,:)   = cs.red.*(1-n0)+[1,1,1].*n0;%cs.magenta.*(1-n0)+[1,1,1].*n0;
    lineColors.simTitinPassive(i,:)   = cs.magenta.*(1-n0)+[1,1,1].*n0;%cs.magenta.*(1-n0)+[1,1,1].*n0;
    lineColors.simTitinActive(i,:)     = cs.red.*(1-n0)+[1,1,1].*n0;%cs.magenta.*(1-n0)+[1,1,1].*n0;
    
    lineColors.simTitinK(i,:)   = lineColors.simTitinF(i,:);
end



  if(trialId==0)
    n=0;
  else
    n = 1-((trialId-1)/2);
  end

  colorTitin            = [1,0,1].*n + ([1,0,1].*0.25).*(1-n);
  colorTitinLight   = colorTitin.*0.1 + [1,1,1].*0.9;  

  colorIgp                = [1,0,0].*n + ([1,0,0].*0.25).*(1-n);
  colorPEVK            = [0.5,0,0].*n + ([0.5,0,0].*0.25).*(1-n);
  colorIgd                = [0.25,0,0].*n + ([0.25,0,0].*0.25).*(1-n);
  colorDistal           = colorTitin;
  colorProximal     = colorDistal.*0.5;


%%
%Plot configuration
%%
  flag_fibrilModel = 0;
  if(muscleModel.sarcomere.extraCellularMatrixPassiveForceFraction==0)
    flag_fibrilModel=1;
  end


  fig_pubTitinCurves=figure;

  numberOfHorizontalPlotColumns = 3;
  numberOfVerticalPlotRows      = 2;

  pageWidth         = 21;
  pageHeight        = 29.7;
  plotWidth         = 3.75;%4.5;
  plotHeight        = 3.75;%4.5;
  plotHorizMarginCm = 1.0;
  plotVertMarginCm  = 2.;

  flag_usingOctave  = 0;
  plotConfigGeneric;  
    
  lineColorTendon = [0,0,0];
  lineWidthTendon = 1;
  
  labelRotationOffset = -10;
  
  plotProps(6) = struct('xlim',[],'ylim',[],'domain',[],...
                     'xticks',[],'xticklabels',{''},...
                     'yticks',[],'yticklabels',{''},...
                     'lineColor',[],'lineWidth',[],...
                     'xlabel','', 'ylabel','','title','');


  %%
  % Data
  %%
  fileTrombitas1998Figure5 = fullfile(projectFolders.experiments_TGFG1998,...
      'Trombitas1998_Figure5.csv');

  dataTrombitas1998Figure5 = loadDigitizedData(fileTrombitas1998Figure5,...
                        'Sarcomere Length','PEVK Width (um)',...
                        {'a','b'},'Trombitas 1998 Figure 5');  
                      
  fileTrombitas1998Figure6 = fullfile(projectFolders.experiments_TGFG1998, ...
      'Trombitas1998_Figure6.csv');
  dataTrombitas1998Figure6 = loadDigitizedData(fileTrombitas1998Figure6,...
                        'Sarcomere Length','PEVK Width (um)',...
                        {'a'},'Trombitas 1998 Figure 6');  

  %%
  % Tendon force-length
  %%
                   
  idxTendon =1;
  idx       =idxTendon;
  
  plotProps(idx).xlim = ...
    [1.0, (1+muscleModel.musculotendon.tendonStrainAtOneNormForce)] ...
    +[-0.001,0.001];
  
  plotProps(idx).ylim = [0,1.0]+[-0.01,0.01];
  plotProps(idx).domain      = plotProps(idx).xlim + [-0.01,0.01];
  plotProps(idx).xticks      = [round(1.0,3), round((1+muscleModel.musculotendon.tendonStrainAtOneNormForce),3)];
  plotProps(idx).xticklabels = {'0','$$e^T_o$$'};
  plotProps(idx).yticks      = [round(0.0,3), round(2/3,3),round(1.0,3)];
  plotProps(idx).yticklabels = {'0','$$f^T_{toe}$$','$$f^M_o$$'};
  plotProps(idx).lineColor   = [0,0,0];
  plotProps(idx).lineWidth   = 0.5;
  plotProps(idx).xlabel = 'Strain ($$\ell/\ell^{T}_{S}-1$$)';
  plotProps(idx).ylabel = 'Normalized Force ($$f/f^{M}_{o}$$)';
  plotProps(idx).title  = 'A. Tendon-force length ($$\mathbf{f}^T$$)';

  %%
  % CE force-length
  %%
  idxCELength = 2;
  idx         =idxCELength;
  normMyosinLength      = muscleModel.sarcomere.normMyosinHalfLength*2;
  normMyosinBareLength  = muscleModel.sarcomere.normMyosinBareHalfLength*2;
  normActinLength       = muscleModel.sarcomere.normActinLength;
  normZLineThickness    = muscleModel.sarcomere.normZLineLength;
  normSarcomereLengthZeroForce = muscleModel.sarcomere.normSarcomereLengthZeroForce;  
  
  %from createFiberActiveForceLength.m
  plotProps(idx).xlim = ...
      (muscleModel.curves.activeForceLengthCurve.xEnd +[-0.02,0.02]).*scaleLength;
  
  plotProps(idx).ylim   = [0,1.0]+[-0.01,0.01];
  plotProps(idx).domain = plotProps(idx).xlim + [-0.01,0.01];
  plotProps(idx).xticks = ...
      [ round(activeForceLengthCurveAnnotationPoints.x(1,1),2),...                                
        1,...
        round(muscleModel.curves.fiberForceLengthCurve.xEnd(1,2),2),...
        round(activeForceLengthCurveAnnotationPoints.x(end,1),2)].*scaleLength;
  
  plotProps(idx).xticks = sort(plotProps(idx).xticks);

  plotProps(idx).xticklabels = {num2str(plotProps(idx).xticks(1,1)),...
                                num2str(plotProps(idx).xticks(1,2)),...
                                num2str(plotProps(idx).xticks(1,3)),...
                                num2str(plotProps(idx).xticks(1,4))};
                              
  plotProps(idx).yticks      = [0,1];
  plotProps(idx).yticklabels = {'0','$$f^M_o$$'};
  plotProps(idx).lineColor   = lineColors.simF(1,:);
  plotProps(idx).lineWidth   = 0.5;

  if(abs(scaleLength-1)<1e-3)
      plotProps(idx).xlabel = 'Normalized Length ($$\ell/\ell^{M}_{o}$$)';
      plotProps(idx).ylabel = 'Normalized Force ($$f/f^{M}_{o}$$)';
  else
      plotProps(idx).xlabel = 'Length ($$\mu$$m)';
      plotProps(idx).ylabel = 'Normalized Force ($$f/f^{M}_{o}$$)';
  end
  plotProps(idx).title  ={ 'A. CE and titin ','force-length relations'};

  %%
  % CE force-velocity
  idxCEVelocity = 3;
  idx           = idxCEVelocity;

  %from createFiberActiveForceLength.m
  plotProps(idx).xlim = [-1.01,1.01].*scaleVelocity;
  
  plotProps(idx).ylim = [0,muscleModel.curves.fiberForceVelocityCurve.yEnd(1,2)] ...
                        +[-0.01,0.01];
                      
  plotProps(idx).domain      = plotProps(idx).xlim + [-0.01,0.01];
  
  plotProps(idx).xticks      = [-1,0,1].*scaleVelocity;

%  plotProps(idx).xticklabels = {'$$-1$$','0','$$1$$'};
                              
  plotProps(idx).yticks      = [0,1,muscleModel.curves.fiberForceVelocityCurve.yEnd(1,2)];
  plotProps(idx).yticklabels = {0,1,num2str(round(plotProps(idx).yticks(1,end),2))};
  plotProps(idx).lineColor   = [0,0,0];
  plotProps(idx).lineWidth   = 0.5;
  plotProps(idx).xlabel = 'Normalized Velocity ($$v/\ell^{M}_{o}$$)';
  plotProps(idx).ylabel = 'Normalized Force ($$f/f^{M}_{o}$$)';
  plotProps(idx).title  = {'B. CE force-velocity relation',''};
  %%

  
  %%
  % PE: ECM + titin force-length
  %%
  idxPeEcmTitinLength = 4;
  idx         =idxPeEcmTitinLength;
  

  %from createFiberActiveForceLength.m
  plotProps(idx).xlim = ...
    (muscleModel.curves.fiberForceLengthCurve.xEnd +[-0.02,0.02]);
  
  plotProps(idx).ylim = [0,1.0]+[-0.01,0.01];
  plotProps(idx).domain      = plotProps(idx).xlim + [-0.01,0.01];
  plotProps(idx).xticks      = ...
      [round(muscleModel.curves.fiberForceLengthCurve.xEnd(1,1),2),...                                
        1,...
        round(muscleModel.curves.fiberForceLengthCurve.xEnd(1,2),2)];

  plotProps(idx).xticklabels = ...
                              {num2str(plotProps(idx).xticks(1,1)),...
                                num2str(plotProps(idx).xticks(1,2)),...
                                num2str(plotProps(idx).xticks(1,3))};
                              
  plotProps(idx).yticks      = [0,1];
  plotProps(idx).yticklabels = {'0','$$f^M_o$$'};
  plotProps(idx).lineColor   = [0,0,0];
  plotProps(idx).lineWidth   = 0.5;
  %if(abs(scaleVelocity-1)<1e-3)  
      plotProps(idx).xlabel = 'Normalized Length ($$\ell/\ell^{M}_{o}$$)';
  %else
  %    plotProps(idx).xlabel = 'Length ($$\mu$$m)';
  %end
  plotProps(idx).ylabel = 'Normalized Force ($$f/f^{M}_{o}$$)';
  plotProps(idx).title  = {'B. Passive force-length','relation'};  
  
  %%
  % PE: IgP, PEVK, Titin
  %%
  idxIgpPevkForceLength = 5;
  idx         =idxIgpPevkForceLength;
  

  %from createFiberActiveForceLength.m
  
  lambdaECM = muscleModel.sarcomere.extraCellularMatrixPassiveForceFraction;  

%Structurally distinct segments of titin
  
  ligpZero  = muscleModel.curves.forceLengthIgPTitinCurve.xEnd(1,1);
  lpevkZero = muscleModel.curves.forceLengthPevkTitinCurve.xEnd(1,1);
  
  
  ligpOpt  = muscleModel.sarcomere.IGPNormLengthAtOptimalFiberLength;
  lpevkOpt = muscleModel.sarcomere.PEVKNormLengthAtOptimalFiberLength;  
  ligdOpt  = muscleModel.sarcomere.IGDFreeNormLengthAtOptimalFiberLength;

  ligpFiso  = calcBezierYFcnXDerivative((1-lambdaECM),...
              muscleModel.curves.('forceLengthIgPTitinInverseCurve'),...
              0);
            
  lpevkFiso = calcBezierYFcnXDerivative((1-lambdaECM),...
              muscleModel.curves.('forceLengthPevkTitinInverseCurve'),...
              0);
  
  normLengthZToT12 = ...
    muscleModel.sarcomere.ZLineToT12NormLengthAtOptimalFiberLength;   
  
  normLengthIgdFixed = ...
    muscleModel.sarcomere.IGDFixedNormLengthAtOptimalFiberLength;    
  
  ltitinOpt = normLengthZToT12+normLengthIgdFixed + ligpOpt+lpevkOpt;
  ltitinFiso = normLengthZToT12+normLengthIgdFixed + ligpFiso+lpevkFiso;
    
%Segments of titin from the Z-line to the PEVK-actin attachment point and
%from the PEVK-actin attachment point to myosin:

normLengthProximalTitinAtOptimalFiberLengthHuman=0;
normLengthDistalTitinAtOptimalFiberLengthHuman=0;

normLengthProximalTitinAtOptimalFiberLengthAnimal=0;
normLengthDistalTitinAtOptimalFiberLengthAnimal=0;

assert(humanSoleusMuscle.sarcomere.titinModelType ...
      == muscleModel.sarcomere.titinModelType,...
      ['Error: the muscleModel and humanSoleusModel need to have' ...
      'the same titinModelType']);


if(muscleModel.sarcomere.titinModelType==0)

    %Input model
    n = muscleModel.sarcomere.normPevkToActinAttachmentPoint;
                          
    normLengthProximalTitinAtOptimalFiberLengthAnimal = ...
          muscleModel.sarcomere.IGPNormLengthAtOptimalFiberLength...
    + (n)*muscleModel.sarcomere.PEVKNormLengthAtOptimalFiberLength;

    normLengthDistalTitinAtOptimalFiberLengthAnimal = ...
        + (1-n)*muscleModel.sarcomere.PEVKNormLengthAtOptimalFiberLength ...
        + muscleModel.sarcomere.IGDFreeNormLengthAtOptimalFiberLength;
 
    %Human
    n = humanSoleusMuscle.sarcomere.normPevkToActinAttachmentPoint;
                          
    normLengthProximalTitinAtOptimalFiberLengthHuman = ...
          humanSoleusMuscle.sarcomere.IGPNormLengthAtOptimalFiberLength...
    + (n)*humanSoleusMuscle.sarcomere.PEVKNormLengthAtOptimalFiberLength;

    normLengthDistalTitinAtOptimalFiberLengthHuman = ...
        + (1-n)*humanSoleusMuscle.sarcomere.PEVKNormLengthAtOptimalFiberLength ...
        + humanSoleusMuscle.sarcomere.IGDFreeNormLengthAtOptimalFiberLength;

end
if(muscleModel.sarcomere.titinModelType==1)
    %Input Model                        
    normLengthProximalTitinAtOptimalFiberLengthAnimal = ...
          muscleModel.sarcomere.IGPNormLengthAtOptimalFiberLength...
        + muscleModel.sarcomere.IGDFreeNormLengthAtOptimalFiberLength;

    normLengthDistalTitinAtOptimalFiberLengthAnimal = ...
        + muscleModel.sarcomere.PEVKNormLengthAtOptimalFiberLength;
    
    %Human
    normLengthProximalTitinAtOptimalFiberLengthHuman = ...
          humanSoleusMuscle.sarcomere.IGPNormLengthAtOptimalFiberLength...
        + humanSoleusMuscle.sarcomere.IGDFreeNormLengthAtOptimalFiberLength;

    normLengthDistalTitinAtOptimalFiberLengthHuman = ...
        + humanSoleusMuscle.sarcomere.PEVKNormLengthAtOptimalFiberLength;

end

  lPZero = muscleModel.curves.forceLengthProximalTitinCurve.xEnd(1,1);
  lDZero = muscleModel.curves.forceLengthDistalTitinCurve.xEnd(1,1);
  
  
  lPOpt = normLengthProximalTitinAtOptimalFiberLengthAnimal;

  lDOpt = normLengthDistalTitinAtOptimalFiberLengthAnimal;

  lPFiso = calcBezierYFcnXDerivative((1-lambdaECM),...
              muscleModel.curves.('forceLengthProximalTitinInverseCurve'),...
              0);
            
  lDFiso = calcBezierYFcnXDerivative((1-lambdaECM),...
              muscleModel.curves.('forceLengthDistalTitinInverseCurve'),...
              0);

  fpeDomain = muscleModel.curves.('fiberForceLengthCurve').xEnd ...
             +[-0.01,0.01]+[0,1];
  curveSampleFPE = calcBezierYFcnXCurveSampleVector(...
                muscleModel.curves.('fiberForceLengthCurve'), ...
                200,fpeDomain);
%%
% Titin segment
%%

  %Now go and evaluate the ECM and the titin elements at all of the equivalent
  %lengths
  z0 = zeros(size(curveSampleFPE.x));
  
  
[ curveSampleECMHalf,...
  curveSampleTitin,...
  curveSampleTitinActive,...
  curveSampleIgp,...
  curveSamplePevk,...
  curveSampleIgd,...
  curveSampleProximalTitin,...
  curveSampleDistalTitin,...
  curveSampleTwoSegmentTitinActive] = ...
        sampleTitinCurves20250612(curveSampleFPE.x.*0.5,...
                                  muscleModel); 

  lengthZ2Igp  = muscleModel.sarcomere.ZLineToT12NormLengthAtOptimalFiberLength ...
                                         .*ones(size(curveSampleFPE.x));
  lengthZ2Pevk = lengthZ2Igp + curveSampleIgp.x;
  lengthZ2Igd  = lengthZ2Pevk + curveSamplePevk.x;     
  lengthZ2Myosin= lengthZ2Igd + curveSampleIgd.x;  
  lengthCE     = curveSampleFPE.x.*0.5;                  
                  
[ curveSampleECMHalfHuman,...
  curveSampleTitinHuman,...
  curveSampleTitinActiveHuman,...
  curveSampleIgpHuman,...
  curveSamplePevkHuman,...
  curveSampleIgdHuman,...
  curveSampleProximalTitinHuman,...
  curveSampleDistalTitinHuman,...  
  curveSampleTwoSegmentTitinActiveHuman] = ...
  sampleTitinCurves20250612(curveSampleFPE.x.*0.5,...
                    humanSoleusMuscle);   
                  
  lengthZ2IgpHuman  = ...
      humanSoleusMuscle.sarcomere.ZLineToT12NormLengthAtOptimalFiberLength...
      .*ones(size(curveSampleFPE.x));
  lengthZ2PevkHuman     = lengthZ2IgpHuman + curveSampleIgpHuman.x;
  lengthZ2IgdHuman      = lengthZ2PevkHuman + curveSamplePevkHuman.x;  
  lengthZ2MyosinHuman   = lengthZ2IgdHuman + curveSampleIgdHuman.x;
  lengthCEHuman     = curveSampleFPE.x.*0.5;    


%%
% Plot configuration
%%

  xposXTickIgpPevkFL = [lPOpt,lPFiso,lDFiso,lDOpt].*2;
  textXTickIgpPevkFL = {'$$\ell^{1}_{o}$$','$$\ell^{2}_{o}$$',...
                   '$$\ell^{1}_{1-\lambda}$$','$$\ell^{2}_{1-\lambda}$$',...
                   '$$\ell^{T12}+\ell^{1}(f)+\ell^{2}(f)$$'}; 

%  xposXTickIgpPevkFL = [ligpOpt,ligpFiso,lpevkFiso,ltitinOpt].*2;
%  textXTickIgpPevkFL = {'$$\ell^{Igp}_{o}$$','$$\ell^{PEVK}_{o}$$',...
%                   '$$\ell^{Igp}_{1-\lambda}$$','$$\ell^{PEVK}_{1-\lambda}$$',...
%                   '$$\ell^{T12}+\ell^{Igp}(f)+\ell^{PEVK}(f)+\ell^{Igd}$$'};  

  plotProps(idx).xlim = ([0,max(lPFiso,lDFiso)] + [-0.01,0.01]).*scaleLength;  
  plotProps(idx).ylim = [0,(1-lambdaECM)] + [-0.01,0.01]*(1-lambdaECM);
  plotProps(idx).domain      = plotProps(idx).xlim + [-0.01,0.01];
                 
                 
  plotProps(idx).xticks      = [0,lDFiso, lPFiso ].*scaleLength;
  plotProps(idx).xticklabels = {num2str(round(plotProps(idx).xticks(1,1),2)),...
                                num2str(round(plotProps(idx).xticks(1,2),2)),...
                                num2str(round(plotProps(idx).xticks(1,3),2))};
                              
  
                              
  plotProps(idx).yticks      = [0,(1-lambdaECM)];
  plotProps(idx).yticklabels = {'0',['$$',num2str(round(1-lambdaECM,2)),'f^M_o$$']};
  plotProps(idx).lineColor   = [0,0,0];
  plotProps(idx).lineWidth   = 0.5;
%  if(abs(scaleLength-1)<1e-3)
    plotProps(idx).xlabel = 'Normalized Length ($$\ell/\ell^{M}_{o}$$)';
%   else
%     plotProps(idx).xlabel = 'Length ($$\mu$$m)';      
%   end
  plotProps(idx).ylabel = 'Normalized Force ($$f/f^{M}_{o}$$)';
  plotProps(idx).title  = {'C. Titin force-length','relations'};  
    
  

  %%
  % PE: ECM + titin force-length
  %%
  idxTrombitasFigure5 = 6;
  idx   =idxTrombitasFigure5;
  
  lopt = muscleModel.sarcomere.optimalSarcomereLength;
  loptHuman=humanSoleusMuscle.sarcomere.optimalSarcomereLength;%um
  %from createFiberActiveForceLength.m
  plotProps(idx).xlim = [2.25,4.75];  
  plotProps(idx).ylim = [0,1.5];
  
  plotProps(idx).domain      = plotProps(idx).xlim + [-0.01,0.01];
  plotProps(idx).xticks      = [2.25:0.5:4.75];
  plotProps(idx).xticklabels = cell(1,length(plotProps(idx).xticks),1);
  for z=1:1:length(plotProps(idx).xticks)
    plotProps(idx).xticklabels{z} = num2str(round(plotProps(idx).xticks(1,z),3));
  end
  plotProps(idx).yticks      = [0:0.2:1.2];
  plotProps(idx).yticklabels = cell(1,length(plotProps(idx).yticks));
  for z=1:1:length(plotProps(idx).yticks)
    plotProps(idx).yticklabels{z} = num2str(round(plotProps(idx).yticks(1,z),3));
  end
  
  plotProps(idx).lineColor   = [0,0,0];
  plotProps(idx).lineWidth   = 0.5;
  plotProps(idx).xlabel = 'Sarcomere Length ($$\mu m$$)';
  plotProps(idx).ylabel = 'epitope to Z-line ($$\mu m$$)';
  plotProps(idx).title  = {'A. Titin segment elongation model','vs. Trombitas et al. 1998 Fig. 5'};  
    
  

  
  figure(fig_pubTitinCurves);
  
  %%
  % Igp + PEVK force-length
  %%    
  figure(fig_pubTitinCurves);  
  subplotTrombitasFigure5 = reshape(subPlotPanel(1,1,:),1,4);
  subplot('Position',subplotTrombitasFigure5);  
  idx=idxTrombitasFigure5; 

  fillx = [(lengthCE.*2).*lopt;...
           ((fliplr(lengthCE')').*2).*lopt;...
           (lengthCE(1).*2).*lopt];
  filly = [(lengthZ2Myosin).*lopt; ...
           ((fliplr(lengthCE')')).*lopt;...
           (lengthZ2Myosin(1)).*lopt];

  fill(fillx,filly,[cs.blue],'EdgeColor','none','FaceAlpha',0.25,...
       'DisplayName','Myosin (rat)');

  hold on;

  fillx = [(lengthCEHuman.*2).*loptHuman;...
           ((fliplr(lengthCEHuman')').*2).*loptHuman;...
           (lengthCEHuman(1).*2).*loptHuman];
  filly = [(lengthZ2MyosinHuman).*loptHuman; ...
           ((fliplr(lengthCEHuman')')).*loptHuman;...
           (lengthZ2MyosinHuman(1)).*loptHuman];

  fill(fillx,filly,[0,0,0],'EdgeColor','none','FaceAlpha',0.25,...
       'DisplayName','Myosin (human)');
  hold on;


  hVis ='off';


  plot([1,1].*loptHuman,[0,1.0],'-k','HandleVisibility','off');
  hold on;
  text(loptHuman,1.0,'$$\ell^M_o$$ (human)',...
            'HorizontalAlignment','left',...
            'VerticalAlignment','bottom',...
            'FontSize',6);
  hold on;

  plot([1,1].*lopt,[0,1.2],'-','Color',cs.blue,'HandleVisibility','off');
  hold on;
  text(lopt,1.2,'$$\ell^M_o$$ (rat)',...
            'HorizontalAlignment','left',...
            'VerticalAlignment','bottom',...
            'FontSize',6,...
            'Color',cs.blue);
  hold on;
  


  plot( dataTrombitas1998Figure5(2).x,...
        dataTrombitas1998Figure5(2).y,...
        'o','MarkerSize',3, ...
        'MarkerFaceColor',[1,1,1].*0.9,...
        'Color',[1,1,1].*0.75,...
        'DisplayName','Exp: Z-line (ZL) to N-end 9D10',...
        'HandleVisibility',hVis);
  hold on;

  plot( dataTrombitas1998Figure5(1).x,...
        dataTrombitas1998Figure5(1).y,...
        'o','MarkerSize',3, ...
        'MarkerFaceColor',[1,1,1].*0.5,...
        'Color',[1,1,1].*0.5,...
        'DisplayName','Exp: ZL to C-end 9D10',...
        'HandleVisibility',hVis);
  hold on;

  plot(muscleModel.sarcomere.extra.lineHumanZToPevkP_raw.x,...
       muscleModel.sarcomere.extra.lineHumanZToPevkP_raw.y,...
       '-',...
       'Color',[1,0,1],...
       'LineWidth',plotProps(idx).lineWidth*2,...
       'HandleVisibility','off');
  hold on;

  plot(muscleModel.sarcomere.extra.lineHumanZToPevkD_raw.x,...
       muscleModel.sarcomere.extra.lineHumanZToPevkD_raw.y,...
       '-',...
       'Color',[1,0,1],...
       'LineWidth',plotProps(idx).lineWidth*2,...
       'HandleVisibility','off');
  hold on;

  
  
  
  plot( (lengthCEHuman.*2).*loptHuman,...
        (lengthZ2PevkHuman).*loptHuman,...
        '-','Color',[1,1,1],...
        'LineWidth',plotProps(idx).lineWidth*3,...
        'HandleVisibility','off');      
  hold on;

  plot( (lengthCEHuman.*2).*loptHuman,...
        (lengthZ2PevkHuman).*loptHuman,...
        '-','Color',[0,0,0],...
        'LineWidth',plotProps(idx).lineWidth*2,...
        'DisplayName','Human Sol. Model: ZL to IgP/PEVK',...
        'HandleVisibility',hVis);      
  hold on;
          
  plot( (lengthCEHuman.*2).*loptHuman,...
        (lengthZ2IgdHuman).*loptHuman,...
        '-','Color',[1,1,1],...
        'LineWidth',plotProps(idx).lineWidth,...
        'HandleVisibility','off');

 hold on;

 plot( (lengthCEHuman.*2).*loptHuman,...
        (lengthZ2IgdHuman).*loptHuman,...
        '-','Color',[0,0,0],...
        'LineWidth',plotProps(idx).lineWidth,...
        'DisplayName','Human Sol. Model: ZL to PEVK/IgD',...
        'HandleVisibility',hVis);

 hold on;



 plot( (lengthCEHuman.*2).*loptHuman,...
        (lengthZ2MyosinHuman).*loptHuman,...
        '--','Color',[0,0,0],...
        'LineWidth',plotProps(idx).lineWidth,...
        'DisplayName','Human Sol. Model: ZL to IgD/Myosin',...
        'HandleVisibility',hVis);

 hold on;

 plot( (lengthCEHuman.*2).*loptHuman,...
        (lengthCEHuman).*loptHuman,...
        '--','Color',[0,0,0],...
        'LineWidth',plotProps(idx).lineWidth,...
        'DisplayName','Human Sol. Model: ZL to M-line',...
        'HandleVisibility',hVis);

 hold on; 

 

  plot( (lengthCE.*2).*lopt,...
        (lengthZ2Pevk).*lopt,...
        '-','Color',[1,1,1],...
        'LineWidth',plotProps(idx).lineWidth*3,...
        'DisplayName','Rat Model: ZL to IgP/PEVK',...
        'HandleVisibility','off');     
  hold on;


  plot( (lengthCE.*2).*lopt,...
        (lengthZ2Pevk).*lopt,...
        '-','Color',cs.blue,...
        'LineWidth',plotProps(idx).lineWidth*2,...
        'DisplayName','Rat Model: ZL to IgP/PEVK',...
        'HandleVisibility',hVis);     
  hold on;

  plot( (lengthCE.*2).*lopt,...
        (lengthZ2Igd).*lopt,...
        '-','Color',[1,1,1],...
        'LineWidth',plotProps(idx).lineWidth*2,...
        'HandleVisibility','off');
      
  hold on;        
  plot( (lengthCE.*2).*lopt,...
        (lengthZ2Igd).*lopt,...
        '-','Color',cs.blue,...
        'LineWidth',plotProps(idx).lineWidth,...
        'DisplayName','Rat Model: ZL to PEVK/IgD',...
        'HandleVisibility',hVis);
      
  hold on;


  plot( (lengthCE.*2).*lopt,...
        (lengthZ2Myosin).*lopt,...
        '--','Color',cs.blue,...
        'LineWidth',plotProps(idx).lineWidth,...
        'DisplayName','Rat Model: ZL to PEVK/IgD',...
        'HandleVisibility',hVis);      
  hold on;      

  plot( (lengthCE.*2).*lopt,...
        (lengthCE).*lopt,...
        '--','Color',cs.blue,...
        'LineWidth',plotProps(idx).lineWidth,...
        'DisplayName','Rat Model: ZL to M-line',...
        'HandleVisibility',hVis);      
    hold on;        

  legend('Position',(subplotTrombitasFigure5 + [0,-0.5,0,0]));
  legend boxoff;

  
  
  %text(3.25,0,'Trombitas et al. 1998',...
  %     'FontSize',8,...
  %     'HorizontalJustification','left',...
  %     'Vertical Justification','top');
  %hold on;
  
 %   lengthZ2Igp(i,1) = normLengthZToT12;
 %   lengthZ2Pevk(i,1) = lengthZ2IgpP(i,1) + ligpH;
 %   lengthZ2Igd(i,1) = lengthZ2Pevk(i,1) + lpevkH;                  
 %   lengthCE(i,1)     = xH;  
  box off

  xticks(plotProps(idx).xticks);
  xticklabels(plotProps(idx).xticklabels);

  yticks(plotProps(idx).yticks);
  yticklabels(plotProps(idx).yticklabels);

  xlim(plotProps(idx).xlim);
  ylim(plotProps(idx).ylim);

  xlabel(plotProps(idx).xlabel);
  ylabel(plotProps(idx).ylabel);
  title(plotProps(idx).title);   
          



  %%
  % Passive curve decomposition: total, ecm, titin
  %%    
  figure(fig_pubTitinCurves);  
  subplotPeEcmTitin = reshape(subPlotPanel(1,2,:),1,4);
  subplot('Position',subplotPeEcmTitin);  
  idx=idxPeEcmTitinLength;
  

%   fill([curveSampleECMHalf.x;fliplr(curveSampleECMHalf.x')'].*(2*scaleLength), ...
%        [zeros(size(curveSampleECMHalf.y)); fliplr(curveSampleECMHalf.y')'],...
%         [0.9,0.9,1],'EdgeColor','none','HandleVisibility','off');
%   hold on;  

%   plot( curveSampleECMHalf.x.*(2*scaleLength), ...
%         curveSampleECMHalf.y,...
%         '-','Color',[0,0,1],...
%         'LineWidth',plotProps(idx).lineWidth*2,...
%         'HandleVisibility','off');
%   hold on;    
 
  

  
%   fill([curveSampleTitin.x;fliplr(curveSampleTitin.x')'].*(2*scaleLength), ...
%          [curveSampleECMHalf.y; fliplr(curveSampleTitin.y'+curveSampleECMHalf.y')'],...
%           colorTitinLight,'EdgeColor','none');
%   hold on;  

  curveSamplePassiveTitin.x = curveSampleTitin.x.*(2);
  curveSamplePassiveTitin.y = curveSampleTitin.y+curveSampleECMHalf.y;

  dydxTitin = calcCentralDifferenceDataSeries(...
                curveSamplePassiveTitin.x,...
                curveSamplePassiveTitin.y);


  plot( curveSampleTitin.x.*(2), ...
        curveSampleTitin.y+curveSampleECMHalf.y,...
        '-','Color',plotProps(idx).lineColor,...
        'LineWidth',plotProps(idx).lineWidth,...
        'DisplayName','Titin model (passive)');
  hold on;     
  
  
%       plot( curveSampleFPE.x, ...
%             curveSampleFPE.y,...
%             '--','Color',[0,0,0],...
%             'LineWidth',plotProps(idx).lineWidth);
%       hold on;   

    for i=1:1:length(passiveForceLengthData)
      displayMod='';
      labelShort = passiveForceLengthData(i).label;
      idxLabel = strfind(labelShort,' ');
      labelShort = labelShort(1:idxLabel(1));
      plotColor = passiveForceLengthData(i).color;
      if(contains(passiveForceLengthData(i).label,'TRSS2017'))
        displayMod=' (fit)';
        plotColor=[0,0,0];
      end

        plot(passiveForceLengthData(i).x./scaleLength,...
             passiveForceLengthData(i).y,...
             passiveForceLengthData(i).mark,...
             'Color',plotColor,...
             'MarkerFaceColor',plotColor,...
             'MarkerSize',2,...
             'LineWidth',plotProps(idx).lineWidth,...
             'DisplayName',[labelShort,displayMod]);
        hold on;

        if(contains(passiveForceLengthData(i).label,'SW1982'))
          displayMod=' (fit)';
          plotColor = [0,0,0];
          idxFit = find(passiveForceLengthData(i).y > 0.3);

          plot(passiveForceLengthData(i).x(idxFit)./scaleLength,...
               passiveForceLengthData(i).y(idxFit),...
               passiveForceLengthData(i).mark,...
               'Color',plotColor,...
               'MarkerFaceColor',plotColor,...
               'MarkerSize',2,...
               'LineWidth',plotProps(idx).lineWidth,...
               'DisplayName',[labelShort,displayMod]);
          hold on;          

        end
        

    end

  flag_plotActiveTitin=0;
  
  if(flag_plotActiveTitin==1)
    plot( curveSampleTwoSegmentTitinActive.x.*(2),...
          curveSampleTwoSegmentTitinActive.y+curveSampleECMHalf.y,...
          '-','Color',colorTitin,...
          'LineWidth',plotProps(idx).lineWidth);
    hold on;   
    
    yTxt = 1;
    xTxt = interp1(curveSampleTwoSegmentTitinActive.y+curveSampleECMHalf.y,...
                                curveSampleTwoSegmentTitinActive.x*(2), yTxt);
    hTa =  text(xTxt+0.1,yTxt,'Titin (active)',...
         'FontSize',6,...
         'HorizontalAlignment','right',...
         'VerticalAlignment','bottom');
    hold on
    angle = atan(curveSampleTwoSegmentTitinActive.dydx(end)...
                 +curveSampleECMHalf.dydx(end))*(180/pi);             
    set(hTa,'Rotation',angle+labelRotationOffset*0.7);
  end


  

  lambdaECM = muscleModel.sarcomere.extraCellularMatrixPassiveForceFraction;
  
  if(lambdaECM < 1 && lambdaECM > 0)
    text(normFiberLengthAtOneNormPassiveForce,...
         0,...
         ['$$\lambda^{ECM}=',sprintf('%1.2f',lambdaECM),'$$'],...
         'FontSize',8,...
         'HorizontalAlignment','right',...
         'VerticalAlignment','bottom');
    hold on;
  end
      
%       plot(passiveForceLengthData(:,1),...
%            passiveForceLengthData(:,2),...
%            '.','Color',[0,0,1],'MarkerSize',5,...
%            'LineWidth',plotProps(idx).lineWidth);
%       hold on;  

  
  
  

  if(flag_fibrilModel==0)
      yTxt = lambdaECM;
      yErr = Inf;
      z=length(curveSampleECMHalf.y);
      while z > 2 && ...
          ((curveSampleECMHalf.y(z,1)-yTxt)*(curveSampleECMHalf.y(z-1,1)-yTxt) > 0)
        z=z-1;
      end      
      xTxt = curveSampleECMHalf.x(z,1);
      yTxt = curveSampleECMHalf.y(z,1);
      angle = atan(curveSampleECMHalf.dydx(z,1))*(180/pi);
      
      hecm= text(xTxt*2, yTxt,...
           'ECM          ',...
           'FontSize',8,...
           'HorizontalAlignment','right',...
           'VerticalAlignment','top');
      hold on;
      set(hecm,'Rotation',angle+labelRotationOffset*2.5);
      hold on;
  end

  yTxt = 0.5;
  yErr = Inf;
  z=length(curveSampleTitin.y);
  while z > 2 && ...
      ((curveSampleTitin.y(z,1)-yTxt)*(curveSampleTitin.y(z-1,1)-yTxt) > 0)
    z=z-1;
  end  
  xTxt = curveSampleTitin.x(z,1);
  yTxt = curveSampleTitin.y(z,1);

  hti = text(xTxt*2, yTxt,...
     '$$\mathbf{f}^{PE}(\ell^M)$$',...
     'FontSize',6,...
     'HorizontalAlignment','right',...
     'VerticalAlignment','bottom');
  hold on;

  yTxt = 1;
  yErr = Inf;
  z=length(curveSampleTitin.y);
  while z > 2 && ...
      ((curveSampleTitin.y(z,1)-yTxt)*(curveSampleTitin.y(z-1,1)-yTxt) > 0)
    z=z-1;
  end  
  xTxt = curveSampleTitin.x(z,1);
  yTxt = curveSampleTitin.y(z,1);

  hti = text(xTxt*2-0.05, yTxt,...
     sprintf('%s=%1.2f%s',...
      '$$\mathbf{k}^{PE}',...
      dydxTitin(end),...
      'f^M_o/\ell^M_o$$'),...
     'FontSize',6,...
     'HorizontalAlignment','right',...
     'VerticalAlignment','top');
  hold on;
   
    
  box off

  legend('Position',(subplotPeEcmTitin + [0,-0.5,0,0]));
  legend boxoff;  
  
  
  xticks(plotProps(idx).xticks);
  xticklabels(plotProps(idx).xticklabels);

  yticks(plotProps(idx).yticks);
  yticklabels(plotProps(idx).yticklabels);

  xlim(plotProps(idx).xlim);
  ylim(plotProps(idx).ylim);

  box off;

  xlabel(plotProps(idx).xlabel);
  ylabel(plotProps(idx).ylabel);
  title(plotProps(idx).title);      
  

  %%
  % Igp + PEVK force-length
  %%    
  figure(fig_pubTitinCurves);  
  subplotIgpPevk = reshape(subPlotPanel(1,3,:),1,4);
  subplot('Position',subplotIgpPevk);  
  idx=idxIgpPevkForceLength;


%   curveSampleECMHalf,...
%   curveSampleTitin,...
%   curveSampleTitinActive,...
%   curveSampleIgp,...
%   curveSamplePevk,...
%   curveSampleIgd,...
%   curveSampleProximalTitin,...
%   curveSampleDistalTitin,...
%   curveSampleTwoSegmentTitinActive  
  colorIgdUpd   = [0,0,1];
  colorPEVKUpd  = [1,0,0];
  colorIgpUpd   = [1,0,1];

  plot(curveSampleIgp.x,...
       curveSampleIgp.y,...
       '-','Color',colorIgpUpd,...
       'LineWidth',plotProps(idx).lineWidth,...
       'DisplayName','$$\mathbf{f}^{IgP}$$');
  hold on;

  xTitinSegTicks = [0];

  yLabel= 1;
  xLabel = interp1(curveSampleIgp.y,curveSampleIgp.x,yLabel);
  xTitinSegTicks = [xLabel;xTitinSegTicks];

  text(xLabel,...
       yLabel,...
       sprintf('%s=%1.2f%s',...
      '$$\mathbf{k}^{IgP}',...
      muscleModel.curves.forceLengthIgPTitinCurve.dydxEnd(2),'$$'),...
       'FontSize',6,...
       'Color',colorIgpUpd,...
       'VerticalAlignment','top',...
       'HorizontalAlignment','right');
  hold on;

  yLabel=0.5;
  xLabel = interp1(curveSampleIgp.y,curveSampleIgp.x,yLabel);
  
  text(xLabel,...
       yLabel,...
       '$$\mathbf{f}^{IgP}$$',...
       'FontSize',6,...
       'Color',colorIgpUpd,...
       'VerticalAlignment','top',...
       'HorizontalAlignment','right');
  hold on;



  plot(curveSamplePevk.x,...
       curveSamplePevk.y,...
       '-','Color',colorPEVKUpd,...
       'LineWidth',plotProps(idx).lineWidth,...
       'DisplayName','$$\mathbf{f}^{PEVK}$$');
  hold on;

  yLabel = 1;
  xLabel = interp1(curveSamplePevk.y,curveSamplePevk.x,yLabel);
  xTitinSegTicks = [xLabel;xTitinSegTicks];  
  xPlotMax = xLabel;

  text(xLabel,...
       yLabel,...
       sprintf('%s=%1.2f%s',...
      '$$\mathbf{k}^{PEVK}',...
      muscleModel.curves.forceLengthPevkTitinCurve.dydxEnd(2),'$$'),...
       'FontSize',6,...
       'Color',colorPEVKUpd,...
       'VerticalAlignment','top',...
       'HorizontalAlignment','left');
  hold on;

  yLabel = 0.5;
  xLabel = interp1(curveSamplePevk.y,curveSamplePevk.x,yLabel);


  text(xLabel,...
       yLabel,...  
       '$$\mathbf{f}^{PEVK}$$',...
       'FontSize',6,...
       'Color',colorPEVKUpd,...
       'VerticalAlignment','top',...
       'HorizontalAlignment','right'); 
  hold on;


  plot(curveSampleIgd.x,...
       curveSampleIgd.y,...
       '-','Color',colorIgdUpd,...
       'LineWidth',plotProps(idx).lineWidth,...
       'DisplayName','$$\mathbf{f}^{IgD}$$');
  hold on;
 
  yLabel = 1;
  xLabel = interp1(curveSampleIgd.y,curveSampleIgd.x,yLabel);  
  xTitinSegTicks = [xLabel;xTitinSegTicks];

  text(xLabel,...
       yLabel,...
       sprintf('%s=%1.2f%s',...
      '$$\mathbf{k}^{IgD}',...
      muscleModel.curves.forceLengthIgDTitinCurve.dydxEnd(2),'$$'),...
       'FontSize',6,...
       'Color',colorIgdUpd,...
       'VerticalAlignment','top',...
       'HorizontalAlignment','right');
  hold on; 

  yLabel = 0.5;
  xLabel = interp1(curveSampleIgd.y,curveSampleIgd.x,yLabel);  

  text(xLabel,...
       yLabel,...
       '$$\mathbf{f}^{IgD}$$',...
       'FontSize',6,...
       'Color',colorIgdUpd,...
       'VerticalAlignment','top',...
       'HorizontalAlignment','right'); 
  hold on;
  box off

  legend('Position',(subplotIgpPevk + [0,-0.5,0,0]));
  legend boxoff;  
  

  xticks( round(sort(xTitinSegTicks),2) ) ;  

  yticks([0,1]);
  yticklabels({'0','$$f^M_o$$'});

  xlim([0,(xPlotMax+0.01)]);
  ylim([0,1.01]);

  box off;

  xlabel(plotProps(idx).xlabel);
  ylabel(plotProps(idx).ylabel);
  title(plotProps(idx).title);   



  
  
  here=1;

  figure(fig_pubTitinCurves);  
  configPlotExporter;
  print('-dpdf', [filePathAndName,'.pdf']); 
  saveas(fig_pubTitinCurves,filePathAndName,'fig');
  success=1;