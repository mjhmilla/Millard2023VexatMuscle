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

function [success] = plotMuscleCurves2025( ...
                      muscleModel,...
                      humanSoleusMuscle,...                      
                      activeForceLengthCurveAnnotationPoints,...
                      activeForceLengthData,...
                      passiveForceLengthData,...   
                      normFiberLengthAtOneNormPassiveForce,...
                      trialId,...
                      updateTitinPlotsOnly,...
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
%%
% Series colors
%%
  n = 1-((trialId-1)/2);
  
  colorTitin = [1,0,1].*n + ([1,0,1].*0.25).*(1-n);
  colorTitinLight = colorTitin.*0.1 + [1,1,1].*0.9;  

  colorIgp   = [1,0,0].*n + ([1,0,0].*0.25).*(1-n);
  colorPEVK  = [0.5,0,0].*n + ([0.5,0,0].*0.25).*(1-n);
  colorIgd   = [0.25,0,0].*n + ([0.25,0,0].*0.25).*(1-n);
  colorDistal = colorTitin;
  colorProximal   = colorDistal.*0.5;


%%
%Plot configuration
%%
  flag_fibrilModel = 0;
  if(muscleModel.sarcomere.extraCellularMatrixPassiveForceFraction==0)
    flag_fibrilModel=1;
  end


  if(updateTitinPlotsOnly==1)
    fig_pubCurves=openfig( [previousFilePathAndName,'.fig'],'new','visible');
  else
    fig_pubCurves=figure;
  end
  
  numberOfHorizontalPlotColumns = 3;
  numberOfVerticalPlotRows      = 2;

  pageWidth         = 21;
  pageHeight        = 29.7;
  plotWidth         = 4.5;
  plotHeight        = 4.5;
  plotHorizMarginCm = 1.5;
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
  fileTrombitas1998Figure5 = fullfile(projectFolders.experiments_TGFG1998, 'Trombitas1998_Figure5.csv');

  dataTrombitas1998Figure5 = loadDigitizedData(fileTrombitas1998Figure5,...
                        'Sarcomere Length','PEVK Width (um)',...
                        {'a','b'},'Trombitas 1998 Figure 5');  
                      
  fileTrombitas1998Figure6 = fullfile(projectFolders.experiments_TGFG1998, 'Trombitas1998_Figure6.csv');;
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
      muscleModel.curves.activeForceLengthCurve.xEnd ...
        +[-0.02,0.02];
  
  plotProps(idx).ylim   = [0,1.0]+[-0.01,0.01];
  plotProps(idx).domain = plotProps(idx).xlim + [-0.01,0.01];
  plotProps(idx).xticks = ...
      [ round(activeForceLengthCurveAnnotationPoints.x(1,1),2),...                                
        1,...
        round(muscleModel.curves.fiberForceLengthCurve.xEnd(1,2),2),...
        round(activeForceLengthCurveAnnotationPoints.x(end,1),2)];
  plotProps(idx).xticks = sort(plotProps(idx).xticks);

  plotProps(idx).xticklabels = {num2str(plotProps(idx).xticks(1,1)),...
                                num2str(plotProps(idx).xticks(1,2)),...
                                num2str(plotProps(idx).xticks(1,3)),...
                                num2str(plotProps(idx).xticks(1,4))};
                              
  plotProps(idx).yticks      = [0,1];
  plotProps(idx).yticklabels = {'0','$$f^M_o$$'};
  plotProps(idx).lineColor   = [0,0,0];
  plotProps(idx).lineWidth   = 0.5;
  plotProps(idx).xlabel = 'Normalized Length ($$\ell/\ell^{M}_{o}$$)';
  plotProps(idx).ylabel = 'Normalized Force ($$f/f^{M}_{o}$$)';
  plotProps(idx).title  = 'B. CE force-length ($$\mathbf{f}^{L}\,\&\,\mathbf{f}^{PE}$$)';

  %%
  % CE force-velocity
  idxCEVelocity = 3;
  idx           = idxCEVelocity;

  %from createFiberActiveForceLength.m
  plotProps(idx).xlim = [-1.01,1.01];
  
  plotProps(idx).ylim = [0,muscleModel.curves.fiberForceVelocityCurve.yEnd(1,2)] ...
                        +[-0.01,0.01];
                      
  plotProps(idx).domain      = plotProps(idx).xlim + [-0.01,0.01];
  
  plotProps(idx).xticks      = [-1,0,1];
  plotProps(idx).xticklabels = {'$$-1$$','0','$$1$$'};
                              
  plotProps(idx).yticks      = [0,1,muscleModel.curves.fiberForceVelocityCurve.yEnd(1,2)];
  plotProps(idx).yticklabels = {0,1,num2str(round(plotProps(idx).yticks(1,end),2))};
  plotProps(idx).lineColor   = [0,0,0];
  plotProps(idx).lineWidth   = 0.5;
  plotProps(idx).xlabel = 'Normalized Velocity ($$v/v^{M}_{max}$$)';
  plotProps(idx).ylabel = 'Normalized Force ($$f/f^{M}_{o}$$)';
  plotProps(idx).title  = 'C. CE force-velocity ($$\mathbf{f}^{V}$$)';
  %%

  
  %%
  % PE: ECM + titin force-length
  %%
  idxPeEcmTitinLength = 4;
  idx         =idxPeEcmTitinLength;
  

  %from createFiberActiveForceLength.m
  plotProps(idx).xlim = ...
    muscleModel.curves.fiberForceLengthCurve.xEnd ...
    +[-0.02,0.02];
  
  plotProps(idx).ylim = [0,1.0]+[-0.01,0.01];
  plotProps(idx).domain      = plotProps(idx).xlim + [-0.01,0.01];
  plotProps(idx).xticks      = [round(muscleModel.curves.fiberForceLengthCurve.xEnd(1,1),2),...                                
                                1,...
                                round(muscleModel.curves.fiberForceLengthCurve.xEnd(1,2),2)];
  plotProps(idx).xticklabels = {num2str(plotProps(idx).xticks(1,1)),...
                                num2str(plotProps(idx).xticks(1,2)),...
                                num2str(plotProps(idx).xticks(1,3))};
                              
  plotProps(idx).yticks      = [0,1];
  plotProps(idx).yticklabels = {'0','$$f^M_o$$'};
  plotProps(idx).lineColor   = [0,0,0];
  plotProps(idx).lineWidth   = 0.5;
  plotProps(idx).xlabel = 'Normalized Length ($$\ell/\ell^{M}_{o}$$)';
  plotProps(idx).ylabel = 'Normalized Force ($$f/f^{M}_{o}$$)';
  plotProps(idx).title  = 'A. ECM ($$\mathbf{f}^{ECM}$$) \& Titin ($$\mathbf{f}^{1}+\mathbf{f}^{2}$$)';  
  
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

normLengthProximalTitinAtOptimalFiberLengthFeline=0;
normLengthDistalTitinAtOptimalFiberLengthFeline=0;

assert(humanSoleusMuscle.sarcomere.titinModelType ...
      == muscleModel.sarcomere.titinModelType,...
      ['Error: the muscleModel and humanSoleusModel need to have' ...
      'the same titinModelType']);


if(muscleModel.sarcomere.titinModelType==0)

    %Feline
    n = muscleModel.sarcomere.normPevkToActinAttachmentPoint;
                          
    normLengthProximalTitinAtOptimalFiberLengthFeline = ...
          muscleModel.sarcomere.IGPNormLengthAtOptimalFiberLength...
    + (n)*muscleModel.sarcomere.PEVKNormLengthAtOptimalFiberLength;

    normLengthDistalTitinAtOptimalFiberLengthFeline = ...
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
    %Feline                          
    normLengthProximalTitinAtOptimalFiberLengthFeline = ...
          muscleModel.sarcomere.IGPNormLengthAtOptimalFiberLength...
        + muscleModel.sarcomere.IGDFreeNormLengthAtOptimalFiberLength;

    normLengthDistalTitinAtOptimalFiberLengthFeline = ...
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
  
  
  lPOpt = normLengthProximalTitinAtOptimalFiberLengthFeline;

  lDOpt = normLengthDistalTitinAtOptimalFiberLengthFeline;

  lPFiso = calcBezierYFcnXDerivative((1-lambdaECM),...
              muscleModel.curves.('forceLengthProximalTitinInverseCurve'),...
              0);
            
  lDFiso = calcBezierYFcnXDerivative((1-lambdaECM),...
              muscleModel.curves.('forceLengthDistalTitinInverseCurve'),...
              0);


  xposXTickIgpPevkFL = [lPOpt,lPFiso,lDFiso,lDOpt].*2;
  textXTickIgpPevkFL = {'$$\ell^{1}_{o}$$','$$\ell^{2}_{o}$$',...
                   '$$\ell^{1}_{1-\lambda}$$','$$\ell^{2}_{1-\lambda}$$',...
                   '$$\ell^{T12}+\ell^{1}(f)+\ell^{2}(f)$$'}; 

%  xposXTickIgpPevkFL = [ligpOpt,ligpFiso,lpevkFiso,ltitinOpt].*2;
%  textXTickIgpPevkFL = {'$$\ell^{Igp}_{o}$$','$$\ell^{PEVK}_{o}$$',...
%                   '$$\ell^{Igp}_{1-\lambda}$$','$$\ell^{PEVK}_{1-\lambda}$$',...
%                   '$$\ell^{T12}+\ell^{Igp}(f)+\ell^{PEVK}(f)+\ell^{Igd}$$'};  

  plotProps(idx).xlim = [0,max(lPFiso,lDFiso)] + [-0.01,0.01];  
  plotProps(idx).ylim = [0,(1-lambdaECM)] + [-0.01,0.01]*(1-lambdaECM);
  plotProps(idx).domain      = plotProps(idx).xlim + [-0.01,0.01];
                 
                 
  plotProps(idx).xticks      = [0,lDFiso, lPFiso ];
  plotProps(idx).xticklabels = {num2str(round(plotProps(idx).xticks(1,1),2)),...
                                num2str(round(plotProps(idx).xticks(1,2),2)),...
                                num2str(round(plotProps(idx).xticks(1,3),2))};
                              
  
                              
  plotProps(idx).yticks      = [0,(1-lambdaECM)];
  plotProps(idx).yticklabels = {'0',['$$',num2str(round(1-lambdaECM,2)),'f^M_o$$']};
  plotProps(idx).lineColor   = [0,0,0];
  plotProps(idx).lineWidth   = 0.5;
  plotProps(idx).xlabel = 'Normalized Length ($$\ell/\ell^{M}_{o}$$)';
  plotProps(idx).ylabel = 'Normalized Force ($$f/f^{M}_{o}$$)';
  plotProps(idx).title  = 'B. Titin segments ($$\mathbf{f}^{1}$$ and $$\mathbf{f}^{2}$$)';  
    
  

  %%
  % PE: ECM + titin force-length
  %%
  idxTrombitasFigure5 = 6;
  idx   =idxTrombitasFigure5;
  
  loptHuman=2.725;%um
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
  plotProps(idx).title  = {'C. Human soleus titin segment elongation','(Trombitas et al. 1998 Fig. 5)'};  
    
  
  
  
  %%
  %Tendon
  %%
  figure(fig_pubCurves);
    if(updateTitinPlotsOnly==0)
      if(isfield(muscleModel.curves,'tendonForceLengthCurve'))
          if(isempty(muscleModel.curves.tendonForceLengthCurve)==0)
              subplotTendon = reshape(subPlotPanel(1,1,:),1,4);
              subplot('Position',subplotTendon);
              
                idx=idxTendon;    
                curveSample = calcBezierYFcnXCurveSampleVector(...
                              muscleModel.curves.('tendonForceLengthCurve'), ...
                              200,plotProps(idx).domain);
            
                plot(curveSample.x, ...
                  curveSample.y,...
                  '-','Color',plotProps(idx).lineColor,...
                      'LineWidth',plotProps(idx).lineWidth);
                hold on;
            
                xticks(plotProps(idx).xticks);
                xticklabels(plotProps(idx).xticklabels);
            
                yticks(plotProps(idx).yticks);
                yticklabels(plotProps(idx).yticklabels);
            
                text( muscleModel.curves.tendonForceLengthCurve.xEnd(1,2),...
                      muscleModel.curves.tendonForceLengthCurve.yEnd(1,2),...
                      '$$k^{T}_o$$',...
                      'VerticalAlignment','top',...
                      'HorizontalAlignment','left',...
                      'FontSize',8);
            
                hold on;
            
                plot(muscleModel.curves.tendonForceLengthCurve.xEnd(1,1),...
                      muscleModel.curves.tendonForceLengthCurve.yEnd(1,1),...
                      '.','Color',plotProps(idx).lineColor);
                hold on;    
                
                plot(muscleModel.curves.tendonForceLengthCurve.xEnd(1,2),...
                      muscleModel.curves.tendonForceLengthCurve.yEnd(1,2),...
                      '.','Color',plotProps(idx).lineColor);
                hold on;
            
                xlim(plotProps(idx).xlim);
                ylim(plotProps(idx).ylim);
            
                box off;
            
                xlabel(plotProps(idx).xlabel);
                ylabel(plotProps(idx).ylabel);
                title(plotProps(idx).title);
          end
      end
    end
    
  %%
  %CE force-length
  %%
    if(updateTitinPlotsOnly==0)  
      figure(fig_pubCurves);  
      subplotCELength = reshape(subPlotPanel(1,2,:),1,4);
      subplot('Position',subplotCELength);  
      
        idx=idxCELength;    
        curveSampleFL = calcBezierYFcnXCurveSampleVector(...
                      muscleModel.curves.('activeForceLengthCurve'), ...
                      200,plotProps(idx).domain);
    
        plot(curveSampleFL.x, ...
          curveSampleFL.y,...
          '-','Color',plotProps(idx).lineColor,...
              'LineWidth',plotProps(idx).lineWidth);
        hold on;
    
        curveSampleFPE = calcBezierYFcnXCurveSampleVector(...
                      muscleModel.curves.('fiberForceLengthCurve'), ...
                      200,plotProps(idx).domain);
    
        plot(curveSampleFPE.x, ...
          curveSampleFPE.y,...
          '-','Color',[1,1,1],...
              'LineWidth',plotProps(idx).lineWidth*2);
        hold on;    
        
        plot(curveSampleFPE.x, ...
          curveSampleFPE.y,...
          '--','Color',plotProps(idx).lineColor,...
              'LineWidth',plotProps(idx).lineWidth);
        hold on;    
            
        plot(activeForceLengthData(:,1),...
             activeForceLengthData(:,2),...
             'x','Color',[0,0,1],'MarkerSize',2,...
             'LineWidth',plotProps(idx).lineWidth);
        hold on;
    
        plot(passiveForceLengthData(:,1),...
             passiveForceLengthData(:,2),...
             '.','Color',[0,0,1],'MarkerSize',5,...
             'LineWidth',plotProps(idx).lineWidth);
        hold on;
           
        
        xticks(plotProps(idx).xticks);
        xticklabels(plotProps(idx).xticklabels);
    
        yticks(plotProps(idx).yticks);
        yticklabels(plotProps(idx).yticklabels);
    
        xlim(plotProps(idx).xlim);
        ylim(plotProps(idx).ylim);
    
    %     text( plotProps(idx).xticks(1,3),0.6,...
    %           'Herzog 2002',...
    %           'FontSize',8,...
    %           'Color',[0,0,1]);
    %     
         box off;
    
        xlabel(plotProps(idx).xlabel);
        ylabel(plotProps(idx).ylabel);
        title(plotProps(idx).title);    
    end
  %%
  % Force-velocity
  %%
    if(updateTitinPlotsOnly==0)  
      figure(fig_pubCurves);    
      subplotCEVelocity = reshape(subPlotPanel(1,3,:),1,4);
      subplot('Position',subplotCEVelocity);  
      
        idx=idxCEVelocity;    
        curveSampleFV = calcBezierYFcnXCurveSampleVector(...
                      muscleModel.curves.('fiberForceVelocityCurve'), ...
                      200,plotProps(idx).domain);
    
        plot(curveSampleFV.x, ...
          curveSampleFV.y,...
          '-','Color',plotProps(idx).lineColor,...
              'LineWidth',plotProps(idx).lineWidth);
        hold on;        
        
        plot(0,...
             1,...
             '.','Color',plotProps(idx).lineColor);
        hold on;    
        
        text(0.01,-0.01, 'Lengthening','FontSize',8,...
             'VerticalAlignment','bottom');
        hold on;
        text(-0.01,-0.01, 'Shortening','FontSize',8,...
             'HorizontalAlignment','Right',...
             'VerticalAlignment','bottom');
        hold on;
        
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
    end
  %%
  % Passive curve decomposition: total, ecm, titin
  %%    
  figure(fig_pubCurves);  
  subplotPeEcmTitin = reshape(subPlotPanel(2,1,:),1,4);
  subplot('Position',subplotPeEcmTitin);  
  idx=idxPeEcmTitinLength;
  
  fpeDomain = muscleModel.curves.('fiberForceLengthCurve').xEnd ...
             +[-0.01,0.01]+[0,1];
  curveSampleFPE = calcBezierYFcnXCurveSampleVector(...
                muscleModel.curves.('fiberForceLengthCurve'), ...
                200,fpeDomain);

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
                    muscleModel);   
                  
  lengthZ2IgpHuman  = ...
      humanSoleusMuscle.sarcomere.ZLineToT12NormLengthAtOptimalFiberLength...
      .*ones(size(curveSampleFPE.x));
  lengthZ2PevkHuman     = lengthZ2IgpHuman + curveSampleIgpHuman.x;
  lengthZ2IgdHuman      = lengthZ2PevkHuman + curveSamplePevkHuman.x;  
  lengthZ2MyosinHuman   = lengthZ2IgdHuman + curveSampleIgdHuman.x;
  lengthCEHuman     = curveSampleFPE.x.*0.5;                  
             
  %Here we use the 'B' version of the curves so that we can compare
  %the kinematics of the Igp/PEVK boundary and the PEVK/Igd boundary as
  %was measured by Trombitas et al.
%   lengthZ2IgpHumanB  = ...
%       humanSoleusMuscle.sarcomere.ZLineToT12NormLengthAtOptimalFiberLength ...
%         .*ones(size(curveSampleFPE.x));
% 
%   lengthZ2PevkHumanB = lengthZ2IgpHumanB + curveSampleIgpHumanB.x;
%   lengthZ2IgdHumanB  = lengthZ2PevkHumanB + curveSamplePevkIgdHumanB.x;                  
%   lengthCEHumanB     = curveSampleFPE.x.*0.5;       
%                 


    if(updateTitinPlotsOnly==0)
      fill([curveSampleECMHalf.x;fliplr(curveSampleECMHalf.x')'].*2, ...
           [zeros(size(curveSampleECMHalf.y)); fliplr(curveSampleECMHalf.y')'],...
            [0.9,0.9,1],'EdgeColor','none');
      hold on;  
    
      plot( curveSampleECMHalf.x.*2, ...
            curveSampleECMHalf.y,...
            '-','Color',[0,0,1],...
            'LineWidth',plotProps(idx).lineWidth*2);
      hold on;    
     
      

      
      fill([curveSampleTitin.x;fliplr(curveSampleTitin.x')'].*2, ...
           [curveSampleECMHalf.y; fliplr(curveSampleTitin.y'+curveSampleECMHalf.y')'],...
            colorTitinLight,'EdgeColor','none');
      hold on;  
    
      plot( curveSampleTitin.x.*2, ...
            curveSampleTitin.y+curveSampleECMHalf.y,...
            '-','Color',colorTitin,...
            'LineWidth',plotProps(idx).lineWidth);
      hold on;     
      
      
      plot( curveSampleFPE.x, ...
            curveSampleFPE.y,...
            '--','Color',[0,0,0],...
            'LineWidth',plotProps(idx).lineWidth);
      hold on;    
    end
  
  plot( curveSampleTwoSegmentTitinActive.x.*2,...
        curveSampleTwoSegmentTitinActive.y+curveSampleECMHalf.y,...
        '-','Color',colorTitin,...
        'LineWidth',plotProps(idx).lineWidth);
  hold on;   
  
  if(updateTitinPlotsOnly==0)
      yTxt = 1;
      xTxt = interp1(curveSampleTwoSegmentTitinActive.y+curveSampleECMHalf.y,...
                     curveSampleTwoSegmentTitinActive.x*2, yTxt);
      hTa =  text(xTxt+0.1,yTxt,'Titin (active)',...
           'FontSize',8,...
           'HorizontalAlignment','right',...
           'VerticalAlignment','bottom');
      hold on
      angle = atan(curveSampleTwoSegmentTitinActive.dydx(end)...
                   +curveSampleECMHalf.dydx(end))*(180/pi);             
      set(hTa,'Rotation',angle+labelRotationOffset*0.7);
  end
  if(flag_fibrilModel==0)
      yTxt = 1;
      xTxt = interp1(curveSampleTitin.y+curveSampleECMHalf.y,...
                     curveSampleTitin.x*2, yTxt);
      hTa =  text(xTxt,yTxt,'ECM + Titin (passive)',...
           'FontSize',8,...
           'HorizontalAlignment','right',...
           'VerticalAlignment','bottom');
      hold on
      angle = atan(curveSampleTitin.dydx(end)...
                   +curveSampleECMHalf.dydx(end))*(180/pi);             
      set(hTa,'Rotation',angle+labelRotationOffset*1.7);
  end

  
  if(updateTitinPlotsOnly==0)
      lambdaECM = muscleModel.sarcomere.extraCellularMatrixPassiveForceFraction;
      
      text(normFiberLengthAtOneNormPassiveForce,...
           0,...
           ['$$\lambda^{ECM}=',sprintf('%1.2f',lambdaECM),'$$'],...
           'FontSize',8,...
           'HorizontalAlignment','right',...
           'VerticalAlignment','bottom');
      hold on;
      
      plot(passiveForceLengthData(:,1),...
           passiveForceLengthData(:,2),...
           '.','Color',[0,0,1],'MarkerSize',5,...
           'LineWidth',plotProps(idx).lineWidth);
      hold on;  
  end
  
  
  

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

  yTxt = lambdaECM;
  if(flag_fibrilModel==1)
    yTxt = 1;
  end
  yErr = Inf;
  z=length(curveSampleTitin.y);
  while z > 2 && ...
      ((curveSampleTitin.y(z,1)-yTxt)*(curveSampleTitin.y(z-1,1)-yTxt) > 0)
    z=z-1;
  end  

  xTxt = curveSampleTitin.x(z,1);
  yTxt = curveSampleTitin.y(z,1);
  angle = atan(curveSampleTitin.dydx(z,1))*(180/pi);

  if(updateTitinPlotsOnly==0)
      hti = text(xTxt*2, yTxt,...
         'Titin (passive)',...
         'FontSize',8,...
         'HorizontalAlignment','right',...
         'VerticalAlignment','bottom');
      hold on;
      
      set(hti,'Rotation',angle+labelRotationOffset*2.5);
      hold on;
  end
    
  box off

  
  
  
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
  figure(fig_pubCurves);  
  subplotIgpPevk = reshape(subPlotPanel(2,2,:),1,4);
  subplot('Position',subplotIgpPevk);  
  idx=idxIgpPevkForceLength;
    
  %textXTickIgpPevkFL
  


  plot( curveSampleProximalTitin.x,...
        curveSampleProximalTitin.y,...
        '-','Color',colorProximal,...
        'LineWidth',plotProps(idx).lineWidth);

%  plot( curveSampleIgp.x,...
%        curveSampleIgp.y,...
%        '-','Color',colorIgp,...
%        'LineWidth',plotProps(idx).lineWidth);
  hold on;
  
  angleStretch  = 1.1;
  
  eIsoIgp = ligpFiso/ligpOpt;

  yTxt =  1-(trialId-1)*0.2;%(1-lambdaECM)*(0.5);
  xTxt = interp1( curveSampleProximalTitin.y, curveSampleProximalTitin.x, yTxt);
%  xTxt = xTxt - 0.025.*(max(curveSampleProximalTitin.x)-min(curveSampleProximalTitin.x));

%   hIgp = text(xTxt, yTxt,...
%         ['$$\mathbf{f}^{1}$$'],'FontSize',8,...
%         'HorizontalAlignment','left',...
%         'VerticalAlignment','bottom',...
%         'Color',colorProximal);     
%   hold on;

%   yTxt = (1-lambdaECM)*0.5;
%   xTxt = interp1( curveSampleProximalTitin.y, curveSampleProximalTitin.x, yTxt);  
%   xTxt = xTxt - 0.025.*(max(curveSampleProximalTitin.x)-min(curveSampleProximalTitin.x));

  l1HNo = calcBezierFcnXGivenY((1-lambdaECM),...
            muscleModel.curves.forceLengthProximalTitinCurve);
  kP = calcBezierYFcnXDerivative(l1HNo,...
           muscleModel.curves.forceLengthProximalTitinCurve,1);



  %kP = muscleModel.curves.forceLengthProximalTitinCurve.dydxEnd(1,2);
  hP = text(xTxt,yTxt,...
       ['$$','\mathbf{f}^{1}: \hat{k}^{1}_o=',num2str(round(kP,2)),'$$'],'FontSize',8,...
       'HorizontalAlignment','right',...
       'VerticalAlignment','top',...
       'Color',colorProximal);%,...
%       'Color',colorProximal);      
  hold on;   
  
  %angleIgp = atan(curveSampleIgp.dydx(end))*(180/pi);
  %set(hIgp,'Rotation',angleIgp+labelRotationOffset);
  %hIgp2 = text(curveSampleIgp.x(end).*2,...
  %      curveSampleIgp.y(end),...
  %      ['$$e_{ISO}=\, ',num2str(round(eIsoIgp,2)),'$$'],'FontSize',8,...
  %      'HorizontalAlignment','right',...
  %      'VerticalAlignment','top');           
  %angleIgp2 = atan(curveSampleIgp.dydx(end))*(180/pi);
  %set(hIgp2,'Rotation',angleIgp2*angleStretch);
  
  
  plot( curveSampleDistalTitin.x,...
        curveSampleDistalTitin.y,...
        '-','Color',colorDistal,...
        'LineWidth',plotProps(idx).lineWidth);
  hold on;  

%   yTxt = (1-lambdaECM)*0.8;
%   xTxt = interp1( curveSampleDistalTitin.y, curveSampleDistalTitin.x, yTxt);  
%   xTxt = xTxt + 0.05.*(max(curveSampleDistalTitin.x)-min(curveSampleDistalTitin.x));
  

  %eIsoPevk = lpevkFiso/lpevkOpt;
%   hpevk = text(xTxt,yTxt,...
%        ['$$\mathbf{f}^{2}$$'],'FontSize',8,...
%        'HorizontalAlignment','left',...
%        'VerticalAlignment','bottom',...
%        'Color',colorDistal);      
%   hold on;  
  
  yTxt = 1-(trialId-1)*0.2;
  xTxt = interp1( curveSampleDistalTitin.y, curveSampleDistalTitin.x, yTxt);  
  xTxt = xTxt + 0.05.*(max(curveSampleDistalTitin.x)-min(curveSampleDistalTitin.x));

  l2HNo = calcBezierFcnXGivenY((1-lambdaECM),...
            muscleModel.curves.forceLengthDistalTitinCurve);
  kD = calcBezierYFcnXDerivative(l2HNo,...
           muscleModel.curves.forceLengthDistalTitinCurve,1);


  %kD = muscleModel.curves.forceLengthDistalTitinCurve.dydxEnd(1,2);
  hD = text(xTxt,yTxt,...
       ['$$\mathbf{f}^2: \hat{k}^{2}_o=',num2str(round(kD,2)),'$$'],'FontSize',8,...
       'HorizontalAlignment','right',...
       'VerticalAlignment','top',...
       'Color',colorDistal);      
  hold on;    

  
  box off

  [valXTick, idxXTick] = sort(plotProps(idx).xticks);

  xticks(valXTick);
  xticklabels(plotProps(idx).xticklabels(idxXTick));

  yticks(plotProps(idx).yticks);
  yticklabels(plotProps(idx).yticklabels);

  xlim(plotProps(idx).xlim);
  ylim(plotProps(idx).ylim);

  xlabel(plotProps(idx).xlabel);
  ylabel(plotProps(idx).ylabel);
  title(plotProps(idx).title);   
  
  %%
  % Igp + PEVK force-length
  %%    
  figure(fig_pubCurves);  
  subplotTrombitasFigure5 = reshape(subPlotPanel(2,3,:),1,4);
  subplot('Position',subplotTrombitasFigure5);  
  idx=idxTrombitasFigure5; 

  hVis = 'on';
  if(updateTitinPlotsOnly)
    hVis ='off';
  end

    if(updateTitinPlotsOnly==0)  
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
    end
  
  plot( (lengthCEHuman.*2).*loptHuman,...
        (lengthZ2PevkHuman).*loptHuman,...
        '-','Color',colorIgp,...
        'LineWidth',plotProps(idx).lineWidth,...
        'DisplayName','Model: ZL to IgD/PEVK',...
        'HandleVisibility',hVis);
      
  hold on;
  

        
  plot( (lengthCEHuman.*2).*loptHuman,...
        (lengthZ2IgdHuman).*loptHuman,...
        '-','Color',colorPEVK,...
        'LineWidth',plotProps(idx).lineWidth,...
        'DisplayName','Model: ZL to PEVK/IgD',...
        'HandleVisibility',hVis);
      
        
  plot( (lengthCEHuman.*2).*loptHuman,...
        (lengthZ2MyosinHuman).*loptHuman,...
        '--','Color',colorIgd,...
        'LineWidth',plotProps(idx).lineWidth,...
        'DisplayName','Model: ZL to IgD/Myosin',...
        'HandleVisibility',hVis);
      
  hold on;

  if(updateTitinPlotsOnly==0)
      legend('Location','NorthWest');
      legend boxoff;
  end
  
  
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
   
  
  
  here=1;

  figure(fig_pubCurves);  
  configPlotExporter;
  print('-dpdf', [filePathAndName,'.pdf']); 
  saveas(fig_pubCurves,filePathAndName,'fig');
  success=1;