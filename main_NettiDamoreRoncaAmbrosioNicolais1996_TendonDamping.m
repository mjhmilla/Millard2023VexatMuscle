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

clc;
close all;
clear all;

flag_usingOctave=0;
flag_NormalizeByKmax=1;

flag_addCurveStorageModulus=0;
%The publication plots are making the point that the loss modulus is
%well approximated as a scaled version of the storage modulus. It is 
%a bit of a distraction to show, on top of this, the storage modulus
%of the tendon-force-lenght curve.

rootDir         = getRootProjectDirectory();
projectFolders  = getProjectFolders(rootDir);

addpath(genpath(projectFolders.curves));
addpath(genpath(projectFolders.experiments));
addpath(genpath(projectFolders.postprocessing));


fileNetti1996Figure5  = fullfile( projectFolders.experiments_NDRAN1996,...
                        'Netti1996Figure5.csv' );
fileNetti1996Figure7  = fullfile( projectFolders.experiments_NDRAN1996,...
                        'Netti1996Figure7.csv' );
fileNetti1996Figure9  = fullfile( projectFolders.experiments_NDRAN1996,...
                        'Netti1996Figure9.csv' );
fileNetti1996Figure11 = fullfile( projectFolders.experiments_NDRAN1996,...
                        'Netti1996Figure11.csv');

dataNetti1996Figure5 = loadDigitizedData(fileNetti1996Figure5,...
                         'Frequency (Hz)',['Storage Modulus, E` MPa'],...
                         {'5N','10N','20N','30N','40N','60N'},'');
                       
                       
dataNetti1996Figure9 = loadDigitizedData(fileNetti1996Figure9,...
                         'Frequency (Hz)',['Loss Modulus, E" MPa'],...
                         {'5N','10N','20N','30N','40N','60N'},'');     
                       

dataNetti1996Figure7 = loadDigitizedData(fileNetti1996Figure7,...
                         'Force (N)',['Storage Modulus, E` MPa'],...
                         {'10Hz','1Hz'},'');                       

dataNetti1996Figure11 = loadDigitizedData(fileNetti1996Figure11,...
                         'Force (N)',['Loss Modulus, E` MPa'],...
                         {'10Hz','1Hz'},'');                       
                       
                 
if(flag_NormalizeByKmax==1)
  kmax = max(dataNetti1996Figure7(1).y);
  disp('Normalizing Storage and Loss coefficients against max. Storage coeff');
  disp([num2str(kmax),' MPa']);
  
  for z=1:1:length(dataNetti1996Figure5)
    dataNetti1996Figure5(z).y = dataNetti1996Figure5(z).y./kmax;
    dataNetti1996Figure5(z).yName = 'Norm. Storage Modulus (MPa/MPa)';
  end
  for z=1:1:length(dataNetti1996Figure7)
    dataNetti1996Figure7(z).y = dataNetti1996Figure7(z).y./kmax;
    dataNetti1996Figure7(z).yName = 'Norm. Storage Modulus (MPa/MPa)';
    
  end
  for z=1:1:length(dataNetti1996Figure9)
    dataNetti1996Figure9(z).y = dataNetti1996Figure9(z).y./kmax;
    dataNetti1996Figure9(z).yName = 'Norm. Storage Modulus (MPa/MPa)';    
  end
  for z=1:1:length(dataNetti1996Figure11)
    dataNetti1996Figure11(z).y = dataNetti1996Figure11(z).y./kmax;
    dataNetti1996Figure11(z).yName = 'Norm. Loss Modulus (MPa/MPa)';
    
  end
  
  
end
                       
%%
% configure the plot
%%

numberOfHorizontalPlotColumns = 3;
numberOfVerticalPlotRows      = 3;

pageWidth         = 21;
pageHeight        = 29.7;
plotWidth         = 5;
plotHeight        = 5;
plotHorizMarginCm = 1;
plotVertMarginCm  = 1;

flag_usingOctave  = 0;
plotConfigGeneric;

%%
% plot the data
%%
figNetti = figure;
figPub = figure;

load(fullfile(projectFolders.output_structs_FittedModels,'defaultFelineSoleus.mat'));



figure(figPub);
    subplot('Position',reshape(subPlotPanel(1,1,:),1,4));
    et = defaultFelineSoleus.musculotendon.tendonStrainAtOneNormForce;

    curveSample = calcBezierYFcnXCurveSampleVector(...
                    defaultFelineSoleus.curves.tendonForceLengthCurve,...
                    200,[0.99,(1+et+0.01)]);    
    plot(curveSample.x, curveSample.y,'Color',[0,0,0]);
    hold on;
    plot(1+et,1,'.k');
    hold on;
    %xticks([round(1.0,3),round(1+et,3)]);
    %xticklabels({'0','$$e^T_o$$'});
    xticks([round(1.0,3),round(1.05,3),round(1.10,4)]);
    xticklabels({'0','0.05','0.10'});
    yticks([round(0.,3),round(2/3,3),round(1.0,3)]);
    yticklabels({'0','$$f^T_{toe}$$','$$f^M_o$$'});

    xlabel('Strain ($$\ell^{T}/\ell^{T}_{s}-1$$)');
    ylabel('Norm. Force ($$f^{T}/f^{M}_{o}$$)');
    title('A. Tendon force-length ($$\mathbf{f}^T$$)');
    box off;
    axis([0.99,1.101,0,1.05]);

    ftIso = 161.3+86.4+24.1;
%    
% From Siebert et al, we know the maximum isometric
% forces that the gastrocs (161.3 +/- 18.2N), plantarius (86.4 +/- 21.3 N), 
% and soleus (24.1 +/- 5.8 N) can apply to the Achilles tendon of a 
% 3.2 kg rabbit is, in total 271.8 +/- 45.3 N. 
%
% If we assume at the rabbits in Siebert et al.'s study are a decent
% approximation of the rabbits in Netti et al.'s study (the species
% differs, and the Netti did not report the mass of the rabbits) then the 
% 60 N load in Netti's study corresponds to a normalized tendon load of
% 60/271.8 = 0.22
% 
% The plot of the storage modulus does not look great. Since it is unclear
% exactly what the maximum isometric force of the triceps surae of the 
%
% Siebert T, Leichsenring K, Rode C, Wick C, Stutzig N, Schubert H, 
% Blickhan R, Böl M. Three-dimensional muscle architecture and 
% comprehensive dynamic properties of rabbit gastrocnemius, plantaris 
% and soleus: input for simulation studies. PLoS one. 2015 Jun 26;10(6):e0130985.

figure(figNetti);
subplotStorageVsFreq = reshape(subPlotPanel(1,1,:),1,4);
subplot('Position',subplotStorageVsFreq);


%Storage modulus
for z=1:1:length(dataNetti1996Figure5)
  n = (z-1)/(length(dataNetti1996Figure5)-1);
  color0 = [0.5,0.25,0.25];
  color1 = [1,0,0];
  lineColor = color0.*(1-n) + color1.*(n);
  loglog(dataNetti1996Figure5(z).x,...
       dataNetti1996Figure5(z).y,...
       '-','Color',lineColor,...
       'DisplayName',[dataNetti1996Figure5(z).seriesName]);
  hold on;
  
  text(dataNetti1996Figure5(z).x(end),...
       dataNetti1996Figure5(z).y(end),...
       [dataNetti1996Figure5(z).seriesName],...
       'HorizontalAlignment','left',...
       'Color',lineColor);
  hold on;
     
  
  if(z==length(dataNetti1996Figure5))
    xlabel(dataNetti1996Figure5(z).xName);
    ylabel(dataNetti1996Figure5(z).yName);           
    box off
  end
  %axis square;
  
end

%xlim([0,100]);


subplotLossVsFreq = reshape(subPlotPanel(1,2,:),1,4);
subplot('Position',subplotLossVsFreq);

%Loss modulus
for z=1:1:length(dataNetti1996Figure9)
  n = (z-1)/(length(dataNetti1996Figure9)-1);
  color0 = [0.5,0.25,0.25];
  color1 = [1,0,0];
  lineColor = color0.*(1-n) + color1.*(n);
  loglog(dataNetti1996Figure9(z).x,...
       dataNetti1996Figure9(z).y,...
       '-','Color',lineColor,...
       'DisplayName',[dataNetti1996Figure9(z).seriesName,'N']);
  hold on;


  text(dataNetti1996Figure9(z).x(end),...
       dataNetti1996Figure9(z).y(end),...
       [dataNetti1996Figure9(z).seriesName],...
       'HorizontalAlignment','left',...
       'Color',lineColor);
  hold on;  
  
  if(z==length(dataNetti1996Figure9))
    xlabel(dataNetti1996Figure9(z).xName);
    ylabel(dataNetti1996Figure9(z).yName);       
   

    box off
    hold on;
  end
  %axis square;
  
end

%%
% 
%%
for idxFigure=1:1:2
    
    if(idxFigure==1)
        figure(figNetti);
        subplotStorageVsForce = reshape(subPlotPanel(2,1,:),1,4);
        subplotLossVsForce = reshape(subPlotPanel(3,1,:),1,4);
    end
    if(idxFigure==2)
        figure(figPub);
        subplotStorageVsForce = reshape(subPlotPanel(1,2,:),1,4);
        subplotLossVsForce = reshape(subPlotPanel(1,3,:),1,4);
    end

    
    subplot('Position',subplotStorageVsForce);
    if(idxFigure==2)
        title('B. Tendon storage modulus');
        hold on;
    end
    %Storage modulus
    for z=1:1:length(dataNetti1996Figure7)
      lineColor = [0,0,0];
      switch z
        case 1
          lineColor = [1,0,0];
        case 2
          lineColor = [1,0.5,0.5];      
        otherwise assert(0)      
      end
      
      plot(dataNetti1996Figure7(z).x,...
           dataNetti1996Figure7(z).y,...
           '-','Color',lineColor,...
           'DisplayName',[dataNetti1996Figure7(z).seriesName]);
      hold on;
      
      text(dataNetti1996Figure7(z).x(end),...
           dataNetti1996Figure7(z).y(end),...
           [dataNetti1996Figure7(z).seriesName, ' Netti'],...
           'HorizontalAlignment','right',...
           'Color',lineColor);
      hold on;
         
          
      if(z==length(dataNetti1996Figure7))
        xlabel(dataNetti1996Figure7(z).xName);
        ylabel(dataNetti1996Figure7(z).yName);       
            
        box off
      end
      %axis square;
      
    end

    if(flag_addCurveStorageModulus==1)
        ftCurve = curveSample.y.*ftIso;
        idx_5N_60N = find(ftCurve <= 60 & ftCurve >= 5);
        dydxAt60N = max(curveSample.dydx(idx_5N_60N));
        curveNormStiffness.f = ftCurve(idx_5N_60N);
        curveNormStiffness.kN = curveSample.dydx(idx_5N_60N)./dydxAt60N;
        plot(curveNormStiffness.f,...
             curveNormStiffness.kN,...
            '-','Color',[0,0,0] );
        hold on;
        text(curveNormStiffness.f(1,1),...
             curveNormStiffness.kN(1,1),'Model');
        hold on;
    end
    xlim([0,70]);
    
    
    
    subplot('Position',subplotLossVsForce);
    if(idxFigure==2)
        title('C. Tendon loss modulus with model fit');
        hold on;
    end

    for z=1:1:length(dataNetti1996Figure11)
      lineColor = [0,0,0];
      switch z
        case 1
          lineColor = [1,0,0];
        case 2
          lineColor = [1,0.5,0.5];      
        otherwise assert(0)      
      end
      
      plot(dataNetti1996Figure11(z).x,...
           dataNetti1996Figure11(z).y,...
           '-','Color',lineColor,...
           'DisplayName',[dataNetti1996Figure11(z).seriesName]);
      hold on;
      
      text(dataNetti1996Figure11(z).x(end),...
           dataNetti1996Figure11(z).y(end),...
           [dataNetti1996Figure11(z).seriesName,' Netti'],...
           'HorizontalAlignment','right',...
           'Color',lineColor);
      hold on;
         
      
      if(z==length(dataNetti1996Figure11))
        xlabel(dataNetti1996Figure11(z).xName);
        ylabel(dataNetti1996Figure11(z).yName);       
    
        box off
      end
      
      %axis square;
      
    end
    xlim([0,70]);
    %%
    % Fit a 2 parameter damping model to the data
    %%
    
    n01 = [0:0.1:1];
    sampleX = (10.^n01 - 1).*((10/9)*(3.96)) + 0.4;
    
    sampleFrequency = [1;10];
    
    sampleDamping  = zeros(length(dataNetti1996Figure9),length(sampleFrequency));
    sampleStiffness= zeros(length(dataNetti1996Figure9),length(sampleFrequency));
    
    
    
    for indexFrequency = 1:1:length(sampleFrequency)
      for indexLoad = 1:1:length(dataNetti1996Figure9)
    
        sampleStiffness(indexLoad,indexFrequency) ...
          = interp1(dataNetti1996Figure5(indexLoad).x,...
                    dataNetti1996Figure5(indexLoad).y,...
                    sampleFrequency(indexFrequency,1));                                  
    
        sampleDamping(indexLoad,indexFrequency) ...
          = interp1(dataNetti1996Figure9(indexLoad).x,...
                    dataNetti1996Figure9(indexLoad).y,...
                    sampleFrequency(indexFrequency,1));
      end
    end
    
    meanStiffness = mean(sampleStiffness,2);
    meanDamping   = mean(sampleDamping,2);
    
    A = [meanStiffness];%[meanStiffness, ones(size(sampleStiffness,1),1)];
    %M = eye(2,2).*max(sampleStiffness);
    b =  meanDamping;
    
    [x,flag_converged,relres,iter] =lsqr(A,b,1e-6,1000);
    
    disp('Damping Model Coefficients')
    disp(['  ',num2str(x(1,1)),' Linear']);
    %disp([num2str(x(2,1)),' Constant']);
    disp(['  ',num2str(relres),' Rel.Error']);
    disp(['  ',num2str(iter),' Iter']);
    disp(['  ',num2str(flag_converged),' Converged']);
    
    modelDamping= A*x;
    
    if(idxFigure==1)    
        subplot('Position',subplotLossVsFreq);
        
        for w = 1:1:length(modelDamping)
          n = (w-1)/(length(modelDamping)-1);
          color0 = [0.5,0.5,0.5];
          color1 = [0,0,0];
          lineColor = color0.*(1-n) + color1.*(n);
          
          loglog( [0.4,40],...
                  [1,1].*modelDamping(w,1),...
                  '-.','Color',lineColor);
          hold on;
          
          text(0.4,...
               modelDamping(w,1),...
               [dataNetti1996Figure9(w).seriesName, ' Model'],...
               'HorizontalAlignment','left',...
               'VerticalAlignment','bottom',...
               'Color',lineColor);
          hold on;  
          
        end
    end
    %xlim([0,100]);
    
    subplot('Position',subplotLossVsForce)
    
    sampleForces  = dataNetti1996Figure11(1).x;
    sampleTendonStiffness = zeros(size(sampleForces));
    
    for z=1:1:length(sampleForces)
      sampleTendonStiffness(z,1) = interp1(dataNetti1996Figure7(1).x,...
                                           dataNetti1996Figure7(1).y,...
                                           sampleForces(z,1),...
                                           'linear','extrap');
                                           
    end
    
    %sampleDamping = [sampleTendonStiffness ones(size(sampleTendonStiffness))]*x;
    sampleDamping = [sampleTendonStiffness]*x;
    plot(sampleForces, sampleDamping,'-.','Color',[0,0,0]);
    hold on;
    
    if(flag_addCurveStorageModulus==1)
        plot(curveNormStiffness.f,...
             curveNormStiffness.kN.*x,...
            '-','Color',[0,0,0] );
        hold on;
        text(curveNormStiffness.f(1,1),...
             curveNormStiffness.kN(1,1).*x,...
             'Model');
        hold on;
    end

    if(flag_addCurveStorageModulus==1)
        text(sampleForces(1),sampleDamping(1),'Scaled Storage Modulus',...
             'HorizontalAlignment','left',...
             'VerticalAlignment','top',...
             'Color',[0,0,0]);
        hold on;     
    else
        text(sampleForces(1),sampleDamping(1),'Model',...
             'HorizontalAlignment','left',...
             'VerticalAlignment','top',...
             'Color',[0,0,0]);
        hold on; 
    end
    
    text(1,max(sampleDamping),'$$\beta = \beta_1 \hat{K}^T$$');
    %text(1,max(sampleDamping),'$$\beta = (\beta_1 \hat{K}^T + \beta_0)\,K^T_{max}$$');
    %text(1,max(sampleDamping)*0.9,['$$\beta_0=',sprintf('%1.3f',x(2,1)),'$$']);
    text(1,max(sampleDamping)*0.8,['$$\beta_1=',sprintf('%1.4f',x(1,1)),'$$']);
    
    xlim([0,70]);
    %axis square;


end

figure(figNetti);
configPlotExporter;
filePath = fullfile(projectFolders.output_plots_NDRAN1996,...
                    'fig_NettiData_TendonDampingModel.pdf');
print('-dpdf', filePath); 

figure(figPub);
configPlotExporter;
filePath = fullfile(projectFolders.output_plots_NDRAN1996,...
                    'fig_Pub_NettiData_TendonDampingModel.pdf');
print('-dpdf', filePath); 
