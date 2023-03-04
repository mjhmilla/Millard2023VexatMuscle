clc;
close all;
clear all;

flag_usingOctave=0;
flag_NormalizeByKmax=1;

experimentsDirectoryTree      = genpath('experiments');
addpath(experimentsDirectoryTree      );

postprocessingDirectoryTree      = genpath('postprocessing');
addpath(postprocessingDirectoryTree   );

parametersDirectoryTreeCurves       = genpath('curves');
addpath(parametersDirectoryTreeCurves);

pubOutputFolder = 'output/plots/NettiDamoreRoncaAmbrosioNicolais1996/';

folderNetti1996 = 'experiments/NettiDamoreRoncaAmbrosioNicolais1996/';

fileNetti1996Figure5  = [folderNetti1996, 'Netti1996Figure5.csv'];
fileNetti1996Figure7  = [folderNetti1996, 'Netti1996Figure7.csv'];
fileNetti1996Figure9  = [folderNetti1996, 'Netti1996Figure9.csv'];
fileNetti1996Figure11 = [folderNetti1996, 'Netti1996Figure11.csv'];

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
    dataNetti1996Figure7(z).yName = 'Norm. Loss Modulus (MPa/MPa)';
    
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

%subPlotSquare(1,1)=subPlotSquare(1,1)*0.5;
%subPlotSquare(1,2)=subPlotSquare(1,2);
load('output/structs/defaultFelineSoleus.mat');

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
    xticks([round(1.0,3),round(1+et,3)]);
    xticklabels({'0','$$e^T_o$$'});
    yticks([round(0.,3),round(2/3,3),round(1.0,3)]);
    yticklabels({'0','$$f^T_{toe}$$','$$f^M_o$$'});

    xlabel('Strain ($$\ell^{T}/\ell^{T}_{s}-1$$)');
    ylabel('Norm. Force ($$f^{T}/f^{M}_{o}$$)');
    title('A. Tendon force-length ($$\mathbf{f}^T$$)');
    box off;
    axis tight;


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
           'HorizontalAlignment','left',...
           'Color',lineColor);
      hold on;
         
      
      if(z==length(dataNetti1996Figure7))
        xlabel(dataNetti1996Figure7(z).xName);
        ylabel(dataNetti1996Figure7(z).yName);       
            
        box off
      end
      %axis square;
      
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
           'HorizontalAlignment','left',...
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
    disp([num2str(x(1,1)),' Linear']);
    %disp([num2str(x(2,1)),' Constant']);
    disp([num2str(relres),' Rel.Error']);
    disp([num2str(iter),' Iter']);
    disp([num2str(flag_converged),' Converged']);
    
    modelDamping= A*x;
    
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
    
    text(sampleForces(1),sampleDamping(1),'Model',...
         'HorizontalAlignment','left',...
         'VerticalAlignment','top',...
         'Color',[0,0,0]);
    hold on;     
    
    text(1,max(sampleDamping),'$$\beta = \beta_1 \hat{K}^T$$');
    %text(1,max(sampleDamping),'$$\beta = (\beta_1 \hat{K}^T + \beta_0)\,K^T_{max}$$');
    %text(1,max(sampleDamping)*0.9,['$$\beta_0=',sprintf('%1.3f',x(2,1)),'$$']);
    text(1,max(sampleDamping)*0.8,['$$\beta_1=',sprintf('%1.4f',x(1,1)),'$$']);
    
    xlim([0,70]);
    %axis square;

end

figure(figNetti);
configPlotExporter;
print('-dpdf', [pubOutputFolder,'fig_NettiData_TendonDampingModel.pdf']); 

figure(figPub);
configPlotExporter;
print('-dpdf', [pubOutputFolder,'fig_Pub_NettiData_TendonDampingModel.pdf']); 