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

%Simulation of
% Hasselman CT, Best TM, Seaber AV, Garrett JR WE. A threshold and 
% continuum of injury during active stretch of rabbit skeletal muscle. 
% The American Journal of Sports Medicine. 1995 Jan;23(1):65-73.

rootDir         = getRootProjectDirectory();
projectFolders  = getProjectFolders(rootDir);

flag_useElasticTendon               = 1;

flag_simulateHillModel              = 1; 
flag_simulateVexatModel             = 1;
    flag_simulateActiveStretch      = 1;
    flag_simulatePassiveStretch     = 1;
    flag_useCalibratedVexatCurves   = 1;

flag_plotData                       = 1;  
    flag_sizePlotsForSlides         = 1;
    flag_savePlotsToFile            = 1;

flag_useOctave = 0;
flag_usingOctave = 0;

addpath( genpath(projectFolders.parameters)     );
addpath( genpath(projectFolders.curves)         );
addpath( genpath(projectFolders.experiments)    );
addpath( genpath(projectFolders.simulation)     );
addpath( genpath(projectFolders.models)         );
addpath( genpath(projectFolders.postprocessing) );


numberOfHorizontalPlotColumns   = 2;
numberOfVerticalPlotRows        = 3;

plotWidth         = 4;
plotHeight        = 4;        
plotHorizMarginCm = 1.5;
plotVertMarginCm  = 2.;  
baseFontSize      = 6;

if(flag_sizePlotsForSlides==1)
    plotWidth         = 4;
    plotHeight        = 4;        
    plotHorizMarginCm = 1.5;
    plotVertMarginCm  = 2.;  
    baseFontSize      = 6;
end
pageWidth   = numberOfHorizontalPlotColumns*(plotWidth+plotHorizMarginCm)...
                +2*plotHorizMarginCm;
pageHeight  = numberOfVerticalPlotRows*(plotHeight+plotVertMarginCm)...
                +2*plotVertMarginCm;

plotConfigGeneric;

figOutputPlot = figure;
%%
% Meta configuration properties: Do not touch.  
%%

dataFolder   = [projectFolders.experiments_HBSG1995,filesep];
structsFolder= [projectFolders.output_structs_HBSG1995,filesep];
plotFolder   = [projectFolders.output_plots_HBSG1995,filesep];

%%
%Load the model parameters
%%
rabbitStudy = 'hasselman';%'siebert';
files = {[rabbitStudy,'RabbitTibialisAnterior.mat'],[rabbitStudy,'RabbitExtensorDigitorumLongus.mat']};
structNames = {[rabbitStudy,'RabbitTA'],[rabbitStudy,'RabbitEDL']};
outputFileEnding = {'TA','EDL'};

%TA
% 53.65N Failure force (avg) of table 4    
% 14.03N Max. Isometric force (avg) from table 5
%  3.82 fo at failure
% Li: Initial length is at point where passive force generates 50g.
% Failure strain L/Li is
% 38.305
%
%EDL
% 124.66N Failure force (avg) of table 4
%  35.21N Max. Isometric force (avg) from table 5
%   3.54 fo at failure
% Li: Initial length is at point where passive force generates 50g.
% Failure strain L/Li is
% 24.025%
%

hasselmanFailureData.initialCELength = [nan,nan];

hasselmanFailureData.strainFailure  = [0.38305,0.24025];
hasselmanFailureData.normForce      = [3.82,3.54]; 
hasselmanFailureData.pathStart      = [0.0889, 0.0948];
hasselmanFailureData.pathOptimal    = [0.0988, 0.1017];

if(flag_simulateVexatModel || flag_simulateHillModel)
    for indexFile = 1:1:length(files)
    
        disp(['Running Hasselman et al. on ',files{indexFile}]);
    
        tmp=load(fullfile(projectFolders.output_structs_FittedModels,...
                         files{indexFile}));
        musculotendonProperties   = tmp.(structNames{indexFile}).musculotendon;
        sarcomereProperties       = tmp.(structNames{indexFile}).sarcomere;
        normMuscleCurves          = tmp.(structNames{indexFile}).curves;
        fitting                   = tmp.(structNames{indexFile}).fitting;    
        
        if(indexFile==1)
            normMuscleCurves.useCalibratedCurves = flag_useCalibratedVexatCurves;
        end

              
        %%
        % Ramp configuration from Hasselman et al.
        %%
        lceOpt      = musculotendonProperties.optimalFiberLength;
        alphaOpt    = musculotendonProperties.pennationAngle;
        ltSlk       = musculotendonProperties.tendonSlackLength;
        etOne       = musculotendonProperties.tendonStrainAtOneNormForce;
        fiso        = musculotendonProperties.fiso;

     

        expKinematics = calcHasselmanBestSeaberGarrett1995ExpKinematics(...
                            musculotendonProperties,...
                            normMuscleCurves);

        hasselmanFailureData.initialCELength(1,indexFile)= expKinematics.lceNStart;
        lengthStart = expKinematics.lpStart;

        
        disp('Hasselman path start vs. model');
        fprintf('\t%1.3f\tvs\t%1.3f\n\n',hasselmanFailureData.pathStart(1,indexFile), expKinematics.lpStart);
        disp('Hasselman path optimal vs. model');
        fprintf('\t%1.3f\tvs\t%1.3f\n\n',hasselmanFailureData.pathOptimal(1,indexFile), expKinematics.lpOpt);

        switch indexFile
            case 1
                lengthEnd   = lengthStart + 2*lceOpt*cos(alphaOpt);
            case 2
                lengthEnd   = lengthStart + 3*lceOpt*cos(alphaOpt);                
            otherwise
                assert(0,'Error: indexFilemust be 1 or 2');
        end
               
        
        timeStretch = (lengthEnd-lengthStart)/0.10;
        
        dt = 0.3;


        lengthRampKeyPoints = [      2*dt, lengthStart;...
                               (timeStretch+2*dt), lengthEnd];
        
        stimulationKeyTimes = [          dt, 1;...
                             (timeStretch+3*dt), 1];
        
        timeSpan  = [0, (timeStretch+dt*4)];
        
        %%
        % Output file endings
        %%
        
        outputFileEndingVexat   = outputFileEnding{indexFile};
        outputFileEndingHill    = outputFileEnding{indexFile};
    
        if(sarcomereProperties.titinModelType==0)
            outputFileEndingVexat = [ outputFileEndingVexat,'_Linear-Titin'];
        end
        if(sarcomereProperties.titinModelType==1)
            outputFileEndingVexat = [outputFileEndingVexat, '_WLC-Titin'];    
        end
        
        nameModification='';
        if(flag_useElasticTendon == 1)
            nameModification = '_ElasticTendon';
        else
            nameModification = '_RigidTendon';          
        end    
        
        outputFileEndingVexat = [outputFileEndingVexat,nameModification];
        outputFileEndingHill  = [outputFileEndingHill,nameModification];    
    
        %%
        % Simulate Hasselman et al. for the Vexat model
        %%
        if(flag_simulateVexatModel==1)
            [success] = runHasselmanBestSeaberGarrett1995SimulationsVexat(...
                                  timeSpan,...
                                  lengthRampKeyPoints,...
                                  stimulationKeyTimes,...
                                  flag_useElasticTendon,...
                                  musculotendonProperties,...
                                  sarcomereProperties,...
                                  normMuscleCurves,...
                                  outputFileEndingVexat, ...
                                  structsFolder,...
                                  flag_simulateActiveStretch,...
                                  flag_simulatePassiveStretch,...
                                  flag_useOctave);  
        end
        
        %%
        % Simulate Hasselman et al. for the Hill model
        %%
        if(flag_simulateHillModel==1)
            flag_useFiberDamping=1;
            fiberDampingCoefficient=0.01;
            [success] = runHasselmanBestSeaberGarrett1995SimulationsDampedEquilibrium( ...
                                      timeSpan,...
                                      lengthRampKeyPoints,...
                                      stimulationKeyTimes,...
                                      flag_useElasticTendon,...
                                      flag_useFiberDamping,...
                                      fiberDampingCoefficient,...
                                      musculotendonProperties,...
                                      sarcomereProperties,...
                                      normMuscleCurves,...
                                      outputFileEndingHill, ...
                                      structsFolder,...
                                      flag_simulateActiveStretch,...
                                      flag_simulatePassiveStretch,...
                                      flag_useOctave);
        end
    
    end
end

%%
% Post process the data
%%
if(flag_plotData==1)
        
    figure(figOutputPlot);



    for indexFile=1:1:length(files)
        %%
        % Load the parameter file
        %%
        tmp=load(fullfile(projectFolders.output_structs_FittedModels,...
                         files{indexFile}));
        musculotendonProperties   = tmp.(structNames{indexFile}).musculotendon;
        sarcomereProperties       = tmp.(structNames{indexFile}).sarcomere;
        normMuscleCurves          = tmp.(structNames{indexFile}).curves;
        fitting                   = tmp.(structNames{indexFile}).fitting;    

        lceOpt = musculotendonProperties.optimalFiberLength;
        fiso   = musculotendonProperties.fiso;
        %%
        % Load the simulation file
        %%
        nameModification = '';

        outputFileEndingVexat   = outputFileEnding{indexFile};
        outputFileEndingHill    = outputFileEnding{indexFile};
    
        if(sarcomereProperties.titinModelType==0)
            outputFileEndingVexat = [ outputFileEndingVexat,'_Linear-Titin'];
        end
        if(sarcomereProperties.titinModelType==1)
            outputFileEndingVexat = [outputFileEndingVexat, '_WLC-Titin'];    
        end

        if(flag_useElasticTendon == 1)
            nameModification = '_ElasticTendon';
        else
            nameModification = '_RigidTendon';          
        end

        vexatSimDataFile = [structsFolder,'benchRecordVexat_',...
                            outputFileEndingVexat,nameModification,'.mat'];
        dampedEqSimDataFile= [structsFolder,'benchRecordHill_',...
                            outputFileEndingHill,nameModification,'.mat'];


        switch indexFile
            case 1
                vexatColor  = [0,0,1];
                vexatLine   = '-';
                hillColor   = [1,0,0];
                hillLine   = '-';                
            case 2
                vexatColor  = [0,0,1];
                vexatLine   = '-';
                hillColor   = [1,0,0];
                hillLine   = '-';                                
            otherwise
                asser(0,'indexFile must be 1 or 2');
        end

        vexatSimData = load(vexatSimDataFile);
        hillSimData = load(dampedEqSimDataFile);

        timeKeyPoints = [vexatSimData.lengthRampKeyPoints(:,1);...
                         vexatSimData.stimulationKeyTimes(:,1)];
        timeKeyPoints=sort(timeKeyPoints);
        timeLengthPoints = vexatSimData.lengthRampKeyPoints(:,1);
        timeStimPoints   = vexatSimData.stimulationKeyTimes(:,1);

        %%
        % Evaluate extras for each plot
        %%
        %Thesholds of injury from Hasselman et al. and Noonan et al.
        tfN = 3.41;
        tfActiveMinorN    = tfN*0.7;
        tfActiveMajorN    = tfN*0.9;
        tfActiveRuptureN  = tfN;
        tfPassiveMinorN   = tfN*0.3;
        tfPassiveMajorN   = tfN*0.8;
        tfPassiveRuptureN = tfN;

        c1 = ones(size(vexatSimData.benchRecord.time(:,1)));

        minorInjuryV = c1.*tfPassiveMinorN ...
            + vexatSimData.benchRecord.activation(:,1).*(tfActiveMinorN-tfPassiveMinorN);
        majorInjuryV = c1.*tfPassiveMajorN ...
            + vexatSimData.benchRecord.activation(:,1).*(tfActiveMajorN-tfPassiveMajorN);

        ruptureV     = c1.*tfN;

        idxA = find(vexatSimData.benchRecord.time >= vexatSimData.lengthRampKeyPoints(1,1),1);
        idxA = idxA-1;
        yA = minorInjuryV(idxA,1);
        idxB=idxA;
        yB = majorInjuryV(idxB,1);
        idxC = idxB;
        yC = ruptureV(idxC,1);

        expKinematics = calcHasselmanBestSeaberGarrett1995ExpKinematics(...
                            musculotendonProperties,...
                            normMuscleCurves);        

        subPlotPanel(1,indexFile,4) = subPlotPanel(1,indexFile,4)*0.33;    

        lceOpt  = musculotendonProperties.optimalFiberLength;
        ltSlk   = musculotendonProperties.tendonSlackLength;
        alphaOpt= musculotendonProperties.pennationAngle;

        pathLength0 = ltSlk+lceOpt*cos(alphaOpt);

        subplot('Position',reshape(subPlotPanel(1,indexFile,:),1,4));
            plot(vexatSimData.benchRecord.time(:,1), ...
                 vexatSimData.benchRecord.pathLength(:,1)./pathLength0,...
                 vexatLine,'Color',vexatColor);
            hold on;
            plot(hillSimData.benchRecord.time(:,1), ...
                 hillSimData.benchRecord.pathLength(:,1)./pathLength0,...
                 hillLine,'Color',hillColor);
            hold on;


            xticks(round(vexatSimData.lengthRampKeyPoints(:,1)',2));
            yticks(vexatSimData.lengthRampKeyPoints(:,2)'./pathLength0);

            box off;
            xlabel('Time (s)');
            ylabel('Norm. Length ($$\ell / \ell_o$$)');
            title('Path Length');

        subplot('Position',reshape(subPlotPanel(2,indexFile,:),1,4));
            plot(vexatSimData.benchRecord.time(:,1), minorInjuryV,'--','Color',[1,1,1].*0.5);
            hold on;
            text(vexatSimData.benchRecord.time(idxA,1), yA, 'Minor Injury',...
                'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',6);
            hold on;
            plot(vexatSimData.benchRecord.time(:,1), majorInjuryV,'--','Color',[1,1,1].*0.25);
            hold on;
            text(vexatSimData.benchRecord.time(idxB,1), yB, 'Major Injury',...
                'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',6);        
            hold on;
            plot(vexatSimData.benchRecord.time(:,1), ruptureV,'-','Color',[0,0,0].*0.5);
            hold on;
            text(vexatSimData.benchRecord.time(idxC,1), yC, 'Rupture',...
                'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',6);        
            hold on;

            plot(hillSimData.benchRecord.time(:,1),...
                 hillSimData.benchRecord.normFiberForceAlongTendon(:,1),...
                 vexatLine,'Color',hillColor);
            hold on;
            
            plot(vexatSimData.benchRecord.time(:,1),...
                 vexatSimData.benchRecord.normFiberForceAlongTendon(:,1),...
                 vexatLine,'Color',vexatColor);
            hold on;

            lpi = expKinematics.lpStart;
            lpf = lpi * (1+hasselmanFailureData.strainFailure(1,indexFile));
            timeF = interp1(vexatSimData.lengthRampKeyPoints(:,2),...
                           vexatSimData.lengthRampKeyPoints(:,1),lpf);

            ffN = hasselmanFailureData.normForce(1,indexFile);

            plot(timeF,ffN,'s','Color',[0,0,0],'MarkerFaceColor',[0,0,0]);
            hold on;

            idxHillFailure  = find(hillSimData.benchRecord.normFiberForceAlongTendon(:,1) > tfN*1.1, 1 );
            timeHillFailure = hillSimData.benchRecord.time(idxHillFailure,1);
            timeMax = hillSimData.benchRecord.time(:,1);

            timeMaxTick = round(timeHillFailure,1);
            
            xticks(round(timeKeyPoints',2));

            yticks(round([0,1,yA,yB,yC],2));
            xlim([0,timeHillFailure]);
            ylim([0,4.1]);
            box off;
            xlabel('Time (s)');
            ylabel('Norm. Force ($$f^M / f_o^M$$)');


        subplot('Position',reshape(subPlotPanel(3,indexFile,:),1,4));
        %%
        % Plot the thresholds of failure
        %%
        plot([0,5], [1,1].*tfN,'-','Color',[1,1,1].*0,'LineWidth',0.5,...
                'DisplayName','none');
        hold on;
        text(0, yC, 'Rupture',...
            'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',6);        
        hold on;        

        plot([0,5], [1,1].*tfActiveMajorN,'-','Color',[1,1,1].*0.25,'LineWidth',0.5,...
             'DisplayName','none');
        hold on;
        text(0, yB, 'Major Injury',...
            'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',6);        
        hold on;
        
        plot([0,5],[1,1].*tfActiveMinorN,'-','Color',[1,1,1].*0.5,'LineWidth',0.5,...
             'DisplayName','none');
        hold on;
        text(0, yA, 'Minor Injury',...
            'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',6);
        hold on;

        %%
        % Plot the force-length envelope
        %%
        lceN0 = 0;
        lceN1 = 1.5*max(vexatSimData.benchRecord.normFiberLength(:,1));

        lceNV = lceN0 + ([0:0.01:1]').*(lceN1-lceN0);
        fpeNV = zeros(length(lceNV),1);
        falNV = zeros(length(lceNV),1);
        fenvNV= zeros(length(lceNV),1);

        for idxPt=1:1:length(lceNV)
            fpeNV(idxPt,1)=calcBezierYFcnXDerivative(lceNV(idxPt,1),...
                            normMuscleCurves.fiberForceLengthCurve, 0);
            falNV(idxPt,1)=calcBezierYFcnXDerivative(lceNV(idxPt,1),...
                            normMuscleCurves.activeForceLengthCurve, 0);
            fenvNV(idxPt,1)= falNV(idxPt,1)+fpeNV(idxPt,1);
        end

        lceNVFill = [lceNV(1,1);lceNV(end,1);fliplr(lceNV')';lceNV(1,1)];
        fpeNVFill = [fpeNV(1,1);fpeNV(1,1);fliplr(fpeNV')';fpeNV(1,1)];

        falNVFill = [falNV(1,1);falNV(end,1);fliplr(falNV')';falNV(1,1)];
        fenvNVFill = [fenvNV(1,1);fenvNV(1,1);fliplr(fenvNV')';fenvNV(1,1)];

        fill(lceNVFill,fenvNVFill,[1,1,1].*0.6,'EdgeColor','none',...
             'DisplayName','none');
        hold on;        

        fill(lceNVFill,falNVFill,[1,1,1].*0.7,'EdgeColor','none',...
             'DisplayName','none');
        hold on;
        fill(lceNVFill,fpeNVFill,[1,1,1].*0.8,'EdgeColor','none',...
             'DisplayName','none');
        hold on;

        plot(lceNV, fenvNV,'-','Color',[1,1,1].*0.5,'LineWidth',1.5,...
            'DisplayName','none');
        hold on;

        plot(hillSimData.benchRecord.normFiberLength(:,1),...
             hillSimData.benchRecord.normFiberForceAlongTendon(:,1),...
             vexatLine,'Color',hillColor,'LineWidth',1,...
             'DisplayName','Hill');
        hold on;
        
        plot(vexatSimData.benchRecord.normFiberLength(:,1),...
             vexatSimData.benchRecord.normFiberForceAlongTendon(:,1),...
             vexatLine,'Color',vexatColor,'LineWidth',1,...
             'DisplayName','VEXAT');
        hold on;

        %%
        %Add annotation
        %%
        fNList = [tfActiveMinorN;tfActiveMajorN;tfActiveRuptureN];
        for i=1:1:length(fNList)
            idx = find(hillSimData.benchRecord.normFiberForceAlongTendon(:,1)>fNList(i,1),1);
            plot(hillSimData.benchRecord.normFiberLength(idx,1),...
                hillSimData.benchRecord.normFiberForceAlongTendon(idx,1),...
                'o','Color',hillColor,'MarkerFaceColor',[1,1,1],'MarkerSize',2,...
                'DisplayName','none');
            hold on;
            text(hillSimData.benchRecord.normFiberLength(idx,1),...
                 hillSimData.benchRecord.normFiberForceAlongTendon(idx,1),...
                 sprintf('%1.2f%s', hillSimData.benchRecord.normFiberLength(idx,1),'$$\ell_o$$'),...
                 'HorizontalAlignment','left','VerticalAlignment','bottom',...
                 'FontSize',6);
            hold on;
        end

        for i=1:1:length(fNList)
            idx = find(vexatSimData.benchRecord.normFiberForceAlongTendon(:,1)>fNList(i,1),1);
            plot(vexatSimData.benchRecord.normFiberLength(idx,1),...
                vexatSimData.benchRecord.normFiberForceAlongTendon(idx,1),...
                'o','Color',vexatColor,'MarkerFaceColor',[1,1,1],'MarkerSize',2,...
                'DisplayName','none');
            hold on;
            text(vexatSimData.benchRecord.normFiberLength(idx,1),...
                 vexatSimData.benchRecord.normFiberForceAlongTendon(idx,1),...
                 sprintf('%1.2f%s', vexatSimData.benchRecord.normFiberLength(idx,1),'$$\ell_o$$'),...
                 'HorizontalAlignment','right','VerticalAlignment','bottom',...
                 'FontSize',6);
            hold on;
        end        


        plot(vexatSimData.benchRecord.normFiberLength(1,1),...
             vexatSimData.benchRecord.normFiberForceAlongTendon(1,1),...
             's','Color',[0,0,0],'MarkerFaceColor',[1,1,1]);
        hold on;
        text(vexatSimData.benchRecord.normFiberLength(1,1)+0.05,...
             vexatSimData.benchRecord.normFiberForceAlongTendon(1,1),...
             '1','FontSize',6, 'HorizontalAlignment','left','VerticalAlignment','bottom');

        timeRampStart = vexatSimData.lengthRampKeyPoints(1,1);
        idxRampStart = find(vexatSimData.benchRecord.time(:,1) >= timeRampStart,1);
        idxRampStart = idxRampStart-1;
        plot(vexatSimData.benchRecord.normFiberLength(idxRampStart,1),...
             vexatSimData.benchRecord.normFiberForceAlongTendon(idxRampStart,1),...
             's','Color',[0,0,0],'MarkerFaceColor',[1,1,1]);
        hold on;
        text(vexatSimData.benchRecord.normFiberLength(idxRampStart,1)-0.05,...
             vexatSimData.benchRecord.normFiberForceAlongTendon(idxRampStart,1),...
             '2','FontSize',6, 'HorizontalAlignment','right','VerticalAlignment','bottom');
        


        %%
        %Evaluate Li
        %%
        %fi = 0.050*9.81; %50g load
        %fiN = fi/fiso;
        %liN = calcBezierFcnXGivenY(fiN,normMuscleCurves.fiberForceLengthCurve,1.0);
        %liN = hasselmanFailureData.initialCELength(1,indexFile);
        lpi = expKinematics.lpStart;
        lpf = lpi * (1+hasselmanFailureData.strainFailure(1,indexFile));
        idx = find(vexatSimData.benchRecord.pathLength >= lpf,1);
        
        lceFN = vexatSimData.benchRecord.normFiberLength(idx,1);
        ffN = hasselmanFailureData.normForce(1,indexFile);

%         plot(lfN,ffN,'s','MarkerFaceColor',[0,0,0],...
%             'DisplayName',['HBSG1995: ',outputFileEnding{indexFile}]);
%         hold on;
        
        plot(lceFN,ffN,'s','Color',[0,0,0],'MarkerFaceColor',[0,0,0],...
            'DisplayName',['HBSG1995: ',outputFileEnding{indexFile}]);
        hold on;

        text(lceFN+0.05,...
             ffN,...
             '3','FontSize',6, 'HorizontalAlignment','left','VerticalAlignment','bottom');
        hold on;
        text(lceFN-0.05,...
             ffN,...
             'Exp.','FontSize',6, 'HorizontalAlignment','right','VerticalAlignment','bottom');
        hold on;

        lceNExp = [vexatSimData.benchRecord.normFiberLength(1,1);...
                   vexatSimData.benchRecord.normFiberLength(idxRampStart,1);...
                   lceFN];
        fceNExp = [vexatSimData.benchRecord.normFiberForceAlongTendon(1,1);...
                   vexatSimData.benchRecord.normFiberForceAlongTendon(idxRampStart,1);...
                   ffN];

        %Manually made legend;
        plot([0.05,0.4],[1.5,1.5],vexatLine,'Color',vexatColor);
        hold on;
        text(0.45,1.5,'VEXAT','FontSize',6);
        hold on;
        plot([0.05,0.4],[1.0,1.0],hillLine,'Color',hillColor);
        hold on;
        text(0.45,1.0,'Hill','FontSize',6);
        hold on;
        plot(0.225,0.5,'s','Color',[0,0,0],'MarkerFaceColor',[0,0,0]);
        hold on;
        text(0.45,0.5,'Exp.','FontSize',6);
        hold on;

        %hasselmanFailureData.strainFailure = [0.38305,0.24025];
        %hasselmanFailureData.normForce  = [3.82,3.54];         

        %lengthHillFailure = hillSimData.benchRecord.normFiberLength(idxHillFailure,1);
        %xtickAnnotation = round(sort([0,1,lceNExp(2,1), lceNExp(end,1)]),2);
        xtickAnnotation = round(sort([0,lceNExp(2:3,1)',2.8]),2);
        xticks(xtickAnnotation);

        ytickAnnotation = round(sort([0,1,yA,yB,yC,fceNExp(end,1)]),2);
        yticks(ytickAnnotation);

        yMax = 4.1;
        idx = find(hillSimData.benchRecord.normFiberForceAlongTendon(:,1) >= yMax,1);
        xMax = 2.27;%hillSimData.benchRecord.normFiberLength(idx,1);
        switch rabbitStudy
            case 'siebert'
                xMax = 2.81;
            case 'hasselman'
                xMax = 2.27;
            otherwise
                assert(0,'Unrecognized rabbit study');
        end

        xlim([0,xMax]);
        ylim([0,yMax]);
        box off;

        xtickangle(0);

        xlabel('Norm. Length ($$\ell/\ell_o^M$$)')
        ylabel('Norm. Force ($$f^T/f_o^M$$)')
        title(outputFileEnding{indexFile});
    end
end

if(flag_savePlotsToFile==1)

    figure(figOutputPlot);
    
    set(figOutputPlot,'Units','centimeters',...
    'PaperUnits','centimeters',...
    'PaperSize',[pageWidth pageHeight],...
    'PaperPositionMode','manual',...
    'PaperPosition',[0 0 pageWidth pageHeight]);     
    %set(findall(figList(i).h,'-property','FontSize'),'FontSize',10);     
    set(figOutputPlot,'renderer','painters');     
    set(gcf,'InvertHardCopy','off')
    
    
    filePlotName = fullfile(projectFolders.output_plots_HBSG1995,...
                   ['fig_Pub_HasselmanBestSeaberGarrett1995_',rabbitStudy,'.pdf']);
    
    
    print('-dpdf', filePlotName);   


end
