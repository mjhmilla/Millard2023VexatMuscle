clc;
close all;
clear all;

rootDir         = getRootProjectDirectory();
projectFolders  = getProjectFolders(rootDir);

addpath( genpath(projectFolders.parameters)     );
addpath( genpath(projectFolders.curves)         );
addpath( genpath(projectFolders.experiments)    );
addpath( genpath(projectFolders.simulation)     );
addpath( genpath(projectFolders.models)         );
addpath( genpath(projectFolders.postprocessing) );

%%
%
% Script configuration variables
%
%%

numberOfHeaderLines = 4;
rowOfColumnLabels   = 3;
rowOfUnitLabels     = 4;

dataFolder = fullfile('experiments','TomalkaWeidnerHahnSieberlSiebert2021','data');

folders = {'passiveStretches',...
           'Figures_4_5'};

plotRowId = {'_passive_','SHO_','SSC_','SHO_Bleb_','SSC_Bleb_'};
plotRowXlim = [28.25,60;...
               nan,nan;...
               nan,nan;...
               nan,nan;...
               nan,nan];

speedId = {'_03_','_06_','_085_'};
speedIdPlotted = zeros(size(speedId));

flagTrimData=1;
%%
%
% Plot configuration
%
%%
plotColors = getPaulTolColourSchemes('vibrant');

plotLayoutSettings = struct('numberOfHorizontalPlotColumns',  3,...
                            'numberOfVerticalPlotRows',       length(plotRowId),...
                            'flag_fixedPlotWidth',            1,...
                            'plotWidth',                      7,...
                            'plotHeight',                     10,...
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
%
% Plot the mat files
%
%%
figMat=figure;


for i=1:1:length(folders)
    fullFolderPath = fullfile(dataFolder,folders{i});
    files = dir(fullFolderPath);    
    isFirstFile=1;
    isDataPlotted=0;
    for j=1:1:length(files)
        if(contains(files(j).name,'.mat')==1)

            fullFilePath=fullfile(dataFolder,folders{i},files(j).name);
            matData=load(fullFilePath);
            matFields = fieldnames(matData);

            isFirstFigure=1;
            isFirstPlot=1;
            for k=1:1:length(matFields)

                lengthData=[];
                forceData=[];
                seriesId = '';
                seriesColor = plotColors.grey;

                if(contains(folders{i},'passiveStretches')==1 ...
                        && contains(matFields{k},'NormForce')==1)                    
                    tmpTime     =  [1:1:length(matData.length_STR_01vmax)]' ./ 1000;
                    lengthData  = [tmpTime, matData.length_STR_01vmax];
                    tmpTime     = [1:1:length(matData.(matFields{k}))]' ./ 1000;
                    forceData   = [tmpTime, matData.(matFields{k})];
                    seriesId    = matFields{k};
                    seriesId    = seriesId(1,11:16);
                end
                if(strcmp(folders{i},'Figures_4_5')==1 ...
                        && contains(matFields{k},'_L_')==1)

                    lengthData=[];
                    forceData=[];
                    seriesId = '';

                    lengthField = matFields{k};
                    forceField  = lengthField;
                    idx=strfind(forceField,'_L_');
                    forceField(1,idx+1)='F';
                    id = lengthField(1,(idx+3):end);

                    tmpTime     =  [1:1:length(matData.(lengthField))]' ./ 1000;
                    lengthData  = [tmpTime, matData.(lengthField)];
                    tmpTime     = [1:1:length(matData.(forceField))]' ./ 1000;
                    forceData   = [tmpTime, matData.(forceField)];
                    seriesId    = id;
                    idx=strfind(seriesId,'_');
                    seriesId(idx)=' ';
                    here=1;

                    for z=1:1:length(speedId)
                        if(contains(forceField,speedId{z}))
                            fieldColors = fields(plotColors);
                            seriesColor = plotColors.(fieldColors{z});
                        end
                    end
                end  

                row=nan;
                for(z=1:1:length(plotRowId))
                    idx = strfind(matFields{k},plotRowId{z});
                    if(~isempty(idx))
                        row=z;
                    end
                end

                if(~isempty(lengthData) && ~isempty(forceData) && ~isnan(row))

                    %Identify the length change times
                    samplingFreq=1000;%Hz
                    filterFrequency=250;
                    [b,a]=butter(2,2/(0.5*samplingFreq),'low');
                    l = filtfilt(b,a,lengthData(:,2));
                    %l=lengthData(:,2);
                    dldt = calcCentralDifferenceDataSeries(lengthData(:,1),l);
                    %dldt = filtfilt(b,a,dldt);
                    dl2dt2 = calcCentralDifferenceDataSeries(lengthData(:,1),dldt);
                    %dl2dt2 = filtfilt(b,a,dl2dt2);
                    threshold = 0.25*max(abs(dl2dt2));
                    
                    timeLengthChange = [];
                    for z=2:1:(length(dl2dt2)-1)
                        zL = dl2dt2(z)-dl2dt2(z-1);
                        zR = dl2dt2(z+1)-dl2dt2(z);
                        
                        if(abs(dl2dt2(z)) > threshold && zL*zR <= 0)
                            timeLengthChange =...
                                [timeLengthChange,lengthData(z,1)]; 
                        end
                    end

                    timeStart = min(timeLengthChange);
                    timeEnd = max(timeLengthChange);
                    timeMin = timeStart - 0.05*(timeEnd-timeStart);
                    timeMax = timeEnd + 0.05*(timeEnd-timeStart);

                    timeLengthInt = find(lengthData(:,1) >= timeMin...
                                        & lengthData(:,1) <= timeMax );
                    timeForceInt = find(forceData(:,1) >= timeMin...
                                        & forceData(:,1) <= timeMax);

                    figure(figMat);
                    isDataPlotted=1;
                    subplot('Position',reshape(subPlotPanel(row,1,:),1,4)); 

                    if(isempty(seriesColor))
                        plot(lengthData(timeLengthInt,1),...
                             lengthData(timeForceInt,2),...
                             'DisplayName',seriesId);
                    else
                        plot(lengthData(timeLengthInt,1),...
                             lengthData(timeForceInt,2),...
                             'Color',seriesColor,...
                             'DisplayName',seriesId,...
                             'HandleVisibility','off');
                    end

                    box off;
                    hold on;
                    subplot('Position',reshape(subPlotPanel(row,2,:),1,4));

                    if(isempty(seriesColor))
                        plot(forceData(timeForceInt,1),...
                             forceData(timeForceInt,2),...
                             'DisplayName',seriesId);
                    else
                        plot(forceData(timeForceInt,1),...
                             forceData(timeForceInt,2),...
                             'Color',seriesColor,...
                             'DisplayName',seriesId,...
                             'HandleVisibility','off');
                    end
                    box off;
                    hold on;

                    timeCommon = ...
                        [timeStart:((timeEnd-timeStart)/999):timeEnd];
                    lengthCommon= ...
                        interp1(lengthData(:,1),lengthData(:,2),timeCommon);
                    forceCommon= ...
                        interp1(forceData(:,1),forceData(:,2),timeCommon);

                    subplot('Position',reshape(subPlotPanel(row,3,:),1,4));

                    if(isempty(seriesColor))
                        plot(lengthCommon,...
                             forceCommon,...
                             'DisplayName',seriesId);
                        box off;
                        hold on;
                    else
                        plot(lengthCommon,...
                             forceCommon,...
                             'Color',seriesColor,...
                             'DisplayName',seriesId);
                        box off;
                        hold on;
                    end
                    
                end
                
            end

        end
    end
end

figure(figMat);
for z=1:1:length(plotRowId)
    row=z;
    subplot('Position',reshape(subPlotPanel(row,1,:),1,4));
    xlabel('Time (s)');
    ylabel('Norm. Length $$\ell / \ell^M_o$$ ');

%     legend('Location','SouthEast');
%     legend("boxoff");

    xLimSettings = plotRowXlim(z,:);
    if(sum(isnan(xLimSettings))==0)
        xlim(xLimSettings);
    end

    titleStr = plotRowId{z};
    idx = strfind(titleStr,'_');
    titleStr(idx)=' ';

    title(titleStr);

    subplot('Position',reshape(subPlotPanel(row,2,:),1,4));

    xLimSettings = plotRowXlim(z,:);
    if(sum(isnan(xLimSettings))==0)
        xlim(xLimSettings);
    end
    
    xlabel('Time (s)');
    ylabel('Norm. Force $$f / f^M_o$$ ');

    subplot('Position',reshape(subPlotPanel(row,2,:),1,4));
    xlabel('Norm. Length ($$\ell / \ell_o^M$$)');
    ylabel('Norm. Force $$f / f^M_o$$ ');
end

figure(figMat);
configPlotExporter;
filePath = fullfile(projectFolders.output_plots_TWHSS2021,...
                    'fig_TWHSS2021_MatFileData.pdf');
print('-dpdf', filePath); 

%%
%
% Plot the dat files
%
%%

plotLayoutSettings = struct('numberOfHorizontalPlotColumns',  2,...
                            'numberOfVerticalPlotRows',       2,...
                            'flag_fixedPlotWidth',            1,...
                            'plotWidth',                      7,...
                            'plotHeight',                     10,...
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

folders = {'passiveStretches',...
           'Fiso'};
plotRow = [1,2];


figDat=figure;

for i=1:1:length(folders)
    fullFolderPath = fullfile(dataFolder,folders{i});
    files = dir(fullFolderPath);    
    isFirstFile=1;
    isDataPlotted=0;
    for j=1:1:length(files)
        if(contains(files(j).name,'.dat')==1)
            filePath = fullfile(dataFolder,folders{i},files(j).name);
            auroraData = loadAuroraDatFile(filePath,...
                                numberOfHeaderLines,...
                                rowOfColumnLabels,...
                                rowOfUnitLabels);
            auroraData.name = files(j).name;

            legendEntry = auroraData.name;
            idxUS = strfind(legendEntry,'_');
            legendEntry(1,idxUS)=' ';

            figure(figDat);
            subplot('Position',reshape(subPlotPanel(plotRow(1,i),1,:),1,4));
            plot(auroraData.data(:,1),auroraData.data(:,2),...
                 'DisplayName',legendEntry);
            hold on;
            box off;

            subplot('Position',reshape(subPlotPanel(plotRow(1,i),2,:),1,4));
            plot(auroraData.data(:,1),auroraData.data(:,3),...
                'DisplayName',legendEntry);
            hold on;
            box off

            if(isFirstFile)
                isDataPlotted=1;
                subplot('Position',reshape(subPlotPanel(plotRow(1,i),1,:),1,4));
                    xlabel([auroraData.columnNames{1,1},' ',auroraData.columnUnits{1,1}]);
                    ylabel([auroraData.columnNames{1,2},' ',auroraData.columnUnits{1,2}]);
                    title([folders{i}, ' col 1 and 2 (.dat)']);
                
                subplot('Position',reshape(subPlotPanel(plotRow(1,i),2,:),1,4));
                    xlabel([auroraData.columnNames{1,1},' ',auroraData.columnUnits{1,1}]);
                    ylabel([auroraData.columnNames{1,3},' ',auroraData.columnUnits{1,3}]);
                    title([folders{i}, ' col 1 and 3 (.dat)']);

                isFirstFile=0;
                
            end
        end
        here=1;
    end
    if(isDataPlotted==1)
        subplot('Position',reshape(subPlotPanel(plotRow(1,i),1,:),1,4));
%         legend('Location','SouthWest');
%         legend("boxoff");
    end
end

figure(figDat);
configPlotExporter;
filePath = fullfile(projectFolders.output_plots_TWHSS2021,...
                    'fig_TWHSS2021_DatFileData.pdf');
print('-dpdf', filePath); 
