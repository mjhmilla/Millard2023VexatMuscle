clc;
close all;
clear all;

filePath=fullfile('passiveStretches','passive_STR_085_13Lo_01vmax_090620_exp.dat');


numberOfHeaderLines = 4;
rowOfColumnLabels   = 3;
rowOfUnitLabels     = 4;


folders = {fullfile('data','passiveStretches'),...
           fullfile('data','Fiso')};



for i=1:1:length(folders)
    files = dir(folders{i});
    figH = figure;
    isFirstFile=1;
    for j=1:1:length(files)
        if(contains(files(j).name,'.dat')==1)
            filePath = fullfile(folders{i},files(j).name);
            auroraData = loadAuroraDatFile(filePath,...
                                numberOfHeaderLines,...
                                rowOfColumnLabels,...
                                rowOfUnitLabels);
            subplot(1,2,1);
            plot(auroraData.data(:,1),auroraData.data(:,2));
            hold on;
            box off;

            subplot(1,2,2);
            plot(auroraData.data(:,1),auroraData.data(:,3));
            hold on;
            box off

            if(isFirstFile)
                subplot(1,2,1);
                xlabel([auroraData.columnNames{1,1},' ',auroraData.columnUnits{1,1}]);
                ylabel([auroraData.columnNames{1,2},' ',auroraData.columnUnits{1,2}]);
                title([folders{i}, ' col 1 and 2']);
                subplot(1,2,2);
                xlabel([auroraData.columnNames{1,1},' ',auroraData.columnUnits{1,1}]);
                ylabel([auroraData.columnNames{1,3},' ',auroraData.columnUnits{1,3}]);
                title([folders{i}, ' col 1 and 3']);

                isFirstFile=0;
                
            end
        end
        here=1;
    end
end