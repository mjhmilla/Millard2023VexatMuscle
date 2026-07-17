function expTRSS2017 = loadTRSS2017Data(projectFolders)

[expDataDigitized, expIndices] = ...
        loadRatSkeletalMuscleData(projectFolders);

expTRSS2017Digitized = expDataDigitized(expIndices.index_TRSS2017);

expTRSS2017RawFields ={'data_TRSS2017_0p70_shortOutput','lN070';...
                   'data_TRSS2017_0p85_shortOutput','lN085';...
                   'data_TRSS2017_1p45_shortOutput','lN145';...
                   'data_TRSS2017_1p45_longOutput','lN145Long'};

expTRSS2017 = struct('lN070',[],'lN085',[],'lN145',[],'lN145Long',[]);

%expTRSS2017.keyIndices = ...
%  struct('lN070',[],'lN085',[],'lN145',[],'lN145Long',[]);
expTRSS2017.trials = {'lN070','lN085','lN145'};
expTRSS2017.trialsLong = {'lN145Long'};


for i=1:1:size(expTRSS2017RawFields,1)
  expTRSS2017.(expTRSS2017RawFields{i,2}) = ...
    readtable(fullfile(projectFolders.experiments_TRSS2017,...
                    [expTRSS2017RawFields{i,1},'.csv']));
  tmp = readtable(fullfile(projectFolders.experiments_TRSS2017,...
                    [expTRSS2017RawFields{i,1},'_keyIndices','.csv'])); 
  expTRSS2017.keyIndices.(expTRSS2017RawFields{i,2})=tmp.keyIndices;
end
expTRSS2017.sampleFrequencyHz = 1000;