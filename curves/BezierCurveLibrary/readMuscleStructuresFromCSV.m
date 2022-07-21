function muscleCurves = readMuscleStructuresFromCSV(folder)



success=0;

workingDir = pwd;
cd(folder);
folderContents = dir;
cd(workingDir)




for indexFile=3:1:length(folderContents)

    if(   contains(folderContents(indexFile).name,'.csv')==1 ...
       && contains(folderContents(indexFile).name,'~')  == 0)
        curveName = folderContents(indexFile).name;
        curveName = curveName(1:1:(end-4));

        muscleCurves.(curveName) = ...
            struct( 'xpts',[],   'ypts',[],...
                    'xEnd',[],   'yEnd',[],...
                    'dydxEnd',[],'d2ydx2End',[],'integral',[]);

        fid = fopen([folder,folderContents(indexFile).name],'r');

        line=fgetl(fid);

        assert(contains(line,'xpts'));
        
        [matrix,line] = readMatrixFromCSV(fid,'ypts');
        assert(contains(line,'ypts'));
        muscleCurves.(curveName).xpts=matrix;
        
        [matrix,line] = readMatrixFromCSV(fid,'xEnd');
        assert(contains(line,'xEnd'));
        muscleCurves.(curveName).ypts=matrix;

        [matrix,line] = readMatrixFromCSV(fid,'yEnd');
        assert(contains(line,'yEnd'));
        muscleCurves.(curveName).xEnd=matrix;

        [matrix,line] = readMatrixFromCSV(fid,'dydxEnd');
        assert(contains(line,'dydxEnd'));
        muscleCurves.(curveName).yEnd=matrix;

        [matrix,line] = readMatrixFromCSV(fid,'d2ydx2End');
        assert(contains(line,'d2ydx2End'));
        muscleCurves.(curveName).dydxEnd=matrix;

        [matrix,line] = readMatrixFromCSV(fid,'dummyEndString');
        assert(~ischar(line));
        muscleCurves.(curveName).d2ydx2End=matrix;

        
        here=1;

    end
end
