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

function success = writeMuscleStructuresToCSV(structOfMuscleCurves,...
    muscleName,folder)

success=0;

curveNames = fieldnames(structOfMuscleCurves);

for indexCurve=1:1:length(curveNames)
    fid = fopen([folder,muscleName,'_',curveNames{indexCurve},'.csv'],'w');

    xpts = structOfMuscleCurves.(curveNames{indexCurve}).xpts;
    ypts = structOfMuscleCurves.(curveNames{indexCurve}).ypts;
    
    xEnd = structOfMuscleCurves.(curveNames{indexCurve}).xEnd;
    yEnd = structOfMuscleCurves.(curveNames{indexCurve}).yEnd;

    dydxEnd = structOfMuscleCurves.(curveNames{indexCurve}).dydxEnd;
    d2ydx2End = structOfMuscleCurves.(curveNames{indexCurve}).d2ydx2End;
    

    fprintf(fid,'xpts\n');
    for i=1:1:size(xpts,1)
        for j=1:1:size(xpts,2)
            if(j==1)
                fprintf(fid,'%1.16e',xpts(i,j));                
            elseif(j==size(xpts,2))
                fprintf(fid,',%1.16e\n',xpts(i,j));
            else
                fprintf(fid,',%1.16e',xpts(i,j));
            end
        end
    end
    fprintf(fid,'ypts\n');
    for i=1:1:size(ypts,1)
        for j=1:1:size(ypts,2)
            if(j==1)
                fprintf(fid,'%1.16e',ypts(i,j));                
            elseif(j==size(ypts,2))
                fprintf(fid,',%1.16e\n',ypts(i,j));
            else
                fprintf(fid,',%1.16e',ypts(i,j));
            end
        end
    end

    fprintf(fid,'xEnd\n');
    fprintf(fid,'%1.16e,%1.16e\n',xEnd(1,1),xEnd(1,2));

    fprintf(fid,'yEnd\n');
    fprintf(fid,'%1.16e,%1.16e\n',yEnd(1,1),yEnd(1,2));

    fprintf(fid,'dydxEnd\n');
    fprintf(fid,'%1.16e,%1.16e\n',dydxEnd(1,1),dydxEnd(1,2));

    fprintf(fid,'d2ydx2End\n');
    fprintf(fid,'%1.16e,%1.16e\n',d2ydx2End(1,1),d2ydx2End(1,2));

    fclose(fid);
end

success=1;