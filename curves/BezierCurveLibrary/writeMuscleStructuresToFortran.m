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

function success = writeMuscleStructuresToFortran(structOfMuscleCurves,...
    muscleName,folder)

success=0;
curveNames = fieldnames(structOfMuscleCurves);

for indexCurve=1:1:length(curveNames)
    fid = fopen([folder,muscleName,'_',curveNames{indexCurve},'.f'],'w');

    xpts = structOfMuscleCurves.(curveNames{indexCurve}).xpts;
    ypts = structOfMuscleCurves.(curveNames{indexCurve}).ypts;
    
    xEnd = structOfMuscleCurves.(curveNames{indexCurve}).xEnd;
    yEnd = structOfMuscleCurves.(curveNames{indexCurve}).yEnd;

    dydxEnd = structOfMuscleCurves.(curveNames{indexCurve}).dydxEnd;
    d2ydx2End = structOfMuscleCurves.(curveNames{indexCurve}).d2ydx2End;
    
%     fprintf(fid,'xpts\n');
%     for i=1:1:size(xpts,1)
%         for j=1:1:size(xpts,2)
%             if(j < size(xpts,2))
%                 fprintf(fid,['%1.16e',delimiter],xpts(i,j));
%             else
%                 fprintf(fid,['%1.16e',delimiter],xpts(i,j));
%             end
%         end
%     end

    fprintf(fid,'      real, dimension(%d,%d)::xpts \n',...
        size(xpts,1),size(xpts,2));
    fprintf(fid,'      real, dimension(%d,%d)::ypts \n',...
        size(ypts,1),size(ypts,2));
    fprintf(fid,'      real, dimension(%d,%d)::xEnd \n',...
        size(xEnd,1),size(xEnd,2));
    fprintf(fid,'      real, dimension(%d,%d)::yEnd \n',...
        size(yEnd,1),size(yEnd,2));
    fprintf(fid,'      real, dimension(%d,%d)::dydxEnd \n',...
        size(dydxEnd,1),size(dydxEnd,2));
    fprintf(fid,'      real, dimension(%d,%d)::d2ydx2End \n',...
        size(d2ydx2End,1),size(d2ydx2End,2));
    
    fid=writeArrayToFortran(fid,xpts,'xpts');
    fid=writeArrayToFortran(fid,ypts,'ypts');
    fid=writeArrayToFortran(fid,xEnd,'xEnd');
    fid=writeArrayToFortran(fid,yEnd,'yEnd');
    fid=writeArrayToFortran(fid,dydxEnd,    'dydxEnd');
    fid=writeArrayToFortran(fid,d2ydx2End,  'd2ydx2End');
    fclose(fid);
end

success=1;