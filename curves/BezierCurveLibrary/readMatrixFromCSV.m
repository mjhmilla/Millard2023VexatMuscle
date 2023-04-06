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

function [matrix, lastLine] = readMatrixFromCSV(fid,terminatingStr)
%%
% This function assumes that different matrices are separated by a row with
% the name of the matrix.
%  
% @param fid: the file handle
% @param terminatingStr: the name of the next matrix/field in the file
% @return matrix: the numerical values of the matrix in the file
% @return lastLine: the last line reached before the routine terminated
%%
line=fgetl(fid);
idxComma = strfind(line,',');
formatStr = '';

if(length(idxComma)==0)
    formatStr = '%e';
elseif(length(idxComma)==1)
    formatStr = '%e,%e\n';
else
    formatStr = '%e,';
    for i=1:1:(length(idxComma)-1)
        formatStr = [formatStr,'%e,'];
    end            
    formatStr = [formatStr,'%e\n'];            
end

matrix = [];
while ischar(line)==1 && contains(line,terminatingStr)==0 
    if(isempty(matrix)==1)                
        matrix = sscanf(line,formatStr)';
    else
        matrix = [matrix;sscanf(line,formatStr)']; 
    end
    line=fgetl(fid);        
end

lastLine=line;