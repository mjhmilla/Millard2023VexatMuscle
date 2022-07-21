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