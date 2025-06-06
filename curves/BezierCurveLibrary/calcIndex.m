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

function idx = calcIndex(x, ptsM,options)
%%
% Given a series of sequential ranges in the 2xRealN matrixReal
% ptsM, this function will calculate the row indexReal that
% has an interval ptsM(n,1)-ptsM(n,2) that contains the
% point xReal. If xReal lies on a border that defines a sub-interval
% the function will return the indexReal such that xReal is at the 
% start of hte interval.
%
% @param xReal: a double 
% @param ptsM: a 2 xReal n matrixReal that defines the range of n
%              sub intervals. For exRealample
%
%
%        ptsM = [0 1 2 3 4; ...
%                1 2 3 4 5]
%
%        defines 5 adjacent sub intervals 
%        [[0,1],[1,2],[2,3],[3,4],[4,5]]
%
% @param tol: how close a value xReal is allowed to be to a
%             sub interval border before it is declared to be
%             on the border.
%
% @returns idx: the column indexReal of ptsM that contains an
%               interval that includes xReal.
%
% ExRealample:
% 1. xReal = .123
%    ptsM = [0 1 2 3 4; ...
%            1 2 3 4 5]
%  then calcIndexReal will return 1
%
% 2. xReal = 1
%    ptsM as before
%    then calcIndexReal will return 2.
%
%%
xReal = real(x);
idx = 0;
flag_found = 0;

rows = size(ptsM,1);

tol = options.tol;

if( isempty( options.moduloRange ) == 0 && ...
           (xReal > options.moduloRange || xReal < 0 ) )
    if(xReal < 0 || xReal > options.moduloRange)        
        xReal = mod(xReal,options.moduloRange);
    end 
    if((xReal >= 0 && xReal <= options.moduloRange) == 0)
       fprintf('xReal: %e xmin: %e xmax: %e\n',...
                xReal, 0, options.moduloRange); 
    end
    assert( (xReal >= 0 && xReal <= options.moduloRange),'xReal is out of range');
end

%xReal is either monotonically increasing or decreasing
dxSign = 1;
if(ptsM(2,1)-ptsM(1,1) < 0)
    dxSign = -1;
end

for i=1:1:size(ptsM,2)
    
    if(dxSign == 1)
        if( (xReal >= ptsM(1,i) && xReal < ptsM(rows,i)) || ... 
            (abs(xReal-ptsM(1,i)) < tol || abs(xReal-ptsM(rows,i)) < tol) )
            idx = i;
            flag_found = 1;
        end
    else
        if( (xReal <= ptsM(1,i) && xReal > ptsM(rows,i)) || ... 
            (abs(xReal-ptsM(1,i)) < tol || abs(xReal-ptsM(rows,i)) < tol) )
            idx = i;
            flag_found = 1;
        end        
    end
    
end

%Check if the value xReal is identically the last point
if(flag_found == 0 && xReal == ptsM(rows,size(ptsM,2)))
    idx = size(ptsM,2);
    flag_found = 1;
end

if(flag_found==0)
   here=1; 
end

assert( (flag_found == 1),... 
    'Error: A value of xReal was used that is not within the Bezier curve set.');
