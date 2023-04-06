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

function rootDir = getRootProjectDirectory()

localDirContents = dir();

rootDirName = 'Millard2023VexatMuscle';

assert(contains(localDirContents(1).folder,rootDirName),...
    sprintf(['Error: Matlab must start in the %s',...
             ' directory for this script to work'],rootDirName) );    

rootDir = localDirContents(1).folder;		
i0 = strfind(rootDir,rootDirName)-1;
i1 = i0+length(rootDirName);
rootDir = rootDir(1,1:i1);

