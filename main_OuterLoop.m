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

%%
% This function only calls 5 other functions, some of which have lengthy
% run times
%
%    36 min: main_CreateModels_OuterLoop;
%   226 min: main_KirschBoskovRymer1994_OuterLoop;
%    52 min: main_HerzogLeonard2002_OuterLoop;
%     6 min: main_LeonardJoumaaHerzog2010_OuterLoop;
%    26 min: main_ClassicForceLengthVelocityExperiments_OuterLoop;
%   345 min: Total (5 hours and 45 minutes)
% 
% *Intel i7-3630QM @ 2.40 GHz, Ubuntu 22
%  8 GB ram, SSD harddrive
%
% All of the extra code that you see beyond these 5 lines is to
% clear the memory between these (sometimes large) main functions, and
% to measure and store timing information.
%%
clc;
close all;
clear all;


 timeStart = tic;
% main_CreateModels_OuterLoop;
 secondsElapsed = toc(timeStart);

rootDir         = getRootProjectDirectory();
cd(rootDir);
fid=fopen('timing_main_OuterLoop.txt','w');
fprintf(fid,'main_CreateModels_OuterLoop, %1.1fmin\n',...
             round(secondsElapsed/60));
fclose(fid);


clc;
close all;
clear all;

timeStart = tic;
main_KirschBoskovRymer1994_OuterLoop;
secondsElapsed = toc(timeStart);

rootDir         = getRootProjectDirectory();
cd(rootDir);
fid=fopen('timing_main_OuterLoop.txt','a');
fprintf(fid,'main_KirschBoskovRymer1994_OuterLoop, %1.1fmin\n',...
        round(secondsElapsed/60));
fclose(fid);



clc;
close all;
clear all;

timeStart = tic;
main_HerzogLeonard2002_OuterLoop;
secondsElapsed = toc(timeStart);

rootDir         = getRootProjectDirectory();
cd(rootDir);
fid=fopen('timing_main_OuterLoop.txt','a');
fprintf(fid,'main_HerzogLeonard2002_OuterLoop, %1.1fmin\n',...
        round(secondsElapsed/60));
fclose(fid);



clc;
close all;
clear all;

timeStart = tic;
main_LeonardJoumaaHerzog2010_OuterLoop;
secondsElapsed = toc(timeStart);

rootDir         = getRootProjectDirectory();
cd(rootDir);
fid=fopen('timing_main_OuterLoop.txt','a');
fprintf(fid,'main_LeonardJoumaaHerzog2010_OuterLoop, %1.1fmin\n',...
        round(secondsElapsed/60));
fclose(fid);



clc;
close all;
clear all;

timeStart = tic;
main_ClassicForceLengthVelocityExperiments_OuterLoop;
secondsElapsed = toc(timeStart);

rootDir         = getRootProjectDirectory();
cd(rootDir);
fid=fopen('timing_main_OuterLoop.txt','a');
fprintf(fid,'main_ClassicForceLengthVelocityExperiments_OuterLoop, %1.1fmin\n',...
        round(secondsElapsed/60));
fclose(fid);


clc;
close all;
clear all;