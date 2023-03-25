%%
% This function only calls 5 other functions:
%
%   main_CreateModels_OuterLoop;
%   main_KirschBoskovRymer1994_OuterLoop;
%   main_HerzogLeonard2002_OuterLoop;
%   main_LeonardJoumaaHerzog2010_OuterLoop;
%   main_ClassicForceLengthVelocityExperiments_OuterLoop;
%
% All of the extra code that you see beyond these 5 lines is to
% clear the memory between these (sometimes large) main functions, and
% to measure and store timing information.
%%
clc;
close all;
clear all;


timeStart = tic;
main_CreateModels_OuterLoop;
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
fid=fopen('timing_main_OuterLoop.txt','w+');
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
fid=fopen('timing_main_OuterLoop.txt','w+');
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
fid=fopen('timing_main_OuterLoop.txt','w+');
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
fid=fopen('timing_main_OuterLoop.txt','w+');
fprintf(fid,'main_ClassicForceLengthVelocityExperiments_OuterLoop, %1.1fmin\n',...
        round(secondsElapsed/60));
fclose(fid);


clc;
close all;
clear all;