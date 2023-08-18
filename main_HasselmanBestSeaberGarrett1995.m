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

%Simulation of
% Hasselman CT, Best TM, Seaber AV, Garrett JR WE. A threshold and 
% continuum of injury during active stretch of rabbit skeletal muscle. 
% The American Journal of Sports Medicine. 1995 Jan;23(1):65-73.

rootDir         = getRootProjectDirectory();
projectFolders  = getProjectFolders(rootDir);

flag_useElasticTendon               = 1;

flag_simulateHillModel              = 0; 
flag_simulateVexatModel             = 1;

flag_useCalibratedVexatCurves       = 1;

flag_simulateActiveStretch          = 1;
flag_simulatePassiveStretch         = 1;

flag_plotData                       = 1;  
flag_savePlotsToFile                = 1;

flag_useOctave = 0;

addpath( genpath(projectFolders.parameters)     );
addpath( genpath(projectFolders.curves)         );
addpath( genpath(projectFolders.experiments)    );
addpath( genpath(projectFolders.simulation)     );
addpath( genpath(projectFolders.models)         );
addpath( genpath(projectFolders.postprocessing) );



plotLayoutSettings = struct('numberOfHorizontalPlotColumns',  1,...
                            'numberOfVerticalPlotRows',       3,...
                            'flag_fixedPlotWidth',            1,...
                            'plotWidth',                      7,...
                            'plotHeight',                     4,...
                            'flag_usingOctave',               0);
% Single column plots of:
% length vs time
% force vs. time
% stiffness vs. time
% damping vs time.
numberOfHorizontalPlotColumns = plotLayoutSettings.numberOfHorizontalPlotColumns;
numberOfVerticalPlotRows      = plotLayoutSettings.numberOfVerticalPlotRows;
flag_fixedPlotWidth           = plotLayoutSettings.flag_fixedPlotWidth;
plotWidth                     = plotLayoutSettings.plotWidth;
plotHeight                    = plotLayoutSettings.plotHeight;
flag_usingOctave              = plotLayoutSettings.flag_usingOctave;
plotConfig;

%Load the model parameters
tmp=load(fullfile(projectFolders.output_structs_FittedModels,...
                 'defaultRabbitTibialisAnterior.mat'));
musculotendonProperties   = tmp.defaultRabbitTA.musculotendon;
sarcomereProperties       = tmp.defaultRabbitTA.sarcomere;
normMuscleCurves          = tmp.defaultRabbitTA.curves;
fitting                   = tmp.defaultRabbitTA.fitting;    
outputFileEndingVexat = 'Linear-Titin';


%%
% Meta configuration properties: Do not touch.  
%%

dataFolder   = [projectFolders.experiments_HBSG1995,filesep];
structsFolder= [projectFolders.output_structs_HBSG1995,filesep];
plotFolder   = [projectFolders.output_plots_HBSG1995,filesep];


%%
% Ramp configuration from Hasselman et al.
%%
lceOpt      = musculotendonProperties.optimalFiberLength;
alphaOpt    = musculotendonProperties.pennationAngle;
ltSlk       = musculotendonProperties.tendonSlackLength;
etOne       = musculotendonProperties.tendonStrainAtOneNormForce;

lengthStart = lceOpt*cos(alphaOpt) + (1+etOne)*ltSlk;
lengthEnd   = lengthStart + lceOpt*cos(alphaOpt);
timeStretch = (lengthEnd-lengthStart)/0.10;

lengthRampKeyPoints = [      0, lengthStart;...
                       timeStretch, lengthEnd];

stimulationKeyTimes = [          0, 1;...
                     (timeStretch), 1];

timeSpan  = [0, (timeStretch+2)];

outputFileEndingVexat = '';
if(sarcomereProperties.titinModelType==0)
    outputFileEndingVexat = 'Linear-Titin';
end
if(sarcomereProperties.titinModelType==1)
    outputFileEndingVexat = 'WLC-Titin';    
end


%%
% Simulate Hasselman et al. for the Vexat model
%%
if(flag_simulateVexatModel==1)
    [success] = runHasselmanBestSeaberGarrett1995SimulationsVexat(...
                          timeSpan,...
                          lengthRampKeyPoints,...
                          stimulationKeyTimes,...
                          flag_useElasticTendon,...
                          musculotendonProperties,...
                          sarcomereProperties,...
                          normMuscleCurves,...
                          outputFileEndingVexat, ...
                          structsFolder,...
                          flag_simulateActiveStretch,...
                          flag_simulatePassiveStretch,...
                          flag_useOctave);  
end

%%
% Simulate Hasselman et al. for the Hill model
%%
if(flag_simulateHillModel==1)

end



