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

function fittedFelineSoleus = fitFelineSoleusPevkActinBondLocation( ...
                                defaultFelineSoleus,...
                                flag_useElasticTendon,...
                                felineSoleusPassiveForceLengthCurveSettings,...
                                projectFolders,...
                                flag_useOctave)

%The parameters updated are the 
%  :normPevkToActinAttachmentPoint (default of 0.5)
%  :normMaxActiveTitinToActinDamping (default of 20)
%
%The hand-tuned default values are quite good, but fitting is required to
%minimize the error. Since the process of both force development and 
%relaxation are nonlinear, there is not an elegant and fast way to find 
%these parameters without simulating the model directly. 

figureNumber       = 7;
subFigureNumber    = 2;
trialNumber        = 3;  

expConfigHerzogLeonard2002 =...
 getHerzogLeonard2002Configuration( figureNumber,...
                                    subFigureNumber, ...
                                    trialNumber,...
                                    projectFolders);

%dataFolder = 'experiments/HerzogLeonard2002/fitting/';
dataFolder = fullfile(projectFolders.output_structs_HL2002,['fitting',filesep]);


fittedFelineSoleus=defaultFelineSoleus;

tendonStr = '';
if(flag_useElasticTendon==1)
    tendonStr = 'ET';
else
    tendonStr = 'RT';
end

fittingStr = sprintf('HL2002_%i%i%i_%s',...
                    figureNumber,subFigureNumber, trialNumber,tendonStr);

fittedFelineSoleus.fitting = [fittedFelineSoleus.fitting;...
                              {fittingStr}];

    
[sarcomerePropertiesUpd,...
 normMuscleCurvesUpd] = ...
    updateActiveTitinParameters(defaultFelineSoleus.musculotendon, ...
                             defaultFelineSoleus.sarcomere,...
                             defaultFelineSoleus.curves,...
                             felineSoleusPassiveForceLengthCurveSettings,...
                             expConfigHerzogLeonard2002,...
                             flag_useElasticTendon,...
                             projectFolders,...
                             flag_useOctave);

fittedFelineSoleus.sarcomere=sarcomerePropertiesUpd;
fittedFelineSoleus.curves=normMuscleCurvesUpd;
