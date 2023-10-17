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

function projectFolders = getProjectFolders(rootProjectDirectoryFullPath)

currDir = pwd();

%check the root directory
cd(rootProjectDirectoryFullPath);
rootDirContents = dir();


flag_rootDirPathValid = 0;
for i=1:1:length(rootDirContents)
    if(strcmp(rootDirContents(i).name,'.rootDirectory') && ...
       rootDirContents(i).isdir == 0)
        flag_rootDirPathValid = 1;
    end
end
cd(currDir);

assert(flag_rootDirPathValid==1, ['Error: the rootProjectDirectoryFullPath ',...
    'does not appear to be the root project directory because it is missing ',...
    'hidden file (.rootDir) that marks it as the root project directory.']);


projectFolders.curves         = fullfile(rootProjectDirectoryFullPath,'curves'        );   
projectFolders.experiments    = fullfile(rootProjectDirectoryFullPath,'experiments'   );  
projectFolders.models         = fullfile(rootProjectDirectoryFullPath,'models'        ); 
projectFolders.output         = fullfile(rootProjectDirectoryFullPath,'output'        );  
projectFolders.parameters     = fullfile(rootProjectDirectoryFullPath,'parameters'    );   
projectFolders.postprocessing = fullfile(rootProjectDirectoryFullPath,'postprocessing');  
projectFolders.simulation     = fullfile(rootProjectDirectoryFullPath,'simulation'    ); 



projectFolders.experiments_HL2002        = fullfile(projectFolders.experiments, 'HerzogLeonard2002'                   );   
projectFolders.experiments_KBR1994       = fullfile(projectFolders.experiments, 'KirschBoskovRymer1994'               ); 
projectFolders.experiments_LJH2010       = fullfile(projectFolders.experiments, 'LeonardJoumaaHerzog2010'             ); 
projectFolders.experiments_NDRAN1996     = fullfile(projectFolders.experiments, 'NettiDamoreRoncaAmbrosioNicolais1996');  
projectFolders.experiments_StandardTests = fullfile(projectFolders.experiments, 'StandardTests'                       );   
projectFolders.experiments_TGFG1998      = fullfile(projectFolders.experiments, 'TrombitasGreaserFrenchGranzier1998'  );
projectFolders.experiments_HBSG1995      = fullfile(projectFolders.experiments, 'HasselmanBestSeaberGarret1995'       ); 


projectFolders.output_plots     = fullfile(projectFolders.output,'plots');
projectFolders.output_structs   = fullfile(projectFolders.output,'structs');
projectFolders.output_tables    = fullfile(projectFolders.output,'tables');


projectFolders.output_plots     = fullfile(projectFolders.output,'plots');

projectFolders.output_plots_HL2002         = fullfile(projectFolders.output_plots, 'HerzogLeonard2002'                   );   
projectFolders.output_plots_KBR1994        = fullfile(projectFolders.output_plots, 'KirschBoskovRymer1994'               ); 
projectFolders.output_plots_LJH2010        = fullfile(projectFolders.output_plots, 'LeonardJoumaaHerzog2010'             ); 
projectFolders.output_plots_NDRAN1996      = fullfile(projectFolders.output_plots, 'NettiDamoreRoncaAmbrosioNicolais1996');  
projectFolders.output_plots_StandardTests  = fullfile(projectFolders.output_plots, 'StandardTests'                       );   
projectFolders.output_plots_MuscleCurves   = fullfile(projectFolders.output_plots, 'MuscleCurves'                       );   
projectFolders.output_plots_Initialization = fullfile(projectFolders.output_plots, 'InitializationBenchmark'            );   
projectFolders.output_plots_HBSG1995       = fullfile(projectFolders.output_plots, 'HasselmanBestSeaberGarret1995'             ); 
projectFolders.output_plots_SystemIdentificationExample ...
                                           = fullfile(projectFolders.output_plots, 'SystemIdentificationExample'             ); 


projectFolders.output_tables_KBR1994      = fullfile(projectFolders.output_tables,'KirschBoskovRymer1994');
projectFolders.output_tables_MuscleCurves = fullfile(projectFolders.output_tables,'MuscleCurves');


projectFolders.output_structs_FittedModels   = fullfile(projectFolders.output_structs, 'FittedModels'                        );
projectFolders.output_structs_HL2002         = fullfile(projectFolders.output_structs, 'HerzogLeonard2002'                   );   
projectFolders.output_structs_KBR1994        = fullfile(projectFolders.output_structs, 'KirschBoskovRymer1994'               ); 
projectFolders.output_structs_LJH2010        = fullfile(projectFolders.output_structs, 'LeonardJoumaaHerzog2010'             ); 
projectFolders.output_structs_NDRAN1996      = fullfile(projectFolders.output_structs, 'NettiDamoreRoncaAmbrosioNicolais1996');  
projectFolders.output_structs_StandardTests  = fullfile(projectFolders.output_structs, 'StandardTests'                       );   
projectFolders.output_structs_MuscleCurves   = fullfile(projectFolders.output_structs, 'MuscleCurves'                       );   
projectFolders.output_structs_Initialization = fullfile(projectFolders.output_structs, 'InitializationBenchmark'            );   

projectFolders.output_structs_HBSG1995  = fullfile(projectFolders.output_structs, 'HasselmanBestSeaberGarret1995'             ); 


