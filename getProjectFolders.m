function projectFolders = getProjectFolders(rootProjectDirectoryFullPath)




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

projectFolders.output_plots     = fullfile(projectFolders.output,'plots');
projectFolders.output_structs   = fullfile(projectFolders.output,'structs');
projectFolders.output_tables    = fullfile(projectFolders.output,'tables');


projectFolders.output_plots     = fullfile(projectFolders.output,'plots');

projectFolders.output_plots_HL2002        = fullfile(projectFolders.output_plots, 'HerzogLeonard2002'                   );   
projectFolders.output_plots_KBR1994       = fullfile(projectFolders.output_plots, 'KirschBoskovRymer1994'               ); 
projectFolders.output_plots_LJH2010       = fullfile(projectFolders.output_plots, 'LeonardJoumaaHerzog2010'             ); 
projectFolders.output_plots_NDRAN1996     = fullfile(projectFolders.output_plots, 'NettiDamoreRoncaAmbrosioNicolais1996');  
projectFolders.output_plots_StandardTests = fullfile(projectFolders.output_plots, 'StandardTests'                       );   

projectFolders.output_plots_MuscleCurves   = fullfile(projectFolders.output_plots, 'MuscleCurves'                       );   
projectFolders.output_plots_Initialization = fullfile(projectFolders.output_plots, 'InitializationBenchmark'            );   


projectFolders.output_tables_KBR1994      = fullfile(projectFolders.output_tables,'KirschBoskovRymer1994');
projectFolders.output_tables_MUscleCurves = fullfile(projectFolders.output_tables,'MuscleCurves');
