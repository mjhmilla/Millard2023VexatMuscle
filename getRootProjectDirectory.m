function rootDir = getRootProjectDirectory()

localDirContents = dir();

rootDirName = 'Millard2021ImpedanceMuscle';

assert(contains(localDirContents(1).folder,rootDirName),...
    sprintf(['Error: Matlab must start in the %s',...
             ' directory for this script to work'],rootDirName) );    

rootDir = localDirContents(1).folder;		
i0 = strfind(rootDir,rootDirName)-1;
i1 = i0+length(rootDirName);
rootDir = rootDir(1,1:i1);

