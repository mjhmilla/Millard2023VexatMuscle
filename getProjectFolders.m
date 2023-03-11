function projectFolders = getProjectFolders(rootProjectDirectoryFullPath)


startingDir = pwd;
cd(rootProjectDirectoryFullPath);

projectFolders.root = rootProjectDirectoryFullPath;

rootDirContents = dir();
experimentNameList = [];
%Automatically build a structure that contains the 
for i=1:1:length(rootDirContents)

    if(rootDirContents(i).isdir==1 ...
        && contains(rootDirContents(i).name,'.')==0) 

                    

        cd( rootDirContents(i).name );
        subDirContents=dir();

        subFolders=[];
        subFolders.('path')=fullfile(rootProjectDirectoryFullPath,rootDirContents(i).name);

        for j=1:1:length(subDirContents)
            if(subDirContents(j).isdir ...
                && contains(subDirContents(j).name,'.')==0)

    
                subFolders.(subDirContents(j).name) = ...
                    fullfile(subFolders.path, ...
                             subDirContents(j).name);

                if(strcmp(rootDirContents(i).name,'experiments'))
                    experimentNameList.(subDirContents(j).name)=...
                        subDirContents(j).name;
                end

            end

        end

        %fieldName = [rootDirContents(i).name,'_'];


        projectFolders.(rootDirContents(i).name) = subFolders;

        if(isempty(experimentNameList)==0)
            projectFolders.experiment_names = experimentNameList;
            experimentNameList=[];
        end

        cd(projectFolders.root);

    end

end




cd(startingDir);