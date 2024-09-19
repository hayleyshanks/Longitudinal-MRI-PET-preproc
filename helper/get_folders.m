function folder_names = get_folders(rootPath)
    allFiles = dir(rootPath);
    %extract only those that are directories.
    allDirFlags = [allFiles.isdir];
    allFolders = allFiles(allDirFlags);
    allFolders = allFolders(~ismember({allFolders(:).name},{'.','..'}));
    folder_names = {allFolders.name};
end 