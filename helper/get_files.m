function files = get_files(direc, filt)
%this function returns a list of files as specificed by the variable input
    %filt = filter string
    %direc = cell array of directory names
%if there is no input or specific input does not exist, 
%throw error letting user know to put an input, or respecify input
if nargin~=2, error('get_files:missing inputs, Please input folder(s) and file filter.');
end
files = [];
%if the directory is already a CHARACTER array
if ischar(direc) % 
    currDir = direc;
    %find all files matching f*.nii
    tmp = dir(fullfile(currDir,filt)); 
    %then build the full path name for these files
    tmp = [repmat([currDir filesep],size(tmp,1),1) char(tmp.name)]; 
    files = char(files,tmp);
%if the directory is a CELL array
else 
    %determine whether the size of the first directory parsed is larger 
    %than the second
    if size(direc,1)>size(direc,2)
        %if it is, use the first directory
        nRuns=size(direc,1);
    else
        %if it is not, use the second
        nRuns=size(direc,2);
    end
    %loop through each EPI session in directory
    for runI=1:nRuns 
        currDir = char(direc{runI});
        %find all files matching f*.nii
        tmp = dir(fullfile(currDir,filt)); 
        %build the full path name for these files
        tmp = [repmat([currDir filesep],size(tmp,1),1) char(tmp.name)]; 
        files = char(files,tmp);
    end
end
files = files(~all(files'==' ')',:);
end