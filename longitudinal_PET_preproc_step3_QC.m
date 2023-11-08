%% PET QC
% code by Hayley R. C. Shanks, Taylor W. Schmitz, Kate M. Onuska
function longitudinal_PET_preproc_step3_QC(N, study_path, pet_to_run, pet_orders,study_name, frames, qc_thresh,mri_ref, withinSubThresh, pet_QC_prefix, varargin)

%%%%%%%%%%%%%%%%%%% lab-specfic INPUTS - may need to be changed %%%%%%%%%%%%%%%%%%%%%%%%%
% add path to spm
addpath '/srv/schmitz/projects/toolboxes/spm12';
% add cat12 path
addpath '/srv/schmitz/projects/toolboxes/spm12/toolbox/cat12';
functionRoot = '/srv/schmitz/projects/toolboxes/incalab_source_code/human';
% catch sight of the functions
addpath(genpath(functionRoot));
%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of specific inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
study_avg_dir = [study_path, filesep, study_name, '_averages'];
% raw participant data should live in a first level subdirectory of the study path
rootPath = [study_path, filesep, 'first_level'];
sheet_path = [study_path, filesep, 'spreadsheets'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parse user inputs
pars = inputParser;
% by default, all modules are set to off.
default_run = 0;
% required inputs
addRequired(pars, 'N', @isnumeric);
addRequired(pars, 'study_path', @ischar);
addRequired(pars, 'pet_to_run', @iscell);
addRequired(pars, 'pet_orders', @iscell);
addRequired(pars, 'study_name', @ischar);
addRequired(pars, 'frames', @iscell); % static or dynamic
addRequired(pars, 'qc_thresh', @isnumeric);
addRequired(pars, 'mri_ref', @iscell);
addRequired(pars, 'withinSubThresh', @isnumeric);
addRequired(pars, 'pet_QC_prefix', @ischar);
%add optional inputs.
addParameter(pars,'getSUVRfails',default_run,@isnumeric); 
addParameter(pars,'withinSubQC',default_run,@isnumeric); 
addParameter(pars,'betweenSubQC',default_run,@isnumeric); 
addParameter(pars,'getMotionIssues',default_run,@isnumeric);  
addParameter(pars,'summarizeQC',default_run,@isnumeric);  
addParameter(pars,'checkReg',default_run,@isnumeric);  
% parse the function inputs
parse(pars, N, study_path, pet_to_run, pet_orders, study_name, frames, qc_thresh, mri_ref, withinSubThresh, pet_QC_prefix, varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get a list of all files and folders in the root path folder.
RIDs = get_folders(rootPath);
mbatch = 0;
sub_to_run = [1:length(RIDs)];
% read in the master info sheet
mastersheet = readtable(fullfile(sheet_path, [study_name, '_PET_mastertable.csv']));
% make the pet_orders variable a 1x2 cell for simplicity if it wasn't
% specified as a 1x2. 
if isequal(size(pet_orders,2),1)
    % if the cell is only 1x1, then the optional argument was not specified
    % so code will default to including cross sectional subjects. specify
    % this here so the code won't throw errors 
    pet_orders = [pet_orders, {'includeCrossSectionalSubjects'}];
end
mastersheet.keys = [1:height(mastersheet)]';    
%% check that the pet_order variables are set correctly 
if ~isnumeric(pet_orders{1}) && ~strcmpi(pet_orders{1}, 'all_times')
    disp(['input to pet_orders is set to ', pet_orders{1}, ' which is not a valid input'])
    disp('please review function documentation for the correct input arguments and try again')
    return
end
% figure out who fails MRI data QC or has none
mri_mastersheet = readtable([sheet_path, filesep, study_name, '_MRI_mastertable.csv']);
% when the summarized QC variable is 1, they have data passing QC. If it's
% 0 it fails, and NaN if not processed
if strcmpi(mri_ref{1}, 'base')
    good_mri = mri_mastersheet.RID(mri_mastersheet.summarizedTimepointFails ==0);
elseif strcmpi(mri_ref{1}, 'avg')
    good_mri = mri_mastersheet.RID(mri_mastersheet.summarizedAvgFails ==0);
else
    disp(['please review mri_ref description and re-define the variable'])
    return
end
for t=1:length(pet_to_run)    
    if pars.Results.getSUVRfails
        if ~contains(mastersheet.Properties.VariableNames, 'hasSUVR')
            mastersheet.hasSUVR = nan(height(mastersheet),1);
        end
        for i=1:height(mastersheet)
            if ~contains(mastersheet.RID(i), good_mri) 
                % won't have an SUVR if they have no mri
                mastersheet.hasSUVR(i) =0;           
            else 
                if ~strcmpi(mastersheet.ScanType(i), pet_to_run{t})
                    % skip
                else
                    current_dir = fullfile(rootPath, char(mastersheet.RID(i)), char(mastersheet.Date(i)), char(mastersheet.Modality(i)), char(mastersheet.ScanType(i)), char(mastersheet.ScanName(i)));                    
                    if strcmpi(frames, {''})
                        dir_mod = '';
                    else
                        possibilities = get_folders(current_dir);
                        to_use = strcmpi(possibilities, frames);
                        dir_mod = [filesep, char(possibilities(to_use))];
                    end
                    current_dir = [current_dir, dir_mod];
                    check = get_files(current_dir, 'suvr*.nii');
                    if isempty(check)
                        mastersheet.hasSUVR(i) =0;
                    else
                        mastersheet.hasSUVR(i) =1;
                    end
                end
            end
        end
        writetable(mastersheet, [sheet_path, filesep, study_name, '_PET_mastertable.csv']);
    end
    if pars.Results.withinSubQC
        if isequal(sum(strcmpi(mastersheet.Properties.VariableNames, 'withinSubWarpQC')), 0)
            mastersheet.withinSubWarpQC = nan(height(mastersheet),1);
        end
        for subI=sub_to_run(N)
            if ~contains(RIDs{subI}, good_mri) 
                % they failed Qc. Skip.
            else
                sub_dir = fullfile(rootPath,RIDs{subI});
                sub_table = get_files(sub_dir, 'scan_info*.csv');
                % this table has all of the dates and types of scans for this P
                sub_table = readtable(sub_table);
                table_to_run = spread_pet_filter(sub_table, pet_to_run{t});
                if height(table_to_run) < 2
                    % skip
                else
                    imgToQC = {};
                    for row=1:height(table_to_run)
                        % collect the images to warp
                        current_dir = fullfile(sub_dir, char(table_to_run.Date(row)), char(table_to_run.Modality(row)), char(table_to_run.ScanType(row)), char(table_to_run.ScanName(row)));                    
                        if strcmpi(frames, {''})
                            dir_mod = '';
                        else
                            possibilities = get_folders(current_dir);
                            to_use = strcmpi(possibilities, frames);
                            dir_mod = [filesep, char(possibilities(to_use))];
                        end
                        current_dir = [current_dir, dir_mod];
                        % if the scan is dynamic, use the mean. 
                        if strcmpi(frames, 'dynamic')
                            imgToQC{row} = get_files(current_dir,[pet_QC_prefix, 'rmean*.nii']);   
                        else
                            imgToQC{row} = get_files(current_dir,[pet_QC_prefix, 'rac*.nii']);     
                        end
                    end
                    sub_rows = strcmp(mastersheet.RID, RIDs{subI});
                    modality_rows = strcmp(mastersheet.ScanType, pet_to_run{t});
                    sub_modality = sum([sub_rows, modality_rows, mastersheet.TimepointReference], 2);% rowsums of these concatenated arrays tell us which scans have the right sub and modality and the ref
                    covStruct = cat_stat_check_cov(struct('data_vol',{{ imgToQC' }'} ,'gap',3,'c',[],'data_xml',{{}}));
                    all_covs = covStruct.covmat; % gives the covariance between each image and all other images
                    % remove the 1's from each row because that's the
                    % covariance with itself which isn't important
                    all_covs(all_covs ==1) = NaN;
                    mean_covs = mean(all_covs,2, 'omitnan'); % get the row means, omitting nan
                    mastersheet.withinSubWarpQC(sub_modality == 3) = mean_covs;
                end
            end
        end
        writetable(mastersheet, [sheet_path, filesep, study_name, '_PET_mastertable.csv']);
    end
    if pars.Results.betweenSubQC
        if isequal(sum(strcmpi(mastersheet.Properties.VariableNames, 'betweenSubWarpQC')),0)
            mastersheet.betweenSubWarpQC = nan(height(mastersheet),1);
        end
        % we're only using one timepoint of data here
        if strcmpi(pet_orders{1}, 'all_times') 
            disp('Only one timepoint of data can be used for between subject QC. Defaulting to time 1');
            time_to_use = 1;
        elseif isnumeric(pet_orders{1}) && length(pet_orders{1}) > 1
            disp('Only one timepoint of data can be used for between subject QC. Defaulting to time 1');
            time_to_use =1;
        else
            time_to_use = pet_orders{1};
        end
        to_qc = {};
        subset = mastersheet(mastersheet.Order == time_to_use & mastersheet.TimepointReference ==1 & mastersheet.hasSUVR ==1 & strcmp(mastersheet.ScanType, pet_to_run{t}), :);
        RID_subset = RIDs';
        RID_subset = RID_subset(N);
        for row=1:height(subset)
            % first check subtable to make sure they're usable
             sub_table = readtable(get_files([rootPath, filesep, char(subset.RID(row))], 'scan_info*.csv'));
             if isempty(sub_table)
                 continue
             else
                 pet_table = spread_pet_filter(sub_table, pet_to_run{t});
                if max(sub_table.Order) < 2 && strcmpi(pet_orders{2}, 'excludeCrossSectionalSubjects')
                    continue
                else
                    current_dir = fullfile(rootPath, char(subset.RID(row)), char(subset.Date(row)), char(subset.Modality(row)), char(subset.ScanType(row)), char(subset.ScanName(row)));
                    if strcmpi(frames, {''})
                        dir_mod = '';
                    else
                        possibilities = get_folders(current_dir);
                        to_use = strcmpi(possibilities, frames);
                        dir_mod = [filesep, char(possibilities(to_use))];
                    end
                    current_dir = [current_dir, dir_mod];
                    if strcmpi(frames, 'dynamic')
                        to_qc{row} = get_files(current_dir, [pet_QC_prefix, 'rmean*.nii']);   
                    else
                        to_qc{row} = get_files(current_dir, [pet_QC_prefix, 'rac*.nii']);   
                    end                
                end
             end
        end
        covStruct = cat_stat_check_cov(struct('data_vol',{{ to_qc' }'} ,'gap',3,'c',[],'data_xml',{{}}));
        all_covs = covStruct.covmat; % gives the covariance between each image and all other images
        % remove the 1's from each row because that's the
        % covariance with itself which isn't important
        all_covs(all_covs ==1) = NaN;
        mean_covs = mean(all_covs,2, 'omitnan'); % get the row means, omitting nan
        % find the mastersheet rows corresponding to the rows in
        % subset
        to_add = ismember(mastersheet.keys, subset.keys);
        mastersheet.betweenSubWarpQC(to_add) = mean_covs;
        writetable(mastersheet, [sheet_path, filesep, study_name, '_PET_mastertable.csv']);
    end
    if pars.Results.getMotionIssues
    end
    if pars.Results.summarizeQC
        % find average corr for the given tracer
        tracer_qc = mastersheet.betweenSubWarpQC(strcmpi(mastersheet.ScanType, pet_to_run{t}));
        mean_corr = mean(tracer_qc, 'omitnan');
        sd = std(tracer_qc, 'omitnan');
        thresh = mean_corr - 2*sd;
        double_check = nan(height(mastersheet),1);
        if isequal(sum(strcmpi(mastersheet.Properties.VariableNames, 'summarizedFails')),0)
            mastersheet.summarizedFails = nan(height(mastersheet),1);
        end
        for i=1:height(mastersheet)
            if ~strcmpi(mastersheet.ScanType(i), pet_to_run{t})
                % skip
            else
                if mastersheet.withinSubWarpQC(i) < withinSubThresh
                    mastersheet.summarizedFails(i) =1;
                    double_check(i) =1;
                else
                    if isnan(mastersheet.betweenSubWarpQC(i))
                        % that means this scan is not the reference scan used for
                        % this subject in between sub QC. find the ref
                        % scan
                        sub = mastersheet(strcmp(mastersheet.RID, mastersheet.RID{i}), :);
                        if isnumeric(pet_orders{1})
                            match = sub(sub.Order == pet_orders{1} & sub.TimepointReference ==1, :);
                        else
                            % order 1 would be used
                            match = sub(sub.Order == 1 & sub.TimepointReference ==1, :);
                        end
                        if match.betweenSubWarpQC >= thresh
                            mastersheet.summarizedFails(i) = 0;
                        else
                            mastersheet.summarizedFails(i) = 1;
                            if ~isnan(match.betweenSubWarpQC)
                                % if this is nan that just means their scan
                                % wasnt processed (e.g. data is cross sec)
                                double_check(i) =1;
                            end
                        end
                    elseif mastersheet.betweenSubWarpQC(i) >= thresh
                        mastersheet.summarizedFails(i) = 0;
                    end
                end
            end
        end
        for j=1:length(double_check)
            if ~isequal(double_check(j), 1)
                % skip
            else
                % display their slice so we can see if it actually
                % looks bad. people 2 std's below the mean corr may
                % still have usable data if the overall sample
                % quality is high
                check_dir = fullfile(rootPath, mastersheet.RID{j}, char(mastersheet.Date(j)), char(mastersheet.Modality(j)), char(mastersheet.ScanType(j)), char(mastersheet.ScanName(j)));
                if strcmpi(frames, {''})
                    dir_mod = '';
                else
                    possibilities = get_folders(check_dir);
                    to_use = strcmpi(possibilities, frames);
                    dir_mod = [filesep, char(possibilities(to_use))];
                end
                check_dir = [check_dir, dir_mod];
                if strcmpi(frames, 'dynamic')
                    to_check = get_files(check_dir, [pet_QC_prefix, 'rmean*.nii']);   
                else
                    to_check = get_files(check_dir, [pet_QC_prefix, 'rac*.nii']);     
                end        
                % display the image 
                mbatch = mbatch +1;
                matlabbatch{mbatch}.spm.util.disp.data = {to_check};
                spm_jobman('run',matlabbatch(mbatch))
                pause
                close 
                disp([mastersheet.RID{j} ' on ' datestr(mastersheet.Date(j)) ' has been flagged as a potential QC fail.'])
                disp(['verify that you would like the previous image to be excluded from further analysis'])
                user_check = input('enter 1 to EXCLUDE the image and 0 to INCLUDE the image \n');
                if user_check
                    mastersheet.summarizedFails(j) =1;
                else
                    mastersheet.summarizedFails(j) =0;
                end
            end
        end
     writetable(mastersheet, [sheet_path, filesep, study_name, '_PET_mastertable.csv']);
    end
    if pars.Results.checkReg
        % randomly display subjects who passed QC 
        num_subs = input('enter the number of subjects you would like to check \n');
        master_subset = mastersheet(mastersheet.summarizedFails ~=1, :);
        rids = RIDs(N); % drop any RIDs the user didn't ask to run
        % also drop any rids who failed QC - they would be checked already
        rids = rids(contains(rids, master_subset.RID));
        % randomly shuffle the rid indices 
        randoms = randperm(length(rids));
        % take the first num_subs from the random indices
        ids_to_check = rids(randoms(1:num_subs));
        for i=1:length(ids_to_check)
            sub_table = master_subset(strcmp(master_subset.RID, ids_to_check(i)), :);
            to_check = {};
            for j=1:height(sub_table)
                check_dir = fullfile(rootPath, sub_table.RID{j}, char(sub_table.Date(j)), char(sub_table.Modality(j)), char(sub_table.ScanType(j)), char(sub_table.ScanName(j)));
                if strcmpi(frames, {''})
                    dir_mod = '';
                else
                    possibilities = get_folders(check_dir);
                    to_use = strcmpi(possibilities, frames);
                    dir_mod = [filesep, char(possibilities(to_use))];
                end
                check_dir = [check_dir, dir_mod];
                if strcmpi(frames, 'dynamic')
                    to_check{end+1} = get_files(check_dir, [pet_QC_prefix, 'rmean*.nii']);   
                else
                    to_check{end+1} = get_files(check_dir, [pet_QC_prefix, 'rac*.nii']);     
                end 
            end
            % check reg
            mbatch = mbatch +1;
            matlabbatch{mbatch}.spm.util.checkreg.data = to_check';
            spm_jobman('run',matlabbatch(mbatch))
            pause
            close
        end
    end
end
end
% get folders function (get list of all folders in directory)
function folder_names = get_folders(rootPath)
allFiles = dir(rootPath);
%extract only those that are directories.
allDirFlags = [allFiles.isdir];
allFolders = allFiles(allDirFlags);
allFolders = allFolders(~ismember({allFolders(:).name},{'.','..'}));
folder_names = {allFolders.name};
end
%GET_FILES FUNCTION%
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
function table_to_run = spread_table_filter(sub_table, type_of_scan, field_str)
% function to do a bunch of the common table subsetting we do repeatedly in
% the SPREAD code. 
% make sure the rows are sorted by earliest date
sub_table = sortrows(sub_table, 'Date');
% drop any data not related to the current MRI modality
sub_table = sub_table(contains(sub_table.ScanType, type_of_scan), :);
% drop any data at the wrong field strength
sub_table = sub_table(sub_table.FieldStrength == field_str, :);
% finally, drop any scans which aren't a reference for
% their timepoint (i.e. a duplicate)
table_to_run = sub_table(sub_table.TimepointReference == 1, :);
end
function table_to_run = spread_pet_filter(sub_table, type_of_scan)
% function to do a bunch of the common table subsetting we do repeatedly in
% the SPREAD code. 
% make sure the rows are sorted by earliest date
sub_table = sortrows(sub_table, 'Date');
% drop any data not related to the current MRI modality
sub_table = sub_table(contains(sub_table.ScanType, type_of_scan), :);
% finally, drop any scans which aren't a reference for
% their timepoint (i.e. a duplicate)
table_to_run = sub_table(sub_table.TimepointReference == 1, :);
end