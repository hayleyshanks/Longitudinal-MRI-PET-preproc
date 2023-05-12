%% INCA Lab MRI/PET Preproc Step 0
% code by Hayley R.C. Shanks (hshanks@uwo.ca)
% edits by Kate Onuska, Taylor Schmitz
% last update: May 12 2023
% Western University

%% Required data organization:
% Organize data as follows:
% Particpant --> dates --> scan modality (e.g. PET) --> scan type (e.g.
% PiB) --> scan Name (e.g. scan1) --> folder with image files. 
% Participant data should live in a study specific folder under
% 'first_level'. See example.
%% note that for humans, date folders should be formtted YYYY-MM-DD
function longitudinal_preproc_step0_v1(N, study_path, study_name, field_str)
% get a list of all folders in the root path (these will be the subjects)
rootPath = [study_path, filesep, 'first_level'];
sheet_path = [study_path, filesep, 'spreadsheets'];
if ~exist(sheet_path, 'dir')
    mkdir(sheet_path)
end
RIDs = get_folders(rootPath);
sub_to_run = [1:length(RIDs)];
    for subI=sub_to_run(N)  
        subj_directory = fullfile(rootPath,RIDs{subI});
        sub_info = {};
        % Each folder will be named by the date the subject was scanned on
        scan_dates = get_folders(subj_directory);
        for i=1:length(scan_dates)
            date_dir = char(fullfile(subj_directory, scan_dates(i)));
            modality = get_folders(date_dir); % scan modalities are MRI or PET
            for j = 1:length(modality) 
                modality_dir = char(fullfile(date_dir, modality(j)));
                % check if the participant has any data for this modality at this time (PET and MRI directories are
                % created automatically for each date)
                image_types = get_folders(modality_dir);% e.g. 
                if isempty(image_types)
                    continue
                end
                % if there are image types in the folder,retreive the scan folder names
                for k=1:length(image_types)
                    image_dir = char(fullfile(modality_dir, image_types(k)));
                    scans = get_folders(image_dir); % list all scans (prob just one per date/modality combo)
                    if length(scans) > 1
                        % the participant was scanned on the same modality
                        % more than once per day. we can use this to create
                        % a within timepoint average (intra average).
                        duplicates = 1; % binary - is it duplicate or no
                    else
                        duplicates = 0;
                    end
                    for l=1:length(scans)
                        if startsWith(scans{l}, {'0','1','2', '3', '4', '5', '6', '7', '8', '9'})
                            % need to add a letter to the beginning of the
                            % folder to make sure the folder name will
                            % always be read in as a string. causes
                            % problems on Linux if not.
                            new_name = ['s', num2str(scans{l})];
                            % rename the folder
                            movefile([image_dir, filesep, scans{l}], [image_dir, filesep, new_name]);
                            scans{l} = new_name; % update the name in the scans list
                        end
                        scan_dir = char(fullfile(image_dir, scans(l)));
                        % get the magnetic field strength of the MRI data
                        dcm_set = get_files(scan_dir, '*.dcm');
                        if isempty(dcm_set)
                            % next step in preproc may have been called
                            % already
                            dcm_set = get_files(fullfile(scan_dir, 'dicom'), '*.dcm');
                        end
                        if ~isempty(dcm_set) 
                            dcm_set = cellstr(dcm_set);
                            hdr = dicominfo(dcm_set{1});
                            if isfield(hdr, 'MagneticFieldStrength')
                                fieldStrength = round(hdr.MagneticFieldStrength,1);
                            else
                                fieldStrength = NaN;
                            end
                        else
                             %% may also be working from nii files 
                            %% assume 3 T for now, since niftiinfo often does not have field strength info
                            if strcmp(modality{j}, 'MRI')
                                fieldStrength =3;
                            else
                                fieldStrength =1.5;
                            end
                        end                  
                        if strcmp(modality{j}, 'PET')
                            % add in what the target is of the PET tracer
                            % so we can separate these out later
                            if sum(strcmpi(image_types{k}, {'AV45', 'FLUT', 'NAV', 'PIB', 'amyloid', 'FBB'})) > 0
                                % these are the ADNI and AIBL amyloid
                                % tracers 
                                target = {'amyloid'};
                            elseif strcmpi(image_types{k}, 'FDG')
                                target = {'glucose'};
                            elseif sum(strcmpi(image_types{k}, {'AV1451', 'tau'})) > 0
                                target = {'tau'};
                            elseif sum(strcmpi(image_types{k}, {'FEOBV', 'FVAT', 'vacht'})) > 0
                                target = {'vacht'};
                            else
                                target = {'unknown'};
                                disp(['unknown tracer target for PET scan in ' scan_dir ' . Double check and manually update tracer targets in this code, line 98-110.'])
                            end
                        else
                            target = {};
                        end 
                        if isequal(duplicates, 0)
                            % don't want to include the duplicate sets
                            % variable as it was because it only applies to
                            % duplicates
                            to_add = {scan_dates(i), modality(j), image_types(k), scans(l), fieldStrength, duplicates,target};
                        else
                            to_add = {scan_dates(i), modality(j), image_types(k), scans(l), fieldStrength, duplicates,target};
                        end
                        sub_info = [sub_info; to_add];
                    end
                end 
            end 
        end
        if isempty(sub_info)
            % create an empty table for them
            names = {'Date' ,'Modality','ScanType' ,'ScanName', 'FieldStrength', 'IsDuplicate', 'TracerTarget'};
            sub_info_table = array2table(zeros(0,length(names)));
            sub_info_table.Properties.VariableNames = names;
            writetable(sub_info_table, fullfile(subj_directory, ['scan_info_',RIDs{subI} '.csv']));
        else
            sub_info_table = cell2table(sub_info, 'VariableNames', {'Date' , 'Modality', 'ScanType', 'ScanName', 'FieldStrength', 'IsDuplicate', 'TracerTarget'});
        end
        %% date time doesn't work for human data - changes the format which no longer matches the folders.
        %% this may need to be updated for some users
        % add the order of the scans for each modality 
        order_vec = zeros(height(sub_info_table),1);
        all_types = unique(sub_info_table.ScanType);
        % also add a variable encoding if a scan is the reference (first scan)
        % for that timepoint
        ref_vec = ones(height(sub_info_table),1);
        for t=1:length(all_types)
            wrong_field =0;
            % logical vec of whether a scan is the current scan type or no 
            log_type = strcmp(sub_info_table.ScanType, all_types{t});
            % need to create a vector which tells us scan order, where
            % order corresponds to the timepoint #. We may have duplicate
            % scans so we need to account for this. 
            % first check how many unique dates we have for the subject and
            % modalty. This will be the max order #
            %% also need to consider if this data is at the correct field str for MRI
            %% do not include in order if it isn't
            type_subset = sub_info_table(log_type, :); % table containing only the scan type we're looking at
            if sum(strcmpi(all_types{t}, {'T1', 'T2', 'T2_FLAIR', 'PD', 'PD_T2'})) > 0
                %MRI variable. Consider field strength
                %drop scans of the wrong field str
                corr_str = type_subset(type_subset.FieldStrength == field_str,:);
                max_order = length(unique(corr_str.Date));
                dates_to_check = unique(corr_str.Date);
                if ~isequal(height(corr_str), height(type_subset))
                    % they have some scans at the wrong field strength
                    wrong_field =1;
                end
            else
                max_order = length(unique(type_subset.Date));
                dates_to_check = unique(type_subset.Date);
            end
            if isequal(sum(type_subset.IsDuplicate), 0) && ~wrong_field
                n_vec = [1:max_order];
                curr_refs = ones(height(type_subset),1);
                % simple case of no duplicates
            elseif isequal(sum(type_subset.IsDuplicate), 0) && wrong_field
                n_vec= NaN(height(type_subset),1);
                n_vec(type_subset.FieldStrength == field_str) = 1:max_order;
                curr_refs = zeros(height(type_subset),1);
                curr_refs(type_subset.FieldStrength == field_str) = 1;
            else 
                % have duplicates. 
                n_vec = nan(height(type_subset),1);
                curr_refs = ones(height(type_subset),1);
                order_counter =0;
                dup_dates_processed = {};
                for row=1:height(type_subset)
                    if ~isequal(type_subset.FieldStrength(row), field_str) && wrong_field
                        % order is 0 because we will not include this scan
                        n_vec(row) =NaN;
                        curr_refs(row) = 0;
                    else
                        if isequal(type_subset.IsDuplicate(row),0) 
                            % this row is its own unique timepoint
                            order_counter = order_counter +1;
                            n_vec(row) = order_counter;                       
                        else                     
                            % it is a duplicate. Need to deal with a few
                            % situations
                            if ~ismember(type_subset.Date(row), dup_dates_processed)
                                % we have not yet encountered this timepoint.
                                % Increment the order by 1
                                order_counter = order_counter +1;
                                n_vec(row) = order_counter;
                                dup_dates_processed = [dup_dates_processed, type_subset.Date(row)]; % mark down that we processed this set
                            else
                                % we've seen a scan from this date before. Do
                                % not increment order
                                n_vec(row) = order_counter;
                                % mark down that this scan IS NOT a reference
                                % for any timepoint
                                curr_refs(row) = 0;
                            end
                        end
                    end
                end
            end
            if ~isequal(length(curr_refs), sum(log_type))
                disp('issue')
            end
            % vec of the max number of times for that scan type
            order_vec(log_type) = deal(n_vec); % deal the order numbers from n_vec back to the larger order vector containing all scan types
            ref_vec(log_type) = deal(curr_refs);
        end
        sub_info_table.TimepointReference = ref_vec;
        % order vec is now populated with each scan type's order. Save this
        % to the table
        sub_info_table.Order = order_vec;
        writetable(sub_info_table, fullfile(subj_directory, ['scan_info_',RIDs{subI} '.csv']));
        sub_col = cell(1,height(sub_info_table));
        sub_col(:) = {RIDs{subI}}; 
        % make RID the first col
        sub_info_table.RID = sub_col';
        sub_info_table = movevars(sub_info_table,'RID','Before','Date');
        sub_info_table = movevars(sub_info_table,'Order','Before','ScanName');
        if ~exist('masterTable', 'var')
            masterTable = sub_info_table;
        else
            masterTable = [masterTable; sub_info_table];
        end
    end 
    PET_table = masterTable(strcmp(masterTable.Modality, 'PET'), :);
    MRI_table = masterTable(strcmp(masterTable.Modality, 'MRI'), :);
    MRI_table = removevars(MRI_table, 'TracerTarget');
    PET_table = removevars(PET_table, 'FieldStrength');
    % before we write out a new table, check to see if either master table
    % exists. If it does, we may want to only add new subjects (not
    % overwrite old subjects who may have additional analyses saved in this
    % table
    old_pet = get_files(sheet_path, [study_name, '_PET_mastertable.csv']);
    old_mri = get_files(sheet_path, [study_name, '_MRI_mastertable.csv']);
    if ~isempty(old_pet)
        % update an existing table
        old_pet = readtable(old_pet);
        % the table may exist but have no data in it
        if ~isempty(old_pet)
            % if a participant had scans in the old pet table, then we just
            % select missing rows from new pet to add to old pet, the order
            % variable may have issues since all scans are not taken into
            % account. to avoid this we will use the newer table as our
            % reference table, and add columns from old_pet to the new
            % table
            
            % find any rows from the new PET table which are not in old pet by
            % checking RID, Date and Scan type combos
            old_pet_check = old_pet(:, {'RID', 'Date', 'ScanType'});
            new_pet_check = PET_table(:, {'RID', 'Date', 'ScanType'});
            % first make sure nothing from old pet is missing from new pet.
            % this would mean something was deleted by mistake probably.
            deleted = ~ismember(old_pet_check, new_pet_check,'rows');
            if sum(deleted > 0)
                disp('scans detected in previous mastersheet version which are missing from the newest iteration.')
                disp('the following scans may have been accidentally deleted:')
            end                
            % if ismember is true, then the rows overlap in the 2 tables
            check = ismember(new_pet_check, old_pet_check, 'rows');
            % for the rows where check is true we want to add any extra
            % values that have later been appended to mastersheet (e.g. qc
            % values
            % find columns not in the newest version of mastersheet
            % (PET_table)
            not_added = old_pet(:, ~contains(old_pet.Properties.VariableNames, PET_table.Properties.VariableNames)); 
            % need to add these columns to the out spreadsheet. these columns
            % need to be padded with blanks of the correct data type before we
            % can do this. 
            % go through each column and check its type
            % initialze an empty array
            to_add = array2table(NaN(height(PET_table), size(not_added,2)));
            to_add.Properties.VariableNames = not_added.Properties.VariableNames;
            for z=1:size(not_added,2)
                col_class = class(not_added{:,z});
                % update the column to be an empty data type of the correct
                % type or concatenation won't work
                if strcmp(col_class, 'cell')
                    to_add.(to_add.Properties.VariableNames{z})= cell(height(to_add),1);
                    % now add in any values from the previous table.
                    % Subjects with check beging true should have values
                    to_add.(to_add.Properties.VariableNames{z})(check) = old_pet.(to_add.Properties.VariableNames{z});
                elseif strcmp(col_class, 'double')
                    % can stay as nan
                    to_add.(to_add.Properties.VariableNames{z})(check) = old_pet.(to_add.Properties.VariableNames{z});
                end
            end
            % now horizontally concatenate empty_data to the PET_table
            % add to the out table 
            out_pet = [PET_table, to_add];
        else
             out_pet = PET_table;
        end
    else
        out_pet = PET_table;
    end
    %% a bit redundant but do the same for MRI
    if ~isempty(old_mri)
        old_mri = readtable(old_mri);
        if ~isempty(old_mri)
            old_mri_check = old_mri(:, {'RID', 'Date', 'ScanType'});
            new_mri_check = MRI_table(:, {'RID', 'Date', 'ScanType'});
            deleted = ~ismember(old_mri_check, new_mri_check,'rows');
            if sum(deleted > 0)
                disp('scans detected in previous mastersheet version which are missing from the newest iteration.')
                disp('the following scans may have been accidentally deleted:')
            end                
            check = ismember(new_mri_check, old_mri_check, 'rows');
            not_added = old_mri(:, ~contains(old_mri.Properties.VariableNames, MRI_table.Properties.VariableNames)); 
            to_add = array2table(NaN(height(MRI_table), size(not_added,2)));
            to_add.Properties.VariableNames = not_added.Properties.VariableNames;
            for z=1:size(not_added,2)
                col_class = class(not_added{:,z});
                % update the column to be an empty data type of the correct
                % type or concatenation won't work
                if strcmp(col_class, 'cell')
                    to_add.(to_add.Properties.VariableNames{z})= cell(height(to_add),1);
                    % now add in any values from the previous table.
                    % Subjects with check beging true should have values
                    to_add.(to_add.Properties.VariableNames{z})(check) = old_mri.(to_add.Properties.VariableNames{z});
                elseif strcmp(col_class, 'double')
                    % can stay as nan
                    to_add.(to_add.Properties.VariableNames{z})(check) = old_mri.(to_add.Properties.VariableNames{z});
                end
            end
            % now horizontally concatenate empty_data to the PET_table
            % add to the out table 
            out_mri = [MRI_table, to_add];
        else
             out_mri = MRI_table;
        end
    else
        out_mri = MRI_table;
    end
    % sort tables by RID and date before writing
    out_pet = sortrows(out_pet, {'RID', 'Date'}, 'ascend');
    out_mri = sortrows(out_mri, {'RID', 'Date'}, 'ascend');
    writetable(out_pet, fullfile(sheet_path, [study_name, '_PET_mastertable.csv']));
    writetable(out_mri, fullfile(sheet_path, [study_name, '_MRI_mastertable.csv']));
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