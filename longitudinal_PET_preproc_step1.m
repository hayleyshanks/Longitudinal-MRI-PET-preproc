%% PET Preproc Step 1
% code by Hayley R. C. Shanks, Taylor W. Schmitz, Kate M. Onuska
%% Not recommended for use with dynamic PET data. 
%% Multiple steps required for dynamic data (e.g. checking motion) are not implemented.
%% still need to explicitly exclude ppl who got 2 PET scans on 1 day
function longitudinal_PET_preproc_step1(N, study_path, pet_to_run, pet_orders, mri_ref, study_name, frames, field_str, qc_thresh, varargin)
%%%%%%%%%%%%%%%%%%% lab-specfic INPUTS - will need to be changed %%%%%%%%%%%%%%%%%%%%%%%%%
% add path to spm
addpath '/srv/schmitz/projects/toolboxes/spm12';
% dicom to nifti path
% available here: https://www.mathworks.com/matlabcentral/fileexchange/42997-xiangruili-dicm2nii
addpath('/srv/schmitz/projects/toolboxes/dicm2nii/2022.03.06');
functionRoot = '/srv/schmitz/projects/toolboxes/incalab_source_code/human';
% catch sight of the functions
addpath(genpath(functionRoot));
%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of specific inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
addRequired(pars, 'mri_ref', @iscell);
addRequired(pars, 'study_name', @ischar);
addRequired(pars, 'frames', @iscell); 
addRequired(pars, 'field_str', @isnumeric); 
addRequired(pars, 'qc_thresh', @isnumeric); 
%add optional inputs.
addParameter(pars,'dicomRecon',default_run,@isnumeric); 
addParameter(pars,'reorientAC',default_run,@isnumeric); 
addParameter(pars,'remove4D',default_run,@isnumeric); % dynamic data only
addParameter(pars,'checkVoxels',default_run,@isnumeric); 
addParameter(pars,'realignFrames',default_run,@isnumeric); %dynamic multiframe or longitudinal 3D/4D
addParameter(pars,'coregMR',default_run,@isnumeric); 
% parse the function inputs
parse(pars, N, study_path, pet_to_run, pet_orders, mri_ref, study_name, frames, field_str,qc_thresh,varargin{:});
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
    
%% check that the pet_order variables are set correctly 
if ~isnumeric(pet_orders{1}) && ~strcmpi(pet_orders{1}, 'all_times')
    disp(['input to pet_orders is set to ', pet_orders{1}, ' which is not a valid input'])
    disp('please review function documentation for the correct input arguments and try again')
    return
end
% figure out who fails MRI data QC or has none
mri_mastersheet = readtable([sheet_path, filesep, study_name, '_MRI_mastertable.csv']);
has_mri = intersect(mri_mastersheet.RID, RIDs);
qc_proportion = qc_thresh/100;
if strcmpi(mri_ref{1}, 'base')
    mri_fails = mri_mastersheet.RID(mri_mastersheet.timepoint_fails ==1 & mri_mastersheet.Order ==1);
elseif strcmpi(mri_ref{1}, 'avg')
    mri_fails = mri_mastersheet.RID(mri_mastersheet.avg_IQR_fails ==1);
else
    disp(['please review mri_ref description and re-define the variable'])
    return
end
all_fails = unique(mri_fails);
if pars.Results.dicomRecon 
    for subI=sub_to_run(N)
        sub_dir = fullfile(rootPath,RIDs{subI});
        sub_table = get_files(sub_dir, 'scan_info*.csv');
        % this table has all of the dates and types of scans for this P
        sub_table = readtable(sub_table);
        if isempty(sub_table)
            continue
        end
        sub_table = sortrows(sub_table, 'Date');
        for a=1:length(pet_to_run)
            % we only want to process one tracer at a time
            table_to_run = sub_table(contains(sub_table.ScanType, pet_to_run{a}), :);
            % check if the data is cross sectional
            if isequal(max(table_to_run.Order),1) && strcmpi(pet_orders{2}, 'excludeCrossSectionalSubjects')
                % skip them. they only have one timepoint of data and we
                % are excluding cross sectional subjects
                continue
            end
            if ~strcmpi(pet_orders{1}, 'all_times')
                % we're only processing some of the times.
                table_to_run = table_to_run(ismember(table_to_run.Order, pet_orders{1}), :);
            end
            if isempty(table_to_run)
                disp(['no applicable ' pet_to_run{a} ' data for ' char(RIDs{subI})])
                continue
            end 
            for i=1:height(table_to_run)
                scan_dir = fullfile(sub_dir, char(table_to_run.Date(i)), char(table_to_run.Modality(i)), char(table_to_run.ScanType(i)),char(table_to_run.ScanName(i)));
                if strcmpi(frames, {''})
                    dir_mod = '';
                else
                    possibilities = get_folders(scan_dir);
                    to_use = strcmpi(possibilities, frames);
                    dir_mod = [filesep, char(possibilities(to_use))];
                end
                scan_dir = [scan_dir, dir_mod];
                if isempty(get_files(scan_dir, '*.dcm'))
                    % has already been run
                    disp(['no dicoms found for ' , char(RIDs{subI}), ' on modality ' pet_to_run{a}])
                else
                    dicm2nii([scan_dir filesep '*.dcm'], scan_dir, '3D nii');
                    if ~exist(fullfile(scan_dir, 'dicom'), 'dir')
                        mkdir(fullfile(scan_dir, 'dicom'))
                    end
                    dicom_dir = fullfile(scan_dir, 'dicom');
                    movefile([scan_dir,filesep, '*.dcm'], dicom_dir);
                end
            end
        end
    end
end
if pars.Results.reorientAC
    %% filter by timepoint reference here if we are only using one PET scan per timepoint (no avging). can take this out if we will somehow use all scans
    % set the origin of the PET data to be around the anterior commissure.
    % Then set the file prefix of the data to be ac_ so we always can find
    % it 
    for subI=sub_to_run(N)
        if contains(RIDs{subI}, all_fails) || ~contains(RIDs{subI}, has_mri)
            %skip. we need mri data for pet preproc
        else
            sub_dir = fullfile(rootPath,RIDs{subI});
            sub_table = get_files(sub_dir, 'scan_info*.csv');
            % this table has all of the dates and types of scans for this P
            sub_table = readtable(sub_table);
            if isempty(sub_table)
                continue
            end
            for a=1:length(pet_to_run)
                table_to_run = spread_pet_filter(sub_table, pet_to_run{a});
                if ~isempty(table_to_run)
                    bl_dir = fullfile(sub_dir, char(table_to_run.Date(1)), char(table_to_run.Modality(1)), char(table_to_run.ScanType(1)),char(table_to_run.ScanName(1)));
                else
                    disp(['no applicable ' pet_to_run{a} ' data for ' char(RIDs{subI})])
                    continue
                end
                if max(table_to_run.Order) < 2 && strcmpi(pet_orders{2}, 'excludeCrossSectionalSubjects')
                    continue
                end
                for i=1:height(table_to_run)
                    scan_dir = fullfile(sub_dir, char(table_to_run.Date(i)), char(table_to_run.Modality(i)), char(table_to_run.ScanType(i)),char(table_to_run.ScanName(i)));
                    % data isn't always in a static or dynamic folder. this is the
                    % case if frames variable is empty               
                    if strcmpi(frames, {''})
                        dir_mod = '';
                    else
                        possibilities = get_folders(scan_dir);
                        to_use = strcmpi(possibilities, frames);
                        dir_mod = [filesep, char(possibilities(to_use))];
                    end
                    scan_dir = [scan_dir, dir_mod];
                    % get the scan file. at the beginning of preproc, there
                    % should be only one scan 
                    scan = cellstr(get_files(scan_dir, '*.nii'));
                    if strcmp(scan, '')
                        disp(['possibly non .nii data or missing data for ' RIDs{subI}, ' in folder ' , scan_dir])
                        continue
                    end
                    % with dynamic data, it may come as a 4d file or
                    % multiple 3d files. either way we want to only
                    % reoirent one frame then apply that matrix to all
                    if strcmpi(frames, 'dynamic') && isequal(length(scan),1)
                        %% data is in 1 4D file
                        is4d =1;
                        % if the scan is dynamic, we need to tell spm to open only
                        % the last frame. to do this we'll read in the hdr 
                        scan_hdr = spm_vol(char(scan));
                        % the size of hdr will tell us how many frames are in
                        % the dynamic scan (one row in the struct for each
                        % volume
                        num_vols = size(scan_hdr, 1);
                        % add this to the scan name so we can specify which
                        % frame to open
                        to_ac = cellstr([char(scan), ',' num2str(num_vols)]);
                        % this should give us the final frame to work from
                    elseif strcmpi(frames, 'dynamic') && length(cellstr(scan)) > 1
                        is4d =0;
                        to_ac = scan(end);
                    else
                        is4d =0;
                        to_ac = scan(end);
                    end
                    %% Display  image for visual inspection and set the origin to the AC
                    %% also make sure you rigidly allign the scan to the lorio (e.g. rotate pi radians on each axis and save the transformation
                    %% matrix. gives better starting estimate for coreg and seg
                    mbatch=mbatch+1;
                    if sum(contains(cellstr(scan),'ac_')) > 0 
                        % origin step has been run before. skip. may need to
                        % modify this if people will want to run multiple times
                        continue
                    end
                    % scan should really only point to one file if you're
                    % starting from scratch. if not, not super sure what to do.
                    if length(cellstr(scan)) > 1 && ~strcmpi(frames, 'dynamic')
                        disp(['multiple scans detected for ' RIDs{subI}, 'when only one file should exist. Skipping']) 
                    else
                        matlabbatch{mbatch}.spm.util.disp.data = to_ac;
                        spm_jobman('run',matlabbatch(mbatch))
                        disp('make sure you save any reorientation matrices. if your data is 4d, only reorient the last frame, and that matrix will be applied to all remaining frames')
                        fprintf('press enter when finished to advance to next scan \n')
                        pause
                        close
                        %% if our data is dynamic, we need to apply the reoirentation matrix to all the other frames
                        if strcmpi(frames, 'dynamic')
                            % load the reorientation matrix 
                            reor_mat = load(get_files(scan_dir, '*reorient*.mat'));
                            % this will load a matrix called M into the
                            % reor_mat struct
                            if is4d
                                % apply the matrix to the remaining frames
                                % of the same file. don't do the last frame
                                % because we did this one manually
                                for vol=1:(num_vols-1)
                                    to_update = [char(scan), ',', num2str(vol)];
                                    mbatch = mbatch +1;
                                    matlabbatch{mbatch}.spm.util.reorient.srcfiles = {to_update};
                                    matlabbatch{mbatch}.spm.util.reorient.transform.transM = reor_mat.M;
                                    matlabbatch{mbatch}.spm.util.reorient.prefix = '';
                                    spm_jobman('run',matlabbatch(mbatch))
                                end
                            else
                                % apply the matrix to the remaining
                                % individual files
                                for vol=1:(length(scan)-1)                                   
                                    mbatch = mbatch +1;
                                    matlabbatch{mbatch}.spm.util.reorient.srcfiles = scan(vol);
                                    matlabbatch{mbatch}.spm.util.reorient.transform.transM = reor_mat.M;
                                    matlabbatch{mbatch}.spm.util.reorient.prefix = '';
                                    spm_jobman('run',matlabbatch(mbatch))
                                end
                                clear reor_mat
                            end
                        end
                        % rename the file to be prefixed with ac --> signifies it's
                        % been oriented and also makes it easier later to grab the
                        % correct file
                        for x=1:(length(scan))                          
                            [path, name, ex] = fileparts(scan{x});
                            movefile(scan{x}, [path, filesep, 'ac_', name, ex]);
                            %create a backup as well, so that we can
                            %always start fresh from the original image if needed (with the
                            %origin set to the AC)
                            zip([path, filesep, 'ac_', name, '.zip'], [path, filesep, 'ac_', name, ex]);
                        end
                    end
                end
            end
        end
    end
end

if pars.Results.checkVoxels 
    % add columns of native voxel sizes to master sheet. check if they
    % exist first
    names_to_check = {'nativeVxX', 'nativeVxY', 'nativeVxZ'}; 
    if isequal(sum(contains(mastersheet.Properties.VariableNames, names_to_check)), 0)
        mastersheet.nativeVxX = nan(height(mastersheet),1);
        mastersheet.nativeVxY = nan(height(mastersheet),1);
        mastersheet.nativeVxZ = nan(height(mastersheet),1);
    end
    mastersheet.hasMRI = zeros(height(mastersheet),1);
    for a=1:length(pet_to_run)
        vx = [];
        for subI=sub_to_run(N)
            sub_dir = fullfile(rootPath,RIDs{subI});   
            if contains(RIDs{subI}, all_fails) || ~contains(RIDs{subI}, has_mri)
                %skip. we need mri data for pet preproc
            else
                sub_table = get_files(sub_dir, 'scan_info*.csv');
                % this table has all of the dates and types of scans for this P
                sub_table = readtable(sub_table);
                if isempty(sub_table)
                    continue
                end
                %make sure the rows are sorted by earliest date
                sub_table = sortrows(sub_table, 'Date');
                % subset only the current data modality
                table_to_run = sub_table(contains(sub_table.ScanType, pet_to_run{a}), :);
                if isempty(table_to_run)
                    disp(['no applicable ' pet_to_run{a} ' data for ' char(RIDs{subI})])
                    continue
                end
                for i=1:height(table_to_run)
                    scan_dir = fullfile(sub_dir, char(table_to_run.Date(i)), char(table_to_run.Modality(i)), char(table_to_run.ScanType(i)),char(table_to_run.ScanName(i)));
                    if strcmpi(frames, {''})
                        dir_mod = '';
                    else
                        possibilities = get_folders(scan_dir);
                        to_use = strcmpi(possibilities, frames);
                        dir_mod = [filesep, char(possibilities(to_use))];
                    end
                    frames_dir = [scan_dir, dir_mod];
                    % get the scan file. should now be prefixed with ac
                    scan = get_files(frames_dir, 'ac_*.nii');
                    if isempty(scan)
                        if isequal(table_to_run.TimepointReference(i),1)
                            disp(['no ac_' pet_to_run{a} ' scan found for ' char(RIDs{subI}), ' make sure you ran reorient AC'])
                        end
                    else
                        [~, vx_size ] = spm_get_bbox(scan);
                        % we want to check if the size of the voxels are different,
                        % not the LAS/RAS orientation
                        vx = [vx; abs(double(vx_size))];
                        % also save this info to mastersheet. find the correct row
                        % of mastersheet
                        check = [strcmp(mastersheet.RID, RIDs{subI}), mastersheet.Date == table_to_run.Date(i), strcmp(mastersheet.Modality,table_to_run.Modality(i)), strcmp(mastersheet.ScanType, table_to_run.ScanType(i)), strcmp(mastersheet.ScanName, table_to_run.ScanName(i))] ;
                        % which ever row of check is all 1's is our row
                        check_sum = sum(check, 2);
                        mastersheet.nativeVxX(check_sum == size(check,2)) = vx_size(1);
                        mastersheet.nativeVxY(check_sum == size(check,2)) = vx_size(2);
                        mastersheet.nativeVxZ(check_sum == size(check,2)) = vx_size(3);
                        % check if they have mri data
                        mri_table = spread_table_filter(sub_table, mri_ref{2}, field_str);
                        if ~isempty(mri_table)
                            mastersheet.hasMRI(check_sum == size(check,2)) = 1;
                        end
                    end
                end
            end
        end
        % all rows of the vx array should be equal across all subjects. use
        % ismember to check this. add a tolerance of 0.01
        check = ismembertol(vx, vx(1,:), 0.01, 'ByRows', true);
        disp('voxel sizes are: ');
        disp(unique(vx))
        if ~isequal(sum(check), size(check,1))
            disp(['UNEQUAL voxel sizes for ' char(pet_to_run{a}) ' across participants'])
        else
            disp(['EQUAL voxel sizes for ' char(pet_to_run{a}) ' across participants.'])
        end
    end
    % save the updated master sheet
    writetable(mastersheet, fullfile(sheet_path, [study_name, '_PET_mastertable.csv']));
end
if pars.Results.remove4D     
    for subI=sub_to_run(N)
        sub_dir = fullfile(rootPath,RIDs{subI});
        if contains(RIDs{subI}, all_fails) || ~contains(RIDs{subI}, has_mri)
            %skip. we need mri data for pet preproc
        else
            sub_table = get_files(sub_dir, 'scan_info*.csv');
            % this table has all of the dates and types of scans for this P
            sub_table = readtable(sub_table);
            if isempty(sub_table)
                continue
            end
            sub_table = sortrows(sub_table, 'Date');
            for a=1:length(pet_to_run)
                % only process the modality we're currently interested in
                table_to_run = spread_pet_filter(sub_table, pet_to_run{a});
                if ~isempty(table_to_run)
                    bl_dir = fullfile(sub_dir, char(table_to_run.Date(1)), char(table_to_run.Modality(1)), char(table_to_run.ScanType(1)),char(table_to_run.ScanName(1)));
                else
                    disp(['no applicable ' pet_to_run{a} ' data for ' char(RIDs{subI})])
                    continue
                end
                for i=1:height(table_to_run)
                    scan_dir = fullfile(sub_dir, char(table_to_run.Date(i)), char(table_to_run.Modality(i)), char(table_to_run.ScanType(i)),char(table_to_run.ScanName(i)));
                    % data isn't always in a static or dynamic folder. this is the
                    % case if frames variable is empty               
                    if strcmpi(frames, {'static'})
                        % skip. doesn't apply to static data
                    else
                        possibilities = get_folders(scan_dir);
                        to_use = strcmpi(possibilities, 'dynamic'); % the only point of doing this is to see if file path is in
                        % caps or not
                        dir_mod = [filesep, char(possibilities(to_use))];
                        scan_dir = [scan_dir, dir_mod];
                        % get the scan file. at the beginning of preproc, there
                        % should be only one scan 
                        to_split = get_files(scan_dir, 'ac_*.nii');
                        y4D = spm_vol(to_split);
                        if length(y4D(:,1)) == 1
                            disp(['no dynamic 4D data for ' char(RIDs{subI})])
                            continue
                        else
                            mbatch = mbatch+1;
                            % split the file
                            matlabbatch{mbatch}.spm.util.split.vol = {to_split};
                            matlabbatch{mbatch}.spm.util.split.outdir = {''};
                            spm_jobman('run',matlabbatch(mbatch))
                            [pt, fn, ex] = fileparts(scan);
                            movefile(scan, [pt, filesep, 'orig4D_', fn, ex]);
                            %create a backup as well
                            zip([pt, filesep, 'orig4D_', fn, '.zip'], [pt, filesep, 'orig4D_', fn, ex]);
                        end
                    end
                end
            end
        end
    end
end
if pars.Results.realignFrames
    for subI=sub_to_run(N)
        if contains(RIDs{subI}, all_fails) || ~contains(RIDs{subI}, has_mri)
            %skip. we need mri data for pet preproc
            continue
        else
            skip =0;
            sub_dir = fullfile(rootPath,RIDs{subI});
            sub_table = get_files(sub_dir, 'scan_info*.csv');
            % this table has all of the dates and types of scans for this P
            sub_table = readtable(sub_table);
            if isempty(sub_table)
                continue
            end
            for a=1:length(pet_to_run)
                table_to_run = spread_pet_filter(sub_table, pet_to_run{a});
                if ~isempty(table_to_run)
                    bl_dir = fullfile(sub_dir, char(table_to_run.Date(1)), char(table_to_run.Modality(1)), char(table_to_run.ScanType(1)),char(table_to_run.ScanName(1)));
                else
                    disp(['no applicable ' pet_to_run{a} ' data for ' char(RIDs{subI})])
                    continue
                end
                % we will only re-align frames if data is dynamic or if static
                % data is longitudinal
                if isequal(max(table_to_run.Order),1) && strcmpi(frames, 'static')
                    % skip. nothing to do here
                    disp(['only 1 timepoint of static PET data for ' RIDs{subI}, ' nothing to re-align'])
                    continue
                elseif max(table_to_run.Order) > 1 && ~strcmpi(frames, 'dynamic')
                   % realign scans across all timepoints 
                    scan_array = {};
                    all_dirs = {};             
                    for i=1:height(table_to_run)
                        scan_dir = fullfile(sub_dir, char(table_to_run.Date(i)), char(table_to_run.Modality(i)), char(table_to_run.ScanType(i)),char(table_to_run.ScanName(i)));
                        if ~strcmp(frames, '')
                            frames_dir = get_folders(scan_dir);
                            to_use = strcmpi(frames_dir, frames);
                            mod = frames_dir(to_use);
                        else
                            mod = '';
                        end
                        final_path = [scan_dir, filesep, char(mod)];
                        all_dirs{end+1} = final_path; % keep track of folder paths we used
                        % check if realign was already run
                        mean_file = get_files(final_path, 'meanac*.nii');
                        if exist(mean_file, 'file')
                            skip =1;
                            disp(['re-align has already been run for ' RIDs{subI} '. They will be skipped'])
                            continue
                        else
                            curr_scans = cellstr(get_files(final_path, 'ac*.nii'));
                            scan_array{i} = curr_scans;
                        end
                    end
                    if ~skip
                        % now that we have all the scans, let's realign
                        mbatch = mbatch+1;
                        %% run realignment. If there are multiple timepoints of PET data, it will treat these
                        %% as separate sessions and perform a two pass realignment.
                        %% pass 1 realign all data within each session to the first image of that session
                        %% pass 2 create a mean of the first images of each session and realign everything to that
                        % note that we won't reslice the idnividual pet images
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.data = scan_array;
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.eoptions.sep = 4;
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.eoptions.fwhm = 5; % recommended fwhm is 7 for PET actually
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.eoptions.interp = 4;
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.eoptions.weight = '';
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.roptions.which = [0 1];
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.roptions.interp = 4;
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.roptions.mask = 1;
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
                        spm_jobman('run',matlabbatch(mbatch))
                        %% motion parameter stuff doesn't make sense for static data at diff timepoints
                    end
                elseif strcmpi(frames, 'dynamic')
                    % if data is dynamic, we need to re-align frames within
                    % a timepoint and also across timepoints. If we put all
                    % dynamic data together across all times, we only get 1
                    % mean image (which reflects many times of data). we
                    % don't want this for SUVR
                    for i=1:height(table_to_run)
                        scan_dir = fullfile(sub_dir, char(table_to_run.Date(i)), char(table_to_run.Modality(i)), char(table_to_run.ScanType(i)),char(table_to_run.ScanName(i)));
                        frames_dir = get_folders(scan_dir);
                        to_use = strcmpi(frames_dir, frames);
                        mod = frames_dir(to_use);
                        final_path = [scan_dir, filesep, char(mod)];
                        mean_file = get_files(final_path, 'meanac*.nii');
                        if exist(mean_file, 'file')
                            disp(['re-align has already been run for ' RIDs{subI} '. They will be skipped'])
                            continue
                        else
                            curr_scans = cellstr(get_files(final_path, 'ac*.nii'));
                        end
                        mbatch = mbatch+1;
                        %% run realignment. If there are multiple timepoints of PET data, it will treat these
                        %% as separate sessions and perform a two pass realignment.
                        %% pass 1 realign all data within each session to the first image of that session
                        %% pass 2 create a mean of the first images of each session and realign everything to that
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.data = {curr_scans};
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.eoptions.sep = 4;
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.eoptions.fwhm = 5; % recommended fwhm is 7 for PET actually
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.eoptions.interp = 4;
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.eoptions.weight = '';
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.roptions.which = [0 1];
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.roptions.interp = 4;
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.roptions.mask = 1;
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
                        spm_jobman('run',matlabbatch(mbatch))
                        rpFiles = get_files(final_path,'rp*.txt');
                        rpMat = importdata(rpFiles); 
                        [mat_path, mat_name, ~] = fileparts(rpFiles);
                        save([mat_path, filesep mat_name '_motionParam.mat'],'rpMat');
                    end
                    if max(table_to_run.Order) > 1
                        % also now run realignment across the timepoints of
                        % dynamic data
                        scan_array = {};
                        all_dirs = {};             
                        for i=1:height(table_to_run)
                            scan_dir = fullfile(sub_dir, char(table_to_run.Date(i)), char(table_to_run.Modality(i)), char(table_to_run.ScanType(i)),char(table_to_run.ScanName(i)));
                            if ~strcmp(frames, '')
                                frames_dir = get_folders(scan_dir);
                                to_use = strcmpi(frames_dir, frames);
                                mod = frames_dir(to_use);
                            else
                                mod = '';
                            end
                            final_path = [scan_dir, filesep, char(mod)];
                            all_dirs{end+1} = final_path; % keep track of folder paths we used
                            % also realign the mean image
                            mean_file = cellstr(get_files(final_path, 'meanac*.nii'));
                            curr_scans = cellstr(get_files(final_path, 'ac*.nii'));
                            scan_array{i} = [curr_scans; mean_file];
                        end
                        % now that we have all the scans, let's realign
                        mbatch = mbatch+1;
                        %% run realignment. If there are multiple timepoints of PET data, it will treat these
                        %% as separate sessions and perform a two pass realignment.
                        %% pass 1 realign all data within each session to the first image of that session
                        %% pass 2 create a mean of the first images of each session and realign everything to that
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.data = scan_array;
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.eoptions.sep = 4;
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.eoptions.fwhm = 5; % recommended fwhm is 7 for PET actually
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.eoptions.interp = 4;
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.eoptions.weight = '';
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.roptions.which = [0 1];
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.roptions.interp = 4;
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.roptions.mask = 1;
                        matlabbatch{mbatch}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
                        spm_jobman('run',matlabbatch(mbatch))
                    end
                end
            end
        end
    end
end

if pars.Results.coregMR
    for subI=sub_to_run(N)
        if contains(RIDs{subI}, all_fails) || ~contains(RIDs{subI}, has_mri)
            %skip. we need mri data for pet preproc
        else
            sub_dir = fullfile(rootPath,RIDs{subI});
            sub_table = get_files(sub_dir, 'scan_info*.csv');
            % this table has all of the dates and types of scans for this P
            sub_table = readtable(sub_table);
            %make sure the rows are sorted by earliest date
            if isempty(sub_table)
                continue
            end
            mri_table = spread_table_filter(sub_table, mri_ref{2}, field_str); % mri_ref tells us if we want T1, t2, etc.   
            if isempty(mri_table)
                continue % some people have MRI data at the wrong field strength
            end
            mri_bl_dir = fullfile(sub_dir, char(mri_table.Date(1)), char(mri_table.Modality(1)), char(mri_table.ScanType(1)), char(mri_table.ScanName(1)));
            if strcmpi(mri_ref{1}, 'avg')
                mri_ref_dir = [mri_bl_dir, filesep 'midpoint_average' filesep, 'mri' ];
                mr_ref_img = get_files(mri_ref_dir, 'mavg*.nii');
            elseif strcmpi(mri_ref{1}, 'base')
                mri_ref_dir = [mri_bl_dir, filesep, 'mri'];
                mr_ref_img = get_files(mri_ref_dir, 'mintra*.nii');
                if isempty(mr_ref_img)
                    mr_ref_img = get_files(mri_ref_dir, 'mrac*.nii');
                end
            end
            if isempty(mr_ref_img)
                disp(['no usable reference MRI image for ' RIDs{subI},'. Their PET preproc will be skipped'])
            else
                for a=1:length(pet_to_run)
                    % only process the modality we're currently interested in
                    pet_table = spread_pet_filter(sub_table, pet_to_run{a});
                    if isempty(pet_table)
                        disp(['no applicable ' pet_to_run{a} ' data for ' char(RIDs{subI})])
                        continue
                    else
                        for i=1:height(pet_table)
                            scan_dir = fullfile(sub_dir, char(pet_table.Date(i)), char(pet_table.Modality(i)), char(pet_table.ScanType(i)),char(pet_table.ScanName(i)));
                            % data isn't always in a static or dynamic folder. this is the
                            % case if frames variable is empty               
                            if strcmpi(frames, {''})
                                dir_mod = '';
                            else
                                possibilities = get_folders(scan_dir);
                                to_use = strcmpi(possibilities, frames);
                                dir_mod = [filesep, char(possibilities(to_use))];
                            end
                            scan_dir = [scan_dir, dir_mod];               
                            % grab the rac file for coreg. this means it's been
                            % realigned
                            if strcmpi(frames, 'dynamic')
                                % we have multiple files to deal with. We want
                                % to coreg the multiple dynamic frames from one
                                % timepoint to the MRI ref using the meanscan
                                % as the PET source file
                                mean_img = get_files(scan_dir, 'meanac*.nii');
                                pet_images = get_files(scan_dir, 'ac*.nii');
                                mbatch = mbatch+1;
                                matlabbatch{mbatch}.spm.spatial.coreg.estwrite.ref = {mr_ref_img};
                                matlabbatch{mbatch}.spm.spatial.coreg.estwrite.source = {mean_img};
                                matlabbatch{mbatch}.spm.spatial.coreg.estwrite.other = cellstr(pet_images);
                                matlabbatch{mbatch}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
                                matlabbatch{mbatch}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
                                matlabbatch{mbatch}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                                matlabbatch{mbatch}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
                                matlabbatch{mbatch}.spm.spatial.coreg.estwrite.roptions.interp = 4;
                                matlabbatch{mbatch}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
                                matlabbatch{mbatch}.spm.spatial.coreg.estwrite.roptions.mask = 0;
                                matlabbatch{mbatch}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
                                spm_jobman('run',matlabbatch(mbatch))
                            else
                                % if their data is static, we will only have
                                % one scan per timepoint. Take the re-aligned
                                % file
                                to_coreg = get_files(scan_dir, 'ac*.nii');
                                mbatch = mbatch+1;
                                matlabbatch{mbatch}.spm.spatial.coreg.estwrite.ref = {mr_ref_img};
                                matlabbatch{mbatch}.spm.spatial.coreg.estwrite.source = {to_coreg};
                                matlabbatch{mbatch}.spm.spatial.coreg.estwrite.other = {''};
                                matlabbatch{mbatch}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
                                matlabbatch{mbatch}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
                                matlabbatch{mbatch}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                                matlabbatch{mbatch}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
                                matlabbatch{mbatch}.spm.spatial.coreg.estwrite.roptions.interp = 4;
                                matlabbatch{mbatch}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
                                matlabbatch{mbatch}.spm.spatial.coreg.estwrite.roptions.mask = 0;
                                matlabbatch{mbatch}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
                                spm_jobman('run',matlabbatch(mbatch))
                            end
                        end
                    end                                    
                end
            end
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