%% MRI Preprocessing Step 1:
%% Co-registration, longitudinal registration, segmentation
% code by Hayley R. C. Shanks, Taylor W. Schmitz, Kate M. Onuska
% report any bugs to hshanks@uwo.ca

% segmentation requires you've downloaded the enhanced TPM files from
% Lorio et al. 2016 and placed them in the SPM TPM directory
% paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4819722/

% CAT12 needs to be run in expert mode --> set cat.extopts.expertgui = 1;
% in line 258 of cat_defaults.m
function longitudinal_MRI_preproc_step1(N, study_path, mri_to_run, study_name, avgMode, indMode, surf_mod, field_str, varargin)

%%%%%%%%%%%%%%%%%%% lab-specfic INPUTS - may need to be changed %%%%%%%%%%%%%%%%%%%%%%%%%
% add path to spm
addpath '/srv/schmitz/projects/toolboxes/spm12';
% dicom to nifti path
% available here: https://www.mathworks.com/matlabcentral/fileexchange/42997-xiangruili-dicm2nii
addpath('/srv/schmitz/projects/toolboxes/dicm2nii');
functionRoot = '/srv/schmitz/projects/toolboxes/incalab_source_code/human';
% catch sight of the functions
addpath(genpath(functionRoot));
% path to the SPM12 tpms (here Lorio 2016 tpm will be used)
templatePath = '/srv/schmitz/projects/toolboxes/spm12/tpm';
% path to the CAT12 MNI geodesic shoot template
CAT12_shoot_path = '/srv/schmitz/projects/toolboxes/spm12/toolbox/cat12/templates_MNI152NLin2009cAsym';
%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of specific inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
study_avg_dir = [study_path, filesep, study_name, '_averages'];
% raw participant data should live in a first level subdirectory of the study path
rootPath = [study_path, filesep, 'first_level'];
sheet_path = [study_path, filesep, 'spreadsheets'];
backup_avg_dir = '/srv/schmitz/projects/toolboxes/spm12/canonical'; % switch this to whatever your spm path is
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parse user inputs
pars = inputParser;
% by default, all modules are set to off.
default_run =0;
% required inputs
addRequired(pars, 'N', @isnumeric);
addRequired(pars, 'study_path', @ischar);
addRequired(pars, 'mri_to_run', @iscell);
addRequired(pars, 'study_name', @ischar);
addRequired(pars, 'avgMode', @iscell);
addRequired(pars, 'indMode', @iscell);
addRequired(pars, 'surf_mod', @isnumeric);
addRequired(pars, 'field_str', @isnumeric);
%add optional inputs.
% dicom to nii reconstruction
addParameter(pars,'dicomRecon',default_run,@isnumeric); 
% set origin of images to the anterior commisure (manually)
addParameter(pars,'reorientAC',default_run,@isnumeric); 
% check whether voxel sizes across participants are the same (helps decide
% if you should perform coregistration (estimate only) or coregistration
% with estimation and reslicing)
addParameter(pars,'checkVoxels',default_run,@isnumeric); 
addParameter(pars,'coreg',default_run,@isnumeric); 
% creation of symmetric midpoint average image. Recommended for
% longitudinal studies
addParameter(pars,'longitudinalRegistration',default_run,@isnumeric); 
addParameter(pars,'segment',default_run,@isnumeric); 
% parse the function inputs
parse(pars, N, study_path, mri_to_run, study_name, avgMode, indMode, surf_mod,field_str, varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get a list of all files and folders in the root path folder.
RIDs = get_folders(rootPath);
mbatch = 0;
sub_to_run = [1:length(RIDs)];
% read in the master info sheet
mastersheet = readtable(fullfile(sheet_path, [study_name, '_MRI_mastertable.csv']));
% check that the mode variables are set correctly 
if isequal(sum(strcmpi(avgMode, {'on', 'off'})),0)
    disp(['input to avgMode is set to ', avgMode{1}, ' which is not a valid input'])
    disp('please review function documentation for the correct input arguments and try again')
    return
elseif ~isnumeric(indMode{1}) && isequal(sum(strcmpi(indMode{1}, {'all_times', 'off'})),0)
    disp(['input to indMode is set to ', indMode{1}, ' which is not a valid input'])
    disp('please review function documentation for the correct input arguments and try again')
    return
end
if isequal(length(indMode),1) && strcmpi(indMode, 'off')
    % add a second cell just to make code easier    
    indMode{1,2} = 'off';
elseif isequal(length(indMode),2) && isequal(sum(strcmpi(indMode{2}, {'excludeCrossSectionalSubjects', 'includeCrossSectionalSubjects'})),0)
    disp(['input to indMode is set to ', indMode{2}, ' which is not a valid input'])
    disp('please review function documentation for the correct input arguments and try again')
    return
end
if pars.Results.dicomRecon || pars.Results.reorientAC
    for subI=sub_to_run(N)
        sub_dir = fullfile(rootPath,RIDs{subI});
        sub_table = mastersheet(strcmp(mastersheet.RID, RIDs(subI)), :);
        if isempty(sub_table)
            continue
        end
        sub_table = sortrows(sub_table, 'Date'); % sort by earliest date
        for a=1:length(mri_to_run)
            % only process the modality we're currently interested in
            table_to_run = sub_table(contains(sub_table.ScanType, mri_to_run{a}), :);
            % drop any data at the wrong field strength
            table_to_run = table_to_run(table_to_run.FieldStrength == field_str, :);
            if isequal(max(table_to_run.Order),1) && strcmpi(indMode{2}, 'excludeCrossSectionalSubjects')
                % skip them. they only have one timepoint of data and we
                % are excluding cross sectional subjects
                continue 
            end
            if isequal(max(table_to_run.Order),1) && strcmpi(avgMode, 'on')
                % skip them. they only have one timepoint of data and we
                % are using longitudinal data in avgMode
                continue 
            end
            if ~strcmpi(indMode{1}, 'all_times') && isnumeric(indMode{1}) && strcmpi(avgMode, 'off')
                % we're only processing some of the times. 
                table_to_run = table_to_run(ismember(table_to_run.Order, indMode{1}), :);
            end
            
            if ~isempty(table_to_run)
                bl_dir = fullfile(sub_dir, char(table_to_run.Date(1)), char(table_to_run.Modality(1)), char(table_to_run.ScanType(1)),char(table_to_run.ScanName(1)));
            else
                disp(['no applicable ' mri_to_run{a} ' data for ' char(RIDs{subI})])
                continue
            end
            for i=1:height(table_to_run)
                scan_dir = fullfile(sub_dir, char(table_to_run.Date(i)), char(table_to_run.Modality(i)), char(table_to_run.ScanType(i)),char(table_to_run.ScanName(i)));
                if pars.Results.dicomRecon
                    if isempty(get_files(scan_dir, '*.dcm'))
                        % has already been run
                        disp(['no dicoms found for ' , char(RIDs{subI}), ' on modality ' mri_to_run{a}])
                    else
                        if contains(scan_dir, 'PD_T2') 
                            %% this needs to be fixed - ideally the images should be moved so each has their own folder
                            % in ADNI and AIBL, they did a combined PD/T2 sMRI
                            % scan. These files live in the same dicom folder but
                            % need to be separated for analyses
                            all_dicoms = cellstr(get_files(scan_dir, '*.dcm'));
                            endings = unique(extractAfter(all_dicoms, ['_' char(table_to_run.ScanName(i)) '_']));
                            for ending=1:size(endings,1)
                                dicm2nii([scan_dir filesep '*' char(endings(ending,:))], scan_dir, '3D nii');
                                %rename the dicom header mat file that gets created
                                %so  it won't be overwritten by the second call
                                [~, mat_name, ext] = fileparts(fullfile(scan_dir, 'dcmHeaders.mat'));
                                other_mat = [mat_name, '_', char(endings(ending, :))];
                                movefile([scan_dir, filesep, mat_name, ext], [scan_dir, filesep, other_mat, ext]);
                                dcm_set = get_files(scan_dir, ['*' char(endings(ending, :))]);
                                hdr = dicominfo(dcm_set(1,:));
                                TE = num2str(hdr.EchoTime);
                                % if we're on the second reconstruction, the nifti
                                % will be the file which does not contain _TE_
                                nii_images = get_files(scan_dir, '.nii');
                                ind = contains(cellstr(nii_images), '_TE_');
                                [~, nii_name, nii_ext] = fileparts( nii_images(~ind,:));
                                new_nii_name = [nii_name, '_TE_', TE];
                                movefile([scan_dir, filesep, nii_name, nii_ext], [scan_dir, filesep, new_nii_name, nii_ext]);
                            end
                        else
                            dcm_set = get_files(scan_dir, '*.dcm');
                            dcm_set = cellstr(dcm_set);
                            hdr = dicominfo(dcm_set{1,:});
                            dicm2nii([scan_dir filesep '*.dcm'], scan_dir, '3D nii');
                        end
                        if ~exist(fullfile(scan_dir, 'dicom'), 'dir')
                            mkdir(fullfile(scan_dir, 'dicom'))
                        end
                        dicom_dir = fullfile(scan_dir, 'dicom');
                        movefile([scan_dir,filesep, '*.dcm'], dicom_dir);
                        % note: dicm2nii flips the z and x of the dicoms
                        % (compared to what is reconstructed in spm)
                    end
                end
                if pars.Results.reorientAC
                    mbatch=mbatch+1;
                    scan=get_files(fullfile(scan_dir),'*.nii');
                    if sum(contains(cellstr(scan),'ac_')) > 0
                        % origin step has been run before. skip
                        continue
                    end
                    for k=1:size(scan,1)
                        disp(subI)
                        %% Display  image for visual inspection and set the origin to the AC
                        %% also make sure you rigidly allign the scan to the lorio (e.g. rotate pi radians on each axis and save the transformation
                        %% matrix. gives better starting estimate for coreg and seg
                        matlabbatch{mbatch}.spm.util.disp.data = {scan(k,:)};
                        spm_jobman('run',matlabbatch(mbatch))
                        pause
                        close
                        % rename the file to be prefixed with ac --> signifies it's
                        % been oriented and also makes it easier later to grab the
                        % correct file
                        [path, name, ex] = fileparts(scan(k,:));
                        movefile(scan(k, :), [path, filesep, 'ac_', name, ex]);
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
if pars.Results.checkVoxels 
    for t=1:length(mri_to_run)
        % add a column of native voxel sizes to master sheet
        names_to_check = {'nativeVxX', 'nativeVxY', 'nativeVxZ'}; 
        if isequal(sum(contains(mastersheet.Properties.VariableNames, names_to_check)), 0)
            mastersheet.nativeVxX = nan(height(mastersheet),1);
            mastersheet.nativeVxY = nan(height(mastersheet),1);
            mastersheet.nativeVxZ = nan(height(mastersheet),1);
        end
        for t=1:length(mri_to_run)
            vx = [];
            for subI=sub_to_run(N)
                sub_dir = fullfile(rootPath,RIDs{subI});
                sub_table = mastersheet(strcmp(mastersheet.RID, RIDs(subI)), :);
                if isempty(sub_table)
                    continue
                end
                %make sure the rows are sorted by earliest date
                sub_table = sortrows(sub_table, 'Date');
                % subset only the current data modality
                anat = sub_table(strcmp(sub_table.ScanType, mri_to_run{t}), :);
                anat = anat(anat.FieldStrength == field_str, :);
                if isequal(max(anat.Order),1) && strcmpi(indMode{2}, 'excludeCrossSectionalSubjects')
                    % skip them. they only have one timepoint of data and we
                    % are excluding cross sectional subjects
                    continue 
                end
                if ~strcmpi(indMode{1}, 'all_times') && isnumeric(indMode{1}) && strcmp(avgMode, 'off')
                    % we're only processing some of the times
                    anat = anat(ismember(anat.Order, indMode{1}), :);
                end
                if isempty(anat)
                    continue
                end
                for j=1:height(anat)
                    current_direc = fullfile(sub_dir, char(anat.Date(j)), char(anat.Modality(j)), char(anat.ScanType(j)),char(anat.ScanName(j)));
                    % read in the image and get the voxel sizes
                    image = get_files(current_direc, 'ac*.nii');
                    if isempty(image) %% added 04/28/2022 KO
                        disp('no image detected')
                    else %% added 04/28/2022 KO
                    [~, vx_size ] = spm_get_bbox(image);
                    % we want to check if the size of the voxels are different,
                    % not the LAS/RAS orientation
                    vx = [vx; abs(double(vx_size))];
                    % also save this info to mastersheet. find the correct row
                   % mastersheet
                    check = [strcmp(mastersheet.RID, RIDs{subI}), mastersheet.Date == anat.Date(j), strcmp(mastersheet.Modality,anat.Modality(j)), strcmp(mastersheet.ScanType, anat.ScanType(j)), strcmp(mastersheet.ScanName, anat.ScanName(j))] ;
                    % which ever row of check is all 1's is our row
                    check_sum = sum(check, 2);
                    mastersheet.nativeVxX(check_sum == size(check,2)) = vx_size(1);
                    mastersheet.nativeVxY(check_sum == size(check,2)) = vx_size(2);
                    mastersheet.nativeVxZ(check_sum == size(check,2)) = vx_size(3);
                    end
                end
            end
            % all rows of the vx array should be equal across all subjects. use
            % ismember to check this. add a tolerance of 0.01
            check = ismembertol(vx, vx(1,:), 0.01, 'ByRows', true);
            disp('voxel sizes are: ');
            disp(unique(vx))
            if ~isequal(sum(check), size(check,1))
                disp('UNEQUAL voxel sizes across participants. We recommend the data undergo co-register: estimate and reslice prior to average creation')
            else
                disp('EQUAL voxel sizes across participants. We recommend the data undergo coregister: estimate only prior to average creation')
            end
            user_ip = input('Enter 1 if you would like to re-slice or 0 if you would like to co-register estimate only \n');
            if isequal(user_ip, 1)
                voxel_check = 'reslice';
            elseif isequal(user_ip, 0)
                voxel_check ='est only';
            else
                disp('invalid input. please run the checkVoxels step again')
            end
            % save the status of the voxel check to a mat file so we know
            % later if data needs to be resliced or just coreged
            save([rootPath, filesep, mri_to_run{t},'_voxel_check.mat'], 'voxel_check');
        end
        % save the updated master sheet
        writetable(mastersheet, fullfile(sheet_path, [study_name, '_MRI_mastertable.csv']));
    end
end

%% now we begin main cross-subject preproc loop
for subI=sub_to_run(N)
    sub_dir = fullfile(rootPath,RIDs{subI});
    sub_table = mastersheet(strcmp(mastersheet.RID, RIDs(subI)), :);
    %make sure the rows are sorted by earliest date
    sub_table = sortrows(sub_table, 'Date');
    if isempty(sub_table)
        continue
    end
    for t=1:length(mri_to_run)
        if pars.Results.coreg
            % use the voxel_check variable to determine if the data needs
            % to be resliced or not
            load([rootPath, filesep, mri_to_run{t},'_voxel_check.mat']);
            anat_table = sub_table(contains(sub_table.ScanType, mri_to_run{t}), :) ;
            anat_table = anat_table(anat_table.FieldStrength == field_str, :);
            if isequal(max(anat_table.Order),1) && strcmpi(indMode{2}, 'excludeCrossSectionalSubjects')
                % skip them. they only have one timepoint of data and we
                % are excluding cross sectional subjects
                continue 
            elseif isequal(max(anat_table.Order),1) && strcmpi(indMode{1}, 'off')
                % no cross sectional data in avgmode
                continue
            end
            if ~strcmpi(indMode{1}, 'all_times') && isnumeric(indMode{1}) && strcmp(avgMode, 'off')
                % we're only processing some of the times
                anat_table = anat_table(ismember(anat_table.Order, indMode{1}), :);
            end
            % during coregistration, we will use some kind of reference average
            % image to coregister to. Ideally it should be an image similar to your
            % study population. If you have run preproc for some subjects
            % before, yoy may already have a study average image
            average_anat = get_files(study_avg_dir, [study_name, '_', mri_to_run{t}, '_wmavg_', '*.nii']);
            % if there's no average image, coreg to a template image 
            if isempty(average_anat)
                average_anat = get_files(backup_avg_dir, ['avg152', mri_to_run{t}, '.nii']);
                if isempty(average_anat)
                    % if it is still empty, display an error message. code
                    % will fail. 
                    disp('no reference/average image found that can be used for co-registration')
                    disp('please update line 32 and 318 to point to a reference image');
                end
            end
            if ~isempty(anat_table)
                bl_dir = fullfile(sub_dir, char(anat_table.Date(1)), char(anat_table.Modality(1)), char(anat_table.ScanType(1)),char(anat_table.ScanName(1)));
            else
                disp(['no ' mri_to_run{t} ' data for ' RIDs{subI}])
                continue
            end
            to_coreg ={};
            for row=1:height(anat_table)
                current_dir = fullfile(sub_dir, char(anat_table.Date(row)), char(anat_table.Modality(row)), char(anat_table.ScanType(row)),char(anat_table.ScanName(row)));
                % unzip to use a fresh nii file
                nii_zip = get_files(current_dir, ['ac_*.zip']);
                if ~isempty(nii_zip)
                    unzip(nii_zip, [current_dir, filesep]);
                end
                to_coreg = get_files(current_dir, ['ac_*.nii']);
                % coreg est reslice of the images to the study average
                % anat image for each timepoint
                if strcmp(voxel_check, 'reslice')
                    % coreg est res because not all images have the same
                    % voxel sizes
                    mbatch = mbatch +1;
                    %% added tws 12/28/2021
                    disp(['performing coreg & reslice for ' RIDs{subI}])
                    matlabbatch{mbatch}.spm.spatial.coreg.estwrite.ref = {average_anat};
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
                elseif strcmp(voxel_check, 'est only')
                    %coreg est only
                    mbatch = mbatch + 1;
                    %% added tws 12/28/2021
                    disp(['performing coreg for ' RIDs{subI}])
                    matlabbatch{mbatch}.spm.spatial.coreg.estimate.ref ={average_anat};
                    matlabbatch{mbatch}.spm.spatial.coreg.estimate.source = {to_coreg};
                    matlabbatch{mbatch}.spm.spatial.coreg.estimate.other = {''};
                    matlabbatch{mbatch}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
                    matlabbatch{mbatch}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
                    matlabbatch{mbatch}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
                    matlabbatch{mbatch}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
                    spm_jobman('run',matlabbatch(mbatch))
                    % prefix the coreged data with an r so that the code
                    % can grab the correct files
                    to_rename = get_files(current_dir, 'ac_*.nii');
                    [path, name, ext] = fileparts(to_rename);
                    movefile(to_rename, [path, filesep, 'r', name, ext]);
                else
                    disp('make sure you ran checkVoxels to determine whether coreg est only or coreg est res is needed')
                end
            end
        end
    end
    for t=1:length(mri_to_run)
        if pars.Results.longitudinalRegistration
            % create a symmetric average of the anatomical images. one
            % average is created per imaging type (e.g. PD)
            anat = sub_table(contains(sub_table.ScanType, mri_to_run{t}), :) ;
            anat = anat(anat.FieldStrength == field_str, :);
            %% longitudinal registration can be run in indMode if the user would like. this will improve cross timepoint alignment within subjects
            % double check that if indMode is on, the times to run are more than 1
            if ~strcmpi(indMode{1}, 'all_times') && isnumeric(indMode{1}) && strcmp(avgMode, 'off')
                if length(indMode{1}) > 1
                    % we're only processing some of the times
                    anat = anat(ismember(anat.Order, indMode{1}), :);
                else
                    disp('you are attempting to run longitudinal registration on only one timepoint of data. This cannot be done')
                    disp('please update your inputs to indMode and avgMode');
                    break
                end
            end
            if ~isempty(anat)
                bl_dir = fullfile(sub_dir, char(anat.Date(1)), char(anat.Modality(1)), char(anat.ScanType(1)),char(anat.ScanName(1)));
            else
                disp(['no ' mri_to_run{t} ' data for ' RIDs{subI}])
                continue
            end
            % should be sorted already by first date but make sure
            anat = sortrows(anat, 'Date');
            if isequal(height(anat),1)
                disp(['only one ' mri_to_run{t} ' scan for participant ' RIDs{subI} '. No average will be created'])
                continue
            elseif height(anat) > 1 && isequal(sum(anat.IsDuplicate), 0)
                % there are no intra-interval averages to be made.
                % subtract each date from the baseline date and convert
                % difference to years
                differences = abs(anat.Date(1) - anat.Date);
                interscan_int = years(differences);
                for j=1:height(anat)
                    current_direc = fullfile(sub_dir, char(anat.Date(j)), char(anat.Modality(j)), char(anat.ScanType(j)),char(anat.ScanName(j)));
                    anat_images{j,1} = get_files(current_direc, 'rac_*.nii');
                end
                mbatch = mbatch +1;
                matlabbatch{mbatch}.spm.tools.longit.series.vols = anat_images;
                matlabbatch{mbatch}.spm.tools.longit.series.times = interscan_int';
                matlabbatch{mbatch}.spm.tools.longit.series.noise = NaN;
                matlabbatch{mbatch}.spm.tools.longit.series.wparam = [0 0 100 25 100];
                matlabbatch{mbatch}.spm.tools.longit.series.bparam = 1000000;
                matlabbatch{mbatch}.spm.tools.longit.series.write_avg = 1;
                matlabbatch{mbatch}.spm.tools.longit.series.write_jac = 1;
                matlabbatch{mbatch}.spm.tools.longit.series.write_div = 0;
                matlabbatch{mbatch}.spm.tools.longit.series.write_def = 0;
                spm_jobman('run',matlabbatch(mbatch))
            else
                % this means multiple scans were taken on one day. they may
                % have only one timepoint or may have multiple timepoints.
                % Create any intraaverages first then do midpoints
                if isequal(length(unique(anat.Date)),1) 
                    if strcmpi(avgMode, 'on') && strcmpi(indMode{1}, 'off')
                        continue % skip them because avgMode is for longitudinal data
                    elseif strcmpi(indMode{2}, 'excludeCrossSectionalSubjects')                   
                        continue 
                    end
                end
                % grab the duplicate timepoints for averaging
                duplicates = anat(anat.IsDuplicate  ==1, :);
                % each duplicate is tethered to the order variable. We know
                % all scans which pertain to a timepoint using this
                % variable
                orders_to_run = unique(duplicates.Order);
                for dup=1:length(orders_to_run)
                    current_timept = duplicates(duplicates.Order == orders_to_run(dup), :);
                    % define the interval between scans (0)
                    intra_int = zeros(height(current_timept),1);
                    % grab all of the images of that modality on that
                    % date
                    for c=1:height(current_timept)
                        current_timept_dir = fullfile(sub_dir, char(current_timept.Date(c)), char(current_timept.Modality(c)), char(current_timept.ScanType(c)),char(current_timept.ScanName(c)));
                        current_timept_images{c} = get_files(current_timept_dir, 'rac_*.nii');
                    end
                    % create the intrascan average for that timepoint
                    mbatch = mbatch +1;
                    matlabbatch{mbatch}.spm.tools.longit.series.vols = current_timept_images';
                    matlabbatch{mbatch}.spm.tools.longit.series.times = intra_int';
                    matlabbatch{mbatch}.spm.tools.longit.series.noise = NaN;
                    matlabbatch{mbatch}.spm.tools.longit.series.wparam = [0 0 100 25 100];
                    matlabbatch{mbatch}.spm.tools.longit.series.bparam = 1000000;
                    matlabbatch{mbatch}.spm.tools.longit.series.write_avg = 1;
                    matlabbatch{mbatch}.spm.tools.longit.series.write_jac = 0;
                    matlabbatch{mbatch}.spm.tools.longit.series.write_div = 0;
                    matlabbatch{mbatch}.spm.tools.longit.series.write_def = 0;
                    spm_jobman('run',matlabbatch(mbatch))
                    % rename the average to intra_avg* so it's
                    % clear what this is. It will be written to the
                    % first directory
                    intra_avg = get_files(fullfile(sub_dir, char(current_timept.Date(1)), char(current_timept.Modality(1)), char(current_timept.ScanType(1)),char(current_timept.ScanName(1))), 'avg_rac*.nii');
                    [avg_path, avg_file, ext] = fileparts(intra_avg);
                    better_name = ['intra_' avg_file];
                    movefile([avg_path filesep avg_file ext], [avg_path filesep better_name ext]);
                    % also zip the file to keep a backup
                    zip([avg_path, filesep, better_name, '.zip'], [avg_path, filesep, better_name, ext]);
                    clear current_timept_images
                end
                % now we have calculated all of the intra-timepoint
                % averages. Create the longitudinal averages
                % incorporating the intra averages
                % keep track of which duplicate dates were
                % processed
                if isequal(length(unique(anat.Date)),1)
                    % don't need to make any long averages because they
                    % have cross sectional data
                    continue
                end
                to_avg ={};
                % drop rows from anat that aren't the timepoint references.
                % All intra avg images are in the timepoint reference
                % folder
                anat = anat(anat.TimepointReference ==1,:);
                for b=1:height(anat)
                    anat_dir = fullfile(sub_dir, char(anat.Date(b)), char(anat.Modality(b)), char(anat.ScanType(b)),char(anat.ScanName(b)));
                    if isequal(anat.IsDuplicate(b), 1)
                        % they have an intra average 
                        % unzip to use the original
                        unzip(get_files(anat_dir, 'intra*rac*.zip'), anat_dir);
                        to_avg{end+1} = get_files(anat_dir, 'intra*rac*.nii');
                    else
                        to_avg{end+1} = get_files(anat_dir, 'rac*.nii');
                    end
                end
                % calculate the interval using only non duplciate
                % dates
                unique_dates = unique(anat.Date, 'stable'); % all dates should be unique now but just in case
                differences = abs(unique_dates(1) - unique_dates);
                interscan_int = years(differences);
                mbatch = mbatch + 1;
                disp(['performing longitudinal registration for ' RIDs{subI}])
                matlabbatch{mbatch}.spm.tools.longit.series.vols = to_avg';
                matlabbatch{mbatch}.spm.tools.longit.series.times = interscan_int';
                matlabbatch{mbatch}.spm.tools.longit.series.noise = NaN;
                matlabbatch{mbatch}.spm.tools.longit.series.wparam = [0 0 100 25 100];
                matlabbatch{mbatch}.spm.tools.longit.series.bparam = 1000000;
                matlabbatch{mbatch}.spm.tools.longit.series.write_avg = 1;
                matlabbatch{mbatch}.spm.tools.longit.series.write_jac = 1;
                matlabbatch{mbatch}.spm.tools.longit.series.write_div = 0;
                matlabbatch{mbatch}.spm.tools.longit.series.write_def = 0;
                spm_jobman('run',matlabbatch(mbatch))
            end
            [path, name, ex] = fileparts(get_files(bl_dir, 'avg*.nii'));
            % move it to it's own directory within the bl_dir
            if ~exist(fullfile(bl_dir, 'midpoint_average'), 'dir')
                mkdir(fullfile(bl_dir, 'midpoint_average'));
            end
            movefile(get_files(bl_dir, 'avg_*.nii'), [bl_dir, filesep, 'midpoint_average']);
            % zip the longitudinal average file so we'll have a nice backup if needed
            zip([path, filesep, 'midpoint_average', filesep, name, '.zip'], [path, filesep, 'midpoint_average',filesep, name, ex]);
            clear anat_images current_timept_images to_avg
        end
        if pars.Results.segment
            % we can start using the spread filter table function now since
            % we no longer need to work with the duplicate scans
            anat_table = spread_table_filter(sub_table, mri_to_run{t}, field_str);
            to_seg = {};
            if isempty(anat_table)
                disp(['no ' mri_to_run{t} ' data for ' RIDs{subI}])
                    continue
            end
            if strcmpi(indMode{1}, 'off') || strcmpi(indMode{2}, 'excludeCrossSectionalSubjects') 
                if isequal(max(anat_table.Order),1) || isequal(height(anat_table),1)
                    % there is only one timepoint for this sub. Skip.
                    % had to add the second check because some people in
                    % AIBL have mostly 1.5T data so when dropping that data
                    % ends up being cross sectional
                    continue
                end
            end
            % if we are in avgMode, we need to segment the midpoint avg image
            if strcmpi(avgMode, 'on')             
                bl_dir = fullfile(sub_dir, char(anat_table.Date(1)), char(anat_table.Modality(1)), char(anat_table.ScanType(1)),char(anat_table.ScanName(1)));
                avg_dir = fullfile(bl_dir, 'midpoint_average');            
                to_seg{end+1} = get_files(avg_dir, ['avg*.nii']);
                if isempty(to_seg)
                    disp(['detected multiple timepoints for ', RIDs{subI}, ' but there is no average image. Make sure you ran the longitudinalRegistration module'])
                    continue
                end
            end
            % if indmode is also on, grab any additional timepoint scans to
            % seg
            if ~strcmpi(indMode{1}, 'off')
                if ~strcmpi(indMode{1}, 'all_times') && isnumeric(indMode{1})
                    % we're only processing some of the times
                    anat_table = anat_table(ismember(anat_table.Order, indMode{1}), :);
                end
                for row=1:height(anat_table)
                    if isequal(anat_table.IsDuplicate(row),1)
                        % grab their intra
                        to_seg{end+1}= get_files(fullfile(sub_dir, char(anat_table.Date(row)), char(anat_table.Modality(row)), char(anat_table.ScanType(row)),char(anat_table.ScanName(row))), ['intra*rac*.nii']);
                    else
                        % just a regular scan
                          to_seg{end+1}= get_files(fullfile(sub_dir, char(anat_table.Date(row)), char(anat_table.Modality(row)), char(anat_table.ScanType(row)),char(anat_table.ScanName(row))), ['rac*.nii']);
                    end
                end
            end
            CAT12_shoot_template = get_files(CAT12_shoot_path, 'Template_0_GS.nii');
            lorio_TPM = get_files(templatePath, 'enhanced_TPM.nii');
            mbatch=mbatch+1;
            %% run SPREAD segmentation routine. Make sure CAT12 default is to run in expert mode (see function below)
            % (cat.extopts.expertgui    = 1;) in cat_default.m
            for s=1:length(to_seg)
                SPREAD_segmentation_routine(mbatch, to_seg{s}, lorio_TPM, CAT12_shoot_template, 0.5, surf_mod); % use processing accuracy of 0.5 to begin (default)
            end
        end
        clear to_seg
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
