%% PET preproc step 2
% code by Hayley R. C. Shanks, Taylor W. Schmitz, Kate M. Onuska
function longitudinal_PET_preproc_step2(N, study_path, pet_to_run, pet_orders, mri_ref, study_name, frames, roi_dir, ref_labels, ref_mask, avgPrefix, qc_thresh,field_str, pvc_opt, varargin)

%%%%%%%%%%%%%%%%%%% lab-specfic INPUTS - may need to be changed %%%%%%%%%%%%%%%%%%%%%%%%%
% add path to spm
addpath '/srv/schmitz/projects/toolboxes/spm12';
functionRoot = '/srv/schmitz/projects/toolboxes/incalab_source_code/human';
addpath(genpath('/srv/schmitz/projects/toolboxes/spm12/toolbox/petpve12')); 
% catch sight of the functions
addpath(genpath(functionRoot));
%%%%%%%%%%%%%%%%%%%%%%%%%%%% End of specific inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
study_avg_dir = [study_path, filesep, study_name, '_averages'];
% raw participant data should live in a first level subdirectory of the study path
rootPath = [study_path, filesep, 'first_level'];
sheet_path = [study_path, filesep, 'spreadsheets'];
annotationTableRoot = '/srv/schmitz/projects/toolboxes/';
shootTemplatePath = [study_path, filesep, study_name, '_shoot_template'];
referenceMaskPath = [study_path, filesep, 'ROIs', filesep, 'wholeBrain'];
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
addRequired(pars, 'frames', @iscell); % static or dynamic
addRequired(pars, 'roi_dir', @ischar);
addRequired(pars, 'ref_labels', @iscell); % static or dynamic
addRequired(pars, 'ref_mask', @iscell); % static or dynamic
addRequired(pars, 'avgPrefix', @ischar);
addRequired(pars, 'qc_thresh', @isnumeric);
addRequired(pars, 'field_str', @isnumeric);
addRequired(pars, 'pvc_opt', @ischar);
%add optional inputs.
addParameter(pars,'warpToShoot',default_run,@isnumeric);  
addParameter(pars,'makeMasks',default_run,@isnumeric);  
addParameter(pars,'computeSUVR_refMask',default_run,@isnumeric); % this should be done if a validated apriori ref mask exists
addParameter(pars,'computeSUVR_noMask',default_run,@isnumeric); 
addParameter(pars,'pvc',default_run,@isnumeric);  %optional, should only be used with certain radiotracers
addParameter(pars,'createAvg',default_run,@isnumeric); 
addParameter(pars,'shootTemplateOverride','',@ischar); 
% parse the function inputs
parse(pars, N, study_path, pet_to_run, pet_orders, mri_ref, study_name, frames, roi_dir, ref_labels, ref_mask, avgPrefix, qc_thresh,field_str,pvc_opt,varargin{:});
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
% when the summarized QC variable is 1, they have data passing QC. If it's
% 0 it fails, and NaN if not processed
if strcmpi(mri_ref{1}, 'base')
    good_mri = mri_mastersheet.RID(mri_mastersheet.timepoint_fails == 0);
elseif strcmpi(mri_ref{1}, 'avg')
    good_mri = mri_mastersheet.RID(mri_mastersheet.summarizedAvgFails ==0);
else
    disp(['please review mri_ref description and re-define the variable'])
    return
end
if pars.Results.warpToShoot
    for subI=sub_to_run(N)
        sub_dir = fullfile(rootPath,RIDs{subI});
        if ~contains(RIDs{subI}, good_mri) 
            %skip. we need mri data for pet preproc
        else
            sub_table = get_files(sub_dir, 'scan_info*.csv');
            % this table has all of the dates and types of scans for this P
            sub_table = readtable(sub_table);
            %make sure the rows are sorted by earliest date
            if isempty(sub_table)
                continue
            end
            mri_table = spread_table_filter(sub_table, mri_ref{2}, field_str); % mri_ref tells us if we want T1, t2, etc.        
            mri_bl_dir = fullfile(sub_dir, char(mri_table.Date(1)), char(mri_table.Modality(1)), char(mri_table.ScanType(1)), char(mri_table.ScanName(1)));
            if strcmpi(mri_ref{1}, 'avg')
                mri_ref_dir = [mri_bl_dir, filesep 'midpoint_average' filesep 'mri'];
            elseif strcmpi(mri_ref{1}, 'base')
                mri_ref_dir = [mri_bl_dir, filesep, 'mri'];
            end  
            y_file = get_files(fullfile(mri_ref_dir),'y_rp1*rigid*Template.nii');
            % now get their pet data 
            for a=1:length(pet_to_run)
                pet_table = spread_pet_filter(sub_table, pet_to_run{a});
                if isempty(pet_table)
                    disp(['no applicable ' pet_to_run{a} ' data for ' char(RIDs{subI})])
                    continue
                else
                    to_warp ={};
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
                        to_warp = [to_warp; cellstr(get_files(scan_dir, 'rac*.nii'))];
                        to_warp = [to_warp; cellstr(get_files(scan_dir, 'rmean*.nii'))];
                    end
                    mbatch = mbatch +1;
                    if ~isempty(pars.Results.shootTemplateOverride)
                        shootTemplate=get_files(pars.Results.shootTemplateOverride,['*Template*.nii']);
                    else
                        shootTemplate=get_files(fullfile(shootTemplatePath),[study_name '*' 'Template*.nii']);
                    end
                    [temp_bb,temp_vx] = spm_get_bbox(shootTemplate(5,:));
                    matlabbatch{mbatch}.spm.tools.cat.tools.defs.field1 = cellstr({y_file});
                    matlabbatch{mbatch}.spm.tools.cat.tools.defs.images = to_warp;
                    matlabbatch{mbatch}.spm.tools.cat.tools.defs.bb = [NaN NaN NaN; NaN NaN NaN]; %% leave as NaN, it will default to the target image
                    matlabbatch{mbatch}.spm.tools.cat.tools.defs.vox = [NaN NaN NaN];
                    matlabbatch{mbatch}.spm.tools.cat.tools.defs.interp = 1; %0 = NN; 1 = trilinear; 4 = 4th degree b-splines
                    % Preserve Concentrations: Smoothed spatially normalised images (sw*) represent weighted averages
                    % of the signal under the smoothing kernel, approximately preserving the intensities of the
                    % original images. This option is currently suggested for eg fMRI.
                    % Preserve Total: Smoothed and spatially normalised images preserve the total amount of signal
                    % from each region in the images (smw*). Areas that are expanded during warping are correspondingly
                    % reduced in intensity. This option is suggested for VBM.
                    matlabbatch{mbatch}.spm.tools.cat.tools.defs.modulate = 0; % hard coded to nomod because PET should never be modulated
                    spm_jobman('run',matlabbatch(mbatch))
                end
            end
        end
    end
end   

if pars.Results.makeMasks
    % we need a grey matter mask from the template, binarized at p >= 0.5,
    % and a whole brain mask for SUVR (to mask out non-brain regions)
    %% first make the GM mask
     if ~isempty(pars.Results.shootTemplateOverride)
        shootTemplate=get_files(pars.Results.shootTemplateOverride,['*Template*.nii']);
        temp_path = pars.Results.shootTemplateOverride;
    else
        shootTemplate=get_files(fullfile(shootTemplatePath),[study_name '*' 'Template*.nii']);
        temp_path = shootTemplatePath;
    end
    template = shootTemplate(5,:);
    mbatch = mbatch +1;
    matlabbatch{mbatch}.spm.util.imcalc.input = {[template, ',1']};% take the first frame of the 4d file
    matlabbatch{mbatch}.spm.util.imcalc.output = 'bin_p_0_05_template4_GM';
    matlabbatch{mbatch}.spm.util.imcalc.outdir = {temp_path};
    matlabbatch{mbatch}.spm.util.imcalc.expression = 'i1 >= 0.5';
    matlabbatch{mbatch}.spm.util.imcalc.var = struct('name', {}, 'value', {});
    matlabbatch{mbatch}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{mbatch}.spm.util.imcalc.options.mask = 0;
    matlabbatch{mbatch}.spm.util.imcalc.options.interp = 1;
    matlabbatch{mbatch}.spm.util.imcalc.options.dtype = 4;
    spm_jobman('run',matlabbatch(mbatch))
    %% now make the whole brain mask - first check if there is a mask from the study average
    mask_check = get_files([study_avg_dir, filesep, 'mri'], 'p0*.nii');
    if ~isempty(mask_check) && ~strcmpi(mask_check, '')
        % we don't need to make anything
    else
        mbatch = mbatch +1;
        images = [[template, ',1']; [template, ',2']; [template, ',3']];
        matlabbatch{mbatch}.spm.util.imcalc.input = cellstr(images);
        matlabbatch{mbatch}.spm.util.imcalc.output = 'bin_gm_wm_csf_mask_Template4';
        matlabbatch{mbatch}.spm.util.imcalc.outdir = {temp_path};
        matlabbatch{mbatch}.spm.util.imcalc.expression = '((i1+i2+(i3*-1)))>=0.5';
        matlabbatch{mbatch}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{mbatch}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{mbatch}.spm.util.imcalc.options.mask = 0;
        matlabbatch{mbatch}.spm.util.imcalc.options.interp = 1;
        matlabbatch{mbatch}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run',matlabbatch(mbatch))
        %fill holes in mask
        wb = spm_vol([temp_path, filesep, 'bin_gm_wm_csf_mask_Template4.nii']);
        wb_mat = spm_read_vols(wb);
        wb_filled = imfill(wb_mat, 'holes');
        spm_write_vol(wb, wb_filled);
    end
end

if pars.Results.computeSUVR_refMask
    if ~isempty(pars.Results.shootTemplateOverride)
        wb_mask=get_files(pars.Results.shootTemplateOverride, 'bin_gm_wm_csf_mask_Template4.nii');
        gm_ref_img = get_files(pars.Results.shootTemplateOverride, 'bin_p_0_05_template4_GM.nii');
    else
        wb_mask=get_files(fullfile(shootTemplatePath),'bin_gm_wm_csf_mask_Template4.nii');
        gm_ref_img = get_files(shootTemplatePath, 'bin_p_0_05_template4_GM.nii');
    end
    mask_check = get_files([study_avg_dir, filesep, 'mri'], 'p0*.nii');
    if ~isempty(mask_check) && ~strcmpi(mask_check, '')
        % use the p0 file from cat12 instead
        wb_mask = mask_check;      
    else
        %check if the shootTemplate Override was specified, and use p0 from
        %there
        if ~isempty(pars.Results.shootTemplateOverride)
            other_avg_path = replace(pars.Results.shootTemplateOverride, '_shoot_template', '_averages');
            wb_mask = get_files([other_avg_path, filesep, 'mri'], 'p0*.nii');
        end
    end
    % read in the gm image and the p0 image
    gm_hdr = spm_vol(gm_ref_img);
    p50_gm = spm_read_vols(gm_hdr);
    p0_hdr = spm_vol(wb_mask);
    p0_mat = spm_read_vols(p0_hdr);
    p50_gm(p50_gm ==0) =NaN; % nan non-GM voxels
    % create columns in mastersheet to store info on SUVR. Check if they
    % exist first
    names_to_check = {'uptake_mean_ref', 'uptake_mean_gm', 'mean_gm_suvr','uptake_stdev_reference','uptake_stdev_gm','stdev_gm_suvr'}; 
    if isequal(sum(contains(mastersheet.Properties.VariableNames, names_to_check)), 0)
        mastersheet.uptake_mean_reference = nan(height(mastersheet),1);
        mastersheet.uptake_mean_gm = nan(height(mastersheet),1);
        mastersheet.mean_gm_suvr = nan(height(mastersheet),1);
        mastersheet.uptake_stdev_reference = nan(height(mastersheet),1);
        mastersheet.uptake_stdev_gm = nan(height(mastersheet),1);
        mastersheet.stdev_gm_suvr = nan(height(mastersheet),1);
    end
    %% get the reference region mask. should be in some pre-defined ROI directory
    % get annotation file for the atlas containing various reference
    % regions
    roi_file = get_files(roi_dir, char(ref_mask)); 
    roi_hdr = spm_vol(roi_file);
    roi_mat = spm_read_vols(roi_hdr);
    roi_mat(roi_mat == 0) =NaN;
    for subI=sub_to_run(N)
        sub_dir = fullfile(rootPath,RIDs{subI});
        if ~contains(RIDs{subI}, good_mri)
            %skip. we need mri data for pet preproc
        else
            sub_table = get_files(sub_dir, 'scan_info*.csv');
            % this table has all of the dates and types of scans for this P
            sub_table = readtable(sub_table);
            %make sure the rows are sorted by earliest date
            if isempty(sub_table)
                continue
            end
            mri_table = spread_table_filter(sub_table, mri_ref{2}, field_str); % mri_ref tells us if we want T1, t2, etc.        
            mri_bl_dir = fullfile(sub_dir, char(mri_table.Date(1)), char(mri_table.Modality(1)), char(mri_table.ScanType(1)), char(mri_table.ScanName(1)));
            if strcmpi(mri_ref{1}, 'avg')
                mri_bl_dir = [mri_bl_dir, filesep, 'midpoint_average'];
                mri_ref_dir = [mri_bl_dir filesep 'mri'];
                mr_ref_img = get_files(mri_ref_dir, 'mavg*.nii');
            elseif strcmpi(mri_ref{1}, 'base')
                mri_ref_dir = [mri_bl_dir, filesep, 'mri'];
                mr_ref_img = get_files(mri_ref_dir, 'mintra*.nii');
                if isempty(mr_ref_img)
                    mr_ref_img = get_files(mri_ref_dir, 'mrac*.nii');
                end
            end
            for a=1:length(pet_to_run)
                pet_table = spread_pet_filter(sub_table, pet_to_run{a});
                if isempty(pet_table)
                    disp(['no applicable ' pet_to_run{a} ' data for ' char(RIDs{subI})])
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
                        % find the index of this scan in
                        % mastersheet
                        check = [strcmp(mastersheet.RID, RIDs{subI}), mastersheet.Date == pet_table.Date(i), strcmp(mastersheet.Modality,pet_table.Modality(i)), strcmp(mastersheet.ScanType, pet_table.ScanType(i)), strcmp(mastersheet.ScanName, pet_table.ScanName(i))] ;
                        % which ever row of check is all 1's is our row
                        check_sum = sum(check, 2);
                        ind = find(check_sum == max(check_sum));
                        % we'll loop through each subject and normalize
                        % their data to a reference uptake value. 
                        if ~strcmpi(frames, 'dynamic') 
                            pet_files = cellstr(get_files(scan_dir, 'wrac*.nii'));
                            ref_pet = char(pet_files);
                        else
                            % for dynamic data, we will have multiple
                            % frames per scan folder. 
                            pet_files = cellstr(get_files(scan_dir, 'wrac*.nii'));
                            % their ref file will be the mean image for
                            % all their frames from that timepoint. We
                            % need to confirm that re-align 
                            ref_pet = get_files(scan_dir, 'wrmean*.nii');
                            % add the ref pet image so it will be
                            % suvr'd too
                            pet_files = [pet_files; cellstr(ref_pet)];
                        end
                        ref_hdr = spm_vol(ref_pet);
                        ref_mat = spm_read_vols(ref_hdr);
                        % first we need to create a brain extracted
                        % version of the ref image. In the p0 file, non-brain is 0
                        ref_mat(p0_mat ==0) = NaN;
                        % multiply brain extracted PET by binary
                        % ref region and get average uptake. 
                        % save the mean value to mastersheet
                        if size(ref_mat) ~= size(roi_mat)
                            disp(['fails: ', num2str(subI)])
                            bb = spm_get_bbox(roi_ref_img);
                        else
%                                 disp(['passes: ', num2str(subI)])
%                                 bb = spm_get_bbox(roi_ref_img)
                            mastersheet.uptake_mean_reference(ind) = mean(ref_mat.*roi_mat,'all', 'omitnan');
                            mastersheet.uptake_stdev_reference(ind) = std(ref_mat.*roi_mat,0, 'all', 'omitnan');                           
                            % find the mean and stdev GM uptake 
                            mastersheet.uptake_mean_gm(ind) = mean(ref_mat.*p50_gm,'all', 'omitnan');
                            mastersheet.uptake_stdev_gm(ind) = std(ref_mat.*p50_gm, 0, 'all','omitnan');
                            for r=1:length(pet_files) % r will be 1 max for static pet
                                % now we work with the individual files to
                                % compute an suvr
                                pet_hdr = spm_vol(pet_files{r});
                                pet_mat = spm_read_vols(pet_hdr);
                                pet_mat(p0_mat ==0) = NaN;
                                % compute the "SUVR"
                                pet_norm = pet_mat./mastersheet.uptake_mean_reference(ind);                                                              
                                % save the finalized file. remove pinfo and
                                % update fname
                                [pet_path, old_name, ex] = fileparts(pet_hdr.fname);
                                pet_hdr.fname = [pet_path, filesep, 'suvr_', old_name, ex] ;
                                pet_hdr = rmfield(pet_hdr, 'pinfo');
                                spm_write_vol(pet_hdr, pet_norm);
                                mastersheet.mean_gm_suvr(ind) = mean(pet_norm.*p50_gm, 'all', 'omitnan');
                                mastersheet.stdev_gm_suvr(ind) = std(pet_norm.*p50_gm,0, 'all', 'omitnan');  
                            end
                        end
                    end
                end
            end
        end
    end
    % write the mastersheet
    writetable(mastersheet, [sheet_path, filesep, study_name, '_PET_mastertable.csv']);
end

if pars.Results.computeSUVR_noMask
    if ~isempty(pars.Results.shootTemplateOverride)
        wb_mask=get_files(pars.Results.shootTemplateOverride, 'bin_gm_wm_csf_mask_Template4.nii');
        gm_ref_img = get_files(pars.Results.shootTemplateOverride, 'bin_p_0_05_template4_GM.nii');
    else
        wb_mask=get_files(fullfile(shootTemplatePath),'bin_gm_wm_csf_mask_Template4.nii');
        gm_ref_img = get_files(shootTemplatePath, 'bin_p_0_05_template4_GM.nii');
    end
    mask_check = get_files([study_avg_dir, filesep, 'mri'], 'p0*.nii');
    if ~isempty(mask_check) && ~strcmpi(mask_check, '')
        % use the p0 file from cat12 instead
        wb_mask = mask_check;
    end
    % read in the gm image and the p0 image
    gm_hdr = spm_vol(gm_ref_img);
    p50_gm = spm_read_vols(gm_hdr);
    p0_hdr = spm_vol(wb_mask);
    p0_mat = spm_read_vols(p0_hdr);
    p50_gm(p50_gm ==0) =NaN; % nan non-GM voxels
    % create columns in mastersheet to store info on SUVR. Check if they
    % exist first
    names_to_check = {'uptake_mean_gm', 'mean_gm_suvr','uptake_stdev_gm','stdev_gm_suvr'};
    if isequal(sum(contains(mastersheet.Properties.VariableNames, names_to_check)), 0)
        mastersheet.uptake_mean_gm = nan(height(mastersheet),1);
        mastersheet.mean_gm_suvr = nan(height(mastersheet),1);
        mastersheet.uptake_stdev_gm = nan(height(mastersheet),1);
        mastersheet.stdev_gm_suvr = nan(height(mastersheet),1);
    end
    for subI=sub_to_run(N)
        sub_dir = fullfile(rootPath,RIDs{subI});
        if ~contains(RIDs{subI}, good_mri)
            %skip. we need mri data for pet preproc
        else
            sub_table = get_files(sub_dir, 'scan_info*.csv');
            % this table has all of the dates and types of scans for this P
            sub_table = readtable(sub_table);
            %make sure the rows are sorted by earliest date
            if isempty(sub_table)
                continue
            end
            mri_table = spread_table_filter(sub_table, mri_ref{2}, field_str); % mri_ref tells us if we want T1, t2, etc.
            mri_bl_dir = fullfile(sub_dir, char(mri_table.Date(1)), char(mri_table.Modality(1)), char(mri_table.ScanType(1)), char(mri_table.ScanName(1)));
            if strcmpi(mri_ref{1}, 'avg')
                mri_bl_dir = [mri_bl_dir, filesep, 'midpoint_average'];
                mri_ref_dir = [mri_bl_dir filesep 'mri'];
                mr_ref_img = get_files(mri_ref_dir, 'mavg*.nii');
            elseif strcmpi(mri_ref{1}, 'base')
                mri_ref_dir = [mri_bl_dir, filesep, 'mri'];
                mr_ref_img = get_files(mri_ref_dir, 'mintra*.nii');
                if isempty(mr_ref_img)
                    mr_ref_img = get_files(mri_ref_dir, 'mrac*.nii');
                end
            end
            
            for a=1:length(pet_to_run)
                pet_table = spread_pet_filter(sub_table, pet_to_run{a});
                if isempty(pet_table)
                    disp(['no applicable ' pet_to_run{a} ' data for ' char(RIDs{subI})])
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
                        % find the index of this scan in
                        % mastersheet
                        check = [strcmp(mastersheet.RID, RIDs{subI}), mastersheet.Date == pet_table.Date(i), strcmp(mastersheet.Modality,pet_table.Modality(i)), strcmp(mastersheet.ScanType, pet_table.ScanType(i)), strcmp(mastersheet.ScanName, pet_table.ScanName(i))] ;
                        % which ever row of check is all 1's is our row
                        check_sum = sum(check, 2);
                        ind = find(check_sum == max(check_sum));
                        % we'll loop through each subject and normalize
                        % their data to a reference uptake value.
                        if ~strcmpi(frames, 'dynamic')
                            pet_files = cellstr(get_files(scan_dir, 'wrac*.nii'));
                            ref_pet = char(pet_files);
                        else
                            ref_pet = get_files(scan_dir, 'wrmean*.nii');
                            % add the ref pet image so it will be
                            % suvr'd too
                            pet_files = [cellstr(ref_pet)];
                        end
                        
                        ref_hdr = spm_vol(ref_pet);
                        ref_mat = spm_read_vols(ref_hdr);
                        % first we need to create a brain extracted
                        % version of the ref image. In the p0 file, non-brain is 0
                        ref_mat(p0_mat ==0) = NaN;
                        % find the mean and stdev GM uptake
                        mastersheet.uptake_mean_gm(ind) = mean(ref_mat.*p50_gm,'all', 'omitnan');
                        mastersheet.uptake_stdev_gm(ind) = std(ref_mat.*p50_gm, 0, 'all','omitnan');
                        for r=1:length(pet_files) % r will be 1 max for static pet
                            % now we work with the individual files to
                            % compute an suvr
                            pet_hdr = spm_vol(pet_files{r});
                            pet_mat = spm_read_vols(pet_hdr);
                            pet_mat(p0_mat ==0) = NaN;
                            % compute the "SUVR"
                            pet_norm = pet_mat./mastersheet.([pet_to_run{a}, '_final_reference_value'])(ind);
                            % save the finalized file. remove pinfo and
                            % update fname
                            [pet_path, old_name, ex] = fileparts(pet_hdr.fname);
                            pet_hdr.fname = [pet_path, filesep, 'suvr_', old_name, ex] ;
                            pet_hdr = rmfield(pet_hdr, 'pinfo');
                            spm_write_vol(pet_hdr, pet_norm);
                            mastersheet.mean_gm_suvr(ind) = mean(pet_norm.*p50_gm, 'all', 'omitnan');
                            mastersheet.stdev_gm_suvr(ind) = std(pet_norm.*p50_gm,0, 'all', 'omitnan');
                        end
                    end
                end
            end
        end
    end
    % write the mastersheet
    writetable(mastersheet, [sheet_path, filesep, study_name, '_PET_mastertable.csv']);
end

if pars.Results.pvc %optional, should only be used with certain radiotracers
    %% partial volume effect correction - may be tracer type specific - so you may wanna use it for AB tracers but not for FDG as an example
    for subI=sub_to_run(N)
        sub_dir = fullfile(rootPath,RIDs{subI});
        if ~contains(RIDs{subI}, good_mri) 
            %skip. we need mri data for pet preproc
        else
            sub_table = get_files(sub_dir, 'scan_info*.csv');
            % this table has all of the dates and types of scans for this P
            sub_table = readtable(sub_table);
            %make sure the rows are sorted by earliest date
            if isempty(sub_table)
                continue
            end
            mri_table = spread_table_filter(sub_table, mri_ref{2}, field_str); % mri_ref tells us if we want T1, t2, etc.        
            mri_bl_dir = fullfile(sub_dir, char(mri_table.Date(1)), char(mri_table.Modality(1)), char(mri_table.ScanType(1)), char(mri_table.ScanName(1)));
            if strcmpi(mri_ref{1}, 'avg')
                mri_ref_dir = [mri_bl_dir, filesep 'midpoint_average' filesep 'mri'];
            elseif strcmpi(mri_ref{1}, 'base')
                mri_ref_dir = [mri_bl_dir, filesep, 'mri'];
            end  
            % get gm, wm and CSF segments
            gm = get_files(mri_ref_dir, 'wp1*.nii');
            wm = get_files(mri_ref_dir, 'wp2*.nii');
            csf = get_files(mri_ref_dir, 'wp3*.nii');
            % now get their pet data 
            for a=1:length(pet_to_run)
                pet_table = spread_pet_filter(sub_table, pet_to_run{a});
                if isempty(pet_table)
                    disp(['no applicable ' pet_to_run{a} ' data for ' char(RIDs{subI})])
                    continue
                else
                    to_pvc ={};
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
                        % get our files. pvc_opt tells you what the file
                        % prefix should be - e.g. suvrw or w (dpeending on
                        % if you suvr or no)
                        to_pvc = [to_pvc; cellstr(get_files(scan_dir, [pvc_opt, '*.nii']))];
                    end
                    %pvc each pet scan
                    for p=1:length(to_pvc)
                        mbatch = mbatch +1;
                        matlabbatch{mbatch}.spm.tools.petpve.spatial.PVEcorrection.PETdata = to_pvc(p);
                        % this all stays the same as for the individual scans
                        matlabbatch{mbatch}.spm.tools.petpve.spatial.PVEcorrection.SegImgs.Tsegs.tiss1 = {gm};
                        matlabbatch{mbatch}.spm.tools.petpve.spatial.PVEcorrection.SegImgs.Tsegs.tiss2 = {wm};
                        matlabbatch{mbatch}.spm.tools.petpve.spatial.PVEcorrection.PVEopts.fwhm_PSF = [6 6 6];
                        % using a slightly more liberal GM threshold (default is 0.5)
                        matlabbatch{mbatch}.spm.tools.petpve.spatial.PVEcorrection.PVEopts.gmthresh = 0.3; 
                        % estimate CSF signal
                        matlabbatch{mbatch}.spm.tools.petpve.spatial.PVEcorrection.PVEopts.CSFsignal.CSFcalc.tiss3 = {csf};
                        matlabbatch{mbatch}.spm.tools.petpve.spatial.PVEcorrection.PVEopts.TissConv = 0;
                        matlabbatch{mbatch}.spm.tools.petpve.spatial.PVEcorrection.PVE_Const_opts.type3.wmcsfthresh = 0.9;
                        spm_jobman('run',matlabbatch(mbatch)) 
                    end
                end
            end
        end
    end
end


if pars.Results.createAvg
    step_size = 5;
    % need to exclude subjects who are not part of our current N
    RID_subset = RIDs';
    RID_subset = RID_subset(N);
    for a=1:length(pet_to_run)
        % check if QC was run. If it has, use this info when averaging
        if isequal(sum(strcmpi(mastersheet.Properties.VariableNames, 'summarizedFails')),0)
            table_to_use = mastersheet(contains(mastersheet.RID , RID_subset) & contains(mastersheet.ScanType, char(pet_to_run{a}),'IgnoreCase',true) & mastersheet.TimepointReference ==1, :);
        else
            % we have QC info 
            table_to_use = mastersheet(contains(mastersheet.RID , RID_subset) & contains(mastersheet.ScanType, char(pet_to_run{a}),'IgnoreCase',true) & mastersheet.TimepointReference ==1 & mastersheet.summarizedFails ==0, :);
        end
        % now collect the images from table_to_use. These have been QC'd and filtered so
        % everything should be good to go
        all_images = {};
        for scan=1:height(table_to_use)
            scan_path = fullfile(rootPath, char(table_to_use.RID(scan)), char(table_to_use.Date(scan)), char(table_to_use.Modality(scan)), char(table_to_use.ScanType(scan)), char(table_to_use.ScanName(scan)));
            if strcmpi(frames, {''})
                dir_mod = '';
            else
                possibilities = get_folders(scan_path);
                to_use = strcmpi(possibilities, frames);
                dir_mod = [filesep, char(possibilities(to_use))];
            end
            dir_path = [scan_path, dir_mod];
            % grab the appropriate scan
            all_images{end+1} = get_files(dir_path, [avgPrefix, '*.nii']);
        end
        if length(all_images) < step_size
            % we only need to create one average
                imgArray = {};
                for c=1:length(all_images)
                    imgVol = spm_vol(all_images{c});
                    [imgCell, ~] = spm_read_vols(imgVol);
                    imgVec = imgCell(:)';
                    imgArray{c} = (imgVec);
                end
                imgSeries=cell2mat(imgArray');
                groupAverageSeries= nanmean(imgSeries);
                %% first dimension
                m = length(imgCell(:,1,1));
                %% second dimension
                n = length(imgCell(1,:,1));
                %% third dimension
                p = length(imgCell(1,1,:));
                groupAverageMap = permute(reshape(groupAverageSeries,m,n,p),[1 2 3]);
                %% update filename
                imgVol.fname = [study_avg_dir, filesep, study_name, '_', char(pet_to_run{a}), '_', avgPrefix, '_', num2str(N(1)), '-', num2str(N(end)), '.nii'];              
                %% force rescale of the output image
                if true(isfield(imgVol,'pinfo'))
                    imgVol=rmfield(imgVol,'pinfo');
                end
                %% write it out
                spm_write_vol(imgVol,groupAverageMap);
        else
            % now create averages in small increments for speed, then average
            % the averages
            %create an output dir to temporarily store the averages
            temp_dir = [study_avg_dir, filesep, 'temp'];
            % make sure there isn't an old directory with files from a previous
            % round
            check = get_files(temp_dir, '*.nii');
            if ~isempty(check)
                rmdir(temp_dir, 's')
            end
            if ~exist(temp_dir, 'dir')
                mkdir(temp_dir)
            end
            num_runs = round(height(table_to_use)/step_size);
            first = 1;
            last = step_size;
            for fl=1:num_runs
                if isequal(fl, num_runs)
                    % on the last run, we may have a different number of scans
                    % to average (too few or too many). just grab from our
                    % starting point to the end of the array
                    to_avg = all_images(first:end);
                else
                    to_avg = all_images(first:last);
                    % now update the values
                    first = last+1;
                    last = (fl+1)*step_size;
                end
                % calculate our average image
                imgArray = {};
                for c=1:length(to_avg)
                    imgVol = spm_vol(to_avg{c});
                    [imgCell, ~] = spm_read_vols(imgVol);
                    imgVec = imgCell(:)';
                    imgArray{c} = (imgVec);
                end
                imgSeries=cell2mat(imgArray');
                groupAverageSeries= nanmean(imgSeries);
                %% first dimension
                m = length(imgCell(:,1,1));
                %% second dimension
                n = length(imgCell(1,:,1));
                %% third dimension
                p = length(imgCell(1,1,:));
                groupAverageMap = permute(reshape(groupAverageSeries,m,n,p),[1 2 3]);
                %% update filename
                imgVol.fname = [temp_dir, filesep, 'temp_average_', num2str(fl), '.nii'];
                %% force rescale of the output image
                if true(isfield(imgVol,'pinfo'))
                    imgVol=rmfield(imgVol,'pinfo');
                end
                %% write it out
                spm_write_vol(imgVol,groupAverageMap);
            end
            % now we have all temp averages. average these.
            final_files = cellstr(get_files(temp_dir, 'temp_average*.nii'));
            avgArray = {};
            for f=1:length(final_files)
                hdr = spm_vol(final_files{f});
                [imgCell, ~] = spm_read_vols(hdr);
                imgVec = imgCell(:)';
                avgArray{f} = (imgVec);
            end
            imgSeries=cell2mat(avgArray');
            groupAverageSeries= nanmean(imgSeries);
            %% first dimension
            m = length(imgCell(:,1,1));
            %% second dimension
            n = length(imgCell(1,:,1));
            %% third dimension
            p = length(imgCell(1,1,:));
            groupAverageMap = permute(reshape(groupAverageSeries,m,n,p),[1 2 3]);
            %% update filename
            imgVol.fname = [study_avg_dir, filesep, study_name, '_', char(pet_to_run{a}), '_', avgPrefix, '_', num2str(N(1)), '-', num2str(N(end)), '.nii'];
            %% force rescale of the output image
            if true(isfield(imgVol,'pinfo'))
                imgVol=rmfield(imgVol,'pinfo');
            end
            %% write it out
            spm_write_vol(imgVol,groupAverageMap);
            % remove temp files
            check = get_files(temp_dir, '*.nii');
            if ~isempty(check)
                rmdir(temp_dir, 's')
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