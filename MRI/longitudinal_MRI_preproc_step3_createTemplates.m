%% MRI Preproc Step 3: template creation
%% function to create a custom population template with geodesic shooting (Ashburner and Friston, 2011)
% code by Hayley R. C. Shanks, edits by Taylor W. Schmitz, Kate M. Onuska
function longitudinal_MRI_preproc_step3_createTemplates(N, study_path, mri_to_run, study_name, field_str, avgMode, indMode,varargin)
rootPath = [study_path, filesep, 'first_level'];
sheet_path = [study_path, filesep, 'spreadsheets'];
shootTemplatePath = [study_path, filesep, study_name, '_shoot_template'];
if ~exist(shootTemplatePath, 'dir')
    mkdir(shootTemplatePath);
end
% parse user inputs
pars = inputParser;
% by default, all modules are set to off.
default_run =0;
% required inputs
addRequired(pars, 'N', @isnumeric);
addRequired(pars, 'study_path', @ischar);
addRequired(pars, 'mri_to_run', @iscell);
addRequired(pars, 'study_name', @ischar);
addRequired(pars, 'field_str', @isnumeric);
addRequired(pars, 'avgMode', @iscell);
addRequired(pars, 'indMode', @iscell);
%add optional inputs.
addParameter(pars,'buildCustomTemplate',default_run,@isnumeric); 
addParameter(pars,'addToShoot',default_run,@isnumeric); 
addParameter(pars,'invertY',default_run,@isnumeric);
addParameter(pars,'shootTemplateOverride','',@ischar); % allows you to specify a path to a shoot template to normalize to. Use if you don't have one for your study
% parse the function inputs
parse(pars, N, study_path, mri_to_run, study_name, field_str, avgMode , indMode, varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
elseif isequal(length(indMode),2) && isequal(sum(strcmpi(indMode{2}, {'includeCrossSectionalSubjects', 'excludeCrossSectionalSubjects'})),0)
    disp(['input to indMode is set to ', indMode{2}, ' which is not a valid input'])
    disp('please review function documentation for the correct input arguments and try again')
    return
end
RIDs = get_folders(rootPath);
mbatch = 0;
sub_to_run = [1:length(RIDs)];
mastersheet = readtable(fullfile(sheet_path, [study_name, '_MRI_mastertable.csv']));
for t=1:length(mri_to_run)
    % create a subtable with only valid times based on QC
    if ~strcmpi(avgMode, 'off') && ~strcmpi(indMode{1}, 'off')
        disp('please select only one mode to use to create the template')
        return
    elseif strcmpi(avgMode, 'on')
        % people with valid midpoint average IQR's will have a value of 0 here
        valid_scans = mastersheet(contains(mastersheet.ScanType, mri_to_run{t}) & mastersheet.avg_IQR_fails == 0 & ~isnan(mastersheet.avg_CAT12_IQR) & mastersheet.FieldStrength == field_str,:);
        dir_mod = 'midpoint_average';
        temp_mod = 'average';
    else
        if strcmpi(indMode{1},'all_times')
            time_to_use =1;
            disp('Using time 1 segmentations to create the template.')
        elseif isnumeric(indMode{1}) && isequal(length(indMode{1}),1)
            time_to_use = indMode{1};
        else
            disp('Please choose only one timepoint to be used for template creation. If you choose all times, we will default to Order1')      
        end
        subset = mastersheet(mastersheet.Order ==time_to_use & mastersheet.TimepointReference ==1, :);
        valid_scans = subset(contains(mastersheet.ScanType, mri_to_run{t}) & subset.timepoint_fails == 0 & mastersheet.FieldStrength == field_str, :);
        dir_mod = '';
        temp_mod = num2str(time_to_use);
    end
    for subI=sub_to_run(N)
        % check if they are in the valid scans table
        sub_ind = strcmp(valid_scans.RID, RIDs{subI});
        if isequal(sum(sub_ind), 0)
            continue
        end
        % use the mastersheet to construct their filepath
        seg_dir = fullfile(rootPath, char(valid_scans.RID(sub_ind)), char(valid_scans.Date(sub_ind)), char(valid_scans.Modality(sub_ind)), char(valid_scans.ScanType(sub_ind)),char(valid_scans.ScanName(sub_ind)), dir_mod, 'mri');      
        if ~exist('output_dir', 'var')
            % this means this is the first person to make it into this
            % loop. Make them the output dir. note that subI may not be 1
            % if the 1st sub has cross sectional data
            output_dir = seg_dir;
        end
        rc1{subI}=get_files(seg_dir,'rp1*_rigid.nii');
        rc2{subI}=get_files(seg_dir,'rp2*_rigid.nii');
    end
    if pars.Results.buildCustomTemplate
        mbatch = mbatch+1;
        %% Remove empty cell array contents
        rc1=rc1(~cellfun('isempty',rc1));
        rc2=rc2(~cellfun('isempty',rc2));
        matlabbatch{mbatch}.spm.tools.shoot.warp.images = {rc1'
            rc2'}';
         spm_jobman('run',matlabbatch(mbatch))
         template = get_files(output_dir, 'Template*.nii');
        for i =1:size(template,1)
            [~, name, ext] = fileparts(template(i,:));
            copyfile(template(i,:), [shootTemplatePath, filesep, study_name, '_', temp_mod, '_', num2str(RIDs{1}) '-'  RIDs{length(N)}, '_',name, ext]);
        end
    end
    if pars.Results.addToShoot
        if ~isempty(pars.Results.shootTemplateOverride)
            shootTemplate=get_files(pars.Results.shootTemplateOverride,['*Template*.nii']);
        else
            shootTemplate=get_files(fullfile(shootTemplatePath),[study_name '*' 'Template*.nii']);
        end
        mbatch = mbatch + 1;
        %% Remove empty cell array contents
        rc1=rc1(~cellfun('isempty',rc1));
        rc2=rc2(~cellfun('isempty',rc2));
        matlabbatch{mbatch}.spm.tools.shoot.warp1.images = {rc1'
            rc2'}';
        matlabbatch{mbatch}.spm.tools.shoot.warp1.templates = {
            shootTemplate(1,:)
            shootTemplate(2,:)
            shootTemplate(3,:)
            shootTemplate(4,:)
            shootTemplate(5,:)
            };        
        spm_jobman('run',matlabbatch(mbatch))
    end
    if pars.Results.invertY
        for subI=sub_to_run(N)
            % check if they are in the valid scans table
            sub_ind = strcmp(valid_scans.RID, RIDs{subI});
            if isequal(sum(sub_ind), 0)
                continue
            end
            % use the mastersheet to construct their filepath
            seg_dir = fullfile(rootPath, char(valid_scans.RID(sub_ind)), char(valid_scans.Date(sub_ind)), char(valid_scans.Modality(sub_ind)), char(valid_scans.ScanType(sub_ind)),char(valid_scans.ScanName(sub_ind)), dir_mod, 'mri');
            if ~exist('output_dir', 'var')
                % this means this is the first person to make it into this
                % loop. Make them the output dir. note that subI may not be 1
                % if the 1st sub has cross sectional data
                output_dir = seg_dir;
            end
            yWarp = get_files(seg_dir,'y_r*Template.nii');
            native_seg = get_files(seg_dir,'p1*.nii');
            [ypt yfn yex] = fileparts(yWarp);
            new_fn = ['i' yfn yex];
            mbatch = mbatch+1;
            matlabbatch{mbatch}.spm.util.defs.comp{1}.inv.comp{1}.def = {yWarp};
            % use the atlas space as the target image bb/vox space
            matlabbatch{mbatch}.spm.util.defs.comp{1}.inv.space = {native_seg};
            matlabbatch{mbatch}.spm.util.defs.out{1}.savedef.ofname = new_fn;
            matlabbatch{mbatch}.spm.util.defs.out{1}.savedef.savedir.saveusr = {[ypt, filesep]};
            spm_jobman('run',matlabbatch(mbatch))
            iyWarp = get_files(fullfile(seg_dir), 'y_iy_r*Template.nii');
            movefile(iyWarp, [ypt, filesep, 'i', yfn, yex]);
        end
    end
    clear output_dir
end
end

