%% MRI Preproc Step 4: longitudinal modulation and warping to template space
% code by Hayley R. C. Shanks, edits by Taylor W. Schmitz, Kate M. Onuska

%% if using warp to shoot in average mode
% default will warp the jac files in the scan dir
% if you need to warp the p1 midpoint average files (for between subject
% QC), set warpPath = {['midpoint_average', filesep, 'mri']} and warpFilt = {'p1'}
function longitudinal_MRI_preproc_step4_modulation(N, study_path, mri_to_run, study_name, field_str, avgMode, indMode, tissue_compartments, warpType, warpPath, warpFilt, avgPrefix, varargin)
shootTemplatePath = [study_path, filesep, study_name, '_shoot_template'];
% raw participant data should live in a first level subdirectory of the study path
rootPath = [study_path, filesep, 'first_level'];
sheet_path = [study_path, filesep, 'spreadsheets'];
study_avg_dir = [study_path, filesep, study_name, '_averages'];
%%%%%%%%%%%%%%%%%%vl%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
addRequired(pars, 'tissue_compartments', @iscell);
addRequired(pars, 'warpType', @iscell);
addRequired(pars, 'warpPath', @iscell);
addRequired(pars, 'warpFilt', @iscell);
addRequired(pars, 'avgPrefix', @ischar);
%add optional inputs.
addParameter(pars,'multiplyJacobian',default_run,@isnumeric); 
addParameter(pars,'warpToShoot',default_run,@isnumeric);
addParameter(pars,'shootTemplateOverride','',@ischar); % allows you to specify a path to a shoot template to normalize to. Use if you don't have one for your study
addParameter(pars,'createAvg',default_run,@isnumeric);
% parse the function inputs
parse(pars, N, study_path, mri_to_run, study_name, field_str, avgMode, indMode, tissue_compartments, warpType, warpPath, warpFilt, avgPrefix, varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RIDs = get_folders(rootPath);
mbatch = 0;
sub_to_run = [1:length(RIDs)];
% check that the mode variables are set correctly 
if isequal(sum(strcmpi(avgMode, {'on', 'off'})),0)
    disp(['input to avgMode is set to ', avgMode{1}, ' which is not a valid input'])
    disp('please review function documentation for the correct input arguments and try again')
elseif ~isnumeric(indMode{1}) && isequal(sum(strcmpi(indMode{1}, {'all_times', 'off'})),0)
    disp(['input to indMode is set to ', indMode{1}, ' which is not a valid input'])
    disp('please review function documentation for the correct input arguments and try again')
end
if isequal(length(indMode),1) && strcmpi(indMode, 'off')
    % add a second cell just to make code easier    
    indMode{1,2} = 'off';
elseif isequal(length(indMode),2) && isequal(sum(strcmpi(indMode{2}, {'includeCrossSectionalSubjects', 'excludeCrossSectionalSubjects'})),0)
    disp(['input to indMode is set to ', indMode{2}, ' which is not a valid input'])
    disp('please review function documentation for the correct input arguments and try again')
end
mastersheet = readtable(fullfile(sheet_path, [study_name, '_MRI_mastertable.csv']));
% add a keys variable so that we can easily keep track
% of rows when subsetting this table
mastersheet.keys = [1:height(mastersheet)]';
if strcmpi(tissue_compartments,'gm')
    seg_match = 'p1';
elseif strcmpi(tissue_compartments,'wm')
    seg_match = 'p2';
elseif strcmpi(tissue_compartments,'csf')
    seg_match = 'p3';
elseif isempty(tissue_compartments)
    seg_match = '';
end
for t=1:length(mri_to_run)
    % define datasets to be used for each mode which considers QC info
    if strcmpi(avgMode, 'on')
        avg_subset = mastersheet(contains(mastersheet.ScanType, mri_to_run{t}) & mastersheet.avg_IQR_fails == 0 & ~isnan(mastersheet.avg_CAT12_IQR) & mastersheet.TimepointReference ==1 & mastersheet.FieldStrength == field_str, :);
    end
    if ~strcmpi(indMode{1}, 'off')
        % ind subset will not contain subjects who only have cross sectional
        % data if indMode has been previously run with
        % 'excludeCrossSectionalSubjects' (they will have a 1 in timepoint
        % fails)
        ind_subset = mastersheet(mastersheet.timepoint_fails == 0 & mastersheet.TimepointReference ==1 , :);
        if isnumeric(indMode{1})
            ind_subset = mastersheet(contains(mastersheet.ScanType, mri_to_run{t}) & mastersheet.timepoint_fails == 0 & mastersheet.TimepointReference ==1 & mastersheet.FieldStrength == field_str, :);
        end
    end
    if pars.Results.multiplyJacobian
        if ~strcmpi(avgMode, 'on')
            disp('avgMode is off, but longitudinal modulation is being attempted. Turn avgMode on to run this module');
        end
        for subI=sub_to_run(N)
            % check if that subject is in the avg_subset (meaning they
            % passed segmentation QC)   
            if isequal(sum(strcmp(avg_subset.RID, RIDs{subI})),0)
                disp(['unable to perform longitudinal modulation for ', RIDs{subI}])
                disp('Data is either not longitudinal or subject failed QC')
            else
                sub_dir = fullfile(rootPath,RIDs{subI});
                sub_table = get_files(sub_dir, 'scan_info*.csv');
                % this table has all of the dates and types of scans for this P
                sub_table = readtable(sub_table);
                if isempty(sub_table)
                    continue
                end
                table_to_run = spread_table_filter(sub_table, mri_to_run{t}, field_str);
                avg_dir = fullfile(sub_dir, char(table_to_run.Date(1)), char(table_to_run.Modality(1)), char(table_to_run.ScanType(1)), char(table_to_run.ScanName(1)), 'midpoint_average');
                if isempty(seg_match)
                    disp('no tissue compartment selected for longitudinal modulation.');
                    return
                end
                seg = get_files(fullfile(avg_dir, 'mri'), [seg_match,  '*.nii']);
                % now we collect the jacobians for each timepoint and
                % multiply with our segmented image
                for row=1:height(table_to_run)
                    current_dir = fullfile(sub_dir, char(table_to_run.Date(row)), char(table_to_run.Modality(row)), char(table_to_run.ScanType(row)), char(table_to_run.ScanName(row)));
                    % check if they have an intra first
                    jacobian = get_files(current_dir, ['j_intra*',  '.nii']);
                    if isempty(jacobian)
                        % check for a regular jacobian
                        jacobian = get_files(current_dir, ['j_rac*',  '.nii']);
                    end
                    [~, jac_name, ex] = fileparts(jacobian);
                    output_name = ['jac' seg_match '_' mri_to_run{t} '_' jac_name ex];
                    % multiply the jacobian from that timepoint by the p file
                    % to get the longitudinally modulated image
                    mbatch = mbatch + 1;
                    matlabbatch{mbatch}.spm.util.imcalc.input = {seg;jacobian};
                    matlabbatch{mbatch}.spm.util.imcalc.output = output_name;
                    matlabbatch{mbatch}.spm.util.imcalc.outdir = {current_dir};
                    matlabbatch{mbatch}.spm.util.imcalc.expression = 'i1.*i2';
                    matlabbatch{mbatch}.spm.util.imcalc.var = struct('name', {}, 'value', {});
                    matlabbatch{mbatch}.spm.util.imcalc.options.dmtx = 0;
                    matlabbatch{mbatch}.spm.util.imcalc.options.mask = 0;
                    % use 4th degree sinc interpolation (0=NN, 1=trilinear)
                    matlabbatch{mbatch}.spm.util.imcalc.options.interp = 4;
                    % data type. int16 = 4 (default), 16 = float32
                    matlabbatch{mbatch}.spm.util.imcalc.options.dtype = 4;
                    spm_jobman('run',matlabbatch(mbatch))
                end
            end
        end
    end      
    if pars.Results.warpToShoot
        no_y = [];
        % define the shoot template
        if ~isempty(pars.Results.shootTemplateOverride)
            shootTemplate=get_files(pars.Results.shootTemplateOverride,['*Template*.nii']);
        else
            shootTemplate=get_files(fullfile(shootTemplatePath),[study_name '*' 'Template*.nii']);
        end
        [temp_BB,temp_vx] = spm_get_bbox(shootTemplate(5,:)); % use the newest template
        % define which type of tissue compartment we are working with
        if strcmpi(avgMode, 'on')
            if isempty(warpFilt) || strcmp(warpFilt, '')
                % default for avg mode is the jacobians
                prefix = ['jac' seg_match];
            else 
                prefix = char(warpFilt);
            end
            for subI=sub_to_run(N)
                % find them in the QC table
                if isequal(sum(strcmp(avg_subset.RID, RIDs{subI})),0)
                    % skip them. They failed QC
                else
                    sub_dir = fullfile(rootPath,RIDs{subI});
                    sub_table = get_files(sub_dir, 'scan_info*.csv');
                    % this table has all of the dates and types of scans for this P
                    sub_table = readtable(sub_table);
                    table_to_run = spread_table_filter(sub_table, mri_to_run{t}, field_str);
                    % table to run should never be empty at this point, since
                    % the subject was in the usable subset
                    % define the directory where their midpoint average lives
                    avg_dir = fullfile(sub_dir, char(table_to_run.Date(1)), char(table_to_run.Modality(1)), char(table_to_run.ScanType(1)), char(table_to_run.ScanName(1)), 'midpoint_average', 'mri');
                    y_file = get_files(avg_dir, 'y_rp1*rigid*Template.nii');
                    if isempty(y_file)
                        no_y = [no_y, subI];
                        continue
                    end
                    to_warp = {};
                    for row=1:height(table_to_run)
                        % collect the images to warp
                        current_dir = fullfile(sub_dir, char(table_to_run.Date(row)), char(table_to_run.Modality(row)), char(table_to_run.ScanType(row)), char(table_to_run.ScanName(row)));
                        if isempty(warpPath)
                            % default location in avgMode would be in the
                            % scan directory, so no modification needed
                        else
                            current_dir = [current_dir, filesep, char(warpPath)];
                        end
                        to_warp{end+1} = get_files(current_dir, [prefix, '*.nii']);
                    end
                    disp(current_dir)
                    %clear out empty cells
                    to_warp=to_warp(~cellfun('isempty',to_warp));
                    % create a cell of y files the same size as the number of images to warp
                    deformations = cell(length(to_warp),1);
                    [deformations{:}] = deal(y_file);  
                    mbatch = mbatch+1;
                    matlabbatch{mbatch}.spm.tools.shoot.norm.template = {''}; % leave this empty so that we don't include a transformation to mni space
                    matlabbatch{mbatch}.spm.tools.shoot.norm.data.subjs.deformations = deformations;
                    matlabbatch{mbatch}.spm.tools.shoot.norm.data.subjs.images = {to_warp'};
                    matlabbatch{mbatch}.spm.tools.shoot.norm.vox = temp_vx;
                    matlabbatch{mbatch}.spm.tools.shoot.norm.bb = temp_BB;
                    % Preserve Concentrations: Smoothed spatially normalised images (sw*) represent weighted averages
                    % of the signal under the smoothing kernel, approximately preserving the intensities of the
                    % original images. This option is currently suggested for eg fMRI.
                    % Preserve Total: Smoothed and spatially normalised images preserve the total amount of signal
                    % from each region in the images (smw*). Areas that are expanded during warping are correspondingly
                    % reduced in intensity. This option is suggested for VBM.
                    if strcmpi(warpType, 'nomod')
                        matlabbatch{mbatch}.spm.tools.shoot.norm.preserve = 0;
                    elseif strcmpi(warpType, 'mod')
                        matlabbatch{mbatch}.spm.tools.shoot.norm.preserve = 1;
                    end
                    matlabbatch{mbatch}.spm.tools.shoot.norm.fwhm = 0;
                    spm_jobman('run',matlabbatch(mbatch))
                end
            end
            save('no_y.mat', 'no_y');
        end
        if ~strcmpi(indMode, 'off')
            if isempty(warpFilt) || strcmp(warpFilt, '')
                % default for ind mode is the segmentation itself
                prefix = seg_match;
            else 
                prefix = char(warpFilt);
            end
            for subI=sub_to_run(N)
                if isequal(sum(strcmp(ind_subset.RID, RIDs{subI})),0)
                    % they failed Qc. Skip.
                else
                    % ind subset already has been filtered to exclude
                    % people who failed QC, and to remove any timepoints we
                    % are not currently interested in. 
                    % work from this table instead of the sub tables
                    table_to_run = ind_subset(strcmp(ind_subset.RID, RIDs{subI}), :);
                    if isequal(max(table_to_run.Order),1) && strcmpi(indMode{2}, 'excludeCrossSectionalSubjects')
                        continue
                        % these people would be removed at QC if exclude
                        % cross sectional subjects was specfied in previous
                        % steps, but check in case
                    end
                    % make sure rows are sorted
                    table_to_run = sortrows(table_to_run, 'Date');
                    y_dir = fullfile(rootPath, RIDs{subI}, char(table_to_run.Date(1)), char(table_to_run.Modality(1)), char(table_to_run.ScanType(1)), char(table_to_run.ScanName(1)), 'mri');              
                    y_file = get_files(y_dir, 'y_rp1*rigid*Template.nii');
                    to_warp = {};
                    for row=1:height(table_to_run)
                        % collect the images to warp
                        current_dir = fullfile(rootPath, RIDs{subI}, char(table_to_run.Date(row)), char(table_to_run.Modality(row)), char(table_to_run.ScanType(row)), char(table_to_run.ScanName(row)));
                        if isempty(warpPath) || strcmpi(warpPath, '')
                            % default location in indMode would be in the
                            % mri subdirectory of the scan directory
                            current_dir = [current_dir, filesep, 'mri'];
                        else
                            current_dir = [current_dir, filesep, char(warpPath)];
                        end
                    to_warp{end+1} = get_files(current_dir, [char(prefix) '*.nii']);                      
                    end
                    %clear out empty cells
                    to_warp=to_warp(~cellfun('isempty',to_warp));
                    % create a cell of y files the same size as the number of images to warp
                    deformations = cell(length(to_warp),1);
                    [deformations{:}] = deal(y_file);  
                    mbatch = mbatch+1;
                    matlabbatch{mbatch}.spm.tools.shoot.norm.template = {''}; % leave this empty so that we don't include a transformation to mni space
                    matlabbatch{mbatch}.spm.tools.shoot.norm.data.subjs.deformations = deformations;
                    matlabbatch{mbatch}.spm.tools.shoot.norm.data.subjs.images = {to_warp'};
                    matlabbatch{mbatch}.spm.tools.shoot.norm.vox = temp_vx;
                    matlabbatch{mbatch}.spm.tools.shoot.norm.bb = temp_BB;
                    % Preserve Concentrations: Smoothed spatially normalised images (sw*) represent weighted averages
                    % of the signal under the smoothing kernel, approximately preserving the intensities of the
                    % original images. This option is currently suggested for eg fMRI.
                    % Preserve Total: Smoothed and spatially normalised images preserve the total amount of signal
                    % from each region in the images (smw*). Areas that are expanded during warping are correspondingly
                    % reduced in intensity. This option is suggested for VBM.
                    if strcmpi(warpType, 'nomod')
                        matlabbatch{mbatch}.spm.tools.shoot.norm.preserve = 0;
                    elseif strcmpi(warpType, 'mod')
                        matlabbatch{mbatch}.spm.tools.shoot.norm.preserve = 1;
                    end
                    matlabbatch{mbatch}.spm.tools.shoot.norm.fwhm = 0;
                    spm_jobman('run',matlabbatch(mbatch))
                end
            end
        end
    end
           
    if pars.Results.createAvg
        step_size = 5;
        % use ind_subset or avg_subset
        if strcmpi(avgMode, 'on') && ~strcmpi(indMode{1}, 'off')
            disp('both indMode and avgMode are on. Please select only one mode to create an average with');
            return
        else
            if strcmpi(avgMode, 'on')
                table_to_use = avg_subset;
                dir_mod = ['midpoint_average', filesep, 'mri'];
            elseif strcmpi(indMode{1}, 'all_times')
                disp('only one timepoint of data will be used for average creation. defaulting to time1')
                table_to_use = ind_subset;
                table_to_use = table_to_use(table_to_use.Order ==1, :);
                dir_mod = ['mri'];
            elseif isnumeric(indMode{1}) && length(indMode{1}) > 1 
                disp('only one timepoint of data will be used for average creation. defaulting to time1')
                table_to_use = ind_subset;
                table_to_use = table_to_use(table_to_use.Order ==1, :);
                dir_mod = ['mri'];
            else
                table_to_use = ind_subset;
                table_to_use = table_to_use(table_to_use.Order ==indMode{1}, :);
                dir_mod = ['mri'];
            end
        end
        % need to exclude subjects who are not part of our current N 
        RID_subset = RIDs';
        RID_subset = RID_subset(N);
        table_to_use = table_to_use(contains(table_to_use.RID , RID_subset) & table_to_use.Order == 1, :);
        % now collect the images from table_to_use. These have been QC'd and filtered so
        % everything should be good to go
        all_images = {};
        for scan=1:height(table_to_use)
            dir_path = fullfile(rootPath, char(table_to_use.RID(scan)), char(table_to_use.Date(scan)), char(table_to_use.Modality(scan)), char(table_to_use.ScanType(scan)), char(table_to_use.ScanName(scan)));
            % add the modifier which will depend on what mode we're in
            dir_path = [dir_path, filesep, dir_mod];
            % grab the appropriate scan
            all_images{end+1} = get_files(dir_path, [avgPrefix, '*.nii']);
        end
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
            % first dimension
            m = length(imgCell(:,1,1));
            % second dimension
            n = length(imgCell(1,:,1));
            % third dimension
            p = length(imgCell(1,1,:));       
            groupAverageMap = permute(reshape(groupAverageSeries,m,n,p),[1 2 3]);
            % update filename
            imgVol.fname = [temp_dir, filesep, 'temp_average_', num2str(fl), '.nii'];
            % force rescale of the output image
            if true(isfield(imgVol,'pinfo'))
                imgVol=rmfield(imgVol,'pinfo');
            end
            % write it out
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
            % first dimension
            m = length(imgCell(:,1,1));
            % second dimension
            n = length(imgCell(1,:,1));
            % third dimension
            p = length(imgCell(1,1,:));       
            groupAverageMap = permute(reshape(groupAverageSeries,m,n,p),[1 2 3]);
            % update filename
            imgVol.fname = [study_avg_dir, filesep, study_name, '_', mri_to_run{t}, '_', avgPrefix, '_', num2str(N(1)), '-', num2str(N(end)), '.nii'];
            % force rescale of the output image
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
