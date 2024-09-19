%%MRI preproc: Final QC step
% code by Hayley R. C. Shanks, with edits by Taylor W. Schmitz, Kate M. Onuska
function longitudinal_MRI_preproc_step5_finalQC(N, study_path, mri_to_run, study_name, field_str, avgMode, indMode, tissue_compartments, warpType, warpPath, qcFiltMRI, withinSubThresh, varargin)
% raw participant data should live in a first level subdirectory of the study path
rootPath = [study_path, filesep, 'first_level'];
sheet_path = [study_path, filesep, 'spreadsheets'];
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
addRequired(pars, 'qcFiltMRI', @iscell);
addRequired(pars, 'withinSubThresh', @isnumeric);
%add optional inputs.
addParameter(pars,'withinSubQC',default_run,@isnumeric);
addParameter(pars,'betweenSubQC',default_run,@isnumeric);
addParameter(pars,'extractICV',default_run,@isnumeric);
addParameter(pars,'summarizeQC',default_run,@isnumeric);
% parse the function inputs
parse(pars, N, study_path, mri_to_run, study_name, field_str, avgMode, indMode, tissue_compartments, warpType, warpPath, qcFiltMRI, withinSubThresh,varargin{:});
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
    % define datasets to be used for each mode 
    if strcmpi(avgMode, 'on')
        avg_subset = mastersheet(contains(mastersheet.ScanType, mri_to_run{t}) & mastersheet.avg_IQR_fails == 0 & mastersheet.TimepointReference ==1 & mastersheet.FieldStrength == field_str & mastersheet.Order ==1, :);
    end
    if ~strcmpi(indMode{1}, 'off')
        % ind subset will not contain subjects who only have cross sectional
        % data if indMode has been previously run with
        % 'excludeCrossSectionalSubjects' (they will have a 1 in timepoint
        % fails)
        ind_subset = mastersheet(contains(mastersheet.ScanType, mri_to_run{t}) & mastersheet.timepoint_fails == 0 & mastersheet.TimepointReference ==1 & mastersheet.FieldStrength == field_str, :);
        if isnumeric(indMode{1})
            ind_subset = mastersheet(contains(mastersheet.ScanType, mri_to_run{t}) & mastersheet.timepoint_fails == 0 & mastersheet.TimepointReference ==1 & mastersheet.FieldStrength == field_str & ismember(mastersheet.Order, indMode{1}), :);
        end
    end
    
    if pars.Results.withinSubQC
        % 2 types of QC to do - a within subject QC if data is
        % longitudinal, and a between subject QC using time 1 or the
        % midpoint average
        if strcmpi(warpType, 'mod')
            mod_type = 'mw';
        elseif strcmpi(warpType, 'nomod')
            mod_type = 'w';
        end
        % start with QC for avgMode
        if strcmpi(avgMode, 'on')
            if ~contains(mastersheet.Properties.VariableNames, 'withinSubAvgWarpQC')
                mastersheet.withinSubAvgWarpQC = nan(height(mastersheet),1);
            end
            if isempty(qcFiltMRI) || strcmp(qcFiltMRI, '')
                % default for avg mode is the jacobians
                prefix = [mod_type 'jac' seg_match];
            else
                prefix = [mod_type, char(qcFiltMRI)];
            end
            % use the QC table of subjects having an average image
            for subI=sub_to_run(N)
                % find them in the QC table
                if isequal(sum(strcmp(avg_subset.RID, RIDs{subI})),0)
                    % skip them. They failed QC
                else
                    sub_dir = fullfile(rootPath,RIDs{subI});
                    sub_table = get_files(sub_dir, 'scan_info*.csv');
                    sub_table = readtable(sub_table);
                    table_to_run = spread_table_filter(sub_table, mri_to_run{t}, field_str);
                    imgToQC = {};
                    for row=1:height(table_to_run)
                        % collect the images to warp
                        current_dir = fullfile(sub_dir, char(table_to_run.Date(row)), char(table_to_run.Modality(row)), char(table_to_run.ScanType(row)), char(table_to_run.ScanName(row)));
                        if isempty(warpPath)
                            % default location in avgMode would be in the
                            % scan directory, so no modification needed
                        else
                            current_dir = [current_dir, filesep, char(warpPath)];
                        end
                        imgToQC{row} = get_files(current_dir, [prefix, '*.nii']);
                    end
                    sub_rows = strcmp(mastersheet.RID, RIDs{subI});
                    modality_rows = strcmp(mastersheet.ScanType, mri_to_run{t});
                    sub_modality = sum([sub_rows, modality_rows, mastersheet.TimepointReference], 2);% rowsums of these concatenated arrays tell us which scans have the right sub and modality and the ref
                    disp(RIDs{subI})
                    covStruct = cat_stat_check_cov(struct('data_vol',{{ imgToQC' }'} ,'gap',3,'c',[],'data_xml',{{}}));
                    all_covs = covStruct.covmat; % gives the covariance between each image and all other images
                    % remove the 1's from each row because that's the
                    % covariance with itself which isn't important
                    all_covs(all_covs ==1) = NaN;
                    mean_covs = mean(all_covs,2, 'omitnan'); % get the row means, omitting nan
                    mastersheet.withinSubAvgWarpQC(sub_modality == 3) = mean_covs;
                end
            end
        end
        if ~strcmpi(indMode{1}, 'off')
            if ~contains(mastersheet.Properties.VariableNames, 'betweenSubTimeWarpQC')
                mastersheet.betweenSubTimeWarpQC = nan(height(mastersheet),1);
            end
            if isempty(qcFiltMRI) || strcmp(qcFiltMRI, '')
                % default for ind mode is the segmentation itself
                prefix = [mod_type, seg_match];
            else 
                prefix = [mod_type, char(qcFiltMRI)];
            end
            % first ensure data is longitudinal
            if isnumeric(indMode{1}) && isequal(length(indMode{1}), 1)
                disp('you have only selected one timepoint of data to analyze. within subject QC will be skipped')
                disp('if you would like to perform this step, update indMode');
            elseif isequal(max(ind_subset.Order),1)
                disp('your data is not longitudinal. within subject QC will be skipped')
            else
                for subI=sub_to_run(N)
                    if isequal(sum(strcmp(ind_subset.RID, RIDs{subI})),0)
                        % they failed Qc. Skip.
                    else
                        % ind subset already has been filtered to exclude
                        % people who failed QC, and to remove any timepoints we
                        % are not currently interested in. 
                        % work from this table instead of the sub tables
                        table_to_run = ind_subset(strcmp(ind_subset.RID, RIDs{subI}), :);
                        if isequal(max(table_to_run.Order),1) 
                            continue
                            % possible some subjects with 1 timepoint are
                            % in ind_subset if
                            % includeCrossSectionalSubjects was on. skip
                            % them
                        end
                        % make sure rows are sorted
                        table_to_run = sortrows(table_to_run, 'Date');
                        imgToQC ={};
                        for row=1:height(table_to_run)
                            % collect the images to warp
                            current_dir = fullfile(rootPath, RIDs{subI}, char(table_to_run.Date(row)), char(table_to_run.Modality(row)), char(table_to_run.ScanType(row)), char(table_to_run.ScanName(row)));
                            if isempty(warpPath)
                                % default location in indMode would be in the
                                % mri subdirectory of the scan directory
                                current_dir = [current_dir, filesep, 'mri'];
                            else
                                current_dir = [current_dir, filesep, char(warpPath)];
                            end
                            imgToQC{end+1} = get_files(current_dir, [prefix, '*.nii']);                      
                        end
                    end
                    covStruct = cat_stat_check_cov(struct('data_vol',{{ imgToQC' }'} ,'gap',3,'c',[],'data_xml',{{}}));
                    all_covs = covStruct.covmat; % gives the covariance between each image and all other images
                    % remove the 1's from each row because that's the
                    % covariance with itself which isn't important
                    all_covs(all_covs ==1) = NaN;
                    mean_covs = mean(all_covs,2, 'omitnan'); % get the row means, omitting nan
                    % find the mastersheet rows corresponding to our sub
                    sub_rows = strcmp(mastersheet.RID, RIDs{subI});
                    modality_rows = strcmp(mastersheet.ScanType, mri_to_run{t});
                    sub_modality = sum([sub_rows, modality_rows], 2);% rowsums of these concatenated arrays tell us which scans have the right sub and modality
                    mastersheet.withinSubTimeWarpQC(sub_modality == 2) = mean_covs;
                end
            end
        end
        writetable(mastersheet, [sheet_path, filesep, study_name, '_MRI_mastertable.csv']);
    end
    % now do between subjects QC. If using avgMode, perform QC on the
    % midpoint average segment. Otherwise use the baseline segment
    if pars.Results.betweenSubQC
        if strcmpi(warpType, 'mod')
            mod_type = 'mw';
        elseif strcmpi(warpType, 'nomod')
            mod_type = 'w';
        end
        % start with QC for avgMode
        if strcmpi(avgMode, 'on')
            if ~contains(mastersheet.Properties.VariableNames,  'betweenSubAvgWarpQC')
                mastersheet.betweenSubAvgWarpQC = nan(height(mastersheet),1);
            end
            if isempty(qcFiltMRI) || strcmp(qcFiltMRI, '')
                % default for avg mode when using between subject QC is the
                % mwp file
                prefix = [mod_type seg_match];
            else
                prefix = [mod_type, char(qcFiltMRI)];
            end
            % use the QC table of subjects having an average image. We only
            % need the midpoint average here so work directly from this
            % table
            to_qc = {};
            keys_used=[];
            % need to exclude subjects who are not part of our current N 
            RID_subset = RIDs';
            RID_subset = RID_subset(N);
            for row=1:height(avg_subset)
                if ~contains(avg_subset.RID(row), RID_subset)
                    % skip
                else
                    sub_bl_dir = fullfile(rootPath, char(avg_subset.RID(row)), char(avg_subset.Date(row)), char(avg_subset.Modality(row)), char(avg_subset.ScanType(row)), char(avg_subset.ScanName(row)));
                    if isempty(warpPath) || strcmp(warpPath, '')
                        % default location in avgMode would be in the
                        % midpoint average mri directory
                        sub_bl_dir = [sub_bl_dir, filesep, 'midpoint_average', filesep, 'mri'];
                    else
                        sub_bl_dir = [sub_bl_dir, filesep, char(warpPath)];
                    end
                    to_qc{end+1} = get_files(sub_bl_dir, [prefix, '*.nii']);     
                    keys_used = [keys_used; avg_subset.keys(row)];
                end
            end
            
            covStruct = cat_stat_check_cov(struct('data_vol',{{ to_qc' }'} ,'gap',3,'c',[],'data_xml',{{}}));
            all_covs = covStruct.covmat; % gives the covariance between each image and all other images
            % remove the 1's from each row because that's the
            % covariance with itself which isn't important
            all_covs(all_covs ==1) = NaN;
            mean_covs = mean(all_covs,2, 'omitnan'); % get the row means, omitting nan
            % find the mastersheet rows corresponding to the rows in
            % avg_subset
            to_add = ismember(mastersheet.keys, keys_used);
            % add covariance ratings only to the scans that were actually
            % included in the check covar
            mastersheet.betweenSubAvgWarpQC(to_add) = mean_covs;
        end
        if ~strcmpi(indMode{1}, 'off')
            if isempty(qcFiltMRI) || strcmp(qcFiltMRI, '')
                % default for ind mode is the segmentation itself
                prefix = [mod_type, seg_match];
            else 
                prefix = [mod_type, char(qcFiltMRI)];
            end
            % we're only using one timepoint of data here
            if strcmpi(indMode{1}, 'all_times') 
                disp('Only one timepoint of data can be used for between subject QC. Defaulting to time 1');
                time_to_use = 1;
            elseif isnumeric(indMode{1}) && length(indMode{1}) > 1
                disp('Only one timepoint of data can be used for between subject QC. Defaulting to time 1');
                time_to_use =1;
            else
                time_to_use = indMode{1};
            end
                to_qc = {};
                for subI=sub_to_run(N)
                    if isequal(sum(strcmp(ind_subset.RID, RIDs{subI})),0)
                        % they failed Qc. Skip.
                    else
                        table_to_run = ind_subset(strcmp(ind_subset.RID, RIDs{subI}), :);
                        if isequal(max(table_to_run.Order),1) && strcmpi(indMode{2}, 'excludeCrossSectionalSubjects')
                            % possible some subjects with 1 timepoint are
                            % in ind_subset if
                            % includeCrossSectionalSubjects was on. skip
                            % them. Also drop them from ind subset so we
                            % have a record of this
                            ind_subset = ind_subset(~strcmpi(ind_subset.RID, RIDs{subI}), :);
                        else
                            table_to_run = spread_table_filter(table_to_run, mri_to_run{t}, field_str);
                            table_to_run = table_to_run(table_to_run.Order == time_to_use, :);
                            % table_to_run should only have 1 row now
                            if isempty(table_to_run)
                                % this could be empty if they have
                                % longitudinal data which at least 1 scan
                                % passes Qc, but the scan that passes QC
                                % isn't the scan we are currently doing QC
                                % on
                                ind_subset = ind_subset(~strcmpi(ind_subset.RID, RIDs{subI}), :);
                            else
                                time_dir = fullfile(rootPath, RIDs{subI}, char(table_to_run.Date(1)), char(table_to_run.Modality(1)), char(table_to_run.ScanType(1)), char(table_to_run.ScanName(1)));
                                if isempty(warpPath) || strcmp(warpPath, '')
                                    % default location in indMode would be in the
                                    % mri subdirectory of the scan directory
                                    time_dir = [time_dir, filesep, 'mri'];
                                else
                                    time_dir = [time_dir, filesep, char(warpPath)];
                                end
                                check = get_files(time_dir, [prefix, '*.nii']);  
                                to_qc{end+1} = get_files(time_dir, [prefix, '*.nii']);    
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
                % ind_subset
                to_add = ismember(mastersheet.keys, ind_subset.keys);
                mastersheet.betweenSubTimeWarpQC(to_add) = mean_covs
            end
        writetable(mastersheet, [sheet_path, filesep, study_name, '_MRI_mastertable.csv']);
    end
          
    if pars.Results.extractICV
        % this should only be run in avgMode or indMode, not both (their
        % TICV as calculated at baseline or on the average should be
        % similar
        if strcmpi(avgMode, 'on') && ~strcmpi(indMode{1}, 'off')
            disp('running extract ICV on the baseline scan and the average is redundant')
            disp('defaulting to using the average image. Set avgMode to off if you want timepoint data used')
            dir_mod = ['midpoint_average', filesep, 'report'];
            qc_subset = avg_subset;
            running = 'avg';
        elseif strcmpi(avgMode, 'on')
            dir_mod = ['midpoint_average', filesep, 'report'];
            qc_subset = avg_subset;
            running = 'avg';
        else
            % using timepoint data
            dir_mod = 'report';
            qc_subset = ind_subset;
            running = 'ind';
        end
        if ~contains(mastersheet.Properties.VariableNames, 'Total')
            % initialize the columns in mastersheet as nans
            mastersheet.Total = nan(height(mastersheet),1);
            mastersheet.GM = nan(height(mastersheet),1);
            mastersheet.WM = nan(height(mastersheet),1);
            mastersheet.CSF = nan(height(mastersheet),1);
            mastersheet.WMH = nan(height(mastersheet),1);
        end
        for subI=sub_to_run(N)
            if isequal(sum(strcmp(qc_subset.RID, RIDs{subI})),0)
                % skip them. They failed QC
            else
                sub_dir = fullfile(rootPath,RIDs{subI});
                sub_table = get_files(sub_dir, 'scan_info*.csv');
                sub_table = readtable(sub_table);
                % automatically do the SPREAD table filtering (see
                % function defined at bottom)
                table_to_run = spread_table_filter(sub_table, mri_to_run{t}, field_str);
                if strcmp(running, 'ind')
                    if isequal(max(table_to_run.Order),1) && strcmpi(indMode{2}, 'excludeCrossSectionalSubjects')
                        continue
                    end
                end
                % define the directory where their cat xml file lives
                cat_dir = fullfile(sub_dir, char(table_to_run.Date(1)), char(table_to_run.Modality(1)), char(table_to_run.ScanType(1)), char(table_to_run.ScanName(1)));
                cat_dir = [cat_dir, filesep, dir_mod];
                xml_file = get_files(cat_dir, 'cat*.xml');
                [xmlpt xmlfl xmlex] = fileparts(xml_file);
                p.data_xml = {xml_file};
                %% this option set to 0 will grab total and individual compartments (GM,WM,CSF,WMH)
                p.calcvol_TIV = 0;
                p.calcvol_name = [xmlpt filesep xmlfl '_extractedICV.txt'];
                icv_outputs = cat_stat_TIV(p);
                % find the row of mastersheet which corresponds to the
                % subject's BL directory. The array below is a bunch of
                % logicals where each column should equal 1 where the row
                % is correct
                sub_check = [strcmp(mastersheet.RID, RIDs{subI}), mastersheet.Date == table_to_run.Date(1), strcmp(mastersheet.Modality,table_to_run.Modality(1)), strcmp(mastersheet.ScanType, table_to_run.ScanType(1)), strcmp(mastersheet.ScanName, table_to_run.ScanName(1))] ;
                % row sums should match the number of cols 
                row_ind = find(sum(sub_check,2) == size(sub_check,2));
                mastersheet.Total(row_ind) = icv_outputs.calcvol(1);
                mastersheet.GM(row_ind) = icv_outputs.calcvol(2);
                mastersheet.WM(row_ind) = icv_outputs.calcvol(3);
                mastersheet.CSF(row_ind) = icv_outputs.calcvol(4);
                mastersheet.WMH(row_ind) = icv_outputs.calcvol(5);
            end
        end
        writetable(mastersheet, [sheet_path, filesep, study_name, '_MRI_mastertable.csv']);
    end
    if pars.Results.summarizeQC
        % make a column for the mastersheet which tells us if people fail
        % any measures of QC so we don't always have to index multiple
        % columns. we will also double check the people who have a mean
        % correlation less than 2 SDs below the mean correlation across
        % subs
        if strcmpi(warpType, 'mod')
            mod_type = 'mw';
        elseif strcmpi(warpType, 'nomod')
            mod_type = 'w';
        end
        % start with QC for avgMode
        if strcmpi(avgMode, 'on')
            mean_corr = mean(mastersheet.betweenSubAvgWarpQC, 'omitnan');
            sd = std(mastersheet.betweenSubAvgWarpQC, 'omitnan');
            thresh = mean_corr - 2*sd;
            if isempty(qcFiltMRI) || strcmp(qcFiltMRI, '')
                % default for avg mode when using between subject QC is the
                % mwp file
                prefix = [mod_type seg_match];
            else
                prefix = [mod_type, char(qcFiltMRI)];
            end
            if sum(contains(mastersheet.Properties.VariableNames, 'summarizedAvgFails')) > 0
                % use the previous column as our qc vec
                qc_vec = mastersheet.summarizedAvgFails;
            else
                qc_vec = nan(height(mastersheet),1);
            end
            for i=1:height(mastersheet)
                % we're only working with the midpoint average which is all
                % saved to order 1
                if isequal(mastersheet.Order(i), 1) && isequal(mastersheet.TimepointReference(i), 1)
                    if isnan(mastersheet.avg_CAT12_IQR(i)) && isnan(mastersheet.withinSubAvgWarpQC(i)) && isnan(mastersheet.betweenSubAvgWarpQC(i)) 
                        % they just haven't been run. stay as nan
                    else
                        if isequal(mastersheet.avg_IQR_fails(i),1)
                            qc_vec(i) =1;
                        else
                            % they passed the initial seg QC
                            if mastersheet.withinSubAvgWarpQC(i) < pars.Results.withinSubThresh
                                qc_vec(i) =1;
                            else
                                % finally check their between sub qc
                                if mastersheet.betweenSubAvgWarpQC(i) >= thresh
                                    % they pass on all accounts
                                    qc_vec(i) =0;
                                elseif isnan(mastersheet.betweenSubAvgWarpQC(i))
                                    qc_vec(i) = NaN;
                                else 
                                    % display their slice so we can see if it actually
                                    % looks bad. people 2 std's below the mean corr may
                                    % still have usable data if the overall sample
                                    % quality is high
                                    check_dir = fullfile(rootPath, mastersheet.RID{i}, char(mastersheet.Date(i)), char(mastersheet.Modality(i)), char(mastersheet.ScanType(i)), char(mastersheet.ScanName(i)));
                                    if isempty(warpPath) || strcmp(warpPath, '')
                                        % default location in avgMode would be in the
                                        % midpoint average mri directory
                                        check_dir = [check_dir, filesep, 'midpoint_average', filesep, 'mri'];
                                    else
                                        check_dir = [check_dir, filesep, char(warpPath)];
                                    end                     
                                    to_check = get_files(check_dir, [prefix, '*.nii']);
                                    % display the image 
                                    mbatch = mbatch +1;
                                    disp([mastersheet.RID{i} ' has relatively low covariance with other subjects.'])
                                    disp('verify that you would like the previous image to be excluded from further analysis')
                                    matlabbatch{mbatch}.spm.util.disp.data = {to_check};
                                    spm_jobman('run',matlabbatch(mbatch))
                                    pause
                                    close
                                    user_check = input('enter 1 to EXCLUDE the subject and 0 to INCLUDE the subject \n');
                                    if user_check
                                        qc_vec(i) =1;
                                    else
                                        qc_vec(i) =0;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            % write the variable to mastersheet
            mastersheet.summarizedAvgFails = qc_vec;
        end
        if ~strcmpi(indMode{1}, 'off')
            mean_corr = mean(mastersheet.betweenSubTimeWarpQC, 'omitnan');
            sd = std(mastersheet.betweenSubTimeWarpQC, 'omitnan');
            thresh = mean_corr - 2*sd;
            if isempty(qcFiltMRI) || strcmp(qcFiltMRI, '')
                prefix = [mod_type seg_match];
            else
                prefix = [mod_type, char(qcFiltMRI)];
            end
            if sum(contains(mastersheet.Properties.VariableNames, 'summarizedTimepointFails')) > 0
                % use the previous column as our qc vec
                qc_vec = mastersheet.summarizedTimepointFails;
            else
                qc_vec = nan(height(mastersheet),1);
            end
            for i=1:height(mastersheet)
                if isequal(mastersheet.timepoint_fails(i),1)
                    qc_vec(i) =1;
                else
                    % they passed the initial seg QC
                    % first ensure data is longitudinal
                    if isnumeric(indMode{1}) && isequal(length(indMode{1}), 1)
                        disp('you have only selected one timepoint of data to analyze. inclusion of within subject QC in summaryQC will be skipped')
                        disp('if you would like to perform this step, update indMode');
                    elseif isequal(max(ind_subset.Order),1)
                        disp('your data is not longitudinal. inclusion of within subject QC in summaryQC will be skipped')
                    elseif mastersheet.withinSubTimeWarpQC(i) < withinSubThresh
                        qc_vec(i) =1;
                    end
                    if isnan(mastersheet.betweenSubTimeWarpQC(i))
                        % that means this scan either failed QC at an
                        % early stage (which we've ruled out by now) or
                        % this scan is not the reference scan used for
                        % this subject in between sub QC. find the ref
                        % scan
                        sub = mastersheet(strcmp(mastersheet.RID, mastersheet.RID{i}), :);
                        if isnumeric(indMode{1})
                            match = sub(sub.Order == indMode{1} & sub.TimepointReference ==1, :);
                        else
                            % order 1 would be used
                            match = sub(sub.Order == 1 & sub.TimepointReference ==1, :);
                        end
                        if match.betweenSubTimeWarpQC >= thresh
                            qc_vec(i) = 0;
                        else
                            qc_vec(i) = 1;
                        end
                    elseif mastersheet.betweenSubTimeWarpQC(i) >= thresh
                        qc_vec(i) = 0;
                    else
                        % display their slice so we can see if it actually
                        % looks bad. people 2 std's below the mean corr may
                        % still have usable data if the overall sample
                        % quality is high
                        check_dir = fullfile(rootPath, mastersheet.RID{i}, char(mastersheet.Date(i)), char(mastersheet.Modality(i)), char(mastersheet.ScanType(i)), char(mastersheet.ScanName(i)));
                        if isempty(warpPath) || strcmp(warpPath, '')
                            % default location in avgMode would be in the
                            % midpoint average mri directory
                            check_dir = [check_dir, 'mri'];
                        else
                            check_dir = [check_dir, filesep, char(warpPath)];
                        end
                        to_check = get_files(check_dir, [prefix, '*.nii']);
                        % display the image
                        mbatch = mbatch +1;
                        matlabbatch{mbatch}.spm.util.disp.data = {to_check};
                        spm_jobman('run',matlabbatch(mbatch))
                        pause
                        close
                        disp([mastersheet.RID{i} ' has relatively low covariance with other subjects.'])
                        disp('verify that you would like the previous image to be excluded from further analysis')
                        user_check = input('enter 1 to EXCLUDE the subject and 0 to INCLUDE the subject \n');
                        if user_check
                            qc_vec(i) =1;
                        else
                            qc_vec(i) =0;
                        end
                    end    
                end
            end
            % write to mastersheet
        mastersheet.summarizedTimepointFails = qc_vec;
     end
     writetable(mastersheet, [sheet_path, filesep, study_name, '_MRI_mastertable.csv']);
    end
end
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
function folder_names = get_folders(rootPath)
allFiles = dir(rootPath);
%extract only those that are directories.
allDirFlags = [allFiles.isdir];
allFolders = allFiles(allDirFlags);
allFolders = allFolders(~ismember({allFolders(:).name},{'.','..'}));
folder_names = {allFolders.name};
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