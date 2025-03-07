%% MRI Preproc Step 2: MRI segmentation QC
% code by Hayley R. C. Shanks, Taylor W. Schmitz, Kate M. Onuska
% Note:
% -1 means data is cross sectional, but was run to exclude cross sectional subjects (or in avgMode)
% -10 means segmentation was missing
% -100 means seg failed 
function longitudinal_MRI_preproc_step2_segQC(N, study_path, mri_to_run, study_name , qc_thresh, avgMode, indMode, field_str, templatePath, CAT12_shoot_path, varargin)
rootPath = [study_path, filesep, 'first_level'];
sheet_path = [study_path, filesep, 'spreadsheets'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF USER INPUTS %%%%%%%%%%%%%%%%%%%%%%%%
% parse user inputs
pars = inputParser;
% by default, all modules are set to off.
default_run =0;
% required inputs
addRequired(pars, 'N', @isnumeric);
addRequired(pars, 'study_path', @ischar);
addRequired(pars, 'mri_to_run', @iscell);
addRequired(pars, 'study_name', @ischar);
addRequired(pars, 'qc_thresh', @isnumeric);
addRequired(pars, 'avgMode', @iscell);
addRequired(pars, 'indMode', @iscell);
addRequired(pars, 'field_str', @isnumeric);
addRequired(pars, 'templatePath', @ischar);
addRequired(pars, 'CAT12_shoot_path', @ischar);
%add optional inputs. 
addParameter(pars,'summarizeSegQC',default_run,@isnumeric); 
addParameter(pars,'redoSegProblems',default_run,@isnumeric);
% parse the function inputs
parse(pars, N, study_path, mri_to_run, study_name, qc_thresh, avgMode, indMode,  field_str, templatePath, CAT12_shoot_path, varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check that the mode variables are set correctly 
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
%get a list of all files and folders in the root path folder.
RIDs = get_folders(rootPath);
mbatch = 0;
% read in the master info sheet
mastersheet = readtable(fullfile(sheet_path, [study_name, '_MRI_mastertable.csv']));
sub_to_run = [1:length(RIDs)];
if pars.Results.summarizeSegQC
    if strcmpi(avgMode, 'on')
        for a=1:length(mri_to_run)
            % add columns to the spreadsheet which hold the IQR and a
            % binary pass/fail column. First check if these already exist.
            % If so, do not overwrite the whole column
            names_to_check = {'avg_CAT12_IQR', 'avg_IQR_fails'}; 
            if isequal(sum(contains(mastersheet.Properties.VariableNames, names_to_check)), 0)
                mastersheet.avg_CAT12_IQR = nan(height(mastersheet),1);
                mastersheet.avg_IQR_fails = zeros(height(mastersheet),1);
            end
            for subI=sub_to_run(N)
                sub_dir = fullfile(rootPath,RIDs{subI});
                sub_table = mastersheet(strcmp(mastersheet.RID, RIDs(subI)), :);
                if isempty(sub_table)
                    disp(['no data of any kind for ' RIDs{subI}])
                else
                    table_to_run = spread_table_filter(sub_table, mri_to_run{a}, field_str);
                    if isempty(table_to_run)
                        disp(['no ' char(mri_to_run) ' data for ' char(RIDs{subI})])
                    else
                        bl_dir = fullfile(sub_dir, char(table_to_run.Date(1)), char(table_to_run.Modality(1)), char(table_to_run.ScanType(1)),char(table_to_run.ScanName(1)));
                        % find their bl directory in the mastersheet
                        check = [strcmp(mastersheet.RID, RIDs{subI}), mastersheet.Date == table_to_run.Date(1), strcmp(mastersheet.Modality,table_to_run.Modality(1)), strcmp(mastersheet.ScanType, table_to_run.ScanType(1)), strcmp(mastersheet.ScanName, table_to_run.ScanName(1))] ;
                        % whichever row of check is all 1's is our row
                        check_sum = sum(check, 2);
                        row_num = find(check_sum == size(check,2));
                        % define the directory where the cat12 report will live 
                        if isequal(length(unique(table_to_run.Date)),1) 
                            % cross sectional data, so they will not have
                            % an average 
                            mastersheet.avg_CAT12_IQR(row_num) = -1;
                            % mark their data as unusable for longitudinal
                            % avg calculations
                            mastersheet.avg_IQR_fails(row_num) = 1;
                        else
                            if length(unique(table_to_run.Date)) > 1
                                % they should have an average 
                                seg_dir = [bl_dir, filesep, 'midpoint_average', filesep, 'report'];
                            end
                            % cat12 outputs a .mat file with useful info, as well as some
                            % text files. The PDF also has a good summary of all this info
                            mat = get_files(seg_dir, 'cat_*.mat');
                            if isempty(mat)
                                % should have data, but their segmentation
                                % may have been accidentally missed
                                mastersheet.avg_CAT12_IQR(row_num) = -10;
                                mastersheet.avg_IQR_fails(row_num) = 1;
                            else
                                mat = load(mat);
                                mat = mat.S;
                                if isfield(mat, 'error')
                                    % this means segmentation failed for this subject. save
                                    % this out
                                    mastersheet.avg_CAT12_IQR(row_num) = -100;
                                    mastersheet.avg_IQR_fails(row_num) = 1;
                                else
                                    % get the quality ratings from the catlog
                                    %% image quality rating
                                    ind = contains(mat.catlog, 'Image Quality Rating');
                                    iqr = mat.catlog{ind};
                                    % strip out the number only
                                    iqr = extractBetween(iqr, ':', '%'); % this could get messed up if cat12 changes their output
                                    mastersheet.avg_CAT12_IQR(row_num) = str2num(cell2mat(iqr));                                   
                                    if str2num(cell2mat(iqr)) < qc_thresh
                                        mastersheet.avg_IQR_fails(row_num) = 1; 
                                    else
                                        mastersheet.avg_IQR_fails(row_num) = 0; 
                                    end
                                    %% could also extract more information from the mat file. The quality ratings or quality measures 
                                    %% look useful, but unsure how they are scored (and why the same variable in both has a different score. Links below are unclear 
                                    % but relevant
                                    % CAT12 QC: http://www.neuro.uni-jena.de/cat12-html/cat_methods_QA.html
                                    % JISC mail question: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=SPM;52923909.1710
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    if ~strcmpi(indMode{1}, 'off')
        for a=1:length(mri_to_run)
            names_to_check = {'timepoint_IQR','timepoint_fails'}; 
            if isequal(sum(contains(mastersheet.Properties.VariableNames, names_to_check)), 0)
                mastersheet.timepoint_IQR = nan(height(mastersheet),1);
                mastersheet.timepoint_fails = nan(height(mastersheet),1); % AW: Is this supposed to read 'mastersheet.timepoint_fails'?
            end
            for subI=sub_to_run(N)
                sub_dir = fullfile(rootPath,RIDs{subI});
                sub_table = get_files(sub_dir, 'scan_info*.csv');
                % this table has all of the dates and types of scans for this P
                sub_table = readtable(sub_table);
                if isempty(sub_table)
                    disp(['no data of any kind for ' RIDs{subI}])
                    % they won't have a row in mastersheet so nothing to do
                else
                    table_to_run = spread_table_filter(sub_table, mri_to_run{a}, field_str);
                    if isempty(table_to_run)
                        disp(['no ' char(mri_to_run) ' data for ' char(RIDs{subI})])
                    elseif isequal(max(table_to_run.Order),1) && strcmpi(indMode{2}, 'excludeCrossSectionalSubjects')
                            % skip them. they only have one timepoint of data and we
                            % are excluding cross sectional subjects
                            check = [strcmp(mastersheet.RID, RIDs{subI}), mastersheet.Date == table_to_run.Date(1), strcmp(mastersheet.Modality,table_to_run.Modality(1)), strcmp(mastersheet.ScanType, table_to_run.ScanType(1)), strcmp(mastersheet.ScanName, table_to_run.ScanName(1))] ;
                            % whichever row of check is all 1's is our row
                            check_sum = sum(check, 2);
                            row_num = find(check_sum == size(check,2));
                            mastersheet.timepoint_IQR(row_num) = -1;
                            mastersheet.timepoint_fails(row_num) = 1;
                    else
                        if ~strcmpi(indMode{1}, 'all_times') && isnumeric(indMode{1})
                            % we're only processing some of the times. drop
                            % times we're not interested in
                            table_to_run = table_to_run(ismember(table_to_run.Order, indMode{1}), :);
                        end
                        % now loop through each row of their sub table, get
                        % the reports for that timepoint
                        for r=1:height(table_to_run)
                            seg_dir = fullfile(sub_dir, char(table_to_run.Date(r)), char(table_to_run.Modality(r)), char(table_to_run.ScanType(r)),char(table_to_run.ScanName(r)), 'report');                      
                            % get the cat12 mat file which contains segmentation QC
                            % info
                            % get their row num from the mastersheet for
                            % this particular time
                            check = [strcmp(mastersheet.RID, RIDs{subI}), mastersheet.Date == table_to_run.Date(r), strcmp(mastersheet.Modality,table_to_run.Modality(r)), strcmp(mastersheet.ScanType, table_to_run.ScanType(r)), strcmp(mastersheet.ScanName, table_to_run.ScanName(r))] ;
                            % whichever row of check is all 1's is our row
                            check_sum = sum(check, 2);
                            row_num = find(check_sum == size(check,2));
                            mat = get_files(seg_dir, 'cat_*.mat');
                            if isempty(mat)
                                mastersheet.timepoint_IQR(row_num) = -10;
                                mastersheet.timepoint_fails(row_num) = 1;
                            else
                                mat = load(mat);
                                mat = mat.S;
                                if isfield(mat, 'error')
                                    % this means segmentation failed for this subject. save
                                    % this out
                                    mastersheet.timepoint_IQR(row_num)= -100;
                                    mastersheet.timepoint_fails(row_num) = 1;
                                else
                                    % get the quality ratings from the catlog
                                    %% image quality rating
                                    ind = contains(mat.catlog, 'Image Quality Rating');
                                    iqr = mat.catlog{ind};
                                    % strip out the number only
                                    iqr = extractBetween(iqr, ':', '%'); % this could get messed up if cat12 changes their output
                                    mastersheet.timepoint_IQR(row_num) = str2num(cell2mat(iqr)); 
                                    if str2num(cell2mat(iqr)) < qc_thresh
                                        mastersheet.timepoint_fails(row_num) = 1;
                                    else
                                        mastersheet.timepoint_fails(row_num) = 0;
                                    end
                                    %% could also extract more information from the mat file. The quality ratings or quality measures 
                                    %% look useful, but unsure how they are scored (and why the same variable in both has a different score. Links below are unclear 
                                    % but relevant
                                    % CAT12 QC: http://www.neuro.uni-jena.de/cat12-html/cat_methods_QA.html
                                    % JISC mail question: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=SPM;52923909.1710
                                end
                            end
                        end
                    end
                end
            end
        end       
    end
    writetable(mastersheet, [sheet_path, filesep, study_name, '_MRI_mastertable.csv']);
end

if pars.Results.redoSegProblems
    % re-segment problem subjects.
    if strcmpi(avgMode, 'on')
        for t =1:length(mri_to_run)
            % define the QC variable we are using from masterSheet
            qc_fails = mastersheet(mastersheet.avg_IQR_fails ==1 & mastersheet.avg_CAT12_IQR ~= -1, :);
            for row=1:height(qc_fails)
                % find their directory
                seg_dir = fullfile(rootPath, char(qc_fails.RID(row)), char(qc_fails.Date(row)), char(qc_fails.Modality(row)), char(qc_fails.ScanType(row)),char(qc_fails.ScanName(row)), 'midpoint_average');
                % define image to segment
                to_seg = get_files(seg_dir, 'avg*.nii');
                if isempty(to_seg)
                    disp(['no ' mri_to_run{t} ' midpoint average found for ', char(qc_fails.RID(row))])
                else                
                    %sometimes cat12 segmentation fails if images have negative
                    %intensities. Replace these with 0s
                    hdr = spm_vol(to_seg);
                    mat = spm_read_vols(hdr);
                    fixed = mat;
                    fixed(fixed < 0) = 0;
                    %scrub the pinfo field of the image header so that the correct
                    %scaling is recalculated given the new data
                    hdr = rmfield(hdr, 'pinfo');
                    spm_write_vol(hdr, fixed);
                    % now that the data is fixed, re-run segmentation at higher
                    % accuracy in case the sub failed QC due to an abnormal brain
                    % shape (accuracy = 0.75 instead of 0.5). 
                    %% run SPREAD segmentation routine. Make sure CAT12 default is to run in expert mode (see function below)
                    % (cat.extopts.expertgui    = 1;) in cat_default.m
                    mbatch = mbatch+1;
                    CAT12_shoot_template = get_files(CAT12_shoot_path, 'Template_0_GS.nii');
                    lorio_TPM = get_files(templatePath, 'enhanced_TPM.nii');
                    SPREAD_segmentation_routine(mbatch, to_seg, lorio_TPM, CAT12_shoot_template, 0.75);
                end
            end 
        end
    end
    if ~strcmpi(indMode{1}, 'off')
        for t =1:length(mri_to_run)
            qc_fails = mastersheet(mastersheet.timepoint_fails == 1, :);
            % if we are not analyzing cross sectional data, drop these p's
            if strcmpi('excludeCrossSectionalSubjects', indMode{2})
                qc_fails = qc_fails(qc_fails.timepoint_IQR ~= -1, :);
                % -1 in the IQR column indicates cross sectional data
            end
            % if we are only processing a certain timepoint(s), drop other
            % rows
            if isnumeric(indMode{1})
                qc_fails = qc_fails(ismember(qc_fails.Order, indMode{1}), :);
            end
            % also drop any rows that aren't a timepoint reference (meaning
            % there is an intra average for that scan date but not in that
            % scan folder
            qc_fails = qc_fails(qc_fails.TimepointReference == 1, :);
            for row=1:height(qc_fails)
                % go through each row and re-segment the data 
                seg_dir = fullfile(rootPath, char(qc_fails.RID(row)), char(qc_fails.Date(row)), char(qc_fails.Modality(row)), char(qc_fails.ScanType(row)),char(qc_fails.ScanName(row)));
                % define image to segment. Need to consider duplicates for
                % a given timepoint
                if isequal(mastersheet.IsDuplicate(row),1)
                    % should have an intra avg
                    to_seg = get_files(seg_dir, 'intra*.nii');   
                else
                    to_seg = get_files(seg_dir, 'rac*.nii');  
                end
                %sometimes cat12 segmentation fails if images have negative
                %intensities. Replace these with 0s
                hdr = spm_vol(to_seg);
                mat = spm_read_vols(hdr);
                fixed = mat;
                fixed(fixed < 0) = 0;
                %scrub the pinfo field of the image header so that the correct
                %scaling is recalculated given the new data
                hdr = rmfield(hdr, 'pinfo');
                spm_write_vol(hdr, fixed);
                % now that the data is fixed, re-run segmentation at higher
                % accuracy in case the sub failed QC due to an abnormal brain
                % shape
                %% run SPREAD segmentation routine. Make sure CAT12 default is to run in expert mode (see function below)
                % (cat.extopts.expertgui    = 1;) in cat_default.m
                mbatch = mbatch+1;
                CAT12_shoot_template = get_files(CAT12_shoot_path, 'Template_0_GS.nii');
                lorio_TPM = get_files(templatePath, 'enhanced_TPM.nii');
                SPREAD_segmentation_routine(mbatch, to_seg, lorio_TPM, CAT12_shoot_template, 0.75);      
            end
        end
    end
end
end
