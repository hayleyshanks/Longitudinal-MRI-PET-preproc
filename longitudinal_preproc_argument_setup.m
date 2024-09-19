function [N, study_path, mri_to_run, pet_to_run, mri_ref, frames,  study_name, qc_thresh, avgMode, ...
    indMode, templatePath, CAT12_shoot_path, coreg_image,  tissue_compartments, field_str, ...
    shootTemplateOverride, warpType, warpPath, warpFilt, ref_mask, avgPrefix, pet_orders, ...
    withinSubThresh, qcFiltMRI, pvc_opt, roi_dir, pet_QC_prefix] = longitudinal_preproc_argument_setup()
%% ------------- Argument set-up for longitudinal MRI and PET pre-processing ----------------------
% please read through the instructions and below and adjust the variables as needed to suit your study

% add path to spm - paste in your spm path below
spmpath =  'C:\Users\SchmitzLab\Documents\spm12';
addpath(spmpath);
% download the dicom to nifti conversion tool from the link below. Then,
% set the path to dicm2nii below. 
% available here: https://www.mathworks.com/matlabcentral/fileexchange/42997-xiangruili-dicm2nii
addpath('D:\HS1_APOE_Tau\toolboxes\xiangruili-dicm2nii-3fe1a27');
% set the path to wherever you've stored this code. 
functionRoot = 'D:\HS1_APOE_Tau\code\MATLAB_code'; 
addpath(genpath(functionRoot));
% specify the path to your data. The raw MRI/PET data should be in a
% subfolder of the study path (see instructions for more details)
study_path = 'D:\HS1_APOE_Tau';
% do not change this unless you did not place the enhanced TPMs in your SPM
% TPM folder as described in the instructions
templatePath = [spmpath, filesep, 'tpm'];
% path to the CAT12 MNI geodesic shoot template. do not change. 
CAT12_shoot_path = fullfile(spmpath, 'toolbox', 'cat12', 'templates_MNI152NLin2009cAsym');
% specify an example image to be used for MRI coregistration. can use a
% template image from SPM if needed. 
coreg_image = 'D:\HS2_risk_factors\SPREAD\BME9709\ADNI_BME9709_averages/ADNI_BME9709_T1_wmavg_1-397.nii'; 
%% ------------------  MRI-specific variables -----------------------------------------------------------
% MRI data you are analyzing (cell).
% e.g. {'T1'}, {'T1', 'T2_FLAIR', 'PD'}
mri_to_run = {'T1'};
% either {'on'} or {'off'}. Set to {'on'} if you want to run longitduinal
% analyses using a midpoint average (Ashburner & Ridgeway 2013)
avgMode = {'on'};
% option to run preprocessing on indiviudal timepoints (e.g. for cross
% sectional data)
% 1x1 cell or 1x2 cell array. 
% first element of the cell is either 'off', 'all_times', or a double array
% of timepoints to be run. 
% e.g. {'off'} to skip; {[1,3,5]} to run processing on MR timepoints 1,3,5.  
% second element of cell array (optional) specifies whether subjects with
% one timepoint of data only should be included in analyses. default
% behaviour will exclude subjects with only cross sectional data. 
% second element is either 'includeCrossSectionalSubjects' or
% 'excludeCrossSectionalSubjects'
% example: {'all_times', 'excludeCrossSectionalSubjects'} will analyze all
% timepoints of data but skip subjects who do not have longitudinal data
indMode = {'off'}; 
% abbreviated study name to be used in template naming etc. (char). 
study_name = 'APOE_Tau';
% numeric threshold for CAT12's image quality rating (IQR). Images with an IQR
% below the threshold are failed and will be checked manually at QC step (double)
qc_thresh = 70;
% determine how any subjects exist in the study specific raw data folder.
% This assumes all subjects have their raw data stored in a subdirectory of
% the study path called 'first_level'
num_subs = length(get_folders([study_path, filesep, 'first_level'])); % get the total number of subjects
RIDs = get_folders([study_path, filesep, 'first_level']);
% N can be over-ridden, but default behaviour will determine how many
% subjects exist and run them all. (double). 
N = 1:num_subs;
% if performing longitudinal modulation or warping, which tissue
% compartment would you like to use? options:
% {'gm'}, {'wm'}, {'csf'} --> note you can't run multiple tissue
% compartments at once.
tissue_compartments = {'gm'};
% What field strength of MRI data do you want to include. Only one field
% strength is accepted (double).
field_str = 3;
% if you are creating a custom geodesic shooting template for tis study,
% leave the variable below empty. if you have an existing shoot template,
% specifiy the path to the folder containing the template files (Templates
% 0-4 as 4D files should be provided)
shootTemplateOverride = ''; 
% either {'mod'} or {'nomod'} to warp to a template space with modulation, or without
% modulation, respectively (cell). 
warpType = {'mod'}; 
% if empty ({''}), the default directory is the subject's scan directory (avgMode)
% or the scan_directory/mri (indMode). Cell. 
warpPath = {''}; 
% file prefix to be selected for warping. Do not include asterisk and file
% extension (e.g. {'mavg'} is acceptable. If empty, ({''}) and in avgMode, defaults to searching for the
% jacobian modulated segment specified by tissue_compartments. If empty and
% in IndMode, defaults to searching for the tissue_compartment. Cell
warpFilt = {''}; 
% prefix of the files you want to use when creating a study-specific
% average of all subject's data. Do not include asterisks. Note: the files should be in 
% warped space. Char. 
avgPrefix = 's6APC_maxInterval';
% minimum within subject covariance to pass withinSubQCWarp
withinSubThresh =0.9;
% image prefix to perform QC on. Leave it blank and will decide which files
% to use based on avgMode and warpType
qcFiltMRI = {''}; 
%% ------------------------- PET-specific variables ---------------------------------------------
% the MRI scan to be used as the anatomical image for pet analyses for coregistration and ROI labelling (1x2 cell)
% first cell is the timepoint of data to use (either 'avg' for the mri
% midpoint average or 'base' to use the first timepoint of mri data) and
% the second cell is the type of mri to use (e.g. 'T1', 'T2')
% example: {'avg', 'T1'} will use the T1 midpoint average
mri_ref = {'avg', 'T1'};
% what time point of PET data should be analyzed? either {'all_times'} or a cell
% array of a double containing orders to be run e.g. {[1,3,5]}. If
% pet_orders is a 1x2 cell array, specify whether to exclude cross
% sectional subjects or not -> 'includeCrossSectionalSubjects' or
% 'excludeCrossSectionalSubjects'. Cross sectional subjects are included by
% default
pet_orders = {'all_times', 'includeCrossSectionalSubjects'}; 
% pet tracers to run. Cell array. Mutliple tracers can be included 
% e.g. {'AV45', 'PIB', 'FDG'} or {'FEOBV'}. Cell
pet_to_run = {'AV1451'};
% are your pet frames dynamic (4D) or static (3D)? cell. 
% {'dynamic'} or {'static'}. If {''}, will assume static frames
frames = {''};
% name of the mask that will be used to create your PET reference region.
% Cell. This needs to be created separately by the user. 
ref_mask = {'inferior_gm_p0_7_wCAT12_whole_cerebellum.nii'};
% path to the folder containing the reference region mask above
roi_dir = 'D:\HS1_APOE_Tau\ROIs\CAT';
% will you run pvc - 'off' or 'on'
pvc_opt = 'off'; 
% the prefix of the PET images you would like to run QC on. Usually will be
% 'suvr_w' (SUVR images that are warped to template space)
pet_QC_prefix = 'suvr_w';
end
% get folders function (get list of all folders in directory). do not edit.
function folder_names = get_folders(rootPath)
    allFiles = dir(rootPath);
    %extract only those that are directories.
    allDirFlags = [allFiles.isdir];
    allFolders = allFiles(allDirFlags);
    allFolders = allFolders(~ismember({allFolders(:).name},{'.','..'}));
    folder_names = {allFolders.name};
end



