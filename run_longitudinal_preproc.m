%% to call this function, you need to have updated longitudinal_preproc_argument_setup.m
function run_longitudinal_preproc(subs, image_modality, step_to_run)
%% ----- function to run longitudinal pre-processing for MRI and PET data -----
%% input arguments:
% subs --> numeric array of subject indices to run. e.g. 1:100, or [1,3,5,67,103]. 
% if subs = 0, function will default to running every subject. 
% modality --> character (case insensitive); either 'MRI' or 'PET'
% step_to_run --> character (case sensitive); the step of pre-processing to run
% options for step_to_run are:
    %% ---------------------------------- General ----------------------------------------
    %     createSheet
    %% ---------------------------------- MRI --------------------------------------------
    %     dicomRecon, reorientAC, checkVoxels, coreg, longitudinalRegistration, segment
    %     buildCustomTemplate, addToShoot, invertY, multiplyJacobian, warpToShoot, 
    %     createAvg, withinSubQC, betweenSubQC, extractICV, sumarizeQC.
    %% ---------------------------------- PET ----------------------------------------------
    %     dicomRecon, reorientAC, remove4D, realignFrames, coregMR,
    %     warpToShoot, makeMasks, computeSUVR, pvc, createAvg,
    %     getSUVRfails, withinSubQC, betweenSubQC, sumarizeQC, checkReg.
%% full descriptions of each of the steps above are in the instruction manual
% -------------------------------------------------------------------------------------------------------  
% call main argument set-up function to define input arguments to the various MRI/PET preproc functions
[N, study_path, mri_to_run, pet_to_run, mri_ref, frames,  study_name, qc_thresh, avgMode, ...
    indMode, templatePath, CAT12_shoot_path, coreg_image,  tissue_compartments, ...
    field_str, shootTemplateOverride, warpType, warpPath, warpFilt, ref_mask, avgPrefix, pet_orders, ...
    withinSubThresh, qcFiltMRI, pvc_opt, roi_dir, pet_QC_prefix] = longitudinal_preproc_argument_setup();
% re-define input as uppercase to make calls case-insensitive.
image_modality = upper(image_modality);
% set subs to N (all subjects in directory) if user defined subs as 0
if subs == 0
    subs = N;
end
% call the appropriate preproc function depending on the step specified by the user
if strcmpi(step_to_run, 'createSheet')
    longitudinal_preproc_step0_v1(subs, study_path, study_name, field_str)
end
% ----------------------------------- MRI preproc  -----------------------------------------
if strcmpi(image_modality, 'MRI')
    if ismember(step_to_run, {'dicomRecon', 'reorientAC', 'checkVoxels', 'coreg', 'longitudinalRegistration', 'segment'})
        longitudinal_MRI_preproc_step1_indiviudal_preproc(subs, study_path, mri_to_run, study_name, avgMode, indMode, ... 
        field_str, templatePath, CAT12_shoot_path,  coreg_image, step_to_run, 1)
    elseif ismember(step_to_run, {'buildCustomTemplate', 'addToShoot', 'invertY'})
        if ~isempty(shootTemplateOverride)
            longitudinal_MRI_preproc_step3_createTemplates(subs, study_path, mri_to_run, study_name, field_str, avgMode, indMode,step_to_run, 1, 'shootTemplateOverride', shootTemplateOverride)
        else
            longitudinal_MRI_preproc_step3_createTemplates(subs, study_path, mri_to_run, study_name, field_str, avgMode, indMode,step_to_run, 1)
        end
    elseif ismember(step_to_run, {'multiplyJacobian', 'warpToShoot', 'createAvg'})
         if ~isempty(shootTemplateOverride)
              longitudinal_MRI_preproc_step4_modulation(subs, study_path, mri_to_run, study_name, field_str, avgMode, indMode, tissue_compartments, warpType, warpPath, warpFilt, avgPrefix, step_to_run,1, 'shootTemplateOverride', shootTemplateOverride)
         else
               longitudinal_MRI_preproc_step4_modulation(subs, study_path, mri_to_run, study_name, field_str, avgMode, indMode, tissue_compartments, warpType, warpPath, warpFilt, avgPrefix, step_to_run,1)
         end
    elseif ismember(step_to_run, {'withinSubQC', 'betweenSubQC', 'extractICV', 'summarizeQC'})
        longitudinal_MRI_preproc_step5_finalQC(subs, study_path, mri_to_run, study_name, field_str, avgMode, indMode, tissue_compartments, warpType, warpPath, qcFiltMRI, withinSubThresh, step_to_run, 1)
    end
% ----------------------------------- PET preproc  -----------------------------------------
elseif strcmpi(image_modality, 'PET')
    if ismember(step_to_run, {'dicomRecon', 'reorientAC', 'remove4D','realignFrames', 'coregMR'})
        longitudinal_PET_preproc_step1(subs, study_path, pet_to_run, pet_orders, mri_ref, study_name, frames, field_str, step_to_run, 1)
    elseif ismember(step_to_run, {'warpToShoot', 'makeMasks', 'computeSUVR', 'pvc', 'createAvg'})
        if ~isempty(shootTemplateOverride)
            longitudinal_PET_preproc_step2(subs, study_path, pet_to_run, pet_orders, mri_ref, study_name, frames, roi_dir, ref_mask, avgPrefix,field_str, pvc_opt, step_to_run,1, 'shootTemplateOverride', shootTemplateOverride)
        else
             longitudinal_PET_preproc_step2(subs, study_path, pet_to_run, pet_orders, mri_ref, study_name, frames, roi_dir, ref_mask, avgPrefix,field_str, pvc_opt, step_to_run,1)
        end
    elseif ismember(step_to_run, {'getSUVRfails', 'withinSubQC', 'betweenSubQC', 'summarizeQC', 'checkReg'})
        longitudinal_PET_preproc_step3_QC(subs, study_path, pet_to_run, pet_orders,study_name, frames, mri_ref, field_str, withinSubThresh, pet_QC_prefix, step_to_run,1)
    end
else
    disp('invalid image modality. Currently, code is only set up for MRI or PET data'); 
end
end