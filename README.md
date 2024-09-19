# Longitudinal MRI/PET Pre-processing Code
Code by Hayley R.C. Shanks, with edits from Taylor W. Schmitz, Kate M. Onuska, Western University <br /> 
This repository contains code to perform longitudinal MRI and PET pre-processing as used in [Shanks et al. (2024), Nature Medicine](https://www.nature.com/articles/s41591-024-02977-w) and [German-Castelan et al. (2023), Alzheimer’s and Dementia](https://alz-journals.onlinelibrary.wiley.com/doi/full/10.1002/alz.13481). 
Code takes raw DICOM or NIfTI images and performs the following pre-processing steps: <br /> 
## MRI:
- DICOM reconstruction (optional)
- manual reorientation
- co-registration to template image
- [longitudinal registration](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3564017/)
- segmentation
- segmentation quality control (QC)
- [custom population-specific template creation with geodesic shooting](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3221052/) (optional)
- longitudinal modulation
- spatial normalization to template
- population average creation
- spatial normalization QC <br />
## PET
- DICOM reconstruction (optional)
- manual reorientation
- co-registration to subject’s anatomical image
- spatial normalization
- SUVr computation with a user-supplied reference mask
- Partial volume correction with the Müller-Gärtner (MG) method (optional) <br />
### Notes
The code was written in MATLAB r2021a using SPM12 v7771 and CAT12.9 (r2577). PET pre-processing relies on MRI data for spatial normalization. PET code can handle multi-frame “dynamic” PET data obtained during the static period, but is not made to perform true kinetic analyses. <br /> This code is meant to help users new to longitudinal MRI/PET analysis get started, and comes with no warranty. It was created with longitudinal T1 MRI data, and has undergone limited testing on cross sectional MRI or other MRI weightings. 
The code is essentially a series of wrappers to automate running a combined version of SPM12 and CAT12’s preprocessing. Brief descriptions of each step are included in the instruction manual. Further details can be found in the [SPM manual](https://www.fil.ion.ucl.ac.uk/spm/doc/spm12_manual.pdf) and [CAT12 manual](https://neuro-jena.github.io/cat12-help/). This code assumes basic knowledge of imaging principles.
## Requirements and Set-up
Download the longitudinal preprocessing code, and unzip the folder. In MATLAB, add it to your path as addpath(genpath(‘path_to_code’))
### Toolboxes:
* Download and install [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) and [CAT12](https://neuro-jena.github.io/cat/), following the instructions on the SPM and CAT websites
  * Once CAT12 is downloaded, set CAT12 to run automatically in expert mode. Open cat_defaults.m and set cat.extopts.expertgui = 1;
* MRI segmentation requires you've downloaded the enhanced TPM files from [Lorio et al. 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4819722/) and placed them in the SPM TPM subdirectory
  * The TPM should be named ‘enhanced_TPM.nii’
* If you will be running DICOM reconstruction, download and install [xiangruili’s dicm2nii](https://www.mathworks.com/matlabcentral/fileexchange/42997-xiangruili-dicm2nii) from MATLAB’s file exchange
* If you are performing partial volume correction, download and follow the SPM installation instructions from the [PETPVE toolbox](https://github.com/GGonEsc/petpve12)
### Data organization:
Raw imaging data should be contained in a subfolder within the main project directory titled ‘first_level’. Within the first level folder, the subfolder organization should be: <br /> **subject --> date --> imaging modality (MRI or PET) --> scan type (e.g. T1, FDG) --> scan name (e.g. s001)** <br /> The raw MRI or PET images should be within the scan directory. The midpoint average subdirectory will be created during pre-processing. 
![image](https://github.com/user-attachments/assets/d944cd9b-8f60-4b9b-bb2f-41a89f0e6245)
### Set input arguments in longitudinal_preproc_argument_setup.m file
Most of the functions across pre-processing will use similar input arguments. Follow the instructions in the longitudinal_preproc_argument_setup.m file to set the appropriate input arguments for your study. 
### Run pre-processing steps with run_longitudinal_preproc function 
Once your arguments have been set in the step above, you can run any step of preprocessing by calling run_longitudinal_preproc:
```
run_longitudinal_preproc(subs, modality, step_to_run)
```
Detailed explanations of each of the pre-processing steps that can be run are in **step_descriptions.pdf**
