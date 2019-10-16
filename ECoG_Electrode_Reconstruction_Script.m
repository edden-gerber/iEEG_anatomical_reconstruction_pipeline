
%% Set global parameters
% Run this once before running any of the sections. Should be fixed for all
% subjects. 

% Global parameters
initial_folder = '/home/hcnl/shared_folders/ECoG_Patient_Data'; % General home folder
BIS_folder = '/usr/local/BIS/bioimagesuite30_64'; % BioImage Suite scripts folder
Fs_folder = '/usr/local/freesurfer'; % FreeSurfer folder
Fs_subjects_folder = '/home/hcnl/FreeSurfer/Subjects'; % FreeSurfer subjects folder. If using a virtual machine, note 
                                                          % that making this a shared folder may not allow Fs to run (due 
                                                          % to short delays in reading and writing to these folders). 

template_brain_data_file = '/home/hcnl/shared_folders/ECoG_Patient_Data/fsaverage/brain_data.mat';
template_vertex_data_file = '/home/hcnl/shared_folders/ECoG_Patient_Data/fsaverage/vertex_data.mat';

% File and folder name formats
recon_folder_generated_images = 'Generated_Images';
recon_folder_transformations = 'Transformations';
recon_folder_mgrids = 'Mgrids';
recon_folder_matlab = 'Matlab';
recon_filename_initial_electrode_data = 'electrode_data_initial.mat';
recon_filename_projected_electrode_data = 'electrode_data_projected.mat';
recon_filename_braindata_fs = 'brain_data_freesurfer.mat';
recon_filename_braindata_suma = 'brain_data_suma.mat';
recon_filename_vertexdata_suma = 'vertex_data_suma.mat';
recon_filename_electrodedata_final = 'electrode_data.mat';
recon_filename_braindata_final = 'brain_data.mat';
recon_filename_vertexdata_final = 'vertex_data.mat';



%% Set subject parameters
% Run this once before running any of the sections. Can be edited for an
% individual subject or left as is which will invoke the GUI for
% subject-specific parameters. 

% Option 1: Set parameters in code

% Enter subject ID
subj_id = [];
% Select reconstruction home directory 
recon_folder_home = [];
% Select original MR scan
path_mr = [];
% Select original MR scan
path_ct = [];

% Option 2: use input dialogs (runs if the variables are empty)
if isempty(subj_id)
    subj_id = inputdlg('Enter subject ID: ');
    subj_id = subj_id{1};
end
if isempty(recon_folder_home)
    recon_folder_home = uigetdir(initial_folder,'Select recon folder');
end
if isempty(path_mr)
    disp('Select original MR scan file. Select ''Cancel'' if not needed.');
    [file_name_mr, directory_mr] = uigetfile({'*.nii*;*.hdr*','Nifti or hdr file (*.nii*, *.hdr*)';'*.*', 'All files (*.*)'},'Select MR file',recon_folder_home,'multiselect','off');
    path_mr = fullfile(directory_mr,file_name_mr);
end
if isempty(path_ct)
    disp('Select original CT scan file. Select ''Cancel'' if not needed.');
    [file_name_ct, directory_ct] = uigetfile({'*.nii*;*.hdr*','Nifti or hdr file (*.nii*, *.hdr*)';'*.*', 'All files (*.*)'},'Select CT file',recon_folder_home,'multiselect','off');
    path_ct = fullfile(directory_ct,file_name_ct);
end



%% Run FreeSurfer on the original MR scan
% This step creates the segmented volume and surface data for the selected
% brain, and applied surface co-registration to the fsaverage template
% brain. 
% NOTE THAT THIS STEP TAKES 10-24 HOURS

fs_shell_initialize_cmd = ['export FREESURFER_HOME=' Fs_folder '; source $FREESURFER_HOME/SetUpFreeSurfer.sh; '];
disp(['Running FreeSurfer for subject ' subj_id '...']);
system([fs_shell_initialize_cmd 'recon-all -all -i ' path_mr ' -s ' subj_id]);

% After the main recon process terminates, also run this part - it produces
% a smooth outer surface mesh based on the pial surface data (originally to
% produce a local gyrification index map, but we will use it to project
% surface electrodes to the smooth outer surface of the cortex). 
% Note: this script uses the image processing toolbox in Matlab. 
disp(['Producing Local Gyrification Index for subject ' subj_id '...']);
system([fs_shell_initialize_cmd 'recon-all -s ' subj_id ' -localGI' ]);



%% Inspect the results using freeview
% Uses FreeSurfer's freeview utility. 
% If there are errors in the reconstruction (most likely regions where Fs
% fails to identify white matter, resulting in an erroneous pial surface
% shape), use freeview to correct the errors and re-run the appropriate
% segments of the FreeSurfer analysis. Consult Fs's wiki and YouTube videos
% for assistance. 

fs_shell_initialize_cmd = ['export FREESURFER_HOME=' Fs_folder '; source $FREESURFER_HOME/SetUpFreeSurfer.sh; '];
system([fs_shell_initialize_cmd 'freeview -v "' Fs_subjects_folder '/' subj_id '/mri/T1.mgz" ' ...
    '-f "' Fs_subjects_folder '/' subj_id '/surf/lh.pial" ' '-f "' Fs_subjects_folder '/' subj_id '/surf/rh.pial"']);



%% Run SUMA to create a standardized mesh based on the FreeSurfer surfaces
% This creates a SUMA folder in the FreeSurfer subject folder. 
% AFNI's SUrface MApping tool translates the pial and inflated surface
% meshes produced by FreeSurfer to a standardized mesh. This means that for
% any two brains, a specific vertex will correspond to homologous points on
% the pial surface (based on the registration to the template brain
% performed by FreeSurfer. It is thus very simple to transfer electrode
% coordinates from one brain to another or to a common template brain (or 
% from a pial to an inflated surface),since the mesh vertices can be mapped 
% directly from surface to surface. 

curr_path = pwd; % We will need to go to the subject FreeSurfer folder 
cd(fullfile(Fs_subjects_folder,subj_id));
system(['@SUMA_Make_Spec_FS -sid ' subj_id]);
cd(curr_path);

disp('Done.');



%% Import FreeSurfer-generated images and convert to nifti, produce stripped cortex volume images
% This section imports volume images produced by FreeSurfer from the FS
% subjects folder into the recon folder, after converting them to Nifti
% format. It also generates new volume images of the left and right
% cortical hemisphere. 
% 
% The following images are created: 
% MR_normalized: This is the original MR scan, after normalization to
%                FreeSurfer's 256x256x256 volume format. 
% MR_T1: this is the same volume after intensity normalization - all voxels
%                of the same tissue type will have the same intensity.
% cortex / cortex_right / cortex_left: These are images of the left, right
%                or both cortical hemispheres after the head, skull and 
%                subcortical structures were stripped away. 

% Create folder in the FreeSurfer subject directory for temporary files
directory_fs_temp= fullfile(Fs_subjects_folder,subj_id,'mri','temp');
if ~exist(directory_fs_temp,'dir')
    mkdir(directory_fs_temp);
end
% Create a results folder for generated images if it doesn't exist
directory_output_images = fullfile(recon_folder_home,recon_folder_generated_images);
if ~exist(directory_output_images,'dir')
    mkdir(directory_output_images);
end

% Convert the resliced and normalized MR image to nifty format (using a Fs 
% function) and save it in the results folder
fs_shell_initialize_cmd = ['export FREESURFER_HOME=' Fs_folder '; source $FREESURFER_HOME/SetUpFreeSurfer.sh; '];
path_mr_norm_mgz = fullfile(Fs_subjects_folder,subj_id,'/mri/nu.mgz');
path_mr_norm = fullfile(directory_output_images,[subj_id '_MR_normalized.nii.gz']);
system([fs_shell_initialize_cmd 'mri_convert ' path_mr_norm_mgz ' ' path_mr_norm]);

% Do the same for the FreeSurfer-generated T1 image (similar to the 
% normalized image but after intensity "standardization", e.g. all white 
% matter voxels are set to 110, etc.). 
path_mr_t1_mgz = fullfile(Fs_subjects_folder,subj_id,'/mri/T1.mgz');
path_mr_t1 = fullfile(directory_output_images,[subj_id '_MR_T1.nii.gz']);
system([fs_shell_initialize_cmd 'mri_convert ' path_mr_t1_mgz ' ' path_mr_t1]);

% Generate cortex image by masking the T1 image with the cortical ribbon mask
path_ribbon = fullfile(Fs_subjects_folder, subj_id, 'mri/ribbon.mgz');
system([fs_shell_initialize_cmd 'mri_mask "' path_mr_t1_mgz '" "' path_ribbon '" "' directory_fs_temp '/cortex.mgz"']);

% Generate left/right hemisphere cortex volumes (first creating masks for them):
system([fs_shell_initialize_cmd 'mri_binarize --i "' path_ribbon '" ' ...
    '--o "' directory_fs_temp '/ribbon_mask_left.mgz" --min 2 --max 3']);
system([fs_shell_initialize_cmd 'mri_binarize --i "' path_ribbon '" ' ...
    '--o "' directory_fs_temp '/ribbon_mask_right.mgz" --min 41 --max 42']);
system([fs_shell_initialize_cmd 'mri_mask "' path_mr_t1_mgz '" ' ...
    '"' directory_fs_temp '/ribbon_mask_left.mgz" "' directory_fs_temp '/cortex_left.mgz"']);
system([fs_shell_initialize_cmd 'mri_mask "' path_mr_t1_mgz '" ' ...
    '"' directory_fs_temp '/ribbon_mask_right.mgz" "' directory_fs_temp '/cortex_right.mgz"']);

% Convert all scans to Nifti and copy them to the output directory 
system([fs_shell_initialize_cmd 'mri_convert "' directory_fs_temp '/cortex.mgz" "' directory_output_images '/' subj_id '_cortex.nii.gz"']);
system([fs_shell_initialize_cmd 'mri_convert "' directory_fs_temp '/cortex_left.mgz" "' directory_output_images '/' subj_id '_cortex_left.nii.gz"']);
system([fs_shell_initialize_cmd 'mri_convert "' directory_fs_temp '/cortex_right.mgz" "' directory_output_images '/' subj_id '_cortex_right.nii.gz"']);

disp('Finished');



%% Register the CT image to the Freesurfer-generated MR scans
% This step is required so that the electrode coordinates, identified based
% on the CT scan, will be aligned to the MR image. An alternative approach
% is to identify the coordinates on the original CT, and then apply the
% necessary transformation matrices to them to align them with the MR
% image. However, this approach is much more prone to errors, as
% transformation matrices generated by registering two images may not
% always be appropriate when used on electrode coordinates: there may be
% additional (implied) transformations which will be applied when
% transforming the image, which will not be reflected in the transformation 
% matrix. Specifically, a single image generally includes two coordinate 
% spaces - an image space and an anatomical/scan space; the transformation 
% matrix produced by the coregistration to a second image may apply to the 
% two images' anatomical space, while the electrode coordinates may be 
% coded in image space, and therefore applying the transfomation matrix to 
% it may not work correctly, depending on various factors such as the
% configuration of the scanners and the software used for image and 
% electrode registration. 
% 
% The registration is done between the original CT image and the
% MR_normalized image produced in the previous section. The output is a
% file with suffix CT_reg_to_mr, which should be CT image used in the next
% step for manually identifying electrodes with BioImage Suite. 
% 
% Occasionally the CT-MR registration fails, which should be clear when
% inspecting the results. In this case try the next "troubleshooting"
% sections. 
% 
% This step uses the registration functionality of FSL. 

% Create a results folder for generated images if it doesn't exist
directory_output_images = fullfile(recon_folder_home,recon_folder_generated_images);
if ~exist(directory_output_images,'dir')
    mkdir(directory_output_images);
end
% Create a results folder for transformation matrices if it doesn't exist
directory_output_matrices = fullfile(recon_folder_home,recon_folder_transformations);
if ~exist(directory_output_matrices,'dir')
    mkdir(directory_output_matrices);
end

% Run CT->MR co-registration 
fs_shell_initialize_cmd = ['export FREESURFER_HOME=' Fs_folder '; source $FREESURFER_HOME/SetUpFreeSurfer.sh; '];
path_mr_norm = fullfile(directory_output_images,[subj_id '_MR_normalized.nii.gz']);
path_reg_output = fullfile(directory_output_images,[subj_id '_CT_reg_to_mr.nii.gz']);
path_transformation_mat = fullfile(directory_output_matrices, [subj_id '_CT_to_MRnorm.matrix']);
CMD = ['flirt -in ' path_ct ' -ref ' path_mr_norm '  -out ' path_reg_output ' -omat ' path_transformation_mat ' -bins 256 -cost mutualinfo -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 6  -interp trilinear'];
disp('Co-registering CT to the FreeSurfer-generated normalized MR image...');
system([fs_shell_initialize_cmd CMD]);

% Remove NaN values from CT image
% This step is usually not necessary. However
% in some cases images generated by the above procedure produce an error
% when read into BioImage Suite. This can be fixed by replacing all NaN
% values in the image by a numerical value (the image's minimal value in
% this case). This does not affect the procedure further since a threshold
% value is used on the CT image in any case to visualize the electrodes. 
CT = MRIread(path_reg_output);
CT.vol(isnan(CT.vol)) = min(CT.vol(:));
MRIwrite(CT,path_reg_output);

disp('Done.');



%% View the registration results in freeview
% This will load the reference MR image and the transformed CT image as two
% layers. Use the "opacity" and "window" sliders to check that the two 
% images are indeed aligned. Note that if a pre-surgical MR scan is used, 
% some discrepancy between the images is expected due to the deformation of 
% the brain after the surgery. 
% 
% If the registration was not successful, try options detailed in the next
% two sections. 

directory_output_images = fullfile(recon_folder_home,recon_folder_generated_images);

path_mr_t1 = fullfile(directory_output_images,[subj_id '_MR_T1.nii.gz']);
path_reg_output = fullfile(directory_output_images,[subj_id '_CT_reg_to_mr.nii.gz']);

fs_shell_initialize_cmd = ['export FREESURFER_HOME=' Fs_folder '; source $FREESURFER_HOME/SetUpFreeSurfer.sh; '];
system([fs_shell_initialize_cmd 'freeview -v "' path_mr_t1 '" -v "' path_reg_output '"']);



%% Registration Troubleshooting 1: Try an inverse transform
% In certain cases while a CT->MR registration may fail, the inverse MR->CT
% registartion may succeed. In this case the resulting transformation matrix
% can simply be inverted and applied to the CT image. 
% After this section finishes running, check the results again by running 
% the previous section. 

% Create a results folder for transformation matrices if it doesn't exist
directory_output_matrices = fullfile(recon_folder_home,recon_folder_transformations);
if ~exist(directory_output_matrices,'dir')
    mkdir(directory_output_matrices);
end

% Run MR->CT co-registration 
path_transformation_mat = fullfile(directory_output_matrices, [subj_id '_MRnorm_to_CT.matrix']);
CMD = ['flirt -in ' path_mr_norm ' -ref ' path_ct ' -omat ' path_transformation_mat ' -bins 256 -cost mutualinfo -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 6  -interp trilinear'];
disp('Co-registering MR to the CT image...');
system([fs_shell_initialize_cmd CMD]);

% Invert the resulting transformation matrix
disp('Inverting transformation matrix...');
path_transformation_mat_inv = fullfile(directory_output_matrices, [subj_id '_CT_to_MRnorm.matrix']);
system(['convert_xfm -omat ' path_transformation_mat_inv ' -inverse ' path_transformation_mat]);

% Apply the inverse transform to the CT image
path_mr_norm = fullfile(directory_output_images,[subj_id '_MR_normalized.nii.gz']);
path_reg_output = fullfile(directory_output_images,[subj_id '_CT_reg_to_mr.nii.gz']);
disp('Applying transformation to the CT image...');
system(['flirt -in ' path_ct ' -applyxfm -init ' path_transformation_mat_inv ' -out ' path_reg_output ' -paddingsize 0.0 -interp trilinear -ref ' path_mr_norm]);

% Remove NaN values from CT image - see within cell "Register the CT image 
% to the Freesurfer-generated MR scans" for explanation.
CT = MRIread(path_reg_output);
CT.vol(isnan(CT.vol)) = min(CT.vol(:));
MRIwrite(CT,path_reg_output);

disp('Done.');



%% Registration Troubleshooting 2: Use the GUI
% If the coregistration still doesn't work, you're on your own - use FSL's
% GUI to try different registration methods. 
% The non-default parameters used in the code above: 
% Model: 3D Rigid Body (6 parameters model)
% Search: Images "Incorrectly Oriented"
% Cost Function: Mutual Information
% 
% When you are satisfied with the results, make sure they are saved with
% the correct file name formats - output CT image should be in the
% Generated_Images folder and be named "<subj_ID>_CT_reg_to_mr.nii"
system('Flirt');



%% Manually identify electrodes on the CT image with BioImage Suite
% Note: this step can also be done on a separate Windows machine if it has 
% BioImage Suite installed. 
%
% Instructions: 
% I. Identify electrodes on CT: 
% When the BioImage Suite menu appears, select the "Editors" tab and click 
% on "Electrode Editor".
% 1. In the viewing window, open the CT image (IMPORTANT: this should be
%    the CT image which was aligned to the MR scan).
% 2. On the upper right pane, choose "3D Only" and "Volume".
% 3. Image Processing->Threshold->Threshold Image - apply a high enough low 
%    threshold so that only the electrodes (and wires) remain.
% 4. In the electrode editor window, select Edit->Full Edit Mode, Display->
%    Show All, and Labels->Label Font <num>.
% 5. In the Patient Info tab, add as many grids as necessary.
% 6. For each grid:
%    * Click on the Electrode Info tab, and enter the grid name, dimensions 
%      and spacing. Click "Update Grid".
%    * Grid->Grid Color, choose a new color.
%    * On the bottom right panel, check "Button Pick".
%    * For each electrode, first check the electrode on the electrode info 
%      window (the first electrode of the grid is the bottom-left!
%      the electrode number is indicated in the top-right panel, starting 
%      from zero) and then shift-click on the corresponding electrode 
%      location on the image. 
%    * If some electrodes merge together, try increasing the image 
%      threshold. If they disappear, decrease the threshold. 
%    * It is recommended to enter the grids in the order in which they are 
%      listed in the recorded channel list (e.g. first the grid with 
%      electrodes 1-64, then the strip with electrodes 65-72, etc.). 
%      Although if not this can be fixed later. 
%    * If only some of the electrodes on the grid were actually recorded, 
%      it is recommended to still identify the full grid (it helps to make
%      the projection of the grid to the surface more accurate). In a later 
%      step, the unnecessary electrodes can be removed from the electrode
%      list. 
% 7. When finished, click File->Save. Save the results as 
%    <recon_folder>/Mgrids/<subject_id>_coordinates.mgrid. 
% 8. Exit BioImage Suite. 

% Create results folder
directory_output_mgrids = fullfile(recon_folder_home,recon_folder_mgrids);
if ~exist(directory_output_mgrids,'dir')
    mkdir(directory_output_mgrids);
end

system([BIS_folder '/start_bioimagesuite']);



%% Read electrode coordinates and grid data from the mgrid to Matlab 
% In this section, the mgrid file created with BioImage Suite is read into
% a data structure in Matlab using the custom read_mgrid_file function. The
% resulting structure has the following fields: 
% subject_id: as defined for the current run. 
% num_grids: number of electrode grids
% num_electrodes: total number of (recorded) electrodes
% grids: a table holding information for each grid, including: 
% - name
% - display color
% - dimensions 
% - spacing
% - depth electrodes (true/false)
% - hemisphere (right/left_
% - coordinates matrix for electrodes in the grid, in IJK and RAS
%   coordinate spaces. 
% electrodes: a table holding information for each electrode, including:
% - grid name
% - electrode name (grid name + electrode number in grid)
% - grid color
% - depth electrode (true/false)
% - hemisphere (right/left)
% - electrode coordinates, in IJK and RAS coordinate spaces. 
% - nearest vertex on the FreeSurfer/SUMA mesh. This information is added
%   in a later step. 
% coordinate_transformations: a structure holding any relevant
% transformations between coordinate spaces (e.g. IJK to RAS) in Matlab's
% "affine3d" format. 
%
% After the mgrid data is read into Matlab, two additional functions - 
% gui_reorder_electrodes and gui_mark_hemisphere_and_depth_electrodes - are
% used to fill-in any information not included in the mgrid and to set the
% correct order of electrode channels (at this stage electrodes can be
% re-ordered, non-recorded electrodes can be removed and non-electrode
% channels can be added). 

DEFAULT_MGRID_FILE_NAME = [subj_id '_coordinates.mgrid'];

% Create output folder if it doesn't exist
folder_matlab_output = fullfile(recon_folder_home,recon_folder_matlab);
if ~exist(folder_matlab_output,'dir')
    mkdir(folder_matlab_output);
end

% Select CT image based on which the electrodes were identified (for
% information on how to transform the image space IJK coordinates to scan
% space RAS coordinates)
path_ct_for_elec = fullfile(recon_folder_home,recon_folder_generated_images,[subj_id '_CT_reg_to_mr.nii.gz']);
if ~exist(path_ct_for_elec,'file');
    disp(['CT image file ' path_ct_for_elec ' not found. Please select the CT image file used to identify electrode coordinates.']);
    [file_name_ct_for_elec, directory_ct_for_elec] = uigetfile({'*.nii*;*.hdr*','Nifti or hdr file (*.nii*, *.hdr*)';'*.*', 'All files (*.*)'},'Select CT file (on which electrodes were identified)',recon_folder_home,'multiselect','off');
    path_ct_for_elec = fullfile(directory_ct_for_elec,file_name_ct_for_elec);
end

% Read mgrid file
path_mgrid = fullfile(recon_folder_home,recon_folder_mgrids,DEFAULT_MGRID_FILE_NAME);
if ~exist(path_mgrid,'file')
    disp(['mgrid file ' path_mgrid ' not found. Please select an mgrid file.']);
    [file_name_mgrid, directory_mgrid] = uigetfile({'*.mgrid*','mgrid file (*.mgrid*)';'*.*', 'All files (*.*)'},'Select mgrid file',[recon_folder_home '/Mgrids'],'multiselect','off');
    path_mgrid = [directory_mgrid file_name_mgrid];
end
electrode_data = read_mgrid_file(subj_id, path_mgrid, path_ct_for_elec);

% Re-order the electrodes so that they correspond to the order of recorded
% channels in the data set. Add or remove electrodes as needed (added
% electrodes will be marked as "External" in the "grid" column). 
disp(['Arrange the electrodes according to the desired order for your data set. You can change electrode order, remove electrodes '...
    'which were not recorded or add additional channels not corresponding to any of the electrode grids.']);
electrode_data = gui_reorder_electrodes(electrode_data);

% Query whether there are any depth electrodes (these will not be
% projected to the surface in the following step), and which hemisphere
% each each grid is mapped to:
disp(['Mark grids as containing either surface or depth electrodes, and identify the cortical hemisphere to which they are mapped. ', ...
    'Initial hemisphere values are suggested based on the mean electrode coordinates of the grid.']);
electrode_data = gui_mark_hemisphere_and_depth_electrodes(electrode_data);

% Save electrode data
path_electrode_data = fullfile(folder_matlab_output,[subj_id '_' recon_filename_initial_electrode_data]);
disp(['Saving electrode information to ' path_electrode_data]);
% save(path_electrode_data,'electrode_data');



%% Read pial, inflated and smoothed outer cortical surfaces from FreeSurfer format to Matlab
% In this section, the pial and inflated surfaces generated by FreeSurfer
% are read into a Matlab data structure (using FS's read_surf function).
% The surfaces are transformed from the native FS coordinate space to the
% standard RAS anatomical space. The result is a data structure with the
% following fields:
% - pial_left/pial_right: pial (non-inflated) cortical surfaces
% inflated_left/inflated_right: inflated cortical surfaces
% smooth_outer_left/smooth_outer_right: the smoothed "envelope" surface
% around the pial surface. Used as a target for subsequent projection of
% electrodes since they are expected to be on the outer surface (i.e. not
% inside a sulcus). 
% 
% Note that since the inflated surfaces are both initially centered around
% the origin, they would overlap if plotted together. Therefore in an
% adidtional step the two inflated hemispheres are pulled apart (this has
% no impact on plotting electrodes on them since for inflated surfaces, the 
% electrodes are mapped to a specific mesh vertex and not to an xyz 
% coordinate. 
% 
% This section is partially based on code by Tal Golan @ Rafael Malach lab.

% Create output folder if it doesn't exist
folder_matlab_output = fullfile(recon_folder_home,recon_folder_matlab);
if ~exist(folder_matlab_output,'dir')
    mkdir(folder_matlab_output);
end

% Get transformation matrices: 
disp('Getting transformation matrices...');
fs_shell_initialize_cmd = ['export FREESURFER_HOME=' Fs_folder '; source $FREESURFER_HOME/SetUpFreeSurfer.sh; '];
[~, cmdout] = system([fs_shell_initialize_cmd 'mri_info --vox2ras ' Fs_subjects_folder '/' subj_id '/mri/orig.mgz']);
transformations.ijk2xyz = affine3d(str2num(cmdout)');

[~, cmdout] = system([fs_shell_initialize_cmd 'mri_info --vox2ras-tkr ' Fs_subjects_folder '/' subj_id '/mri/orig.mgz']);
transformations.ijk2xyz_FsMesh = affine3d(str2num(cmdout)');

brain_data = struct;
% Read pial surfaces to Matlab, transform and save: 
disp('Reading pial surfaces...');
% Right hemisphere
[brain_data.pial_right.vertices,brain_data.pial_right.faces] = read_surf(fullfile(Fs_subjects_folder,subj_id,'surf','rh.pial'));
brain_data.pial_right.faces = brain_data.pial_right.faces + 1;
brain_data.pial_right.vertices = transformations.ijk2xyz_FsMesh.transformPointsInverse(brain_data.pial_right.vertices);
brain_data.pial_right.vertices = transformations.ijk2xyz.transformPointsForward(brain_data.pial_right.vertices);
% Left hemisphere
[brain_data.pial_left.vertices,brain_data.pial_left.faces] = read_surf(fullfile(Fs_subjects_folder,subj_id,'surf','lh.pial'));
brain_data.pial_left.faces = brain_data.pial_left.faces + 1;
brain_data.pial_left.vertices = transformations.ijk2xyz_FsMesh.transformPointsInverse(brain_data.pial_left.vertices);
brain_data.pial_left.vertices = transformations.ijk2xyz.transformPointsForward(brain_data.pial_left.vertices);

% Read inflated surfaces: 
disp('Reading inflated surfaces...');
% Right hemisphere
[brain_data.inflated_right.vertices,brain_data.inflated_right.faces] = read_surf(fullfile(Fs_subjects_folder,subj_id,'surf','rh.inflated'));
brain_data.inflated_right.faces = brain_data.inflated_right.faces + 1;
brain_data.inflated_right.vertices = transformations.ijk2xyz_FsMesh.transformPointsInverse(brain_data.inflated_right.vertices);
brain_data.inflated_right.vertices = transformations.ijk2xyz.transformPointsForward(brain_data.inflated_right.vertices);
% Left hemisphere
[brain_data.inflated_left.vertices,brain_data.inflated_left.faces] = read_surf(fullfile(Fs_subjects_folder,subj_id,'surf','lh.inflated'));
brain_data.inflated_left.faces = brain_data.inflated_left.faces + 1;
brain_data.inflated_left.vertices = transformations.ijk2xyz_FsMesh.transformPointsInverse(brain_data.inflated_left.vertices);
brain_data.inflated_left.vertices = transformations.ijk2xyz.transformPointsForward(brain_data.inflated_left.vertices);
% Pull apart the two inflated hemispheres so they do not overlap if plotted together:
[brain_data.inflated_left, brain_data.inflated_right] = pull_apart_inflated_hemispheres(brain_data.inflated_left, brain_data.inflated_right, 10);


% Read smoothed outer surface mesh for left and right cortices (for
% subsequent projection of surface electrode grids to these surfaces):
disp('Reading smoothed outer surfaces...');
% Right 
surface_file_name_r = fullfile(Fs_subjects_folder,subj_id,'surf','rh.pial-outer-smoothed');
[brain_data.smooth_outer_right.vertices, brain_data.smooth_outer_right.faces] = read_surf(surface_file_name_r);
brain_data.smooth_outer_right.faces = brain_data.smooth_outer_right.faces + 1;
brain_data.smooth_outer_right.vertices = transformations.ijk2xyz_FsMesh.transformPointsInverse(brain_data.smooth_outer_right.vertices);
brain_data.smooth_outer_right.vertices = transformations.ijk2xyz.transformPointsForward(brain_data.smooth_outer_right.vertices);
% Left
surface_file_name_l = fullfile(Fs_subjects_folder,subj_id,'surf','lh.pial-outer-smoothed');
[brain_data.smooth_outer_left.vertices, brain_data.smooth_outer_left.faces] = read_surf(surface_file_name_l);
brain_data.smooth_outer_left.faces = brain_data.smooth_outer_left.faces + 1;
brain_data.smooth_outer_left.vertices = transformations.ijk2xyz_FsMesh.transformPointsInverse(brain_data.smooth_outer_left.vertices);
brain_data.smooth_outer_left.vertices = transformations.ijk2xyz.transformPointsForward(brain_data.smooth_outer_left.vertices);

% Save
output_path = fullfile(recon_folder_home,recon_folder_matlab,[subj_id '_' recon_filename_braindata_fs]); % freesurfer mesh data
disp(['Saving brain mesh data to ' output_path ]);
save(output_path,'brain_data');

disp('Done.');



%% Reduce errors in electrode localization by projecting them to the cortical surface
% Errors in electrode localization can occur in the CT->MR coregistration
% phase, especially if the MR image is pre-operation. These errors often
% result in electrodes being localized to coordinates within or outside of
% the cortical volume rather than on its surface (for surface electrodes of
% course). In this step, the electrode locations are projected to the
% smoothed corical surface, using a custom gradient-descent method which 
% iteratively moves and rotates the entire grid until the mean square
% distance of its electrodes from the surface is minimized. 
% 
% NOTE: After finishing this function an already-published similar solution
% was found. See:
% Dykstra, A.R. et al., 2012. Individualized localization and cortical 
% surface-based registration of intracranial electrodes. NeuroImage, 59(4),
% pp.3563�3570.
%
% To test the results, you can compare the old vs. the new mgrid file in 
% BioImage Suite using the appropriate cortex hemisphere image. 
%
% NOTE 2: There may be a case where a surface (non-depth) grid should not 
% be projected to the outer surface (e.g. when it is placed within a 
% sulcus). In this case you can initially set the grid to "depth" (with the 
% gui_mark_hemisphere_and_depth_electrodes function), run this section, and
% then run gui_mark_hemisphere_and_depth_electrodes again to set it to
% non-depth so it will be mapped to a pial/inflated surface vertex (run the
% code under the comment "Find and record the nearest vertex on the surface
% for each electrode" again). 

show_graphical_output = true; % Watch the projection process graphically. 

% load electrode data if it is not loaded (or a wrong subject's data is
% loaded)
if ~exist('electrode_data','var') || ~strcmp(subj_id,electrode_data.subject_id)
    % first check default path
    path_electrode_data = fullfile(recon_folder_home,recon_folder_matlab,[subj_id '_' recon_filename_initial_electrode_data]);
    if ~exist(path_electrode_data,'file')
        % if it is not found, select file 
        disp(['Electrode data file not found at ' path_electrode_data '. Select electrode data .mat file:']);
        [file_name_electrodes, directory_electrodes] = uigetfile({'*.mat*','mat file (*.mat*)';'*.*', 'All files (*.*)'},'Select electrodes file',[recon_folder_home '/Matlab'],'multiselect','off');
        path_electrode_data = [directory_electrodes file_name_electrodes];
    end
    disp(['Loading electrode data from: ' path_electrode_data ]);
    load(path_electrode_data);
else 
    disp(['Electrode data for ' subj_id ' already loaded']);
end

% load surface data if it is not loaded
if ~exist('brain_data','var')
    % first check default path
    path_brain_data = fullfile(recon_folder_home,recon_folder_matlab,[subj_id '_' recon_filename_braindata_fs]);
    if ~exist(path_brain_data,'file')
        % if it is not found, select file 
        disp(['Brain surface mesh data file not found at ' path_brain_data '. Select brain data .mat file:']);
        [file_name_brain_data, directory_brain_data] = uigetfile({'*.mat*','mat file (*.mat*)';'*.*', 'All files (*.*)'},'Select electrodes file',[recon_folder_home '/Matlab'],'multiselect','off');
        path_brain_data = [directory_brain_data file_name_brain_data];
    end
    disp(['Loading brain surface data from: ' path_brain_data ]);
    load(path_brain_data);
else 
    disp('Brain surface data already loaded');
end

% Go over each grid/strip and project it to the surface using a gradient 
% descent algorithm. At each iteration, the entire grid is both displaced a 
% fixed step size in the x-y-z space, and rotated a fixed degree step size 
% around the x, y, and z axes. The displacement and rotation vectors are 
% determined by the gradient of the error signal, computed as the sum of 
% square distances between each electrodes and the nearest vertex on the 
% surface.
% New figure (docked) 
if show_graphical_output
    h_fig = figure;
    drawnow;
    set(h_fig,'Name','Electrode Projection','WindowStyle','Docked');
end

disp('Starting electrode projection...');
for g = 1:electrode_data.num_grids % 8:9 for ST45 demo
    % get grid information
    grid_name = electrode_data.grids.name{g};
    grid_hemis = electrode_data.grids.hemisphere{g};
    depth = electrode_data.grids.depth(g);
    coordinates = electrode_data.grids.coordinates_ras{g};
    
    if depth % depth electrodes - do not project 
        disp(['Skipping projection for ' grid_name]);
    else % not depth electrodes - project
        disp(['Projecting electrodes in ' grid_name '...']);
       if strcmp(grid_hemis,'Right')
           surf = brain_data.smooth_outer_right;
       else
           surf = brain_data.smooth_outer_left;
       end
       
       % Project electrodes to the surface
        proj_electrodes = project_electrodes_to_surf(coordinates, surf, [], [], [], show_graphical_output);
        
        % Update grid data table
        electrode_data.grids.coordinates_ras{g} = proj_electrodes;
        electrode_data.grids.coordinates_ijk{g} = electrode_data.coordinate_transformations.ijk2ras.transformPointsInverse(proj_electrodes);
        % For each electrode in the grid, update it in the electrode table 
        % if it is there:
        for e = 1:size(proj_electrodes,1)
            elec_name = [grid_name '_' num2str(e)];
            elec_idx = find(strcmp(electrode_data.electrodes.name,elec_name));
            if ~isempty(elec_idx)
                electrode_data.electrodes.coordinates_ras(elec_idx,:) = electrode_data.grids.coordinates_ras{g}(e,:);
                electrode_data.electrodes.coordinates_ijk(elec_idx,:) = electrode_data.grids.coordinates_ijk{g}(e,:);
            end
        end
    end
    
    if show_graphical_output
        % Inspect the results and continue
        str = input('Press return to continue... ','s');
    end
end
disp('Projection finished');

% Find and record the nearest vertex on the surface for each electrode 
for e = 1:electrode_data.num_electrodes
    coord = electrode_data.electrodes.coordinates_ras(e,:);
    if isnan(coord(1)) || electrode_data.electrodes.depth(e)
        vert_idx(e) = nan;
    else
        if strcmp(electrode_data.electrodes.hemisphere{e},'Right');
            surf = brain_data.pial_right;
        else
            surf = brain_data.pial_left;
        end
        [vert_idx(e), err(e)] = map_electrode_to_vertex(coord,surf);
    end
end
% Add a new column for nearest mesh vertex:
new_table_column = table(vert_idx','VariableNames',{'nearest_vertex'});
% if it already exists from a previous run, remove it:
if any(strcmp(fieldnames(electrode_data.electrodes),'nearest_vertex'))
    field_idx = find(strcmp(fieldnames(electrode_data.electrodes),'nearest_vertex'));
    electrode_data.electrodes = electrode_data.electrodes(:,[1:(field_idx-1) (field_idx+1):end]);
end
electrode_data.electrodes = [electrode_data.electrodes new_table_column];


% Save results
out_file_name = fullfile(recon_folder_home,recon_folder_matlab,[subj_id '_' recon_filename_projected_electrode_data]);
disp(['Saving results to ' out_file_name]);
save(out_file_name,'electrode_data');

% Save new mgrid file
% Get file names
out_mgrid_name = fullfile(recon_folder_home,recon_folder_mgrids,[subj_id '_coordinates_projected.mgrid']);
in_mgrid_name = fullfile(recon_folder_home,recon_folder_mgrids,[subj_id '_coordinates.mgrid']);
if ~exist(in_mgrid_name,'file');
    disp(['mgrid file at ' in_mgrid_name ' not found. Please select original mgrid file']);
    [file_name_mgrid, directory_mgrid] = uigetfile({'*.mgrid','mgrid file (*.mgrid)';'*.*', 'All files (*.*)'},'Select original mgrid file (on which electrodes were identified)',directory_home,'multiselect','off');
    in_mgrid_name = fullfile(directory_mgrid,file_name_mgrid);
end
% Save mgrid
disp(['Saving new mgrid file to ' out_mgrid_name]);
edit_mgrid_file(subj_id, in_mgrid_name,out_mgrid_name, electrode_data);

% Close figure
if show_graphical_output
    close(h_fig)
end



%% Visualize 3D cortex with electrodes (freesurfer mesh) 

ELECTRODE_PLOT_COLOR = 'k';
ELECTRODE_PLOT_SIZE = 20;

% load electrode data if it is not loaded (or a wrong subject's data is loaded)
if ~exist('electrode_data','var') || ~strcmp(subj_id,electrode_data.subject_id)
    % first check default path
    path_electrode_data = fullfile(recon_folder_home,recon_folder_matlab,[subj_id '_' recon_filename_projected_electrode_data]);
    if ~exist(path_electrode_data,'file')
        % if it is not found, select file 
        disp(['Electrode data file not found at ' path_electrode_data '. Select electrode data .mat file:']);
        [file_name_electrodes, directory_electrodes] = uigetfile({'*.mat*','mat file (*.mat*)';'*.*', 'All files (*.*)'},'Select electrodes file',[recon_folder_home '/Matlab'],'multiselect','off');
        path_electrode_data = [directory_electrodes file_name_electrodes];
    end
    disp(['Loading electrode data from: ' path_electrode_data ]);
    load(path_electrode_data);
else 
    disp(['Electrode data for ' subj_id ' already loaded']);
end
% load surface data if it is not loaded
if ~exist('brain_data','var')
    % first check default path
    path_brain_data = fullfile(recon_folder_home,recon_folder_matlab,[subj_id '_' recon_filename_braindata_fs]);
    if ~exist(path_brain_data,'file')
        % if it is not found, select file 
        disp(['Brain surface mesh data file not found at ' path_brain_data '. Select brain data .mat file:']);
        [file_name_brain_data, directory_brain_data] = uigetfile({'*.mat*','mat file (*.mat*)';'*.*', 'All files (*.*)'},'Select electrodes file',[recon_folder_home '/Matlab'],'multiselect','off');
        path_brain_data = [directory_brain_data file_name_brain_data];
    end
    disp(['Loading brain surface data from: ' path_brain_data ]);
    load(path_brain_data);
else 
    disp('Brain surface data already loaded');
end


hemisphere = {'right','left'};
% Plot pial surfaces with electrodes
for h = 1:2
    if any(strcmpi(hemisphere{h},electrode_data.electrodes.hemisphere)) % Plot hemisphere only if there are any electrodes on it
        figure;
        coordinates = electrode_data.electrodes.coordinates_ras(strcmpi(hemisphere{h},electrode_data.electrodes.hemisphere),:);
        % plot brain
        plot_mesh_brain(brain_data.(['pial_' hemisphere{h}])); 
        % plot electrodes - scatter
        plot_data_on_mesh(coordinates, ELECTRODE_PLOT_COLOR, 'markersize', ELECTRODE_PLOT_SIZE); % Add a dot at each electrode location
        % plot electrodes - surface
        plot_data_on_mesh(coordinates,ones(length(coordinates),1), 'surface'); % plot electrodes as colored patches on the surface instead of dots
        % fix lighting and set title
        view([(h-1)*180+90 0]);fix_lighting; 
        title(['Subject ' subj_id ': ' hemisphere{h} ' pial surface']);
    end
end

% Plot inflated surfaces with electrodes
for h = 1:2
    if any(strcmpi(hemisphere{h},electrode_data.electrodes.hemisphere)) % Plot hemisphere only if there are any electrodes on it
        figure;
        coordinates = electrode_data.electrodes.nearest_vertex(strcmpi(hemisphere{h},electrode_data.electrodes.hemisphere),:);
        % plot brain
        plot_mesh_brain(brain_data.(['inflated_' hemisphere{h}])); 
        % plot electrodes - scatter
        plot_data_on_mesh(coordinates, ELECTRODE_PLOT_COLOR, 'markersize', ELECTRODE_PLOT_SIZE); % Add a dot at each electrode location
        % plot electrodes - surface
        plot_data_on_mesh(coordinates,ones(length(coordinates),1), 'surface'); % plot electrodes as colored patches on the surface instead of dots
        % fix lighting and set title
        view([(h-1)*180+90 0]);fix_lighting; 
        title(['Subject ' subj_id ': ' hemisphere{h} ' inflated surface']);
    end
end



%% Load the SUMA surface mesh and anatomical metadata, map electrode coordinates to standard mesh vertices
% Based on code by Tal Golan @ Rafael Malach Lab

% load electrode data if it is not loaded (or a wrong subject's data is loaded)
if ~exist('electrode_data','var') || ~strcmp(subj_id,electrode_data.subject_id)
    % first check default path
    path_electrode_data = fullfile(recon_folder_home,recon_folder_matlab,[subj_id '_' recon_filename_projected_electrode_data]);
    if ~exist(path_electrode_data,'file')
        % if it is not found, select file 
        disp(['Electrode data file not found at ' path_electrode_data '. Select electrode data .mat file:']);
        [file_name_electrodes, directory_electrodes] = uigetfile({'*.mat*','mat file (*.mat*)';'*.*', 'All files (*.*)'},'Select electrodes file',[recon_folder_home '/Matlab'],'multiselect','off');
        path_electrode_data = [directory_electrodes file_name_electrodes];
    end
    disp(['Loading electrode data from: ' path_electrode_data ]);
    load(path_electrode_data);
else 
    disp(['Electrode data for ' subj_id ' already loaded']);
end


% Get transformation matrix from the AFNI coordinate space to RAS
afni_shift1d_file = fullfile(Fs_subjects_folder,subj_id,'SUMA','brainmask_shft.1D');
[err, ras2afni_matrix] = Read_1D(afni_shift1d_file); assert(err==0); % Read_1D is an AFNI function
ras2afni_matrix = [ reshape(ras2afni_matrix,4,3) [0;0;0;1] ];
ras2afni_matrix(4,3) = -ras2afni_matrix(4,3); % flip Z
transformations.ras2afni = affine3d(ras2afni_matrix);


% Get pial and inflated meshes
% Define paths
mesh_type = 'std.141';
surface_names = {'pial_left','pial_right','inflated_left','inflated_right'};
map_names = {'left','right'};
path_suma_surfaces =   {fullfile(Fs_subjects_folder,subj_id,'SUMA',[mesh_type '.lh.pial.asc']), ...
                        fullfile(Fs_subjects_folder,subj_id,'SUMA',[mesh_type '.rh.pial.asc']), ...
                        fullfile(Fs_subjects_folder,subj_id,'SUMA',[mesh_type '.lh.inflated.asc']),...
                        fullfile(Fs_subjects_folder,subj_id,'SUMA',[mesh_type '.rh.inflated.asc'])};
path_suma_curvatures = {fullfile(Fs_subjects_folder,subj_id,'SUMA',[mesh_type '.lh.curv.niml.dset']), ...
                        fullfile(Fs_subjects_folder,subj_id,'SUMA',[mesh_type '.rh.curv.niml.dset'])};
path_suma_thickness =  {fullfile(Fs_subjects_folder,subj_id,'SUMA',[mesh_type '.lh.thickness.niml.dset']), ...
                        fullfile(Fs_subjects_folder,subj_id,'SUMA',[mesh_type '.rh.thickness.niml.dset'])};
                    
path_fs_aparc =        {fullfile(Fs_subjects_folder,subj_id,'SUMA',[mesh_type '.lh.aparc.a2009s.annot.niml.dset']),...
                        fullfile(Fs_subjects_folder,subj_id,'SUMA',[mesh_type '.rh.aparc.a2009s.annot.niml.dset'])};

% Generate left and right pial and inflated surface structures
brain_data = struct;
for i=1:4
    % Read ASCII files into mesh structures
    [brain_data.(surface_names{i}).vertices brain_data.(surface_names{i}).faces] = read_suma_asc(path_suma_surfaces{i});
    brain_data.(surface_names{i}).faces = brain_data.(surface_names{i}).faces(:,[1 3 2]) + 1;
    % Transform coordinates to RAS
    brain_data.(surface_names{i}).vertices = transformations.ras2afni.transformPointsInverse(brain_data.(surface_names{i}).vertices);
end
% Pull apart the two inflated hemispheres so they do not overlap if plotted together:
[brain_data.inflated_left, brain_data.inflated_right] = pull_apart_inflated_hemispheres(brain_data.inflated_left, brain_data.inflated_right, 10);
    

% Load FreeSurfer annotation color lookup table:
fs_lookup=read_fs_color_lookup_table;

% Generate maps for these surfaces
vertex_data = struct;
for i=1:2
    % Curvature (sucli/gyri) maps 
    curvmap = afni_niml_read(path_suma_curvatures{i});
    vertex_data.(map_names{i}).vertex_curvature_index = curvmap.nodes{1}.data;
    % Thickness maps
    thickmap = afni_niml_read(path_suma_thickness{i});
    vertex_data.(map_names{i}).vertex_thickness_index = thickmap.nodes{1}.data;
    % FreeSurfer Auto-Parcellation
    aparc = afni_niml_read(path_fs_aparc{mod(i+1,2)+1});
    node_label_data = aparc.nodes{1}.data;
    label_idx = zeros(length(node_label_data),1);
    for j = 1:length(node_label_data)
        label_idx(j) = find(floor(fs_lookup.code/1000)==10+i & fs_lookup.rgbCode == node_label_data(j));
    end
    vertex_data.(map_names{i}).auto_parcellation.code = double(fs_lookup.code(label_idx));
    vertex_data.(map_names{i}).auto_parcellation.name = fs_lookup.name(label_idx);
    vertex_data.(map_names{i}).auto_parcellation.rgb = double(fs_lookup.rgbv(label_idx,1:3))/255;
    vertex_data.(map_names{i}).auto_parcellation.lookup_table = fs_lookup;
    
    % add Kastner map data
end


% Find nearest vertex for each electrode
for e = 1:electrode_data.num_electrodes
    coord = electrode_data.electrodes.coordinates_ras(e,:);
    if isnan(coord(1)) || electrode_data.electrodes.depth(e)
        vert_idx(e) = nan;
    else
        if strcmp(electrode_data.electrodes.hemisphere{e},'Right');
            surf = brain_data.pial_right;
        else
            surf = brain_data.pial_left;
        end
        % The "err" vector holds the distance from each electrode to the
        % nearest point on the surface which it was projected to (in mm):
        [vert_idx(e), err(e)] = map_electrode_to_vertex(coord,surf);
    end
end
% Add a new column for nearest SUMA mesh vertex:
new_table_column = table(vert_idx','VariableNames',{'suma_vertex'});
% if it already exists from a previous run, remove it:
if any(strcmp(fieldnames(electrode_data.electrodes),'suma_vertex'))
    field_idx = find(strcmp(fieldnames(electrode_data.electrodes),'suma_vertex'));
    electrode_data.electrodes = electrode_data.electrodes(:,[1:(field_idx-1) (field_idx+1):end]);
end
electrode_data.electrodes = [electrode_data.electrodes new_table_column];

% Save brain data
output_path = fullfile(recon_folder_home,recon_folder_matlab,[subj_id '_' recon_filename_braindata_suma]); % freesurfer mesh data
disp(['Saving brain mesh data to ' output_path ]);
save(output_path,'brain_data');
% Save vertex data
path_map_data = fullfile(recon_folder_home,recon_folder_matlab,[subj_id '_' recon_filename_vertexdata_suma]);
disp(['Saving vertex data to ' path_map_data]);
save(path_map_data,'vertex_data');
% Save electrode data
path_electrode_data = fullfile(recon_folder_home,recon_folder_matlab,[subj_id '_' recon_filename_projected_electrode_data]);
disp(['Saving electrode data to ' path_electrode_data]);
save(path_electrode_data,'electrode_data');

disp('Done.');



%% Visualize subject vs. template brains with electrodes (SUMA mesh) 

ELECTRODE_PLOT_COLOR = 'w';
ELECTRODE_PLOT_SIZE = 30;

% load subject electrode data if it is not loaded (or a wrong subject's data is loaded)
if ~exist('electrode_data','var') || ~strcmp(subj_id,electrode_data.subject_id)
    % first check default path
    path_electrode_data = fullfile(recon_folder_home,recon_folder_matlab,[subj_id '_' recon_filename_projected_electrode_data]);
    if ~exist(path_electrode_data,'file')
        % if it is not found, select file 
        disp(['Electrode data file not found at ' path_electrode_data '. Select electrode data .mat file:']);
        [file_name_electrodes, directory_electrodes] = uigetfile({'*.mat*','mat file (*.mat*)';'*.*', 'All files (*.*)'},'Select electrodes file',[recon_folder_home '/Matlab'],'multiselect','off');
        path_electrode_data = [directory_electrodes file_name_electrodes];
    end
    disp(['Loading electrode data from: ' path_electrode_data ]);
    load(path_electrode_data);
else 
    disp(['Electrode data for ' subj_id ' already loaded']);
end

% load subject surface data if it is not loaded
if ~exist('brain_data','var')
    % first check default path
    path_brain_data = fullfile(recon_folder_home,recon_folder_matlab,[subj_id '_' recon_filename_braindata_suma]);
    if ~exist(path_brain_data,'file')
        % if it is not found, select file 
        disp(['Brain surface mesh data file not found at ' path_brain_data '. Select (SUMA) brain data .mat file:']);
        [file_name_brain_data, directory_brain_data] = uigetfile({'*.mat*','mat file (*.mat*)';'*.*', 'All files (*.*)'},'Select brain data file',[recon_folder_home '/Matlab'],'multiselect','off');
        path_brain_data = [directory_brain_data file_name_brain_data];
    end
    disp(['Loading subject brain surface data from: ' path_brain_data ]);
    load(path_brain_data);
    subject_brain_data = brain_data;
else 
    disp('Brain surface data already loaded');
    subject_brain_data = brain_data;
end

% load subject vertex data if it is not loaded
if ~exist('vertex_data','var')
    % first check default path
    path_vertex_data = fullfile(recon_folder_home,recon_folder_matlab,[subj_id '_' recon_filename_vertexdata_suma]);
    if ~exist(path_vertex_data,'file')
        % if it is not found, select file 
        disp(['Brain vertex data file not found at ' path_vertex_data '. Select vertex data .mat file:']);
        [file_name_vertex_data, directory_vertex_data] = uigetfile({'*.mat*','mat file (*.mat*)';'*.*', 'All files (*.*)'},'Select vertex data file',[recon_folder_home '/Matlab'],'multiselect','off');
        path_vertex_data = [directory_vertex_data file_name_vertex_data];
    end
    disp(['Loading subject vertex data from: ' path_vertex_data ]);
    load(path_vertex_data);
    subject_vertex_data = vertex_data;
else 
    disp('Subject vertex data already loaded');
    subject_vertex_data = vertex_data;
end

% Load template brain mesh if it is not loaded
if ~exist('template_brain_data','var')
    if ~exist(template_brain_data_file,'file')
        disp(['Template brain surface mesh data file not found at ' template_brain_data_file '. Select template brain data .mat file:'])
        [file_name_template_brain, directory_template_brain] = uigetfile({'*.mat','mat file (*.mat)'},'Select template brain data file',initial_folder,'multiselect','off');
        template_brain_data_file = fullfile(directory_template_brain,file_name_template_brain);    
    end
    disp(['Loading template brain surface data from: ' template_brain_data_file ]);
    Load = load(template_brain_data_file);
    template_brain_data = Load.brain_data;
else 
    disp('Template surface data already loaded');
end

% load template vertex data if it is not loaded
if ~exist('template_vertex_data','var')
    if ~exist(template_vertex_data_file,'file')
        disp(['Template vertex data file not found at ' template_vertex_data_file '. Select template vertex data .mat file:'])
        [file_name_template_vertex, directory_template_vertex] = uigetfile({'*.mat','mat file (*.mat)'},'Select template vertex data file',initial_folder,'multiselect','off');
        template_vertex_data_file = fullfile(directory_template_vertex,file_name_template_vertex);    
    end
    disp(['Loading template brain vertex data from: ' template_vertex_data_file ]);
    Load = load(template_vertex_data_file);
    template_vertex_data = Load.vertex_data;
else
    disp('Template vertex data already loaded');
end


hemisphere = {'right','left'};
% Plot pial surfaces with electrodes
% subject brain
for h = 1:2
    if any(strcmpi(hemisphere{h},electrode_data.electrodes.hemisphere)) % Plot hemisphere only if there are any electrodes on it
        figure;
        coordinates = electrode_data.electrodes.coordinates_ras(strcmpi(hemisphere{h},electrode_data.electrodes.hemisphere),:);
        % plot brain
        plot_mesh_brain(subject_brain_data.(['pial_' hemisphere{h}]), [0 0], subject_vertex_data.(hemisphere{h}).auto_parcellation.rgb); 
        % plot electrodes - scatter
        plot_data_on_mesh(coordinates, ELECTRODE_PLOT_COLOR, 'markersize', ELECTRODE_PLOT_SIZE); % Add a dot at each electrode location
        % plot electrodes - surface
%         plot_electrodes(coordinates,ones(length(coordinates),1), 'surface'); % plot electrodes as colored patches on the surface instead of dots
        % fix lighting and set title
        view([(h-1)*180+90 0]);fix_lighting; 
        title(['Subject ' subj_id ': ' hemisphere{h} ' pial surface']);
        pull_3d_scatter_dots(1)
    end
end
% Template brain
for h = 1:2
    if any(strcmpi(hemisphere{h},electrode_data.electrodes.hemisphere)) % Plot hemisphere only if there are any electrodes on it
        figure;
        coordinates = electrode_data.electrodes.suma_vertex(strcmpi(hemisphere{h},electrode_data.electrodes.hemisphere)); 
        % plot brain
        plot_mesh_brain(template_brain_data.(['pial_' hemisphere{h}]), [0 0], template_vertex_data.(hemisphere{h}).auto_parcellation.rgb); 
        % plot electrodes - scatter
        plot_data_on_mesh(coordinates, ELECTRODE_PLOT_COLOR, 'markersize', ELECTRODE_PLOT_SIZE); % Add a dot at each electrode location
        % plot electrodes - surface
%         plot_electrodes(coordinates,ones(length(coordinates),1), 'surface'); % plot electrodes as colored patches on the surface instead of dots
        % fix lighting and set title
        view([(h-1)*180+90 0]);fix_lighting; 
        title(['Template brain: ' hemisphere{h} ' pial surface']);
        pull_3d_scatter_dots(1)
    end
end

% Plot inflated surfaces with auto-parcellation colormap and electrodes
% subject brain
for h = 1:2
    if any(strcmpi(hemisphere{h},electrode_data.electrodes.hemisphere)) % Plot hemisphere only if there are any electrodes on it
        figure;
        coordinates = electrode_data.electrodes.suma_vertex(strcmpi(hemisphere{h},electrode_data.electrodes.hemisphere),:);
        % plot brain
        plot_mesh_brain(subject_brain_data.(['inflated_' hemisphere{h}]), [0 0], subject_vertex_data.(hemisphere{h}).auto_parcellation.rgb); 
        % plot electrodes - scatter
        plot_data_on_mesh(coordinates, ELECTRODE_PLOT_COLOR, 'markersize', ELECTRODE_PLOT_SIZE); % Add a dot at each electrode location
        % fix lighting and set title
        view([(h-1)*180+90 0]);fix_lighting; 
        title(['Subject ' subj_id ': ' hemisphere{h} ' inflated surface']);
        pull_3d_scatter_dots(1)
    end
end
% Template brain
for h = 1:2
    if any(strcmpi(hemisphere{h},electrode_data.electrodes.hemisphere)) % Plot hemisphere only if there are any electrodes on it
        figure;
        coordinates = electrode_data.electrodes.suma_vertex(strcmpi(hemisphere{h},electrode_data.electrodes.hemisphere),:);
        % plot brain
        plot_mesh_brain(template_brain_data.(['inflated_' hemisphere{h}]), [0 0], template_vertex_data.(hemisphere{h}).auto_parcellation.rgb); 
        % plot electrodes - scatter
        plot_data_on_mesh(coordinates, ELECTRODE_PLOT_COLOR, 'markersize', ELECTRODE_PLOT_SIZE); % Add a dot at each electrode location
        % fix lighting and set title
        view([(h-1)*180+90 0]);fix_lighting; 
        title(['Template brain: ' hemisphere{h} ' inflated surface']);
        pull_3d_scatter_dots(1)
    end
end

% Show auto-parcellation color lookup table: 
show_fs_lookup_table(subject_vertex_data.right.auto_parcellation.lookup_table, subject_vertex_data.right.auto_parcellation.code);



%% Export final reconstruction results 

% Load subject electrode data
path_electrode_data = fullfile(recon_folder_home,recon_folder_matlab,[subj_id '_' recon_filename_projected_electrode_data]);
if ~exist(path_electrode_data,'file')
    disp('Select electrodes Matlab file');
    [file_name_elect, directory_elect] = uigetfile({'*.mat','mat file (*.mat)'},'Select electrodes Matlab file',recon_folder_home,'multiselect','off');
    path_electrode_data = fullfile(directory_elect,file_name_elect);    
end
disp(['Loading electrode data from: ' path_electrode_data ]);
load(path_electrode_data);

% Load subject brain mesh data
path_brain_data = fullfile(recon_folder_home,recon_folder_matlab,[subj_id '_' recon_filename_braindata_suma]);
if ~exist(path_brain_data,'file')
    disp('Select subject brain data .mat file');
    [file_name_brain_data, directory_brain_data] = uigetfile({'*.mat','mat file (*.mat)'},'Select (SUMA) brain data .mat file',initial_folder,'multiselect','off');
    path_brain_data = fullfile(directory_brain_data,file_name_brain_data);    
end
disp(['Loading brain mesh data from: ' path_brain_data ]);
load(path_brain_data);

% Load subject vertex data
path_vertex_data = fullfile(recon_folder_home,recon_folder_matlab,[subj_id '_' recon_filename_vertexdata_suma]);
if ~exist(path_brain_data,'file')
    disp('Select vertex data .mat file');
    [file_name_vertex_data, directory_vertex_data] = uigetfile({'*.mat','mat file (*.mat)'},'Select vertex data .mat file',initial_folder,'multiselect','off');
    path_vertex_data = fullfile(directory_vertex_data,file_name_vertex_data);    
end
disp(['Loading vertex map data from: ' path_vertex_data ]);
load(path_vertex_data);

% Save electrode data
output_path = fullfile(recon_folder_home, recon_filename_electrodedata_final); % freesurfer mesh data
disp(['Saving electrode data to ' output_path ]);
save(output_path,'electrode_data');
% Save brain data
output_path = fullfile(recon_folder_home, recon_filename_braindata_final); % freesurfer mesh data
disp(['Saving brain mesh data to ' output_path ]);
save(output_path,'brain_data');
% Save vertex data
path_vertex_data = fullfile(recon_folder_home, recon_filename_vertexdata_final);
disp(['Saving vertex data to ' path_vertex_data]);
save(path_vertex_data,'vertex_data');

disp('Done.');

