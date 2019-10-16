%% Set parameters

% Global parameters
initial_folder = '/home/hcnl/shared_folders/ECoG_Patient_Data'; % General home folder
BIS_folder = '/usr/local/BIS/bioimagesuite30_64'; % BioImage Suite scripts folder
Fs_folder = '/usr/local/freesurfer'; % FreeSurfer folder
Fs_subjects_folder = '/home/hcnl/FreeSurfer/Subjects'; % FreeSurfer subjects folder. If using a virtual machine, note 
                                                          % that making this a shared folder may not allow Fs to run (due 
                                                          % to short delays in reading and writing to these folders).                                                           % to short delays in reading and writing to these folders).                                                    

recon_folder_generated_images = 'Generated_Images';
recon_folder_matlab = 'Matlab';
recon_filename_braindata_suma = 'brain_data_suma.mat';
recon_filename_vertexdata_suma = 'vertex_data_suma.mat';
recon_filename_braindata_final = 'brain_data.mat';
recon_filename_vertexdata_final= 'vertex_data.mat';

% Enter "subject ID" for the template brain (e.g. MNI152)
subj_id = [];
if isempty(subj_id)
    subj_id = inputdlg('Enter template brain ID: ','Template brain ID');
    subj_id = subj_id{1};
end

% Select reconstruction home directory 
recon_folder_home = [];
if isempty(recon_folder_home)
    recon_folder_home = uigetdir(initial_folder,'Select recon folder');
end

% Select MR image file
path_mr = [];
if isempty(path_mr)
    disp('Select original MR scan file');
    [file_name_mr, directory_mr] = uigetfile({'*.nii*','Nifti file (*.nii*)';'*.*', 'All files (*.*)'},'Select MR file',recon_folder_home,'multiselect','off');
    path_mr = fullfile(directory_mr,file_name_mr);
end

%% Run FreeSurfer on the MR scan
% This step creates the segmented volume and surface data for the selected
% brain. 
% NOTE THAT THIS STEP CAN TAKE UP TO 24 HOURS

fs_shell_initialize_cmd = ['export FREESURFER_HOME=' Fs_folder '; source $FREESURFER_HOME/SetUpFreeSurfer.sh; '];
disp(['Running FreeSurfer for subject ' subj_id '...']);
system([fs_shell_initialize_cmd 'recon-all -all -i ' path_mr ' -s ' subj_id]);

%% Look at the results using freeview
% Using FreeSurfer's freeview utility. 
% If there are errors in the reconstruction (most likely regions where Fs
% fails to identify white matter, resulting in an erroneous pial surface
% shape), use freeview to correct the errors and re-run the appropriate
% segments of the FreeSurfer analysis. Consult Fs's wiki and YouTube videos
% for assistance. 

fs_shell_initialize_cmd = ['export FREESURFER_HOME=' Fs_folder '; source $FREESURFER_HOME/SetUpFreeSurfer.sh; '];
system([fs_shell_initialize_cmd 'freeview -v "' Fs_subjects_folder '/' subj_id '/mri/T1.mgz" ' ...
    '-f "' Fs_subjects_folder '/' subj_id '/surf/lh.pial" ' '-f "' Fs_subjects_folder '/' subj_id '/surf/rh.pial"']);

%% Import FreeSurfer-generated images and convert to nifti, produce stripped cortex volumes

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

%% Run SUMA to create a standardized mesh based on the FreeSurfer surfaces
% This creates a SUMA folder in the FreeSurfer subject folder. 

curr_path = pwd; % We will need to go to the subject FreeSurfer folder, so this is to get back to the current folder. 
cd(fullfile(Fs_subjects_directory,subj_id));
system(['@SUMA_Make_Spec_FS -sid ' subj_id]);
cd(curr_path);

%% Load the SUMA surface mesh and anatomical metadata
% Based on code by Tal Golan @ Rafael Malach Lab

% Get transformation matrix from the AFNI coordinate space to RAS
afni_shift1d_file = fullfile(Fs_subjects_folder,subj_id,'SUMA','brainmask_shft.1D');
[err, ras2afni_matrix] = Read_1D(afni_shift1d_file); assert(err==0);
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

% Generate vertex maps for these surfaces
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
end


% Save brain data
output_path = fullfile(recon_folder_home,recon_folder_matlab,[subj_id '_' recon_filename_braindata_suma]); % freesurfer mesh data
disp(['Saving brain mesh data to ' output_path ]);
save(output_path,'brain_data');
% Save vertex data
path_vertex_data = fullfile(recon_folder_home,recon_folder_matlab,[subj_id '_' recon_filename_vertexdata_suma]);
disp(['Saving vertex data to ' path_vertex_data]);
save(path_vertex_data,'vertex_data');

disp('Done.');

%% Plot 3D cortex

% load surface data if it is not loaded
if ~exist('brain_data','var')
    % first check default path
    path_brain_data = fullfile(recon_folder_home,recon_folder_matlab,[subj_id '_' recon_filename_braindata_suma]);
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
% load vertex data if it is not loaded
if ~exist('vertex_data','var')
    % first check default path
    path_vertex_data = fullfile(recon_folder_home,recon_folder_matlab,[subj_id '_' recon_filename_vertexdata_suma]);
    if ~exist(path_vertex_data,'file')
        % if it is not found, select file 
        disp(['Brain vertex data file not found at ' path_vertex_data '. Select vertex data .mat file:']);
        [file_name_vertex_data, directory_vertex_data] = uigetfile({'*.mat*','mat file (*.mat*)';'*.*', 'All files (*.*)'},'Select electrodes file',[recon_folder_home '/Matlab'],'multiselect','off');
        path_vertex_data = [directory_vertex_data file_name_vertex_data];
    end
    disp(['Loading brain vertex data from: ' path_vertex_data ]);
    load(path_vertex_data);
else 
    disp('Brain vertex data already loaded');
end

% Plot pial surfaces
figure;
plot_mesh_brain(brain_data.pial_right, vertex_data.right.auto_parcellation.rgb); % Plot right cortical hemisphere with auto-parcellation
plot_mesh_brain(brain_data.pial_left, vertex_data.left.auto_parcellation.rgb); % Plot left cortical hemisphere with auto-parcellation
view([90 0]);fix_lighting;
title(['Subject ' subj_id ': Pial surface']);

% Plot inflated surfaces
figure;
plot_mesh_brain(brain_data.inflated_right, vertex_data.right.auto_parcellation.rgb); % Plot right cortical hemisphere with auto-parcellation
plot_mesh_brain(brain_data.inflated_left, vertex_data.left.auto_parcellation.rgb); % Plot left cortical hemisphere with auto-parcellation
view([90 0]);fix_lighting; 
title(['Subject ' subj_id ': Inflated surface']);

% Show auto-parcellation color lookup table
show_fs_lookup_table(vertex_data.right.auto_parcellation.lookup_table, vertex_data.right.auto_parcellation.code);

%% Export final results to recon home folder 

% load surface data if it is not loaded
if ~exist('brain_data','var')
    % first check default path
    path_brain_data = fullfile(recon_folder_home,recon_folder_matlab,[subj_id '_' recon_filename_braindata_suma]);
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
% load vertex data if it is not loaded
if ~exist('vertex_data','var')
    % first check default path
    path_vertex_data = fullfile(recon_folder_home,recon_folder_matlab,[subj_id '_' recon_filename_vertexdata_suma]);
    if ~exist(path_vertex_data,'file')
        % if it is not found, select file 
        disp(['Brain vertex data file not found at ' path_vertex_data '. Select vertex data .mat file:']);
        [file_name_vertex_data, directory_vertex_data] = uigetfile({'*.mat*','mat file (*.mat*)';'*.*', 'All files (*.*)'},'Select electrodes file',[recon_folder_home '/Matlab'],'multiselect','off');
        path_vertex_data = [directory_vertex_data file_name_vertex_data];
    end
    disp(['Loading brain vertex data from: ' path_vertex_data ]);
    load(path_vertex_data);
else 
    disp('Brain vertex data already loaded');
end

% Save brain data
output_path = fullfile(recon_folder_home, recon_filename_braindata_final); % freesurfer mesh data
disp(['Saving brain mesh data to ' output_path ]);
save(output_path,'brain_data');
% Save vertex data
path_vertex_data = fullfile(recon_folder_home, recon_filename_vertexdata_final);
disp(['Saving vertex data to ' path_vertex_data]);
save(path_vertex_data,'vertex_data');

disp('Done.');

