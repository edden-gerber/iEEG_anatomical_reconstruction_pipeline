%% Extra and obsolete
%% Convert .mgz to .nii
fs_shell_initialize_cmd = ['export FREESURFER_HOME=' Fs_directory '; source $FREESURFER_HOME/SetUpFreeSurfer.sh; '];

% Select original file
disp('Select file to convert');
[file_name_conv, directory_conv] = uigetfile({'*.*', 'All files (*.*)'},'Select MR file',pwd,'multiselect','off');
path_conv = [directory_conv file_name_conv];

% Select output file name (with extension)
disp('Enter new file name');
[file_name_conv2, directory_conv2] = uiputfile({'*.*', 'All files (*.*)'},'Save mgrid file',directory_conv);
path_conv2 = [directory_conv2 file_name_conv2];

% Convert
system([fs_shell_initialize_cmd 'mri_convert ' path_conv ' ' path_conv2]);

%% Co-register CT to native MR scan - Automatic
% Note: Ideally, this should be done with a post-op MRI scan, as the
% operation leads to significant changes in the brain's shape and thus to
% larger errors in the co-registration of post-op CT to pre-op MR. 

% Define parameters
% linear_coreg_alg = {'rigid','affine','similarity','affine2d','rigid2d','similarity2d'};
% nonlin_coreg_alg = {'none','rigid','affine','similarity'};
%
% The suitable method for this step (CT->MR) is usually the rigid linear transform 
% (since it's the same brain): 
linear_coreg_alg = {'rigid'};
nonlin_coreg_alg = {};

% Create results folder
directory_output_images = [directory_home '/Registrations'];
if ~exist(directory_output_images,'dir')
    mkdir(directory_output_images);
end

% Run BioImage Suite's co-registration functions
% Linear: 
for i = 1:length(linear_coreg_alg)
    alg = linear_coreg_alg{i};
    out_file_name = [directory_output_images filesep 'Transform_CT2MR_Linear_' alg '.grd'];
    BIS_console_command = ['bis_linearintensityregister --dogui 0 --inp "' path_mr '" --inp2 "' path_ct '" --out "' out_file_name '" '...
        '--mode ' alg ' --numberoflevels 3 --resolution 1.5 --useinitial 0 --reslimage Resliced --metric NMI --useweightimage 0 '...
        '--iterations 15 --resolutionrate 2 --autonormalize 1 --optimization default --numberofbins 64 --numberofsteps 1 '...
        '--stepsize 1.0 --autooptscalefactor 1 --optscalefactor 1.0'];
    
    disp(['Running linear coregistration algorithm: ' alg '...']);
    system([BIS_directory filesep  'bin' filesep BIS_console_command]);
end
% Nonlinear: 
for i = 1:length(nonlin_coreg_alg)
    alg = nonlin_coreg_alg{i};
    out_file_name = [directory_output_images filesep 'Transform_CT2MR_Nonlin_' alg '.grd'];
%     out_image_file_name = [directory_output filesep 'Transformed_MR_Nonlin_' alg '.' ref_file((strfind(ref_file,'.')+1):end)];
    BIS_console_command = ['bis_nonlinearintensityregister -dogui 0 --inp "' path_mr '" --inp2 "' path_ct '" --out "' out_file_name '" '...
        '--initialmode ' alg ' --spacing 15.0 --smoothness 0.001 --numberoflevels 3 --resolution 1.5 --useinitial 0 --reslimage Resliced '...
        '--metric NMI --useweightimage 0 --iterations 15 --resolutionrate 2 --autonormalize 1 --optimization default --numberofbins 64 '...
        '--spacingrate 2.0 --extralevels 0 --windowsize 1.0 --numberofsteps 1 --stepsize 1.0'];

    disp(['Running nonlinear coregistration algorithm: ' alg '...']);
    system([BIS_directory filesep 'bin' filesep BIS_console_command]);
end

%% Co-register CT to native MR scan - Manual (in case there is a problem with the automatic script)
% Note: this step should eventually be replaced by a better co-registration
% procedure, based on inflation of the brain and co-registration of the
% resulting surfaces. 
%
% For each function (linear and non-linear), load the reference image (MR or MNI) file into the
% Reference Image" slot (under the Input tab), and the transform image (CT or native MR) file 
% into the "Transform Image" slot. Then go to the Options tab and select the Transformation Mode. 
% Click "Run Algorithm" and when it''s done, save the result (it's recommended to do this in a 
% designated "Coregistration Results" folder). In the file name, include "linear" or "nonlinear" 
% and the transformation mode ("rigid", "affine" etc.). Continue selecting each transformation 
% mode you want to run and saving the results with the appropriate file name. You will need to 
% repeat the file selection procedure for the non-linear function.

BIS_console_command = 'bis_linearintensityregister -dogui 1';
system([BIS_directory filesep 'bin' filesep BIS_console_command]);

BIS_console_command = 'bis_nonlinearintensityregister -dogui 1';
system([BIS_directory filesep 'bin' filesep BIS_console_command]);

%% Create brain-extracted and cropped MRI scan (manually using BioImage Suite)
% NOTE: THIS SECTION IS NOW OBSOLETE AS IT IS DONE AUTOMATICALLY BY
% FREESURFER. 
%
% This phase takes the coriginal MRI scan, applies a "brain extraction"
% algorithm to remove the skull, and then an additional manual step (can be
% somewhat tedious) to manually crop out any remaining non-brain tissue, as
% well as any parts of the brain you don't want to see (usually cerebellum
% and brainstem, which may ventral electrode strips). This step (or just the
% manual part) may be skipped if you don't need to present native-space 3D 
% models of your patient (e.g. if you're only using a single MNI model), but 
% you will not be able to use the "3D volume" display to check and correct the
% electrode locations co-registered from the CT image (as they will be
% hidden inside the head), only the standard "3-Slice" display - which may
% be less convenient. 
% 
% Instructions: 
% When the BioImage Suite menu appears, select the "Editors" tab and click on
% the "Orthogonal Object Editor"
% 1. Load the original MR image.
% 2. Munu-select Segmentation->Brain Extracion Tool (if it is not there,
% you have not successfully installed the FSL package support). 
% 3. Click "Exectute Brain Extraction Tool". No need to enter the output
% file name, it will be set automatically to <original_file>_stripped. Wait
% for a message to appear indicating the completion of the process. You
% should look at the result to make sure it's satisfactory (e.g. using
% MRIcon, with the original image as background and stripped brain as
% overlay). 
% 4. Next we manually clean any remaining non-brain tissue and remove
% unwanted brain parts (usually brainstem+cerebellum). Close the brain
% extraction tool and load the new stripped brain image. 
% 5. Menu-select Objectmap->Show Paint Controls to open the painting
% window.
% 6. In the paint window, choose a brush color and size (starts at "off").
% Check "3D", "Use Threshold", and "Connect" (It's recommended to change the
% lower threhold above the checkboxes to 1). 
% 7. Start "painting" over the parts to remove. Shift-click switches
% the paint mode on/off. Painting with the black bruch erases the object
% map. Remember that you are painting with a 3D brush, so keep checking the
% result along all relevant dimensions (it takes some getting used to). 
% 8. When finished menu-select Objectmap->Save Objectmap (in the viewer
% window, not the paint control window). 
% 
% 9. The script will now ask you for the file names of the MR file (choose
% the scan you used for creating the object map, i.e. the skull-stripped
% scan), and the object map file. To apply the object map such that
% selected volumes will be removed, the map will need to be inversed. In
% the final step the map will be applied to the scan, resulting in a new
% image file with the suffix "_clean". 
% You should look at this image as well using MRIcron as an overlay over the
% original scan. 

% Run BioImage Suite:
system([BIS_directory '/start_bioimagesuite']);

% Apply the generated objectmap to the image:
disp('Select MR scan file');
[file_name_mr, directory_mr] = uigetfile({'*.nii*','Nifti file (*.nii*)';'*.*', 'All files (*.*)'},'Select MR file',initial_directory,'multiselect','off');
path_mr = [directory_mr file_name_mr];

disp('Select objectmap file');
[file_name_om, directory_om] = uigetfile({'*.nii*','Nifti file (*.nii*)';'*.*', 'All files (*.*)'},'Select MR file',directory_mr,'multiselect','off');
om_file = [directory_om file_name_om];

% Invert object map
inv = questdlg('Create inverted objectmap (selected volume will be deleted)?','Apply Objectmap','Yes','No','Yes');
if strcmp(inv,'yes')
    disp('Inverting object map...');
    inv_om_file = [om_file(1:(strfind(om_file,'.')-1)) '_inv' om_file(strfind(om_file,'.'):end)]; % just add the string "_inv" before the file extension. 
    system(['sh ' BIS_directory '/bin/bis_piecewiseimagemap --inp ' om_file ' --out ' inv_om_file ' --on1 1 --x1 0 --y1 1 --on2 1 --x2 1 --y2 0']);
    disp(['Saved inverted objectmap as ' inv_om_file '.']);
else
    inv_om_file = om_file;
end

% Apply object map
disp('Applying objectmap...');
out_file = [path_mr(1:(strfind(path_mr,'.')-1)) '_clean' path_mr(strfind(path_mr,'.'):end)]; % just add the string "_clean" before the file extension. 
system(['fslmaths "' path_mr '" -mul "' inv_om_file '" "' out_file '"' ]);
disp(['Saved cleaned scan as ' out_file'.']); 

%% Transform mgrid electrode coordinates for scans generated by FreeSurfer
% Recommended to save as <subj_id>_coordinates_MR_Fs.mgrid
% When this step is finished, check the results with BioImage Suite by
% loading the new mgrid file with the skull-stripped cortex images. 

% Select input mgrid file
disp('Select original mgrid file');
[file_name_mgrid, directory_mgrid] = uigetfile({'*.mgrid*','mgrid file (*.mgrid*)';'*.*', 'All files (*.*)'},'Select mgrid file',[directory_home '/Mgrids'],'multiselect','off');
path_mgrid = [directory_mgrid file_name_mgrid];

% Read grid into Matlab
Grid = read_mgrid_file(path_mgrid);
elec = [Grid.coordinates ones(length(Grid.coordinates),1)];

% Get transformation matrices
fs_shell_initialize_cmd = ['export FREESURFER_HOME=' Fs_directory '; source $FREESURFER_HOME/SetUpFreeSurfer.sh; '];
[~, cmdout] = system([fs_shell_initialize_cmd 'mri_info --vox2ras-tkr ' Fs_subjects_directory '/' subj_id '/mri/orig/001.mgz']);
transform_matrix1 = str2num(cmdout);
[~, cmdout] = system([fs_shell_initialize_cmd 'mri_info --ras2vox-tkr ' Fs_subjects_directory '/' subj_id '/mri/orig.mgz']);
transform_matrix2 = str2num(cmdout);

% Transform electrode coordinates
new_elec = (transform_matrix2 * transform_matrix1 * elec')';
new_elec = new_elec(:,1:3);

% Select output mgrid file
disp('Enter new mgrid file name');
[file_name_mgrid2, directory_mgrid2] = uiputfile({'*.mgrid*','mgrid file (*.mgrid*)';'*.*', 'All files (*.*)'},'Save mgrid file',directory_mgrid);
path_mgrid2 = [directory_mgrid2 file_name_mgrid2];

% Save coordinates in a new mgrid file
edit_mgrid_file(path_mgrid,path_mgrid2,new_elec);

disp('done.');

%% Correct errors in electrode localization by projecting to the cortical surface - using generate_volume_surface
% Errors in electrode localization can occur in the CT->MR coregistration
% phase, especially if the MR image is pre-operation. These errors often
% result in electrodes being localized to coordinates within or outside of
% the cortical volume rather than on its surface (for surface electrodes of
% course). In this step, the electrode locations are projected to the
% nearest point on the smoothed surface surrounding the cortical
% hemisphere. 
% To check that the generated surfaces are appropriate, load them as
% overlays over the original MR scan (e.g. using MRICron). 
% 
% To test the results, compare the old vs. the new mgrid file in BioImage
% Suite using the appropriate cortex hemisphere image. 
%
% The procedure uses a modified version of code from the CTMR toolbox by 
% Dora Hermes et al. - make sure you cite them for the reconstruction 
% procedure if using this tool: 
% 
% Hermes, D., Miller, K. J., Noordmans, H. J., Vansteensel, M. J., & Ramsey, 
% N. F. (2010). Automated electrocorticographic electrode localization on 
% individually rendered brain surfaces. Journal of Neuroscience Methods, 
% 185(2), 293–8.

% Generate surface images for left and right cortical hemispheres
smooth_level = 9; % This is the degree of smoothing applied to the surface 
    % (lower values mean closer tracking of the cortical surface - a low
    % enough value is wanted so that electrodes will not float above the
    % surface due to exessive smoothing, but very low values will mean that
    % the surface will follow the pial surface into sulci, which may lead
    % to erroneous projections as electrodes are not expected to be placed
    % inside sulci. NOTE THAT EVEN VALUES INDUCE AN ERROR IN THE SMOOTHING
    % FUNCTION. 
brain_cutoff = 0.2; % This is the intensity cutoff value where the surface 
    % is generated around the brain. The volume is initially set to 1's and
    % 0's and then smoothed to create a [0 1] gradient around the edges of
    % the volume. A smaller value means a "tighter" surface. 

lateral_axis_certainty_boundary = 10; % grids whose center is less than or 
    % equal to this number of voxels from the center of the image (on the
    % lateral axis) will be automatically attributed to the left or right 
    % hemisphere - the user will be prompted for the information instead. 

% Generate surfaces
l_cortex_image = [ directory_home '/Generated_Images/' subj_id '_cortex_left.nii.gz'];
r_cortex_image = [ directory_home '/Generated_Images/' subj_id '_cortex_right.nii.gz'];
surface_file_name_l = [ directory_home '/Generated_Images/' subj_id '_outer_surface_left_sl' num2str(smooth_level) '_bc' num2str(brain_cutoff) '.nii.gz'];
r_surface_file_name = [ directory_home '/Generated_Images/' subj_id '_outer_surface_right_sl' num2str(smooth_level) '_bc' num2str(brain_cutoff) '.nii.gz'];
disp('Generating cortical surface images...');
generate_volume_surface(l_cortex_image, surface_file_name_l, smooth_level, brain_cutoff);
generate_volume_surface(r_cortex_image, r_surface_file_name, smooth_level, brain_cutoff);

% Select input mgrid file
disp('Select electrode coordinates mgrid file');
[file_name_mgrid, directory_mgrid] = uigetfile({'*.mgrid*','mgrid file (*.mgrid*)';'*.*', 'All files (*.*)'},'Select mgrid file',[directory_home '/Mgrids'],'multiselect','off');
path_mgrid = [directory_mgrid file_name_mgrid];

% Read the electrode coordinates into Matlab 
electrodes = read_mgrid_file(path_mgrid);

% Read the surface images into Matlab
surfl = load_nifti(surface_file_name_l);
surfr = load_nifti(r_surface_file_name);

% Transform surface volume to vertex list, and transform surface to the same coordinate system as the electrodes
[x,y,z]=ind2sub(size(surfl.vol),find(surfl.vol>0)); 
% surface_vert_l = [256-x 256-y z];
surface_vert_l = [x y z];
[x,y,z]=ind2sub(size(surfr.vol),find(surfr.vol>0)); 
% surface_vert_r = [256-x 256-y z];
surface_vert_r = [x y z];

% Go over each individual electrode grid/strip and project it to the
% surface. For 1-dimensional strips, electrodes are projected to the
% nearest point on the surface. For NxM (N=2,M>2) grids, electrodes are 
% projected to the nearest point on the surface along the vector 
% perpendicular to the surface determined by the electrode and its 3 
% nearest neighbors. For NxM (N>2,M>2) grids, the surface is determined by
% the electrode and its 4 nearest neighbors. 

% First, query whether there are any depth electrodes - these will not be
% projected to the surface. 
choice = questdlg('Are there any depth electrode strips (or any strips that should not be projected to the cortical surface)?','Depth electrodes','Yes','No','No');
if strcmp(choice,'Yes')
    name_list = [];
    for i=1:electrodes.num_grids
        name_list = [name_list ' ' num2str(i) '. ' electrodes.grid{i}.name];
    end
    user_input_str = inputdlg(['Electrode grid list: ' name_list '. Enter depth electrode strip indexes:']);
    depth_elect_strip_indexes = str2num(user_input_str{1});
else
    depth_elect_strip_indexes = [];
end

corrected_coordinates = [];
for g = 1:electrodes.num_grids
    
    % get grid information
    coordinates = electrodes.grid{g}.coordinates;
    dimensions = electrodes.grid{g}.dimensions;
    name = electrodes.grid{g}.name;
    
    if ismember(g,depth_elect_strip_indexes) % depth electrodes - do not project 
        proj_electrodes = coordinates;
    else % not depth electrodes - project
        % determine projection method
        if min(dimensions) == 1
            num_elec_on_plane = 1;
        elseif min(dimensions) == 2
            num_elec_on_plane = 4;
        else
            num_elec_on_plane = 5;
        end

        % choose hemisphere to project to based on the mean location of the
        % electrodes on the lateral axis. If no clear decision can be make, ask 
        % the user. 
        mean_x_loc = mean(coordinates(:,1)); 
        if mean_x_loc <= 128 - lateral_axis_certainty_boundary
            surf = surface_vert_r;
            disp(['Electrode grid ' name ' is projected to the right hemisphere surface.']);
        elseif mean_x_loc > 128 + lateral_axis_certainty_boundary
            surf = surface_vert_l;
            disp(['Electrode grid ' name ' is projected to the left hemisphere surface.']);
        else % too close to center - prompt the user
            if mean_x_loc <= 128
                suggested_side = 'right';
            else
                suggested_side = 'left';
            end
            choice = questdlg(['Electrode grid ' name ' is too close to the middle to be attributed to the left or right hemisphere. '...
                'Which one is it? More probably option is the ' suggested_side ' side.'],'Pick hemisphere for grid','left','right',suggested_side);
            if strcmp(choice,'right')
                surf = surface_vert_r;
                disp(['Electrode grid ' name ' is projected to the right hemisphere surface.']);
            elseif strcmp(choice,'left')
                surf = surface_vert_l;
                disp(['Electrode grid ' name ' is projected to the left hemisphere surface.']);
            else
                error('No hemisphere selected.');
            end
        end

        % Project electrodes to the surface
%         proj_electrodes = project_electrodes_to_surf(coordinates, surf);
        proj_electrodes = project_electrodes_to_surf_debug(coordinates, surf, cortex);
    end
    corrected_coordinates = [corrected_coordinates ; proj_electrodes];
end

% Save mgrid file
out_mgrid_name = [directory_mgrid file_name_mgrid(1:strfind(file_name_mgrid,'.')-1) '_Projected' file_name_mgrid(strfind(file_name_mgrid,'.'):end)];
edit_mgrid_file(path_mgrid,out_mgrid_name,corrected_coordinates);
disp(['Saved projected electrode coordinates in mgrid file: ' out_mgrid_name ]);

%% Map electrode coordinates to MNI space using an MNI template

DEFAULT_MNI_BRAIN_DATA_PATH = fullfile(initial_folder,'fsaverage','brain_data.mat');
DEFAULT_ELECTRODE_DATA_MAT_FILE_NAME = 'electrode_data_projected.mat';

% Load MNI brain mesh
disp('Loading MNI template brain');
path_mni_brain_data = DEFAULT_MNI_BRAIN_DATA_PATH;
if ~exist(path_mni_brain_data,'file')
    disp('Select MNI brain data .mat file');
    [file_name_mni, directory_mni] = uigetfile({'*.mat','mat file (*.mat)'},'Select electrodes Matlab file',initial_folder,'multiselect','off');
    path_mni_brain_data = fullfile(directory_mni,file_name_mni);    
end
mni_brain = load(path_mni_brain_data);

disp('Loading electrode data');
% Load electrode data
path_electrode_data = fullfile(recon_folder_home,recon_folder_matlab,DEFAULT_ELECTRODE_DATA_MAT_FILE_NAME);
if ~exist(path_electrode_data,'file')
    disp('Select electrodes Matlab file');
    [file_name_elect, directory_elect] = uigetfile({'*.mat','mat file (*.mat)'},'Select electrodes Matlab file',recon_folder_home,'multiselect','off');
    path_electrode_data = fullfile(directory_elect,file_name_elect);    
end
Load = load(path_electrode_data);
electrode_data = Load.electrode_data;

% Get MNI coordinates for each electrode by taking the coordinates of
% corresponding vertices on the MNI brain
elec_vertices = electrode_data.electrodes.nearest_vertex;
right_hemisphere_electrodes = find(strcmp(electrode_data.electrodes.hemisphere,'Right'));
left_hemisphere_electrodes = find(strcmp(electrode_data.electrodes.hemisphere,'Left'));
mni_coord(right_hemisphere_electrodes,:) = mni_brain.mesh_suma.pial_right.vertices(elec_vertices(right_hemisphere_electrodes),:);
mni_coord(left_hemisphere_electrodes,:) = mni_brain.mesh_suma.pial_left.vertices(elec_vertices(left_hemisphere_electrodes),:);

% Add MNI coordinates column to the electordes table
new_table_column = table(mni_coord,'VariableNames',{'coordinates_mni'});
% if it already exists from a previous run, remove it:
if any(strcmp(fieldnames(electrode_data.electrodes),'coordinates_mni'))
    field_idx = find(strcmp(fieldnames(electrode_data.electrodes),'coordinates_mni'));
    electrode_data.electrodes = electrode_data.electrodes(:,[1:(field_idx-1) (field_idx+1):end]);
end
electrode_data.electrodes = [electrode_data.electrodes(:,1:(end-1)) new_table_column electrode_data.electrodes(:,end)];

% path_electrode_data = 'home/hcnl/shared_folders/ECoG_Patient_Data/ST18/Recon_HCNL/Matlab/electrode_data_projected_mni.mat';
% Save electrode data
disp(['Saving electrode data to ' path_electrode_data]);
save(path_electrode_data,'electrode_data');
