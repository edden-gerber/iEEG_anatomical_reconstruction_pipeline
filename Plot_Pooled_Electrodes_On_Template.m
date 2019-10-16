%% Set parameters

initial_folder = 'L:\Experiments\ECoG Patient Data'; % General home folder
template_brain_data_file = 'L:\Experiments\ECoG Patient Data\fsaverage\brain_data.mat';
template_vertex_data_file = 'L:\Experiments\ECoG Patient Data\fsaverage\vertex_data.mat';

%% Select subjects and load data

% Option 1: Provide subject path list in code
electrode_data_path_list = {};

% Option 2: Use GUI 
if isempty(electrode_data_path_list)
    disp('Select subjects'' electrode data files (press Cancel to close the list).');
    cancel = 0;
    while ~cancel
        [file_name_electrodes, directory_electrodes] = uigetfile({'*.mat*','mat file (*.mat*)';'*.*', 'All files (*.*)'},'Select electrodes file',initial_folder,'multiselect','off');
        if file_name_electrodes == 0
            cancel = true;
        else
            path = [directory_electrodes file_name_electrodes];
            electrode_data_path_list{end+1,1} = path;
        end
    end
end
num_subjects = length(electrode_data_path_list);

% Load subjects' electrode data
clear('pooled_electrodes');
for s = 1:num_subjects
    Load = load(electrode_data_path_list{s});
    pooled_electrodes(s) = Load.electrode_data;
end

%% Create pooled electrodes table
pooled_electrode_table = table;
for s = 1:num_subjects
    subj_id = pooled_electrodes(s).subject_id; 
    subj_id_column = cell(pooled_electrodes(s).num_electrodes,1);
    subj_id_column(:) = {subj_id};
    subj_id_table = table(subj_id_column,'VariableNames',{'Subject_ID'});
    subject_electrode_table = [subj_id_table pooled_electrodes(s).electrodes];
%     subject_electrode_table = [subj_id_table pooled_electrodes(s).electrodes(:,[1:5 9])];
    pooled_electrode_table = [pooled_electrode_table ; subject_electrode_table];
end
num_electrodes = size(pooled_electrode_table,1);

%% Load template brain

if ~exist(template_brain_data_file,'file')
    % if it is not found, select file 
    disp(['Template brain data file not found at ' template_brain_data_file '. Select brain data .mat file:']);
    [file_name, directory] = uigetfile({'*.mat*','mat file (*.mat*)';'*.*', 'All files (*.*)'},'Select electrodes file',initial_folder,'multiselect','off');
    template_brain_data_file = [directory file_name];
end
disp(['Loading template brain mesh data from ' template_brain_data_file]);
load(template_brain_data_file);

if ~exist(template_vertex_data_file,'file')
    % if it is not found, select file 
    disp(['Template vertex data file not found at ' template_vertex_data_file '. Select vertex data .mat file:']);
    [file_name, directory] = uigetfile({'*.mat*','mat file (*.mat*)';'*.*', 'All files (*.*)'},'Select electrodes file',initial_folder,'multiselect','off');
    template_vertex_data_file = [directory file_name];
end
disp(['Loading template vertex map data from ' template_vertex_data_file]);
load(template_vertex_data_file);

%% Plot electrodes on template surface

% Surface type: options are 'pial' or 'inflated'
surface_type = 'inflated';
% Surface coloring: options are 'none', 'electrodes', 'parcellation', 'curvature', or 'subjects'
surface_color_mode = 'subjects';

ELECTRODE_PLOT_COLOR = 'b';
ELECTRODE_PLOT_SIZE = 200;
ELECTRODE_PATCH_RADIUS = 1;

electrode_color_values = ones(num_electrodes,1);
% electrode_color_values = (1:size(pooled_electrode_table,1))';

subj_colors_table = table(unique(pooled_electrode_table.Subject_ID),'variableNames',{'subj_id'});
subj_colors_table = [subj_colors_table table(distinguishable_colors(length(subj_colors_table.subj_id),[0.5 0.5 0.5]),'variableNames',{'rgb'})];

color_list = zeros(length(pooled_electrode_table.name),3);
for i=1:length(color_list)
    subj_id = pooled_electrode_table.Subject_ID{i};
    color_list(i,:) = subj_colors_table.rgb(find(strcmp(subj_id,subj_colors_table.subj_id)),:);
end

hemisphere = {'right','left'};
% Template brain
for h = 1:2
    if any(strcmpi(hemisphere{h},pooled_electrode_table.hemisphere)) % Plot hemisphere only if there are any electrodes on it
%         figure;

        hold all
        for i=1:num_subjects
            scatter(0,0,[],subj_colors_table.rgb(i,:),'filled')
        end
        legend(subj_colors_table.subj_id');

        coordinates = pooled_electrode_table.suma_vertex(strcmpi(hemisphere{h},pooled_electrode_table.hemisphere)); 
        
        switch surface_color_mode
            case 'none'
                plot_mesh_brain(brain_data.([surface_type '_' hemisphere{h}]), ones(length(brain_data.([surface_type '_' hemisphere{h}]).vertices),1)*1.2); 
                plot_electrodes(coordinates, ELECTRODE_PLOT_COLOR, 'markersize', ELECTRODE_PLOT_SIZE); % Add a dot at each electrode location
            case 'electrodes'
                plot_mesh_brain(brain_data.([surface_type '_' hemisphere{h}])); 
                plot_electrodes(coordinates, electrode_color_values, 'surface', 'patchsize', ELECTRODE_PATCH_RADIUS); % plot electrodes as colored patches on the surface instead of dots
            case 'parcellation'
                plot_mesh_brain(brain_data.([surface_type '_' hemisphere{h}]), vertex_data.(hemisphere{h}).auto_parcellation.rgb);
                plot_electrodes(coordinates, ELECTRODE_PLOT_COLOR, 'markersize', ELECTRODE_PLOT_SIZE); % Add a dot at each electrode location
            case 'curvature'
                plot_mesh_brain(brain_data.([surface_type '_' hemisphere{h}]), -vertex_data.(hemisphere{h}).vertex_curvature_index);
                plot_electrodes(coordinates, ELECTRODE_PLOT_COLOR, 'markersize', ELECTRODE_PLOT_SIZE); % Add a dot at each electrode location
                colormap gray
            case 'subjects'
                colors = color_list(strcmpi(hemisphere{h},pooled_electrode_table.hemisphere),:);
                curv = -vertex_data.(hemisphere{h}).vertex_curvature_index;
                if strcmp(surface_type,'inflated')
                    curv = sign(curv);
                    caxis([-3 1]);
                end
                plot_mesh_brain(brain_data.([surface_type '_' hemisphere{h}]), curv);
                plot_electrodes(coordinates, colors, 'markersize', ELECTRODE_PLOT_SIZE); % Add a dot at each electrode location
                colormap gray
            otherwise
        end
        % fix lighting and set title
        view([(h-1)*180+90 0]);fix_lighting; 
        title(['Template brain: ' hemisphere{h} ' ' surface_type ' surface']);
    end
end

if strcmp(surface_color_mode,'parcellation')
    show_fs_lookup_table(vertex_data.right.auto_parcellation.lookup_table, vertex_data.right.auto_parcellation.code);
end
