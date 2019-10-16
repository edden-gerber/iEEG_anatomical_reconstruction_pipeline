function [electrode_data] = read_mgrid_file(subject_id, in_file, reference_image_file)
% This function read a text file of the .mgrid format created by BioImage
% Suite and containing intracranial electrode coordinates. 
% The only input is the mgrid file path. The output is a structure holding
% grid data and electrode coordinates, both separately for each grid as
% well as for all electrodes together. 
% 
% The raw coordinates in the mgrid file are in the voxel space (i.e. IJK
% coordnates = where 0,0,0 is the corner of the cube). To also generate
% coordinates in scanner/anatomical space (i.e. XYZ/RAS coordinates), 
% include the file name for the (CT) image based on which the electrodes 
% were localized (or any image with the same coordinate system). 
%
% Written by Edden M. Gerber 2015

fid = fopen(in_file);

curr_grid = 0;
% initialize output structure
electrode_data.subject_id = subject_id;
electrode_data.num_grids = 0;
electrode_data.num_electrodes = 0;
electrode_data.electrodes = table();
electrode_data.grids = table();
electrode_data.coordinate_transformations = struct;

elec_grid_names = {};
elec_names = {};
elec_colors = [];
elec_depth = [];
elec_hemis = [];
elec_ijk = [];
elec_ras = [];
grid_names = {};
grid_colors = [];
grid_dimensions = [];
grid_spacing = [];
grid_depth = [];
grid_hemis = [];
grid_coordinates_ijk = {};
grid_coordinates_ras = {};

% read file
tline = fgetl(fid);
while ischar(tline)
    if strcmp(tline,'#Number of Grids')
        tline = fgetl(fid);
        electrode_data.num_grids = str2num(tline);
    end
    
    if strfind(tline,'# Electrode Grid')
        curr_grid = str2num(tline(17:end)) + 1; % original format is zero-based
        elec_count = 0;
        for g = 1:(curr_grid-1)
            elec_count = elec_count + grid_dimensions(g,1) * grid_dimensions(g,2);
        end
        grid_depth(curr_grid) = nan; % can't determine from mgrid
        grid_hemis{curr_grid} = nan; % can't determine from mgrid
    end
    
    if strfind(tline,'#Description') & curr_grid > 0
        tline = fgetl(fid);
        grid_names{curr_grid} = tline;
        grid_names{curr_grid}(grid_names{curr_grid}==' ') = '_'; % replace spaces with underscores
        if isempty(grid_names{curr_grid})
            grid_names{curr_grid} = ['Grid' num2str(curr_grid)];
        end
    end
    
    if strfind(tline,'#Dimensions')
        tline = fgetl(fid);
        grid_dimensions(curr_grid,:) = str2num(tline);
        
%         grid_coordinates_ijk{curr_grid} = nan(grid_dimensions(1)*grid_dimensions(2),3);
%         grid_coordinates_ras{curr_grid} = nan(grid_dimensions(1)*grid_dimensions(2),3);
    end
    
    if strfind(tline,'#Color')
        tline = fgetl(fid);
        grid_colors(curr_grid,:) = str2num(tline);
    end
    
    if strfind(tline,'#Electrode Spacing')
        tline = fgetl(fid);
        grid_spacing(curr_grid,:) = str2num(tline);
    end
    
    if strfind(tline,'# Electrode') & ~isempty(str2num(tline(12:end)))
        el_loc = str2num(tline(12:end));
        curr_elec_in_grid = el_loc(1) * grid_dimensions(curr_grid,2) + grid_dimensions(curr_grid,2) - el_loc(2); % Yes, the numbering format is very confusing. But this seems to work. 
        curr_elec = elec_count + curr_elec_in_grid;
    end
    
    if strfind(tline,'#Position')
        tline = fgetl(fid);
        pos = str2num(tline);
        
        % update grid table 
        grid_coordinates_ijk{curr_grid}(curr_elec_in_grid,:) = pos;
        grid_coordinates_ras{curr_grid}(curr_elec_in_grid,:) = [nan nan nan]; % to be filled in later if reference image file is available. 
        
        % update electrodes table
        elec_grid_names{curr_elec} = grid_names{curr_grid};
        elec_names{curr_elec} = [grid_names{curr_grid} '_' num2str(curr_elec_in_grid)];
        elec_colors(curr_elec,:) = grid_colors(curr_grid,:);
        elec_depth(curr_elec) = grid_depth(curr_grid);
        elec_hemis{curr_elec} = {grid_hemis{curr_grid}};
        elec_ijk(curr_elec,:) = pos;
        elec_ras(curr_elec,:) = [nan nan nan]; % to be filled in later if reference image file is available. 
    end
    
    % advance to next line:
    tline = fgetl(fid);
end

fclose(fid);
% Arrange data in tables
electrode_data.num_electrodes = length(elec_names);
electrode_data.electrodes = table(elec_grid_names',elec_names',elec_colors,elec_depth',elec_hemis',elec_ijk,elec_ras,...
    'VariableNames',{'grid','name','color','depth','hemisphere','coordinates_ijk','coordinates_ras'});
electrode_data.grids = table(grid_names',grid_colors,grid_dimensions,grid_spacing,grid_depth',grid_hemis',...
    grid_coordinates_ijk', grid_coordinates_ras', 'VariableNames',{'name','color','dimensions','spacing',...
    'depth','hemisphere','coordinates_ijk','coordinates_ras'});

% Read image file if given
if nargin > 2
    % read nifti file
    img = load_nifti(reference_image_file);
    % get transformation matrix
    nk = size(img.vol,3);
    invert_k_coordinate = [1 0 0 0 ; 0 1 0 0 ; 0 0 -1 nk ; 0 0 0 1]; % the inversion of the K coordinate appears to be particular to BioImage Suite
    vox2ras_transform = affine3d((img.vox2ras * invert_k_coordinate)');
    % apply transformation
    electrode_data.electrodes.coordinates_ras = vox2ras_transform.transformPointsForward(electrode_data.electrodes.coordinates_ijk);
    for g = 1:electrode_data.num_grids
        electrode_data.grids.coordinates_ras{g} = vox2ras_transform.transformPointsForward(electrode_data.grids.coordinates_ijk{g});
    end
    % record transformation matrix
    electrode_data.coordinate_transformations.ijk2ras = vox2ras_transform;
end

end