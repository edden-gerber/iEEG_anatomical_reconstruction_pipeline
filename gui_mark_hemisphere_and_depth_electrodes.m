function [ electrode_data_out ] = gui_mark_hemisphere_and_depth_electrodes( electrode_data )
%MARK_DEPTH_ELECTRODES_GUI Summary of this function goes here
%   Detailed explanation goes here


% Initialize data
grid_list = electrode_data.grids.name;
depth_list = electrode_data.grids.depth;
depth_list(isnan(depth_list)) = false;
depth_list = logical(depth_list);
hemisphere_list = electrode_data.grids.hemisphere;

table_data = [grid_list hemisphere_list mat2cell(depth_list,ones(length(depth_list),1),1)];
electrode_data_out = electrode_data;

% Initialize GUI parameters
fig_width = 600;
fig_height = 200;
fig_margin = 20;
horiz_control_size_ratio = 0.7;
button_height = 40;
screen_size = get(0,'ScreenSize');
selected_cells = [];

% Get initial guesses for hemisphere value
for g = 1:electrode_data.num_grids
    if isnan(electrode_data.grids.hemisphere{g})
        grid_elec_idxs = find(strcmp(electrode_data.grids.name{g},electrode_data.electrodes.grid));
        grid_coordinates = electrode_data.electrodes.coordinates_ras(grid_elec_idxs,:);
        mean_x_loc = mean(grid_coordinates(:,1)); 
        if mean_x_loc > 0
            table_data{g,2} = 'Right';
        else
            table_data{g,2} = 'Left';
        end
    end
end

% Generate GUI
h_figure = figure('Position',[(screen_size(3)-fig_width)/2 (screen_size(4)-fig_height)/2 fig_width fig_height],'Name','Mark Depth Electrodes');
h_table = uitable(h_figure,'Data',table_data,'ColumnName',{'Grid ID','Hemisphere','Depth Electrodes'},...
    'ColumnEditable',[false true true],'ColumnFormat',{'char',{'Left','Right'},'logical'},...
    'Position',[fig_margin fig_margin fig_width*horiz_control_size_ratio-fig_margin*1.5 (fig_height-2*fig_margin)]);
table_pos = get(h_table,'Position');
set(h_table,'ColumnWidth',{table_pos(3)*0.3 table_pos(3)*0.3 table_pos(3)*0.3});
h_button_save = uicontrol('Style','PushButton','String','Save And Exit',...
    'Position',[fig_width*horiz_control_size_ratio+fig_margin*0.5 fig_height-(fig_margin+button_height)*1 fig_width*(1-horiz_control_size_ratio)-fig_margin*1.5 button_height],...
    'Callback',@callback_button_save);

% Run
waitfor(h_figure);

% GUI Object functions
    function callback_button_save(hObject, eventdata)
        table_data = get(h_table,'Data');;
        % Update grid data
        electrode_data_out.grids.hemisphere = table_data(:,2);
        electrode_data_out.grids.depth = cell2mat(table_data(:,3));
        % Update electrode data
        for g = 1:electrode_data.num_grids
            grid_elec_idxs = find(strcmp(electrode_data.grids.name{g},electrode_data.electrodes.grid));
            electrode_data_out.electrodes.depth(grid_elec_idxs) = electrode_data_out.grids.depth(g);
            electrode_data_out.electrodes.hemisphere(grid_elec_idxs) = {electrode_data_out.grids.hemisphere{g}};
        end
        close(h_figure);
        drawnow;
        
    end
end

