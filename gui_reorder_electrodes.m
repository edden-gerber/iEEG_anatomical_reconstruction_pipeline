function [ electrode_data_out ] = gui_reorder_electrodes( electrode_data )
%REORDER_ELECTRODES_GUI Summary of this function goes here
%   Detailed explanation goes here

% Initialize data
elec_list = electrode_data.electrodes.name;
electrode_data_out = electrode_data;

% Initialize GUI parameters
fig_width = 600;
fig_height = 550;
fig_margin = 20;
horiz_control_size_ratio = 0.6;
button_height = 40;
screen_size = get(0,'ScreenSize');
selected_cells = [];

% Generate GUI
h_figure = figure('Position',[(screen_size(3)-fig_width)/2 (screen_size(4)-fig_height)/2 fig_width fig_height],'Name','Set Electrode Numbering');
h_table = uitable(h_figure,'Data',elec_list,'ColumnName','Electrode ID','ColumnWidth',{fig_width*horiz_control_size_ratio},...
    'Position',[fig_margin fig_margin fig_width*horiz_control_size_ratio-fig_margin*1.5 (fig_height-2*fig_margin)],...
    'CellSelectionCallback',@call_back_table_select);
h_button_add = uicontrol('Style','PushButton','String','Add Electrode Index',...
    'Position',[fig_width*horiz_control_size_ratio+fig_margin*0.5 fig_height-fig_margin-(button_height+fig_margin)*1 fig_width*(1-horiz_control_size_ratio)-fig_margin*1.5 button_height],...
    'Callback',@callback_button_add);
h_button_remove = uicontrol('Style','PushButton','String','Remove Electrode',...
    'Position',[fig_width*horiz_control_size_ratio+fig_margin*0.5 fig_height-fig_margin-(button_height+fig_margin)*2 fig_width*(1-horiz_control_size_ratio)-fig_margin*1.5 button_height],...
    'Callback',@callback_button_remove);
h_button_movedn = uicontrol('Style','PushButton','String','Move Electrodes Down',...
    'Position',[fig_width*horiz_control_size_ratio+fig_margin*0.5 fig_height-fig_margin-(button_height+fig_margin)*3 fig_width*(1-horiz_control_size_ratio)-fig_margin*1.5 button_height],...
    'Callback',@callback_button_movedn);
h_button_moveup = uicontrol('Style','PushButton','String','Move Electrodes Up',...
    'Position',[fig_width*horiz_control_size_ratio+fig_margin*0.5 fig_height-fig_margin-(button_height+fig_margin)*4 fig_width*(1-horiz_control_size_ratio)-fig_margin*1.5 button_height],...
    'Callback',@callback_button_moveup);
h_button_movebt = uicontrol('Style','PushButton','String','Move Electrodes to Bottom',...
    'Position',[fig_width*horiz_control_size_ratio+fig_margin*0.5 fig_height-fig_margin-(button_height+fig_margin)*5 fig_width*(1-horiz_control_size_ratio)-fig_margin*1.5 button_height],...
    'Callback',@callback_button_movebt);
h_button_movetp = uicontrol('Style','PushButton','String','Move Electrodes to Top',...
    'Position',[fig_width*horiz_control_size_ratio+fig_margin*0.5 fig_height-fig_margin-(button_height+fig_margin)*6 fig_width*(1-horiz_control_size_ratio)-fig_margin*1.5 button_height],...
    'Callback',@callback_button_movetp);
h_button_undo = uicontrol('Style','PushButton','String','Undo All Changes',...
    'Position',[fig_width*horiz_control_size_ratio+fig_margin*0.5 fig_height-fig_margin-(button_height+fig_margin)*7 fig_width*(1-horiz_control_size_ratio)-fig_margin*1.5 button_height],...
    'Callback',@callback_button_undo);
h_button_save = uicontrol('Style','PushButton','String','Save And Exit',...
    'Position',[fig_width*horiz_control_size_ratio+fig_margin*0.5 fig_height-fig_margin-(button_height+fig_margin)*8 fig_width*(1-horiz_control_size_ratio)-fig_margin*1.5 button_height],...
    'Callback',@callback_button_save);

% Run
waitfor(h_figure);

% GUI Object functions
    function call_back_table_select(hObject, eventdata)
        selected_cells = eventdata.Indices(:,1);
    end

    function callback_button_add(hObject, eventdata)
        str = inputdlg('Enter electrode name. If it does not match any of the original electrode IDs, this electrode will have no coordinates. Make sure the name is unique: ','New electrode');
        if ~isempty(str) && ~strcmp(str{1},'')
            if isempty(selected_cells)
                new_row_index = 1;
            else
                new_row_index = selected_cells(end)+1;
            end
            d = get(h_table,'Data');
            d((new_row_index+1):(end+1)) = d(new_row_index:end);
            d(new_row_index) = str;
            if size(d,1)==1
                d = d';
            end
            set(h_table,'Data',d);
        end
    end

    function callback_button_remove(hObject, eventdata)
        d = get(h_table,'Data');
        d(selected_cells) = [];
        set(h_table,'Data',d);
    end

    function callback_button_movedn(hObject, eventdata)
        d = get(h_table,'Data');
        nrows = size(d,1);
        if selected_cells(end) < nrows
            d2(1:selected_cells(1)-1) = d(1:selected_cells(1)-1);
            d2(selected_cells(1)) = d(selected_cells(end)+1);
            d2(selected_cells+1) = d(selected_cells);
            d2((selected_cells(end)+2):nrows) = d((selected_cells(end)+2):nrows);
            if size(d2,1)==1
                d2 = d2';
            end
            set(h_table,'Data',d2);
        end
    end

    function callback_button_moveup(hObject, eventdata)
        d = get(h_table,'Data');
        nrows = size(d,1);
        if selected_cells(end) > 1
            d2(1:selected_cells(1)-2) = d(1:selected_cells(1)-2);
            d2(selected_cells-1) = d(selected_cells);
            d2(selected_cells(end)) = d(selected_cells(1)-1);
            d2((selected_cells(end)+1):nrows) = d((selected_cells(end)+1):nrows);
            if size(d2,1)==1
                d2 = d2';
            end
            set(h_table,'Data',d2);
        end
    end

    function callback_button_movebt(hObject, eventdata)
        d = get(h_table,'Data');
        nrows = size(d,1);
        if selected_cells(end) < nrows
            d2(1:selected_cells(1)-1) = d(1:selected_cells(1)-1);
            d2(selected_cells(1):nrows-length(selected_cells)) = d((selected_cells(end)+1):nrows);
            d2(nrows-length(selected_cells)+1:nrows) = d(selected_cells);
            if size(d2,1)==1
                d2 = d2';
            end
            set(h_table,'Data',d2);
        end
    end

    function callback_button_movetp(hObject, eventdata)
        d = get(h_table,'Data');
        nrows = size(d,1);
        if selected_cells(1) > 1
            d2(1:length(selected_cells)) = d(selected_cells);
            d2((length(selected_cells)+1):(selected_cells(end))) = d(1:selected_cells(1)-1);
            d2(selected_cells(end)+1:nrows) = d(selected_cells(end)+1:nrows);
            if size(d2,1)==1
                d2 = d2';
            end
            set(h_table,'Data',d2);
        end
    end

    function callback_button_undo(hObject, eventdata)
        d = elec_list;
        set(h_table,'Data',d);
    end

    function callback_button_save(hObject, eventdata)
        
        if length(unique(get(h_table,'Data'))) < length(get(h_table,'Data'))
            errordlg('There are non-unique electrode names in the list! Fix before saving.');
            return;
        end
        
        elec_list = get(h_table,'Data');
        num_electrodes = length(elec_list);
        
        % Recreate the electrodes table
        cell_data = table2cell(electrode_data.electrodes);
        new_cell_data = {};
        num_variables_in_table = size(cell_data, 2);
        num_coord_systems = num_variables_in_table - 4;
        for i = 1:num_electrodes
            curr_elec_name = elec_list{i};
            idx = find(ismember(electrode_data.electrodes.name,curr_elec_name));
            if isempty(idx)
                new_cell_data(i,:) = [{'External',curr_elec_name,[0 0 0], NaN} mat2cell(repmat([nan nan nan],1,num_coord_systems),1,repmat(3,1,num_coord_systems))];
            else
                new_cell_data(i,:) = cell_data(idx,:);
            end
        end
        electrode_data_out.electrodes = cell2table(new_cell_data,'VariableNames',electrode_data.electrodes.Properties.VariableNames);
        electrode_data_out.num_electrodes = size(electrode_data_out.electrodes,1);
        close(h_figure);
    end

end

