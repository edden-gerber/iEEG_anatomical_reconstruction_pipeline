function edit_mgrid_file( subj_id, in_file, out_file, electrode_data )
%EDIT_MGRID_FILE Summary of this function goes here
%   Detailed explanation goes here
% 
% Written by Edden M. Gerber 2015

% initialize
curr_point_in_file = 1;
curr_grid = 0;

% first read mgrid file
electrode_data_old = read_mgrid_file(subj_id,in_file);

% read file (as binary, because of the wacky and inconsistent way mgrid 
% files treat end-of-lines): 
fid_i = fopen(in_file);
T = fread(fid_i)';

tline = char(get_next_line);
while ~isempty(tline)
    
    if strfind(tline,'# Electrode Grid')
        curr_grid = str2num(tline(17:end)) + 1; % original format is zero-based
        grid_dimensions = electrode_data_old.grids.dimensions(curr_grid,:);
    end
    
    if strfind(tline,'# Electrode') & ~isempty(str2num(tline(12:end)))
        el_loc = str2num(tline(12:end));
        curr_elec_in_grid = el_loc(1)*grid_dimensions(2) + grid_dimensions(2) - el_loc(2); % Yes, the numbering format is very confusing. But this seems to work. 
%         elec_name = [electrode_data_old.grids.name{curr_grid} '_' num2str(curr_elec_in_grid)];
%         elec_idx = find(strcmp(elec_name,electrode_data.electrodes.name));
    end
    
    if strfind(tline,'#Position')
        tline = char(get_next_line);
        new_pos = electrode_data.grids.coordinates_ijk{curr_grid}(curr_elec_in_grid,:);
        new_line = sprintf(' %.4f %.4f %.4f',new_pos(1),new_pos(2),new_pos(3));
        prev_point_in_file = curr_point_in_file - find(T(curr_point_in_file-2:-1:1)==10,1,'first');
        % insert the new line instead of the old:
        T = [T(1:prev_point_in_file) double(new_line) T(curr_point_in_file-1:end)];
        curr_point_in_file = prev_point_in_file + length(new_line) + 1;
    end
    
    % advance to next line:
    tline = char(get_next_line);
end

fclose(fid_i);

fid_o = fopen(out_file,'wb');
fwrite(fid_o,T);
fclose(fid_o);

% Nested functions: 
function line = get_next_line
    eol = find(T(curr_point_in_file:end)==10,1,'first') + curr_point_in_file - 1;
    if isempty(eol)
        line = [];
    else
        line = T(curr_point_in_file:eol);
        curr_point_in_file = eol + 1;
    end
end

end

