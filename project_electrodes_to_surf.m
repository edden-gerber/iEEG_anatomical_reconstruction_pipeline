function out_electrodes = project_electrodes_to_surf(electrode_coordinates, surface_mesh, ...
    initial_step_size_displacement, initial_step_size_rotation, error_reduction_threshold, graphical_output)

% Default parameters
INIT_STEP_SIZE_DISPLACEMENT = 0.1;     % Initial step size for movement along the x-y-z axes
INIT_STEP_SIZE_ROTATION = 0.1;         % Initial step size in degrees for rotation around the x-y-z axes (radians)
ERROR_REDUCTION_THRESHOLD = 0.005;     % Minimal amount of error reduction for a step before reducing step size
marker_size = 10;                      % Marker size in electrode scatter plot

% Handle input
if nargin < 2
    error('ERROR: Not enough input arguments');
end
if size(electrode_coordinates,2) ~= 3
    error('ERROR: electrode coordinate matrix size should be Nx3.');
end
if nargin < 3 || isempty(initial_step_size_displacement)
    initial_step_size_displacement = INIT_STEP_SIZE_DISPLACEMENT;
end
if nargin < 4 || isempty(initial_step_size_rotation)
    initial_step_size_rotation = INIT_STEP_SIZE_ROTATION;
end
initial_step_size_rotation = initial_step_size_rotation * pi / 180; % convert to radians
if nargin < 5 || isempty(error_reduction_threshold)
    error_reduction_threshold = ERROR_REDUCTION_THRESHOLD;
end
if nargin < 6 || isempty(graphical_output)
    graphical_output = true;
end

if graphical_output
    % plot cortical outer surface
    subplot(1,2,1);
    hold off
    plot_mesh_brain(surface_mesh,[320 30],[],0.3);
    hold all;
    scatter3(electrode_coordinates(:,1),electrode_coordinates(:,2),electrode_coordinates(:,3),marker_size,'filled','b');
    h_scat_new = [];
    drawnow;
    set(gcf,'WindowStyle','Docked');
    rotate3d on
end

% Set parameters
num_electrodes = size(electrode_coordinates,1);

% Inline functions for rotation matrices
Rx = @(theta) [1 0 0 ; 0 cos(theta) -sin(theta) ; 0 sin(theta) cos(theta)];
Ry = @(theta) [cos(theta) 0 sin(theta) ; 0 1 0 ; -sin(theta) 0 cos(theta)];
Rz = @(theta) [cos(theta) -sin(theta) 0 ; sin(theta) cos(theta) 0 ; 0 0 1];
% Inline functions for error signal (square of distance from nearest
% surface vertex, averaged across electrodes). 
err = @(elecs,targets) sum(sum((elecs - targets).^2,2)) / num_electrodes;

% Initialize
curr_step = 0;
prj_electrodes = electrode_coordinates;
stop_disp = false;
surface_vert = surface_mesh.vertices;
if num_electrodes == 1 % a single electrode cannot be rotated
    stop_rot = true;
else
    stop_rot = false;
end
step_size_displacement = initial_step_size_displacement;
step_size_rotation = initial_step_size_rotation;
disp_gradient = zeros(3,1); 
rot_gradient = zeros(3,1); 
prev_error = inf;

% find closest vertex for each electrode
nearest_vert = zeros(num_electrodes,3);
for k=1:num_electrodes
    npls=[surface_vert(:,1)-prj_electrodes(k,1) surface_vert(:,2)-prj_electrodes(k,2) surface_vert(:,3)-prj_electrodes(k,3)]; %x,y,z distaces
    npls_dist=sqrt(sum((npls).^2,2));
    [~,idx] = min(npls_dist);
    nearest_vert(k,:) = surface_vert(idx,:);
end

% Run gradient descent
while ~stop_rot || ~stop_disp % Stop only when both rotation and displacement have reached their stopping condition. 
    
    curr_step = curr_step + 1;
    
    % Get current error signal value
    curr_error = err(prj_electrodes,nearest_vert);
    
    % Run one iteration of displacement along the x-y-z axes
    if ~stop_disp
        % displace on x axis
        step_pos = bsxfun(@plus,prj_electrodes,[step_size_displacement 0 0]);
        disp_gradient(1) = curr_error - err(step_pos,nearest_vert);
        % displace on y axis
        step_pos = bsxfun(@plus,prj_electrodes,[0 step_size_displacement 0]);
        disp_gradient(2) = curr_error - err(step_pos,nearest_vert);
        % displace on z axis
        step_pos = bsxfun(@plus,prj_electrodes,[0 0 step_size_displacement]);
        disp_gradient(3) = curr_error - err(step_pos,nearest_vert);
        % normalize the gradient vector
        disp_gradient = disp_gradient / norm(disp_gradient);
        % displace grid one step size along the gradient vector
        prj_electrodes = bsxfun(@plus,prj_electrodes,disp_gradient'*step_size_displacement);
        
        % compute new error
        curr_error = err(prj_electrodes,nearest_vert);
        error_decrease = prev_error - curr_error;
        prev_error = curr_error;
        % If the error has not decreased by more than min_error_difference,
        % (or it has increased, meaning the local minimum has been
        % crossed), decrease the step size. 
        if error_decrease < error_reduction_threshold
            step_size_displacement = step_size_displacement / 2;
            if step_size_displacement < (initial_step_size_displacement / 256)
                stop_disp = true;
            end
        end
    end
    
    % Run one iteration of rotation around x-y-z axes
    if ~stop_rot
        % For rotation, compute the center coordiate of the grid and the
        % electrodes' relative locations
        center = [mean(prj_electrodes(:,1)) mean(prj_electrodes(:,2)) mean(prj_electrodes(:,3))];
        relative_pos = bsxfun(@minus,prj_electrodes,center);
        % rotation on x axis
        step_pos = bsxfun(@plus,center,(double(Rx(step_size_rotation))*relative_pos')');
        rot_gradient(1) = curr_error - err(step_pos,nearest_vert);
        % rotation on y axis
        step_pos = bsxfun(@plus,center,(double(Ry(step_size_rotation))*relative_pos')');
        rot_gradient(2) = curr_error - err(step_pos,nearest_vert);
        % rotation on x axis
        step_pos = bsxfun(@plus,center,(double(Rz(step_size_rotation))*relative_pos')');
        rot_gradient(3) = curr_error - err(step_pos,nearest_vert);
        % normalize the gradient vector
        rot_gradient = rot_gradient / norm(rot_gradient);

        % rotate grid one step size along the gradient vector
        rot_x = double(Rx(step_size_rotation*rot_gradient(1)));
        rot_y = double(Ry(step_size_rotation*rot_gradient(2)));
        rot_z = double(Rz(step_size_rotation*rot_gradient(3)));
        rot_matrix = rot_x*rot_y*rot_z;
        % to rotate, center the coordinates around zero, left-multiply by
        % the rotation matrix, and finally move the center back to where it
        % was: 
        prj_electrodes = bsxfun(@plus,center,(rot_matrix*relative_pos')');
        
        % compute new error
        curr_error = err(prj_electrodes,nearest_vert);
        error_decrease = prev_error - curr_error;
        prev_error = curr_error;
        % If the error has not decreased by more than min_error_difference,
        % (or it has increased, meaning the local minimum has been
        % crossed), decrease the step size. 
        if error_decrease < error_reduction_threshold
            step_size_rotation = step_size_rotation / 2;
            if step_size_rotation < (initial_step_size_rotation / 256)
                stop_rot = true;
            end
        end
    end
    
    % find closest vertex for each electrode at each step
    nearest_vert = zeros(num_electrodes,3);
    for k=1:num_electrodes
        npls=[surface_vert(:,1)-prj_electrodes(k,1) surface_vert(:,2)-prj_electrodes(k,2) surface_vert(:,3)-prj_electrodes(k,3)]; %x,y,z distaces
        npls_dist=sqrt(sum((npls).^2,2));
        [~,idx] = min(npls_dist);
        nearest_vert(k,:) = surface_vert(idx,:);
    end
        
    if graphical_output
    % plot electrodes
        step_res = 1;
        if mod(curr_step,step_res) == 0
            subplot(1,2,1);
            delete(h_scat_new);
            h_scat_new = scatter3(prj_electrodes(:,1),prj_electrodes(:,2),prj_electrodes(:,3),marker_size,'filled','r');
        end
        track_error(curr_step) = curr_error;
        track_step_size_r(curr_step) = step_size_rotation;
        track_step_size_d(curr_step) = step_size_displacement;
        subplot(1,2,2);
        hold off
        plot(track_error/track_error(1));
        hold all
        plot(track_step_size_r/track_step_size_r(1));
        plot(track_step_size_d/track_step_size_d(1));
        ylim([0 1]);
        set(gca,'ytick',{});
        legend('Error signal %','Rotation step size %','Displacement step size %');
        xlabel(['Rotation step size: ' num2str(step_size_rotation*180/pi) 'deg. Displacement step size: ' num2str(step_size_displacement) 'mm. Error: ' num2str(curr_error)]);
        drawnow;
    end
end

% Return output
if graphical_output
    subplot(1,2,2);
    xlabel('Finished projecting electrodes.');
end
out_electrodes = prj_electrodes;

end

