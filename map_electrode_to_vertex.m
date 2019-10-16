function [ vertex_idx, error ] = map_electrode_to_vertex( coordinates, surface_mesh )
%MAP_ELECTRODES_TO_VERTICES Summary of this function goes here
%   Detailed explanation goes here

v = surface_mesh.vertices;

distances = sqrt((coordinates(1)-v(:,1)).^2 + (coordinates(2)-v(:,2)).^2 + (coordinates(3)-v(:,3)).^2);
[error, vertex_idx] = min(distances);

end


