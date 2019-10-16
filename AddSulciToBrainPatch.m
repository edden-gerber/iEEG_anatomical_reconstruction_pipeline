function AddSulciToBrainPatch( vertex_map, double_sided_colormap )
%ADDSULCITOBRAINPATCH Summary of this function goes here
%   Detailed explanation goes here

GY_SUL_INTENSITY_DIFF = 0.3;

h_mesh = findobj(gcf,'Type','Patch');
h_mesh = h_mesh(1);

if isempty(h_mesh)
    error('No patch object found in the current figure (i.e. no brain plotted).');
end

vertex_map(vertex_map == min(vertex_map)) = 0;
vertex_map(vertex_map == max(vertex_map)) = 1;

color_map = colormap;

if double_sided_colormap
    gyrus_color = color_map(32,:);
    gyrus_colormap_idx = 33;
    sulcus_colormap_idx = 32;
else
    gyrus_color = color_map(1,:);
    gyrus_colormap_idx = 2;
    sulcus_colormap_idx = 1;
end
color_map(gyrus_colormap_idx,:) = gyrus_color;
color_map(sulcus_colormap_idx,:) = gyrus_color - GY_SUL_INTENSITY_DIFF;
color_map(color_map<0) = 0; % just in case
colormap(color_map);

color_scale = caxis;

exposed_surface = h_mesh.FaceVertexCData == 0;

if double_sided_colormap
    sulcus_color_val = -(color_scale(2)-color_scale(1)) / 64;
    gyrus_color_val = (color_scale(2)-color_scale(1)) / 65;
else
    sulcus_color_val = 0;
    gyrus_color_val = (color_scale(2)-color_scale(1)) / 64;
end

h_mesh.FaceVertexCData(exposed_surface & vertex_map == 0) = sulcus_color_val;
h_mesh.FaceVertexCData(exposed_surface & vertex_map == 1) = gyrus_color_val;

end

