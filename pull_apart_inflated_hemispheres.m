function [ left_mesh, right_mesh ] = pull_apart_inflated_hemispheres(left_mesh, right_mesh, distance )
% by default, inflated hemispheres are centered so they are
% overlapping. this script pulls them apart.
% Based on code by Tal Golan @ Malach Lab

lh_rightmost_point = max(left_mesh.vertices(:,1));
rh_leftmost_point = min(right_mesh.vertices(:,1));

left_mesh.vertices(:,1) = left_mesh.vertices(:,1) - lh_rightmost_point - distance/2;
right_mesh.vertices(:,1) = right_mesh.vertices(:,1) - rh_leftmost_point + distance/2;

end