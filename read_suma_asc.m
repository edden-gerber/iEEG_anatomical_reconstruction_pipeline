function [v, f] = read_suma_asc(fname)
%
% Reads asc-file with patch-information
% modified from FreeSurfers's read_asc to read SUMA ascii file
%
% Written by Tal Golan @ Malach Lab

fp = fopen(fname, 'r');

% Dump first line
fgets(fp);

% Nr of vertices and faces
S = zeros(1, 3);
S(1) = 1;
[S(2:3)] = fscanf(fp, '%d', 2);

% Read vertices and its indices
v = fscanf(fp, '%f', S(2)*4);
v = reshape(v, [4 S(2)])';
% remove the forth column
v(:,4)=[]; 

% Read faces and its indices
f = fscanf(fp, '%d', S(3)*4);
f = reshape(f, [4 S(3)])';
f(:,4)=[]; % remove forth column

fclose(fp);
