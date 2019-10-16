function FSLUT = read_fs_color_lookup_table
% written by Tal Golan @ Malach Lab

FSLUT=struct;
[FSLUT.code, FSLUT.name, FSLUT.rgbv] = read_fscolorlut;
FSLUT.name=cellstr(FSLUT.name);
FSLUT.rgbCode=sum(bsxfun(@times,FSLUT.rgbv,256.^(0:(size(FSLUT.rgbv,2)-1))),2);
end