function omdTable = read_data_tbl(filename)
% This function reads a .tbl file and stores it in the Matlab workspace as
% an array.
%
% Inputs
% filename - the name of the .tbl file to be read
%
% Outputs
% omdTable - the name of the array to store the data in
omdTable = dread(filename);
end