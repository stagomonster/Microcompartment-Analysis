function omdTable = read_data_tbl(filename)
% This function reads a .tbl file and stores it in the Matlab workspace as
% an array.
%
% Inputs
% filename - the name of the .tbl file to be read
%
% Outputs
% omdTable - the name of the array to store the data in
%
% read_data_tbl.m © 2025 is licensed under CC BY-NC-SA 4.0

omdTable = dread(filename);
end