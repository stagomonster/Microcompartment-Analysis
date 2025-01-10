function preprocess_datasets(varargin)
% USAGE: odd and even files: preprocess_datasets(odd_file, even_file, outfile, clean_value, is_bin2)
%        single file: preprocess_datasets(infile, outfile, clean_value, is_bin2) infile must be in bin2
% if input separated into even and odd files, provide odd first
% and even second. otherwise provide input file name as first
% parameter. Then provide the output file name as the next parameter. 
% The next parameter is the size of the clean to perform on the data. 
% Lastly you may provide a true or false value as the last parameter to
% indicate if input data is bin2. If not provided default is bin2=true.
% If bin2 flag is set as false, we assume data is in bin1 and will
% convert it to bin2.

% Preprocess given data sets: if given separated in even and odd, we will
% concatenate into a single table. if there is a nonzero value in columns
% 4-6 (shift), those values need to be added to columns 24-26,
% respectively, and then 4-6 need to be zeroed out. We will also
% transform all values in the table to real values to ensure that no
% imaginary values exist in the data.
%
% preprocess_datasets.m © 2025 is licensed under CC BY-NC-SA 4.0
    
    % concatenate even and odd data if 2 input files are given
    if nargin < 2
        error("You must provide at least 2 parameters: input file name and output file name");
    else
        if nargin == 4
            new_table = dread(varargin{1});
            clean_value = varargin{3};
            % Run a clean to ensure that there are no duplicate particles
            new_table = dpktbl.exclusionPerVolume(new_table, clean_value);
            outfile = varargin{2};
            is_bin2 = varargin{4};
        else
            odd_data = dread(varargin{1});
            clean_value = varargin{4};
            odd_data = dpktbl.exclusionPerVolume(odd_data, clean_value);
            even_data = dread(varargin{2});
            even_data = dpktbl.exclusionPerVolume(even_data, clean_value);
            new_table = [odd_data; even_data];
            outfile = varargin{3};
            is_bin2 = varargin{5};
        end
    end
    
    % Eliminate imaginary values from the table
    new_table = real(new_table);
    
    % adding shifts from 4-6 to 24-26
    new_table(:, 24:26) = new_table(:, 24:26) + new_table(:, 4:6); 
    
    % zero out columns 4-6 after adding shifts
    new_table(:,4:6)=0;
    
    % Transform from bin1 to bin2. Only do this step if data is bin1 and you
    % wish to conver to bin2.
    if ~is_bin2
        new_table(:,24) = new_table(:,24) / 2;
        new_table(:,25) = new_table(:,25) / 2;
        new_table(:,26) = new_table(:,26) / 2;
    end
    
    %% Test for Invalid Entries
    has_nan = any(isnan(new_table(:)));  % Check if any NaN values exist
    has_inf = any(isinf(new_table(:)));  % Check if any Inf values exist
    
    % Display results
    if has_nan
        disp('The data contains NaN values.');
    else
        disp('No NaN values found in the data.');
    end
    
    if has_inf
        disp('The data contains Inf values.');
    else
        disp('No Inf values found in the data.');
    end
    
    
    % write processed table
    fprintf('\nWriting bin2 file: %s\n', outfile);
    dwrite(new_table, outfile);
end