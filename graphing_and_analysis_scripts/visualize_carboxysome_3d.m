function [] = visualize_carboxysome_3d(filename)
% This function conveniently helps the user run the Visualize Carboxysome
% 3D javascript. First it ensures the user has a .tbl file in the directory
% Visualize_Carboxysome_3D/files for them to visualize. If the input
% filename does not exist in that folder, it is written there. Then it runs
% the npx vite command and tells the user to open the webpage using one of
% two methods. First, they can copy and paste the url next to "Local:".
% Second, they could press the "o" key and then "enter". When finished the
% user should press "q" and then "enter".
%
% Inputs:
% filename - the path to the file that should make the carboxysome
%            visualization. Must be a .tbl file.
%
% visualize_carboxysome_3d.m © 2025 is licensed under CC BY-NC-SA 4.0

    %% Make sure the desired file can be loaded by the visualizer
    file_data = dread(filename); % read in the file to a matlab array

    % set the name of the file to remove the path leading to it and keep 
    % just the name of the file
    k = strfind(filename, '/');
    filename = filename(max(k)+1:end);

    % Write the file to the directory Visualize_Carboxysome_3D/files if it
    % does not already exist there
    if ~isfile(filename)
        path = strcat('Visualize_Carboxysome_3D/files/', filename);
        dwrite(file_data, path);
        fprintf('\nThe file has been copied to %s where the Visualization script can find it.\n\n', path)
    end

    % Inform the user of how to open the local webpage
    fprintf('When VITE is ready either copy and paste the link into your browser or hit "o" and then "enter".\n');
    fprintf('Press any key to continue.\n');
    pause;

    % Run the Visualize_Carboxysome_3D script
    system('cd Visualize_Carboxysome_3D && npx vite');
end