function display_chain_tags(carb_chains, carb_index)
% This function prints to the command window the tags of every rubisco 
% in chains in a given carboxysome. The formatting is designed to be 
% copied and pasted into the 3D carboxysome visualizer prompt about 
% which chains to highlight.
%
% Inputs
% carb_chains - an array of carboxysome objects containing the desired
%               carboxysome
% carb_index - the index of the carboxysome you want the chains of
%
% display_chain_tags.m © 2025 is licensed under CC BY-NC-SA 4.0

    % Access the specified Carboxysome in the data array
    carboxysome = carb_chains([carb_chains.carb_index] == carb_index);
    
    % Initialize an empty array to store all tags
    all_tags = cell(1, length(carboxysome.chains));
    
    % Loop over each chain in the first Carboxysome
    for i = 1:length(carboxysome.chains)
        % Get the tags for the current Rubisco_Chain and concatenate to all_tags
        all_tags{i} = carboxysome.chains(i).tags;
    end
    
    % Manually format and print the output as list of lists
    fprintf('[');
    for i = 1:length(all_tags)
        fprintf('[');
        fprintf('%d, ', all_tags{i});  % Print tags in each chain
        fprintf('\b\b], ');  % Remove last space and close the array
    end
    fprintf('\b\b]\n');  % Remove last comma and space, then close the outer array
end