function display_lattice_tags(carb_chains, carb_index)
% This function prints to the command window the tags of every rubisco 
% in a lattice in a given carboxysome. The formatting is designed to be 
% copied and pasted into the 3D carboxysome visualizer prompt about 
% which chains to highlight.
%
% Inputs
% carb_chains - an array of carboxysome objects with data populated through
%               at least lattice_gen.m
% carb_index - the index of the carboxysome you want the lattice of

    % Access the specified Carboxysome in the data array
    carboxysome = carb_chains(carb_index);

    % tell the user what type the lattice is
    disp(carboxysome.lattice_type);
    
    % get the number of lattices that need to be printed (usually just one)
    num_lattices = length(carboxysome.lattice);

    for i = 1:num_lattices % for each lattice
        this_lattice = carboxysome.lattice{i}; % a list of the chains in the lattice
        num_chains = length(this_lattice); % the number of chains in the lattice

        % Initialize an empty array to store all tags
        all_tags = cell(1, num_chains);        

        for j = 1:num_chains % for each chain
            % Get the tags for the current Rubisco_Chain and concatenate to all_tags
            all_tags{j} = this_lattice(j).tags;
        end
    
        % Manually format and print the output as list of lists
        fprintf('Lattice %i:\n[', i);
        for j = 1:length(all_tags)
            fprintf('[');
            fprintf('%d, ', all_tags{j});  % Print tags in each chain
            fprintf('\b\b], ');  % Remove last space and close the array
        end
        fprintf('\b\b]\n');  % Remove last comma and space, then close the outer array
    end
end