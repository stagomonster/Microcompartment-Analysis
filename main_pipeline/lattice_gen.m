function carboxysome_data = lattice_gen(filename, carboxysome_data, threshold)
% This function takes in the chain linkages and turns them into viable 
% lattices based on the same principles as the chain_maker (using a 
% disjoint set object). Sets whose sizes clear the threshold are called a
% lattice. Each carboxysome is assigned a lattice type and the chains which
% make up the lattice are stored in the carboxysome object.
%
% Inputs
% filename - the name of the file the data originates from. If the only
%            input is filename then it will run every script in the main 
%            pipeline up to here automatically
% carboxysome_data - an array of Carboxysome objects filled with rubisco
%                    data
% threshold - the minimum number of chains that must be connected to form a
%             lattice, defaults to 5
%
% Outputs
% carboxysome_data - an array of Carboxysome objects with the linkages and
%                    max_connections properties populated
%
% lattice_gen.m © 2025 is licensed under CC BY-NC-SA 4.0
    
    % Allow the user to run this script as a starting point, will call
    % the previous scripts in the pipeline and obtain the necessary data
    if nargin == 1
        carboxysome_data = chain_linkages(filename);
    end

    % used to store data about full-type lattices that will be printed to 
    % the command window
    g = zeros(length(carboxysome_data), 3);
    g_last_index = 0;
    % for each carboxysomes
    for carb = carboxysome_data
        chains = carb.chains; % all the chains in the carboxysome

        % creates a disjoint set to find how all the chains are linked to
        % one another
        carb_set = DisjointSet(length(chains));
        % a graph that is later used to map the lattice connections
        adjacencies = zeros(length(chains));
        for chain_link = carb.chain_links
            % for each chain link, connects those chains together
            carb_set.unite(chain_link.I_index, chain_link.J_index);
            adjacencies(chain_link.I_index, chain_link.J_index) = chain_link.distance;
            adjacencies(chain_link.J_index, chain_link.I_index) = chain_link.distance;
        end
        
        % initializes arrays to contain lattice information, leaving room 
        % for there to be multiple lattices in each carboxysome
        largest = zeros(length(carb_set.all)); % number of chains in lattice
        centers = zeros(length(carb_set.all)); % parent element of lattice
        largest_centers_last_index = 0; % index for preallocation
        
        % this is the default threshold that each group of chains must meet
        % to be considered a lattice (user can change this)
        if nargin < 3
            threshold = 5;
        end

        % goes through every chain and sees if it is part of a set that
        % meets the threshold requirements and if so, adds its whole set to
        % a list

        for element = carb_set.all % for each chain
            root = carb_set.locate(element); % find its root
            rank = carb_set.rank(root); % find how many chains the are connected by this root
            if rank >= threshold && ~ismember(root, centers) % if this is a new root and has 5+ chains
                
                largest(largest_centers_last_index+1) = rank; % add the number of chains it has to a list
                centers(largest_centers_last_index+1) = root; % add it to a list of unique lattices
            end
        end

        if largest_centers_last_index == 0
            largest = []; 
            centers = []; 
        else
            largest = largest(1:largest_centers_last_index);
            centers = centers(1:largest_centers_last_index);
        end

        max = 0; % variable to store the maximum number of linkages a chain in the lattice makes
        lattice = cell(1, length(centers)); % cell array to contain lists of chains in each lattice

        % finds the graphical representation of the set and calculates the
        % maximum number of linkages any one chain has in each lattice
        for i = 1:length(centers)
            center = centers(i);
            lattice_indices = carb_set.find_children(center);
            % finds the maximum number of connections in the graph
            max = find_graph(adjacencies, lattice_indices);
            
            % add each Rubisco Chain object to the list of chains in the
            % lattice
            for index = lattice_indices
                lattice{1, i}(end+1) = chains(index);
            end

            % record the maximum number of connections as a property of the
            % Carboxysome object
            if max > carb.max_connections
                carb.max_connections = max;
            end
        end

        % counts how many of the chains are Edge-type
        edges = 0;
        for chain = carb.chains
            if chain.type == ChainType.edge
                edges = edges+1;
            end
        end

        % Determine what type of lattice the carboxysome has

        if length(largest) > 1 % if there is more than one lattice in the carboxysome
            carb.lattice_type = LatticeType.multiple; % multiple-type lattice

        elseif isempty(largest) % if there is no lattice in the carboxysome
            carb.lattice_type = LatticeType.none; % none-type lattice
        
        elseif edges > length(carb.chains) / 2 % if more than half the chains are along the shell
            carb.lattice_type = LatticeType.shell; % shell-type lattice
        
        elseif max >= 6 % if at least one chain in the lattice is connected to at least six others
            carb.lattice_type = LatticeType.full; % full-type lattice
        
        elseif largest >= 4 % if lattice size is >= 4
            carb.lattice_type = LatticeType.incomplete; % incomplete-type lattice
        end

        % displays which carboxysomes have full lattices
        if carb.lattice_type == LatticeType.full
            g(g_last_index+1, 1:3) = [max, carb.carb_index, carb.tomo];
            g_last_index = g_last_index + 1;
        end
        carb.lattice = lattice; % stores the list of chains in the lattice in the carboxysome's lattice property

    end
    if g_last_index == 0
        g = [];
    else
        g = g(1:g_last_index, 1:3);
    end
    % Print information about the full-type lattices
    fprintf('Max Number of Connections in Full-Type Lattices\n');
    max_connections = g(:, 1);
    CB = g(:, 2);
    Tomo = g(:, 3);
    T = table(max_connections, CB, Tomo);
    disp(T);
end

    %% Helper Functions
function max = find_graph(graph, nodes)
% generates the graphical representation and finds the maximum number 
% of connections to any one node
    max = 0;

    % goes through every node in the graph that are connected to each other
    for node = nodes
        % gets the row of the adjacency matrix (graph) that each node represents
        connections = graph(node, :);
        % sees how many of the other nodes in the graph are connected
        % to this node
        connections = length(find(connections));
        % if this is a new record, sets it as the max
        if connections > max
            max = connections;
        end
    end
end