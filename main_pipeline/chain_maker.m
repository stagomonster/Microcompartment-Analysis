function carboxysome_data = chain_maker(filename, carboxysome_data, min_chain_length)
% this is the function following linkage formation in the pipeline. it
% takes all viable linkages and uses a data structure called a disjoint set
% to find all rubiscos that are connected to one another, which are defined
% to be rubisco chains. the min_chain_length variable represents the
% minimum length required for a set of linked rubiscos to be considered a
% chain and is set to a default of 3. After rubiscos are grouped into
% chains, they are put in order in the chains. Which rubiscos are above and
% below them are stored in the rubisco object. Consecutive rubiscos without
% a linkage will return an error and will skip above/below designation for
% the chain. These are likely bad data, so the particles with the 
% lowest cc involved in each error are deleted from the dataset. The 
% original input table is unaltered, but a new table with bad particles 
% deleted is created in the working directory.
%
% Inputs
% filename - the name of data file to be read if pre-computed 
%            carboxysome_data is not provided.
% carboxysome_data - an array of Carboxysome objects already containing
%                    computed data and pairs from linkages.m.
% min_chain_length - the minimum length of linked rubiscos to be considered
%                    a chain
%
% Outputs
% carboxysome_data - an array of Carboxysome objects with chains populated
%
% chain_maker.m © 2025 is licensed under CC BY-NC-SA 4.0

    tag = 1; % the tag of the chain
    
    % Allow the user to run this script as a starting point, will call
    % the previous scripts in the pipeline and obtain the necessary data
    if nargin == 1
        min_chain_length = 3;
        carboxysome_data = linkages(filename);
    elseif nargin < 3
        min_chain_length = 3;
    end

    original_table = dread(filename);
    table_to_edit = original_table;
    
    for carb = carboxysome_data % for each carboxysome
        % uses a DisjointSet to turn the links data into longer chain data
        carb_set = DisjointSet(length(carb.rubisco));
        for link = carb.links
            carb_set.unite(link.I_index, link.J_index);
        end
    
        % translates the disjoint set into a set of arrays, each of which
        % contains a chain
        chains = carb_set.chain_builder();
        rubisco_chains = Rubisco_Chain.empty;
        rubiscos = Rubisco.empty;
        id = 1;
        
        % generates the chains using the output from the previous section
        for i = 1:length(chains) % for each chain in the carboxysome
            tags = uint32([]);
            chain = chains{i};
            % does not allow for chains of length smaller than the user specifies
            if length(chain) < min_chain_length; continue; end
            % finds whether every rubisco on the chain is floating in the
            % carboxysome, attached to the edge, or growing along it
            type = ChainType.free;
            edges = 0;
            % loops through all the indices that appear in the chain
            for index = chain
                rubisco = carb.rubisco(index);
                tags(end+1) = rubisco.tag;
                rubiscos(end+1) = rubisco;
    
                % this will likely require fine tuning in the future, but it is 
                % meant to define chains as edge chains, just attached to the 
                % shell, or free floating within the carboxysome
                if ~rubisco.inside && type == ChainType.free
                    type = ChainType.attached;
                    edges = 1;
                elseif ~rubisco.inside && type == ChainType.attached && edges < length(chain)/3
                    edges = edges+1;
                elseif ~rubisco.inside && type == ChainType.attached && edges >= length(chain)/3
                    type = ChainType.edge;
                end
            end

            % create an instance of a rubisco chain object
            rubisco_chain = Rubisco_Chain(carb.carb_index, id, tag, tags, uint32(chain), type);
            table_to_edit = fill_data(carb, rubisco_chain, table_to_edit);
            rubisco_chains(end+1) = rubisco_chain;
            id = id+1;
            tag = tag+1;
        end
    
        % Store all the rubisco chain objects from this carboxysome
        carb.chains = rubisco_chains;
    end

    if size(table_to_edit, 1) < size(original_table, 1) % if particles were deleted
        dwrite(table_to_edit, '../data/problems_deleted.tbl'); % write the new file without bad particles

        % tell the user about particle deletion and recommend they start over
        fprintf('%i particles were deleted that caused a rubisco to have 2 linkages above (or below) it.\n', size(original_table, 1) - size(table_to_edit, 1));
        fprintf('Particle deletion was done to the particle with the lowest cc value among those involved in the error.\n')
        fprintf('A file called "problems_deleted.tbl" was created in the "data" directory without these problem particles.\n');
        fprintf('The original table is unaltered, but it is not recommended to proceed with the current data.\n');
        fprintf('Future analysis depends upon knowing which rubiscos are above and below each other in a chain.\n');
        fprintf('It is recommended to start the pipeline over with the new .tbl file.\n');
    end
end


    %% Helper Functions
function [new_table] = fill_data(carb, chain, table)
% fills in the missing data after chains are formed
    chain.S = op(carb, chain); % calculates S tensor
    chain.eigenvalues = eigs(chain.S); % calculates eigenvalues of S tensor
    chain.S_val = (3/2)*max(chain.eigenvalues); % calculates scalar S value
    chain.centroid = calc_center(carb, chain); % calculates chain centroid
    chain.average_vector = calc_ave_vector(carb, chain); % calculates average vector of rubisco in chain
    [chain.xy_angle, chain.z_angle] = calc_angles(chain); % calculates angles chain makes to xyz directions
    % organizes the chains in the order the tags/indices occur in
    if chain.length > 1
        furthest = find_furthest(carb, chain);
        [chain.tags, chain.indices] = find_order(carb, chain, furthest);
    end
    new_table = align_rubisco_axes_in_chains(carb, chain, table); % points all rubisco axes in the same direction
end

function S = op(carb, chain)
% compute average S tensor for a rubisco chain
    len = length(carb.rubisco(1).vector);
    rubiscos = carb.rubisco;
    
    S = zeros(len, len);
    for j = chain.indices
        S = S + calc_S(rubiscos(j).vector);
    end
    S = S/length(chain.length);
end

function S = calc_S(vector)
% compute the S alignment tensor for a rubisco
    S = vector' * vector;
    S = S - (1/3)*[1 0 0; 0 1 0; 0 0 1];
end

function center = calc_center(carb, chain)
% compute the location of the center of a rubisco chain
    center = [0 0 0];
    rubiscos = carb.rubisco;
    for j = chain.indices
        center = center + [rubiscos(j).x, rubiscos(j).y, rubiscos(j).z];
    end
    center = center/length(chain.indices);
end

function ave_vec = calc_ave_vector(carb, chain)
% compute the average vector for a rubisco chain using the average of its
% rubisco orientations

    %% calculate mean of orientations method
    vec = [0 0 0];
    rubiscos = carb.rubisco;
    for j = chain.indices % for each rubisco in the chain
        % add its orientation vector to the sum of all orientations in the
        % chain, making sure they all point the same way
        tent_vec1 = vec + rubiscos(j).vector;
        tent_vec2 = vec - rubiscos(j).vector;
        if norm(tent_vec1) >= norm(tent_vec2); vec = tent_vec1; else; vec = tent_vec2; end
    end

    % the average vector of all the rubisco orienation vectors in the chain
    ave_vec = vec / length(chain.indices);
end

function [xy_angle, z_angle] = calc_angles(chain)
    % finds the angle that the chain's average vector forms in the xy plane
    % and with the xy plane in degrees
    chain_vec = chain.average_vector;

    x = chain_vec(1);
    y = chain_vec(2);
    z = chain_vec(3);

    % Calculate the angle it makes in the xy plane
    xy_angle = atan(y/x) * 180/pi;

    % Ensure the angle is between 0 and 360 degrees
    if x < 0
        xy_angle = xy_angle + 180;
    elseif y < 0
        xy_angle = xy_angle + 360;
    end

    % Calculate the angle it makes with the xy plane
    xy_dist = sqrt(x^2+y^2);
    z_angle = atan(z/xy_dist) * 180/pi; % angle is between -90 and 90 degrees

    % Store the angles in the rubisco chain object
    chain.xy_angle = xy_angle;
    chain.z_angle = z_angle;
end

function furthest = find_furthest(carb, chain)
% locates the rubisco that is the furthest from the centroid
    centroid = chain.centroid;
    furthest = centroid;
    max_dist = 0;

    for rubisco_index = chain.indices
        rubisco = carb.rubisco(rubisco_index);
        dist = norm(centroid - [rubisco.x rubisco.y rubisco.z]); % calculates distance from centroid to rubisco
        if dist > max_dist
            furthest = rubisco;
            max_dist = dist;
        end
    end
end

function [ordered_tags, ordered_indices] = find_order(carb, chain, point)
% orders the rubiscos based on how far they are from the furthest
% rubisco from the centroid, since that should be an end rubisco
    distances = [];
    rubiscos = carb.rubisco;
    for rubisco_index = chain.indices
        rubisco = rubiscos(rubisco_index);
        distance = norm([point.x point.y point.z] - [rubisco.x rubisco.y rubisco.z]);
        % stores all the distances of the rubiscos from the end point
        distances(end+1, :) = [distance, rubisco.tag, rubisco.index];
    end
    
    % uses bubble sorting to order the rubiscos from least to greatest
    % distance from one of the endpoints
    n = length(distances(:, 1));
    for k = 1:n-1
        for j = 1:n-k
            if distances(j, 1) > distances(j+1, 1)
                temp = distances(j, :);
                distances(j, :) = distances(j+1, :);
                distances(j+1, :) = temp;
            end
        end
    end

    ordered_tags = [];
    ordered_indices = [];

    for distance = distances(:, 2:3)'
        ordered_tags(end+1) = distance(1);
        ordered_indices(end+1) = distance(2);
    end
end

function [orig_table] = align_rubisco_axes_in_chains(carb, chain, orig_table)
% This function makes every rubisco axis in a chain point the same way
% and then calculates which rubiscos are above and below each other.
    problem_particles = [];

    for rubisco = 1:length(chain.indices) - 1 % for every linkage in the chain
    
        % Get the position and vector of a rubisco in the chain
        source_x = carb.rubisco(chain.indices(rubisco)).x;
        source_y = carb.rubisco(chain.indices(rubisco)).y;
        source_z = carb.rubisco(chain.indices(rubisco)).z;
        source_vector = carb.rubisco(chain.indices(rubisco)).vector;

        % Get the position and vector of the next rubisco in the chain
        target_x = carb.rubisco(chain.indices(rubisco + 1)).x;
        target_y = carb.rubisco(chain.indices(rubisco + 1)).y;
        target_z = carb.rubisco(chain.indices(rubisco + 1)).z;
        target_vector = carb.rubisco(chain.indices(rubisco + 1)).vector;

        % calculate the vector pointing from the first rubisco to the next
        displacement = [target_x target_y target_z] - [source_x source_y source_z];

        % If first rubisco vector and displacement vector point opposite
        % directions, permanently flip the rubisco vector
        if dot(source_vector, displacement) < 0
            carb.rubisco(chain.indices(rubisco)).vector = -1 * source_vector;
            carb.rubisco(chain.indices(rubisco)).narot = -1 * carb.rubisco(chain.indices(rubisco)).narot;
        end

        % If secnd rubisco vector and displacement vector point opposite
        % directions, permanently flip the rubisco vector
        if dot(target_vector, displacement) < 0
            carb.rubisco(chain.indices(rubisco + 1)).vector = -1 * target_vector;
            carb.rubisco(chain.indices(rubisco + 1)).narot = -1 * carb.rubisco(chain.indices(rubisco + 1)).narot;
        end

        % Find the linkage object that contains this linkage. A forward
        % pair has the first rubisco as the I_index while a reverse pair
        % has the first rubisco as the J_index
        forward_pair = find([carb.links.I_index] == chain.indices(rubisco) & [carb.links.J_index] == chain.indices(rubisco + 1));
        reverse_pair = find([carb.links.J_index] == chain.indices(rubisco) & [carb.links.I_index] == chain.indices(rubisco + 1));
        
        if any([forward_pair reverse_pair]) % if the rubiscos have a linkage
            carb.rubisco(chain.indices(rubisco)).rubisco_above_me = chain.tags(rubisco + 1); % this rubisco is below a rubisco
            carb.rubisco(chain.indices(rubisco + 1)).rubisco_below_me = chain.tags(rubisco); % this rubisco is above a rubisco
        else % else print an error message
            fprintf('Error on carboxysome %i chain %i. ', carb.carb_index, chain.index);
            fprintf('%i or %i may be a bad data point. Skipping above/below designation.\n', chain.tags(rubisco), chain.tags(rubisco + 1));
            fprintf('Chain Contains Rubisco Tags: [[%i', chain.tags(1));
            for i = 2:length(chain.tags)
                fprintf(',%i', chain.tags(i));
            end
            fprintf(']]\n');
            problem_particles = [problem_particles; chain.indices(rubisco), chain.indices(rubisco + 1)];
        end
    end
    
    % delete any problem particles from data table if there are any
    if ~isempty(problem_particles)
        orig_table = delete_problems(problem_particles, carb, orig_table);
    end
end

function orig_table = delete_problems(problem_particles, carb, orig_table)
% This function deletes the particle with the worst cross correlation among
% particles that may be false positives and returns the table now without
% any bad data points.
    different_problems = cell(1, size(problem_particles, 1)); % make a cell array to hold individual issue locations
    
    % assume each issue is independent at first
    for i = 1:size(problem_particles, 1)
        different_problems{1, i} = problem_particles(i, :);
    end
   
    % find issues that share a particle and combine them into one issue
    for i = 1:size(problem_particles, 1)
        if sum(problem_particles(1:i-1, :) == problem_particles(i, 1) | problem_particles(1:i-1, :) == problem_particles(i, 2), 'all') > 0
            [row, ~] = find(problem_particles(1:i-1, :) == problem_particles(i, 1) | problem_particles(1:i-1, :) == problem_particles(i, 2));
            different_problems{row} = unique([different_problems{row}, different_problems{i}]);
            different_problems{:, i} = nan;
        end
    end
    
    % find the particle with minimum cross correlation for each issue and
    % mark it for deletion
    for i = 1:size(different_problems, 2) % for each issue
        if ~isnan(different_problems{1, i})
            min_cc = inf;
            for j = 1:length(different_problems{1, i}) % for each possible culprit rubisco
                cross_correlation = carb.rubisco(different_problems{1,i}(j)).cc; % get the rubisco's cc
                if cross_correlation < min_cc % if this new cc is smaller than the previous minimum
                    min_cc = cross_correlation; % keep track of the minimum cc for each issue
                    particle_to_delete = different_problems{1,i}(j); % mark this particle for deletion
                end
            end
        end
        % tell the user what particle is getting deleted
        fprintf('Deleting rubisco with lowest cc value: Tag is %i\n\n', carb.rubisco(particle_to_delete).tag);
        % delete the particle from the data table
        orig_table = orig_table(orig_table(:, 1) ~= carb.rubisco(particle_to_delete).tag, :);
    end
end