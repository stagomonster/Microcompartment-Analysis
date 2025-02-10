function carboxysome_data = convex_hull_and_volume(filename)
% Calculate the volume of carboxysomes from a data set based on a convex
% hull calculation. filename is a .tbl file containing the necessary data,
% carboxysome_data is an array of Carboxysome objects containing calculated
% information.
%
% convex_hull_and_volume.m © 2025 is licensed under CC BY-NC-SA 4.0
    
    % add paths to required helper functions, classes, and data
    addpath('../helper_functions/');
    addpath('../classes/');
    addpath('../data/');

    %% important constants
    % load constants from external file. change constants in file depending on 
    % the dataset being used.
    CONSTANTS = constants();
    rubisco_diameter = CONSTANTS.RUBISCO_DIAMETER_M / CONSTANTS.PIXEL_SIZE; % in pixels

    fprintf("Running with pixel size = %E. If this is not the correct pixel size, please adjust it in the constants file.\n", CONSTANTS.PIXEL_SIZE);
    
    %% Read data from input file
    data = read_data_tbl(filename);
    carboxysomes = get_carboxysomes_indices_from_tomogram(data);
    
    % containerize data in Rubisco objects for easy reference
    data = read_rubisco_objects_from_tomogram(data);
    
    %% Create Carboxysome objects and compute their properties
    carbs = [Carboxysome.empty];
    for carb_index = carboxysomes
        mask = arrayfun(@(x)x.reg == carb_index, data); % This filters out all rubisco except the ones in this carboxysome
        curr_data = data(mask);
        num_rubisco = length(curr_data);
        
        %% get coordinates and do convex hull
        x = [curr_data.x]';
        y = [curr_data.y]';
        z = [curr_data.z]';
    
        [convex_hull, ~] = convhull(x, y, z);
    
        %% find inner and outer rubisco by checking if each rubisco in curr_data is one of the
        % convhull vertices. If it is remove it, if its not keep it.

        verts_x = [x(convex_hull(:, 1)); x(convex_hull(:, 2)); x(convex_hull(:, 3))];
        verts_y = [y(convex_hull(:, 1)); y(convex_hull(:, 2)); y(convex_hull(:, 3))];
        verts_z = [z(convex_hull(:, 1)); z(convex_hull(:, 2)); z(convex_hull(:, 3))];
    
        % filter by unique vertices
        verts = unique([verts_x verts_y verts_z], 'rows');
        inside_rubisco(curr_data) = Rubisco();
        last_index = 0; % index of last valid rubisco in inside_rubiscos
        for i = 1:length(curr_data)
            rubisco = curr_data(i);
            rubisco_coordinates = [rubisco.x rubisco.y rubisco.z];
            if ~ismember(rubisco_coordinates, verts, 'rows')
                rubisco.inside = true; % Set the inside property in the rubisco to be true
                inside_rubisco(last_index + 1) = rubisco;
                last_index = last_index + 1;
            end
        end
        if last_index == 0 % edge case handling
            inside_rubisco = [];
        else
            inside_rubisco = inside_rubisco(1:last_index);
        end
        if ~isempty(inside_rubisco)
            tags_inside = [inside_rubisco.tag];
            num_rubisco_inner = length(tags_inside);
            num_rubisco_outer = length(curr_data) - num_rubisco_inner;
        else
            tags_inside = [];
            num_rubisco_inner = 0;
            num_rubisco_outer = length(curr_data);
        end
    
        %% Compute Normals and Centroids
        normals = zeros(size(convex_hull));
        centroids = zeros(size(convex_hull));
        k = 1;
        for row = convex_hull'
            rubisco1 = curr_data(row(1));
            rubisco2 = curr_data(row(2));
            rubisco3 = curr_data(row(3));
            vec1 = [rubisco2.x - rubisco1.x, rubisco2.y - rubisco1.y, rubisco2.z - rubisco1.z];
            vec2 = [rubisco3.x - rubisco1.x, rubisco3.y - rubisco1.y, rubisco3.z - rubisco1.z];
            vec3 = cross(vec1, vec2); % cross product gives normal for each facet
            normal = vec3/norm(vec3);
            
            % calc centroid for visualization
            centroid = (1/3)*[rubisco1.x + rubisco2.x + rubisco3.x, rubisco1.y + rubisco2.y + rubisco3.y, rubisco1.z + rubisco2.z + rubisco3.z];
            
            normals(k,:) = normal;
            centroids(k, :) = centroid;
            k = k+1;
        end
    
        %% Find average normal vector for each vertex
        for rubisco_idx = 1:length([curr_data])
            if ~ismember(curr_data(rubisco_idx).tag, tags_inside)
                % find all facets with vertex = carb.rubisco(rubisco_idx)
                
                idxs = zeros(length(size(convex_hull, 1))); % preallocate
                last_index = 0; % last used index of idxs
                for k = 1:size(convex_hull, 1)
                    if ismember(rubisco_idx, convex_hull(k, :))
                        idxs(last_index + 1) = k;
                        last_index = last_index + 1;
                    end
                end
                if last_index == 0
                    idxs = [];
                else
                    idxs = idxs(1:last_index);
                end
                
    
                norms = normals(idxs, :);
                ave_vec = sum(norms, 1)/size(norms, 1);
        
                ave_norm = ave_vec/norm(ave_vec);
                curr_data(rubisco_idx).ave_normal = ave_norm;
            end
        end

        %% Compute carboxysome volume
        carb_volume = 0;
        reference = [mean([curr_data.x]), mean([curr_data.y]), mean([curr_data.z])];
        % adjust volume by rubisco diameter and calculate volume using triple product
        for row = convex_hull'
            a = [curr_data(row(1)).x, curr_data(row(1)).y, curr_data(row(1)).z];
            b = [curr_data(row(2)).x, curr_data(row(2)).y, curr_data(row(2)).z];
            c = [curr_data(row(3)).x, curr_data(row(3)).y, curr_data(row(3)).z];
            a_ave_normal = curr_data(row(1)).ave_normal;
            b_ave_normal = curr_data(row(2)).ave_normal;
            c_ave_normal = curr_data(row(3)).ave_normal;
            
            % expansion normal to the surface of the carboxysome
            a = a + sqrt(3) * rubisco_diameter/2 * a_ave_normal - reference;
            b = b + sqrt(3) * rubisco_diameter/2 * b_ave_normal - reference;
            c = c + sqrt(3) * rubisco_diameter/2 * c_ave_normal - reference;
         
            % triple product for volume
            carb_volume = carb_volume + abs(dot(a, cross(b, c)));
        end
        carb_volume = (carb_volume/6)*(CONSTANTS.PIXEL_SIZE^3); % volume in m^3
    
        %% Compute carboxysome concentration
        concentration = 1000*(num_rubisco/carb_volume)/CONSTANTS.AVOGADRO_NUMBER;
    
        %% Create Carboxysome object with all computed properties
        next_index = length(carbs) + 1;
        carbs(next_index) = Carboxysome(carb_index, curr_data, num_rubisco, num_rubisco_inner, ...
            num_rubisco_outer, tags_inside, convex_hull, carb_volume, verts, [x y z], concentration, nan, nan, normals, centroids);
    
        % make carboxysome part of same tomogram as its rubiscos
        carbs(next_index).tomo = curr_data(1).tomo;
    end
    
    %% Compute inner concentration
    for carb = carbs
        x = zeros(length(carb.rubisco));
        y = zeros(length(carb.rubisco));
        z = zeros(length(carb.rubisco));
        ave_normals = [];
        last_index = 0;
        reference = [mean([carb.rubisco.x]), mean([carb.rubisco.y]), mean([carb.rubisco.z])];

        % shrink vertices normal to the surface of the carboxysome
        for rubisco = [carb.rubisco]
            if ~ismember(rubisco.tag, carb.tags_inside)
                new_x = rubisco.x - rubisco_diameter*rubisco.ave_normal(1);
                new_y = rubisco.y - rubisco_diameter*rubisco.ave_normal(2);
                new_z = rubisco.z - rubisco_diameter*rubisco.ave_normal(3);
                x(last_index + 1) = new_x;
                y(last_index + 1) = new_y;
                z(last_index + 1) = new_z;

                ave_normals = [ave_normals; rubisco.ave_normal];
            end
        end
        if last_index == 0
            x = [];
            y = [];
            z = [];
        else
            x = x(1:last_index);
            y = y(1:last_index);
            z = z(1:last_index);
        end
    
        % new hull and volume calculation
        [new_hull, ~] = convhull(x, y, z);
        new_vol = 0;
        for row = new_hull'
            % get the positions of the rubiscos in the new convex hull
            a = [x(row(1)), y(row(1)), z(row(1))];
            b = [x(row(2)), y(row(2)), z(row(2))];
            c = [x(row(3)), y(row(3)), z(row(3))];

            % get the average normal vectors of the rubiscos in the new 
            % convex hull
            a_ave_norm = ave_normals(row(1), :);
            b_ave_norm = ave_normals(row(2), :);
            c_ave_norm = ave_normals(row(3), :);
            
            % expansion normal to the surface of the carboxysome
            a = a + sqrt(3) * rubisco_diameter/2 * a_ave_norm - reference;
            b = b + sqrt(3) * rubisco_diameter/2 * b_ave_norm - reference;
            c = c + sqrt(3) * rubisco_diameter/2 * c_ave_norm - reference;
        
            % triple product for volume
            new_vol = new_vol + abs(dot(a, cross(b, c)));
        end
    
        new_vol = (new_vol/6)*(CONSTANTS.PIXEL_SIZE^3); % new volume in cubic meters
        carb.inner_concentration = 1000*(carb.num_rubisco_inner/new_vol)/CONSTANTS.AVOGADRO_NUMBER;
    end

    % return the carboxysome data to the user
    carboxysome_data = carbs;
end