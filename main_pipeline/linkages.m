function carboxysome_data = linkages(filename, carboxysome_data, min_distance_factor, max_distance_factor, radius_factor, max_angle, mutual)
% Compute the linkages between rubisco in each carboxysome. We define a
% linkage as a valid Pair of Rubisco. All possible pairs are created
% beforehand by local_global_alignment, and this script will filter
% through them to select only valid bindings. There are 3 (optionally
% 4) filters. To create a linkage, target rubiscos must be within a
% cylindrical space at one end of the source rubisco and aligned with 
% the vector of the source rubisco. We define that cylindrical space 
% with "distance" (the height of the cylinder), "angle" (the angle 
% between the orientation vectors of the rubiscos), and "radius" (the
% target rubisco's distance from the axis of the source rubisco). The
% optional 4th parameter, "mutual" requires that all 3 filters be
% passed in both directions, that is, it doesn't matter which rubisco
% is the "source" and which is the "target", they both see the other as
% a linkage. This script is not built to handle cases where multiple
% rubiscos fit into the same cylinder. For this reason we recommend
% using the helper function calc_cylinder_radius.m to calculate what
% maximum radius_factor you should use with this program to prevent
% multiple particles in the same cylinder.
%
% Inputs
% filename - the name of data file to be read if pre-computed 
%            carboxysome_data is not provided.
% carboxysome_data - an array of Carboxysome objects already containing
%                    computed data and pairs from 
%                    local_global_alignment.m. 
% min_distance_factor - a multiple of the rubisco diameter to search 
%                       along the source rubisco axis, so if 
%                       min_distance_factor = 1.3, then the bottom of
%                       the search cylinder will be 1.3 rubisco
%                       diameters from the center of the source rubisco
% max_distance_factor - a multiple of the rubisco diameter to search 
%                       along the source rubisco axis, so if 
%                       max_distance_factor = 1.3, then the top of
%                       the search cylinder will be 1.3 rubisco
%                       diameters from the center of the source rubisco
% radius_factor - the multiple of the rubisco diameter to search 
%                 outward from the axis. If radius_factor = 0.7 then 
%                 the search cylinder will have a radius of 0.7 rubisco
%                 diameters
% max_angle - the maximum allowable angle between the rubiscos in a
%             linkage in degrees.
% mutual - a boolean that requires the linkages to reference one 
%          another if it is true (encouraged if radius_factor is larger
%          than 0.45).
%
% Outputs
% carboxysome_data - an array of Carboxysome objects with rubisco linkages
%                    populated
%
% linkages.m © 2025 is licensed under CC BY-NC-SA 4.0
    
    % Allow the user to run this script as a starting point, will call
    % the previous scripts in the pipeline and obtain the necessary data
    if nargin == 1
        carboxysome_data = local_global_alignment(filename);
    end

    % predetermined parameters to define the cylinder:
    % max dist factor 1.4933
    % min dist factor 0.7947
    % "Tight" parameters: radius factor = 0.31, max angle = 25deg
    % "Pivot" parameters: radius_factor = 0.45, max_angle = 50deg

    % NOTE: For Rubisco, 0.45 should be the maximum allowed radius_factor
    % since this was calculated to be the maximum that would not permit to
    % have multiple bindings at the same binding site
    
    if nargin < 3 % Set values to use for search cylinder if not specified in function call
        max_distance_factor = 1.4933;
        min_distance_factor = 0.7947;
        mutual = 1;

        % Ask user which preset parameter combination to use
        userInput = input('Enter 0 to use "Tight" binding parameters, or 1 to use "Pivot" binding parameters: ');
        if userInput == 0
            disp("Using Tight parameters radius_factor = 0.31, max_angle = 25deg.");
            radius_factor = 0.31;
            max_angle = 25;
        elseif userInput == 1
            disp("Using Pivot parameters radius_factor = 0.45, max_angle = 50deg.");
            radius_factor = 0.45;
            max_angle = 50;
        else
            error('Invalid input. Please enter 0 or 1. Exiting program...');
        end
    end
    
    %% import useful data
    % load constants from external file. change in file depending on dataset.
    CONSTANTS = constants();
    rubisco_diameter = CONSTANTS.RUBISCO_DIAMETER_M / CONSTANTS.PIXEL_SIZE;
    max_distance = rubisco_diameter*max_distance_factor;
    min_distance = rubisco_diameter*min_distance_factor;
    search_radius = rubisco_diameter*radius_factor;
    
    % search through each carboxysome
    for carb = carboxysome_data
        disp(carb.carb_index);
        
        % Initialize an array that will store all the valid linkages
        links = Rubisco_Pair.empty;

        % goes through all pairs in the carboxysome
        for pair = carb.rubisco_pairs
            valid = true;

            % filter by distance along the source axis, upper and lower 
            % limits determined from histogram analysis. forms a flat space
            % between two planes that are perpendicular to the source axis
            if abs(pair.projection) > max_distance || abs(pair.projection) < min_distance
                valid = false;
            end

            % Filter by angle between the axes of both rubisco
            if abs(pair.angle) >= max_angle && abs(pair.angle - 180) >= max_angle
                continue;
            end

            % radius is the linear distance from the target rubisco's center 
            % to the closest point in the central axis of the source rubisco
            radius = sqrt(pair.distance^2 - pair.projection^2);
            if radius >= search_radius
                valid = false;
            end

            source = carb.rubisco(pair.I_index);
            target = carb.rubisco(pair.J_index);

            if (~valid && ~mutual) || (valid && mutual)
                % ensures that each target also references the source as a viable link

                % Calculate the projection in the opposite direction
                inv_projection = calc_projection(target, source);

                % Filter by distance along target axis (differs from source)
                if abs(inv_projection) > max_distance || abs(inv_projection) < min_distance
                    continue;
                end

                % Angle will be the same

                % Calculate radius in opposite direction
                inv_radius = sqrt(pair.distance^2 - inv_projection^2);
                
                % Filter by inv radius (as defined above; differs from source)
                if inv_radius >= search_radius
                    continue;
                end
            elseif ~valid
                continue;
            end

            % All checks passed, linkage is valid
            
            source_cords = [source.x source.y source.z];
            target_cords = [target.x target.y target.z];
            % vector pointing from source to target Rubisco
            displacement = target_cords - source_cords;

            links(end+1) = pair;
        end
    
        carb.links = links;
    end
end

%% helper functions
function projection = calc_projection(rubisco1, rubisco2)
    % calculate the projection between two rubiscos
    dist_vec = [rubisco1.x rubisco1.y rubisco1.z] - [rubisco2.x rubisco2.y rubisco2.z];
    vec1 = rubisco1.vector;
    projection = dot(dist_vec, vec1)/norm(vec1);
end