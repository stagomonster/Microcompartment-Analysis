function [] = radius_factor_calculator()
% This function calculates what size the radius_factor input to linkages.m
% should be to maximize the search cylinder size while disallowing multiple
% multiple monomers in the search cylinder. Two non-overlapping monomers
% will be too large for both to be detected inside the resulting search 
% cylinder. The user inputs three other parameters of the search cylinder, 
% then this function chooses the radius_factor. User should use the output 
% plot of link_arc_analysis.m to choose the other 3 values.
%
% linkage_limit_calculator.m © 2025 is licensed under CC BY-NC-SA 4.0
    
    %% Retrieve search cylinder parameters from user
    angle_limit = input('Enter the maximum angle (degrees): ');
    angle_limit = angle_limit * (pi / 180); % convert to radians
    cylinder_far_end = input('Enter the location of the far end of the cylinder (particle diameters): ');
    cylinder_near_end = input('Enter the location of the near end of the cylinder (particle diameters): ');

    %% Prepare for loop
    % Force one target monomer to the top of the search cylinder. We will
    % force its x value to the edge of the cylinder later, ensuring it is
    % in the top corner.
    target_center_y = cylinder_far_end;

    % if the target monomer is far enough away from the source monomer that
    % overlap is not possible
    if target_center_y > 0.5 + sqrt(10 / 4)
        can_overlap = false; % if too far away, can't overlap
    else
        can_overlap = true; % if close enough, could overlap
    end

    % Calculate a starting search cylinder width using the user-input angle
    % limit. At this width, the only way to fit multiple monomers in the
    % search cylinder is for one to exceed the angle limit.
    cylinder_width = cos(angle_limit);
    target_center_x = cylinder_width / 2; % force target monomer to the edge of the cylinder
    
    % if overlap is possible
    if can_overlap
        if angle_limit > atan(3) % if angle limit is greater than ~71 degrees
            % Coming up, we are going to check if there is overlap between
            % monomers. 71 degrees is the most extreme angle possible for
            % the configuration of two monomers side by side in the search
            % cylinder. If there is no overlap at 71 degrees, there will be
            % no overlap at any angle.

            % Calculate the locations of points on the edges of the target 
            % monomers. These will be used to calculate if a target monomer
            % overlaps with the source monomer.
            target_1_edge_x = target_center_x + 0.5 * sin(atan(3));
            target_2_edge_x = target_1_edge_x - 1.5 * cos(atan(3));
            target_1_edge_y = target_center_y - 0.5 * cos(atan(3));
            target_2_edge_y = target_1_edge_y - 1.5 * sin(atan(3));
        else % if angle limit is not greater than ~71 degrees
            % Calculate the locations of points on the edges of the target
            % monomers when angle is maximal. These will be used to
            % calculate if a target monomer overlaps with the source
            % monomer.
            target_1_edge_x = target_center_x + 0.5 * sin(angle_limit);
            target_2_edge_x = target_1_edge_x - 1.5 * cos(angle_limit);
            target_1_edge_y = target_center_y - 0.5 * cos(angle_limit);
            target_2_edge_y = target_1_edge_y - 1.5 * sin(angle_limit);
        end
    
        % Calculate the x coordinate of the point of intersection between
        % the top of the source monomer and a line passing between
        % target_1_edge and target_2_edge
        m = (target_1_edge_y - target_2_edge_y) / (target_1_edge_x - target_2_edge_x);
        intersect_x = (0.5 - target_2_edge_y) / m + target_2_edge_x;
        
        % if the point of intersection is right of the source monomer's top
        % left corner and target_2_edge
        if intersect_x > max(target_2_edge_x, -0.5)
            % if true then there must be overlap
            is_overlapping = true;
        else
            % if false there is no overlap and there cannot be overlap at
            % any angle
            is_overlapping = false;
        end
        
        % if there is overlap
        if is_overlapping

            % while there is still overlap
            while is_overlapping
                % calculate x position of target_2_edge when
                % target_2_edge_y = 0.5
                target_2_edge_x = target_center_x - sqrt(10/4 - (target_center_y - 0.5) ^ 2);

                % if target_2_edge_x is left of the corner of the source
                % monomer
                if target_2_edge_x < -0.5
                    % express the location of an edge of the source monomer 
                    % relative to the target monomer that was forced into a 
                    % corner
                    p = -0.5 - target_center_x; % x coordinate
                    q = 0.5 - target_center_y; % y coordinate
            
                    % Calculate the maximum possible angle that does not cause 
                    % a target monomer to overlap with the source monomer
                    m = (-1*p*q + 0.5*sqrt(p^2 + q^2 - 1/4)) / (1/4 - p^2);
                    no_overlap_angle = atan(m);
                else
                    % otherwise, this is the maximum angle that does not
                    % cause a target monomer to overlap with the source
                    % monomer
                    no_overlap_angle = asin((target_center_y - 0.5) / sqrt(10/4)) - atan(1/3);
                end

                % if the calculated angle exceeds the user-specified angle 
                % limit, use the angle limit instead
                if no_overlap_angle > angle_limit
                    no_overlap_angle = angle_limit;
                end
    
                % calculate the location of the center of the second target
                % monomer based on the location of the target monomer in
                % the corner and the angle which prevents overlaps
                second_target_x = target_center_x - cos(no_overlap_angle);
    
                % if the cylinder is skinnier than the distance between the
                % target centers
                if target_center_x - second_target_x > cylinder_width
                    % If true, the cylinder can be made wider. Its new 
                    % width is the distance between the x positions of the 
                    % two target monomers.
                    cylinder_width = target_center_x - second_target_x;
                    target_center_x = cylinder_width / 2; % update the location of the monomer in the corner
                else
                    % If false, there is no more overlap and the loop can
                    % be exited
                    is_overlapping = false;
                end

                % Updating the location of the center of the monomer in the
                % corner reintoduces overlap, so we must loop this process
                % until the cylinder width stops growing
            end
        end
    end
    
    % calculate the cylinder height
    cylinder_height = cylinder_far_end - cylinder_near_end;

    % if the diagonal of the search cylinder is too small to fit
    % two monomers in
    if cylinder_width ^ 2 + cylinder_height ^ 2 < 1
        % expand its width until it is just barely too small to fit
        % two monomers
        cylinder_width = sqrt(1 - cylinder_height ^ 2);
    end

    radius_factor = cylinder_width / 2; % calculate the radius factor

    % print a statement to the user about what value to use for the radius
    % factor
    fprintf(['Use a radius factor value less than %.4f to safely eliminate' ...
        ' multiple real particles at the same binding site.\n'], radius_factor);
end