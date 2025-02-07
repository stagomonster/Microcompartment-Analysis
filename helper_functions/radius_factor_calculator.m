function [] = radius_factor_calculator()
% This function calculates what size the radius_factor input to linkages.m
% should be to maximize the search cylinder size while disallowing multiple
% multiple monomers in the search cylinder. Two non-overlapping monomers
% will be too large for both to be detected inside the resulting search 
% cylinder. The user inputs two other parameters of the search cylinder, 
% then this function chooses the radius_factor. User should use the output 
% plot of link_arc_analysis.m to choose the other 2 values.
%
% linkage_limit_calculator.m © 2025 is licensed under CC BY-NC-SA 4.0
    
    %% Retrieve search cylinder parameters from user
    angle_limit = input('Enter the maximum angle (degrees): ');
    angle_limit = angle_limit * (pi / 180); % convert to radians
    cylinder_far_end = input('Enter the location of the far end of the cylinder (particle diameters): ');

    %% Prepare for loop
    % Force one target monomer to the top of the search cylinder. We will
    % force its x value to the edge of the cylinder later, ensuring it is
    % in the top corner.
    target_center_y = cylinder_far_end;

    % Calculate a starting search cylinder width using the user-input angle
    % limit. At this width, the only way to fit multiple monomers in the
    % search cylinder is for one to exceed the angle limit.
    cylinder_width = cos(angle_limit);

    % initialize variables to control looping
    too_small = true;
    i=1;

    % while the cylinder can still be made wider
    while too_small

        % force one target monomer into the corner of the search cylinder
        target_center_x = cylinder_width / 2;

        % express the location of an edge of the source monomer relative to
        % the target monomer that was forced into a corner
        p = -0.5 - target_center_x;
        q = 0.5 - target_center_y;
    
        % Calculate the maximum possible angle that does not cause a target
        % monomer to overlap with the source monomer
        m = (-1*p*q + 0.5*sqrt(p^2 + q^2 - 1/4)) / (1/4 - p^2);
        no_overlap_angle = atan(m);
    
        % if the calculated angle exceeds the user-specified angle limit,
        % use the angle limit instead
        if no_overlap_angle > angle_limit
            no_overlap_angle = angle_limit;
        end
    
        % Based on the calculated angle, find the x position of the second
        % target monomer. It is assumed that it is exactly 1 particle
        % diameter away from the first target monomer (the one forced into
        % the corner).
        second_target_x = target_center_x - cos(no_overlap_angle);

        %% Determine whether the source monomer and second target monomer 
        %% overlap when oriented at the angle limit
        % Calculate the location of a point on the line that connects the
        % edge of the source monomer with the edge of the first target
        % monomer
        r = target_center_x + 0.5 * sin(angle_limit); % x position
        s = target_center_y - 0.5 * cos(angle_limit); % y position

        % Calculate where that line intersects the x=-0.5 line and the
        % y=0.5 line
        x_intersect = r + (0.5 - s) / tan(angle_limit); % x value when y=0.5
        y_intersect = s + tan(angle_limit) * (-1 * 0.5 - r); % y value when x=-0.5

        % if the line cuts through the source monomer
        if (x_intersect > -0.5 && x_intersect < 0.5) || (y_intersect > -0.5 && y_intersect < 0.5)
            is_overlapping = true; % then two monomers are overlapping at the angle limit
        else
            is_overlapping = false; % otherwise there is no overlap at the angle limit
        end
        
        %% Determine if the cylinder can be made wider
        % if the x positions of the two target monomers are further apart
        % than the cylinder is wide and there is monomer overlap
        if target_center_x - second_target_x > cylinder_width && is_overlapping
            % The cylinder can be made wider. Its new width is the distance
            % between the x positions of the two target monomers.
            cylinder_width = target_center_x - second_target_x;
            d = cylinder_width / 2;
        else
            % The cylinder cannot be made wider. Set the loop condition to
            % false and do not increase the width.
            too_small = false;
            d = cylinder_width / 2;
        end

        %% Prevent the while loop from running forever
        % In theory, this loop would run forever under the right
        % conditions. In practice, improvements get smaller with each
        % iteration until they are too small for Matlab to track. Just to
        % be safe, the loop forcefully exits after 20 iterations.
        i = i + 1; % increment an iteration counter
        if i > 20 % after 20 iterations
            too_small = false; % set the loop condition to false
        end
    end

    % print a statement to the user about what value to use for the radius
    % factor
    fprintf(['With the given input parameters, use a radius factor value less ' ...
        'than %.4f to safely eliminate multiple real particles at the same binding site.\n'], d);
end