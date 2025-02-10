function link_arc_analysis(carboxysome_data, label, pointSize, max_angle)
% This function generates a scatter plot of the positions of target rubisco
% relative to source rubisco. On the x-axis is the distance between the
% rubisco centers perpendicular to the symmetry axis of the source rubisco.
% On the y-axis is the distance between the rubisco centers along the
% direction of the symmetry axis of the source rubisco. Both distances are
% represented as a multiple of the rubisco diameter. The color of each
% point on the plot represents the angle between the axes of symmetry of
% the two rubisco. The user can specify a label to be included in the
% figure title, the size of the points on the plot, and the maximum angle
% between rubisco that are still plotted.
%
% Inputs
% carboxysome_data - an array of Carboxysome objects already containing
%                    computed data and pairs from local_global_alignment.m.
% label - a string to be printed in the title of the figure
% pointSize - the size of the data points in the figure (default is 2)
% max_angle - the maximum angle between the symmetry axes of rubisco to
%             allow to be plotted
%
% link_arc_analysis.m © 2025 is licensed under CC BY-NC-SA 4.0
    
    % Important Constants
    CONSTANTS = constants();
    rubisco_diameter = CONSTANTS.RUBISCO_DIAMETER_M /CONSTANTS.PIXEL_SIZE; % rubisco diameter in pixels

    % default angle
    if nargin < 4
        max_angle = 90;
    end
    % default size of the points on the figure
    if nargin < 3
        pointSize = 2;
    end
    % defaults to no label in figure title
    if nargin < 2
        label = '';
    end

    % initialize arrays used to store data for the plot
    xData = [];
    yData = [];
    zData = [];

    for carb = carboxysome_data % for each carboxysome
        xData = [xData, zeros(length(carb.rubisco_pairs))];
        yData = [yData, zeros(length(carb.rubisco_pairs))];
        zData = [zData, zeros(length(carb.rubisco_pairs))];
        last_index = 0;

        for link = carb.rubisco_pairs % for each rubisco pair in the carboxysome
            % target's distance along the source's axis
            projection = link.projection;
            % target's distance away from the source's axis
            d = sqrt(link.distance^2 - projection^2);
            % angle between axes of source and target
            z = link.angle;
            
            % store data in arrays to be plotted
            if z <= max_angle % filter out data by angle between axes
                xData(last_index + 1) = d / rubisco_diameter; % convert units to rubisco diameters
                if projection < 0
                    yData(last_index + 1) = -1*projection / rubisco_diameter; % convert units to rubisco diameters
                else
                    yData(last_index + 1) = projection / rubisco_diameter; % convert units to rubisco diameters
                end
                zData(last_index + 1) = z;
                last_index = last_index + 1;
            end
        end
    end
    if last_index == 0
        xData = [];
        yData = [];
        zData = [];
    else
        xData = xData(1:last_index);
        yData = yData(1:last_index);
        zData = zData(1:last_index);
    end

    visualize(xData, yData, zData, pointSize, label); % create the figure
end

function visualize(xData, yData, zData, pointSize, label)
% This function makes the scatterplot, colored by angle between source and
% target axes. On the y axis is the distance between rubisco along the
% source's axis. On the x axis is the distance between rubisco away from
% the source's axis.

    % if the user specifies a label, displays it in the figure title
    if ~isempty(label)
        label = ['(' label ')'];
    end

    % displays the data with all the coloration that the angles specify
    figure;
    hold on;
    scatter(xData, yData, pointSize, zData, 'filled');

    title(['All Rubisco Pair Targets Relative to Source Rubisco ' label]); % figure title
    xlabel('Deviation from Central Axis (Monomer Diameters)'); % x axis label
    ylabel('Projection along Central Axis (Monomer Diameters)'); % y axis label
    xlim([0, 3.5]); % x axis limits
    ylim([0, 3.5]); % y axis limits
    c = colorbar('fontsize', 14); % create a colorbar for color data
    c.Label.String = 'Angle (Degrees)'; % colorbar title
end