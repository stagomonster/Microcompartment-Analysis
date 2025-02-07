function link_arc_analysis(carboxysome_data, label, pointSize, max_angle)
% This function generates a scatter plot of the positions of target rubisco
% relative to source rubisco. On the y-axis is the distance between the 
% rubisco centers along the direction of the symmetry axis of the source 
% rubisco. On the x-axis is the distance between the rubisco centers 
% perpendicular to the symmetry axis of the source rubisco. Both distances 
% are represented as a multiple of the rubisco diameter. The color of each
% point on the plot represents the angle between the axes of symmetry of
% the two rubisco. Then the function plots the x data and y data as 
% histograms on separate figures. The data can be filtered before plotting
% the histograms for clarity. The user can specify a label to be included 
% in the scatter plot title, the size of the points on the scatter plot, 
% and the maximum angle between rubisco that are still plotted.
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
        for link = carb.rubisco_pairs % for each rubisco pair in the carboxysome
            % target's distance along the source's axis
            projection = link.projection;
            % target's distance away from the source's axis
            d = sqrt(link.distance^2 - projection^2);
            % angle between axes of source and target
            z = link.angle;
            
            % store data in arrays to be plotted
            if z <= max_angle % filter out data by angle between axes
                xData(end+1) = d / rubisco_diameter; % convert units to rubisco diameters
                if projection < 0
                    yData(end+1) = -1*projection / rubisco_diameter; % convert units to rubisco diameters
                else
                    yData(end+1) = projection / rubisco_diameter; % convert units to rubisco diameters
                end
                zData(end+1) = z;
            end
        end
    end

    % create the three plots
    visualize_scatter(xData, yData, zData, pointSize, label); % create the scatter plot
    visualize_histogram(xData, yData, 'from'); % create histogram of distances from axis
    visualize_histogram(yData, xData, 'along'); % create histogram of distances along axis
end

function visualize_scatter(xData, yData, zData, pointSize, label)
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

function visualize_histogram(x_data, filter_data, label)
% This function makes a histogram. On the y axis is the number of instances
% of an x axis value. On the x axis is the distance between rubisco either 
% along the source's axis or from the source's axis (radially). The data
% depends on what data was passed into x_data.

    % Ask the user if they want to filter the data
    if strcmp(label, 'from')
        % filter radial data to exclude points with "along axis" values
        % greater than a user-defined maximum
        filter = input(['Plotting Deviation from Central Axis:\nEnter the maximum ' ...
            'deviation along the central axis a data point can have in this plot (use "inf" to include all): ']);
    else
        % filter "along axis" data to exclude points with radial values
        % greater than a user-defined maximum
        filter = input(['Plotting Deviation along Central Axis:\nEnter the maximum ' ...
            'deviation from the central axis a data point can have in this plot (use "inf" to include all): ']);
    end

    x_data = x_data(filter_data <= filter); % filter the data

    figure;
    histogram(x_data); % create the histogram

    % Make some labels for the plot
    title(['Distances Between Monomers ', label, ' Central Axis']);
    xlabel(['Deviation ', label, ' Central Axis (Monomer Diameters)']);
    ylabel('Number of Monomers');
end