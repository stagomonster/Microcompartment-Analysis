function [] = chain_spacing_analysis(carboxysome_data, min_chain_length, max_distance, bin_width)
% This function creates two plots regarding the spatial relationships 
% between chains. First, it makes a histogram of the angles between 
% neighboring chains. The user can select the width of the histogram bins, 
% the maximum distance between neighbors, and the minimum chain length.
% Second, it makes a normal Kernel probability density function 
% (bandwidth=0.5) of the distances between chains. The user can select the 
% minimum chain length. Both plots'data are colored based on the inner 
% concentration of the chains' parent carboxysome.
%
% Inputs
% carboxysome_data - an array of carboxysome objects filled with data at
%                    least through chain_maker.m
% min_chain_length - the minimum length a chain can be and still be
%                    included in this analysis
% max_distance - the maximum distance (in nanometers) between two chains 
%                that can still be neighbors
% bin_width - the width to make bins in the histogram (in degrees)
%
% chain_spacing_analysis.m © 2025 is licensed under CC BY-NC-SA 4.0

    %% important constants
    % load constants from external file
    CONSTANTS = constants();
    rubisco_diameter = CONSTANTS.RUBISCO_DIAMETER_M / CONSTANTS.PIXEL_SIZE;

    %% Calculate Data from Carboxysome Objects
    % create arrays to hold the data that will be plotted
    angles = [];
    angle_inner_concs = [];
    distances = [];
    carb_ids = [];
    dist_inner_concs = [];
    
    % for each carboxysome in the dataset
    for carb = carboxysome_data
        % for each chain with length >= min_chain_length
        valid_chains = carb.chains([carb.chains.length] >= min_chain_length);
        for i = 1:length(valid_chains)
            % for each other unchecked valid chain
            for j = i+1:length(valid_chains)
                % calculate distance and angle between chains
                [distance, angle] = chain_distance_closest(valid_chains(i), valid_chains(j), carb, rubisco_diameter);

                % Get data for distances plot
                distances(end+1) = distance * 1e9 * CONSTANTS.PIXEL_SIZE; % distance in nanometers
                carb_ids(end+1) = carb.carb_index;
                dist_inner_concs(end+1) = carb.inner_concentration;
                
                if distance < (max_distance*(10^-9)/CONSTANTS.PIXEL_SIZE) % comparison done in pixels
                    % Get data for angles plot
                    angles(end+1) = angle;
                    angle_inner_concs(end+1) = carb.inner_concentration;
                end
            end
        end
    end

    %% Angles Between Chains Plot
    custom_bins = 0:bin_width:90; % the edges of the custom bins
    [~, ~, bin_relation] = histcounts(angles, custom_bins); % finds to which bin each bend data point went
    plotz = NaN(length(custom_bins) - 1, length(angles)); % COLORBAR DATA IN SPECIAL FORMAT
    
    % Group the values in z by the vertical bar their corresponding y
    % values belong to. Each row in plotz represents a vertical bar. Each
    % column in plotz represents a layer in the bar.
    for i = 1:length(angles)
        plotz(bin_relation(i), i) = angle_inner_concs(i);
    end
    
    % Reorder the data in plotz so each column holds only one value that is
    % not NaN and each row is sorted from largest to smallest
    bookmark = 1; % keep track of the end of the data from the last row
    for i = 1:size(plotz, 1) % for each row in plotz
        data_to_move = sort(plotz(i, ~isnan(plotz(i, :))), 'descend')'; % extract all the data from the row
        plotz(i, :) = nan; % clear the row
        plotz(i, bookmark:bookmark + length(data_to_move) - 1) = data_to_move; % paste in the extracted data
        bookmark = bookmark + length(data_to_move); % update the bookmark
    end

    z_values = unique(sort(angle_inner_concs, 'descend'), 'stable'); % all the unique z values from large to small
    plotdata = zeros(length(custom_bins) - 1, length(z_values)); % the data that will be bar heights

    % calculates how many identical values exist in each row of plotz and
    % condenses them into one column. If there are n copies of a value, 
    % then num_repeats will have value n for that row and store n in some
    % position in plotdata
    for z_value = z_values
        num_repeats = sum(plotz == z_value, 2); % the number of times z_value appears in each row of plotz
        plotdata(:, z_values == z_value) = num_repeats; % each column of plotdata contains instances of the same z_value
    end

    figure;
    b = bar(1:length(custom_bins)-1, plotdata, 'stacked', 'FaceColor', 'flat'); % make a stacked bar graph

    % link the z data to a colormap
    cmap = colormap('winter');
    zmap = linspace(min(z_values), max(z_values), length(cmap));
    
    % Color each data point based on where it is between the min and max
    for i = 1:length(b)
        % make the bar's color proportional to its z value's distance between z_min and z_max
        b(i).CData = interp1(zmap, cmap, z_values(i));
        b(i).EdgeColor = 'none'; % remove the edges of the bars
    end

    % Make some lables for the plot
    title(['Angles Between Neighboring Chains of length >=', num2str(min_chain_length), '']);
    xlabel('Angle (deg)');
    ylabel('Number of Chain Pairs');

    % Make x axis labels that reflect the actual bin edges
    x_axis_labels = {};

    % Make labels for the x axis of the format [-45,-40), for example.
    % Matlab puts values on the edge of two bins in the larger bin, so we
    % use an open parenthesis on the right
    for i = 1:length(custom_bins) - 1
        if i == length(custom_bins) - 1
            % the largest bin gets a closed parenthesis on the right
            this_label = ['[',num2str(custom_bins(i)),',',num2str(custom_bins(i+1)),']'];
        else
            this_label = ['[',num2str(custom_bins(i)),',',num2str(custom_bins(i+1)),')'];
        end

        x_axis_labels{end+1} = this_label; % add label to list of labels
    end
    xticks(1:length(custom_bins) - 1); % make enough ticks for each bin
    xticklabels(x_axis_labels); % load the x axis labels
    
    % create and edit the colorbar
    caxis([min(z_values), max(z_values)]);
    c = colorbar('Location', 'southoutside');
    c.Label.String = 'Inner Rubisco Concentration (\muM)'; % colorbar title

    %% Distances Between Chains Plot
    % get bandwidth from user
    bandwidth = input('Please enter the bandwidth for the plot (recommended value is 0.5): ');

    % create a kernel probability density function for each carboxysome
    [pdf_array, ~, group_names] = fitdist(distances', 'Kernel', 'By', carb_ids, 'Kernel', 'Normal', 'Width', bandwidth);

    % start a new plot
    figure;
    hold on;

    % Link the inner concentration data to a colormap
    cmap = colormap('winter');
    zmap = linspace(min(dist_inner_concs), max(dist_inner_concs), length(cmap));

    % for each carboxysome, make a plot
    for i = 1:length(group_names)
        % filter out distances and concentrations in other carboxysomes
        x_data = distances(carb_ids == str2double(group_names{i}));
        z_data = dist_inner_concs(carb_ids == str2double(group_names{i}));

        % create 200 evenly spaced data points to sample the density plot
        density_x = linspace(min(x_data), max(x_data), 200);

        this_pdf = pdf(pdf_array{i}, density_x); % calculate the y values from the pdf
        plot_color = interp1(zmap, cmap, unique(z_data)); % set the color of the line
        plot(density_x, this_pdf, 'LineWidth', 2, 'color', plot_color); % plot the density
    end

    % create an all-inclusive kernel probability density function and plot
    overall_pdf = fitdist(distances', 'Kernel', 'Kernel', 'Normal', 'Width', bandwidth);

    % create 100 evenly spaced data points to sample the density plot
    density_x = linspace(min(distances), max(distances));

    overall_y = pdf(overall_pdf, density_x); % calculate the y values from the pdf
    plot(density_x, overall_y, 'LineWidth', 2, 'color', 'r'); % plot the overall density in red

    % Make some labels for the plot
    title('Probability Density of Lateral Distance Between Chains');
    xlabel('Distance Between Chains (nm)');
    ylabel('Probability Density');

    % create the colorbar
    caxis([min(dist_inner_concs), max(dist_inner_concs)]); % set colorbar tick labels
    c = colorbar('Location', 'southoutside');
    c.Label.String = 'Inner Rubisco Concentration (\muM)'; % colorbar title
end

 %% Helper Functions

function [shortest_distance, angle] = chain_distance_closest(chain_i, chain_j, carb, rubisco_diameter)
% Determine the distance between two chains. Each chain is broken up into 
% straight segments which connect adjacent rubisco. The shortest distance 
% between two segments is calculated for every combination of one segment 
% from chain i and one segment from chain j. The minimum of all these 
% segment-to-segment distances is returned as the distance between the 
% chains.
% 
% Inputs
% chain_i and chain_j are Rubisco_Chain objects whose distance is
% calculated. carb is a Carboxysome object which contains chains i and j.
% rubisco_diameter is a number which represents the diameter of a rubisco
% in pixels.
%
% Outputs
% This function returns distance, the distance between the two chains, and 
% angle, the angle between the two chains.

    % Fetch the lists of rubisco that are in these chains
    rubiscos_i = chain_i.tags;
    rubiscos_j = chain_j.tags;

    % Fetch the xyz positions of the rubiscos in each chain
    rubiscos_i_positions = get_rubisco_positions(rubiscos_i, carb);
    rubiscos_j_positions = get_rubisco_positions(rubiscos_j, carb);
    
    % Create arrays to store segments. Each row represents one
    % segment. Each segment has a xyz point in space (columns 1-3), three
    % vector components describing its direction (columns 5-7), and a
    % parameter lambda (column 4) which is used to calculate the distance
    % between the two endpoints of the segment.
    chain_i_segments = generate_segments(rubiscos_i_positions, carb, rubiscos_i, rubisco_diameter);
    chain_j_segments = generate_segments(rubiscos_j_positions, carb, rubiscos_j, rubisco_diameter);

    % Initialize a variable to store the minimum distance between the
    % chains.
    shortest_distance = inf;

    % Loop through the segments in chain i and calculate each one's minimum
    % distance to each segment in chain j
    for k = 1:size(chain_i_segments, 1)    
        segment_i = chain_i_segments(k, :);

        % Calculate where the endpoints of segment i are.
        i_endpoints = [segment_i(1:3); segment_i(1:3) + segment_i(4) * segment_i(5:7)];

        % For each segment in chain j
        for l = 1:size(chain_j_segments, 1)
            segment_distance = inf; % the least distance between these two segments
            segment_j = chain_j_segments(l, :);

            % Calculate where the endpoints of segment j are
            j_endpoints = [segment_j(1:3); segment_j(1:3) + segment_j(4) * segment_j(5:7)];

            % Calculate dot products used to find the location of the
            % minimum distance between the chains
            db = dot(segment_j(5:7), segment_i(5:7));
            bca = dot(segment_i(5:7), segment_j(1:3)-segment_i(1:3));
            bb = dot(segment_i(5:7), segment_i(5:7));
            dca = dot(segment_j(5:7), segment_j(1:3)-segment_i(1:3));
            dd = dot(segment_j(5:7), segment_j(5:7));

            % These equations represent the solutions to a minimization
            % problem. They were derived by setting the partial derivatives
            % of the distance between two points on opposite lines equal to
            % zero. s and t are scalars used to find how far along the line
            % segments the points that minimize distance are.
            s = (db * bca - bb * dca) / (bb * dd - db * db); % for segment j
            t = (dca + dd * s) / db; % for segment i

            % if the points that minimize distance both lay on the segments
            if s >= 0 && t >= 0 && s <= segment_j(4) && t <= segment_i(4)
                % find the points that minimize distance
                point_on_j = segment_j(1:3) + s * segment_j(5:7);
                point_on_i = segment_i(1:3) + t * segment_i(5:7);

                % calculate the distance between them
                segment_distance = sqrt(sum((point_on_i - point_on_j) .^ 2));
            end

            % find the shortest distance between segment i's first endpoint
            % and segment j
            j_point_1 = closest_point_on_segment_to_point(i_endpoints(1, :), j_endpoints);
            endpoints_1_dist = sqrt(sum((i_endpoints(1, :) - j_point_1) .^ 2));

            % find the shortest distance between segment i's second 
            % endpoint and segment j
            j_point_2 = closest_point_on_segment_to_point(i_endpoints(2, :), j_endpoints);
            endpoints_2_dist = sqrt(sum((i_endpoints(2, :) - j_point_2) .^ 2));

            % find the shortest distance between segment j's first endpoint
            % and segment i
            i_point_1 = closest_point_on_segment_to_point(j_endpoints(1, :), i_endpoints);
            endpoints_3_dist = sqrt(sum((j_endpoints(1, :) - i_point_1) .^ 2));

            % find the shortest distance between segment j's second 
            % endpoint and segment i
            i_point_2 = closest_point_on_segment_to_point(j_endpoints(2, :), i_endpoints);
            endpoints_4_dist = sqrt(sum((j_endpoints(2, :) - i_point_2) .^ 2));

            % of all 4 or 5 distances calculated, keep the minimum
            segment_distance = min([segment_distance, endpoints_1_dist, endpoints_2_dist, endpoints_3_dist, endpoints_4_dist]);
            
            % if the distance between these segments is the new shortest,
            % keep it
            if segment_distance < shortest_distance
                shortest_distance = segment_distance;
            end
        end
    end

    % Calculate the angle between the two chains
    angle = acos(dot(chain_i.average_vector, chain_j.average_vector)/(norm(chain_i.average_vector) * norm(chain_j.average_vector)))*180/pi;

    if angle > 90
        angle = 180 - angle;
    end
end

function rubisco_positions = get_rubisco_positions(rubisco_tags, carbox)
% Fetches the xyz positions of each rubisco in a chain and stores them.
% Each row represents one rubisco.
    rubisco_positions = zeros(length(rubisco_tags), 3); % array to hold positions

    for rubisco_tag = rubisco_tags % for each rubisco in the chain
        x_position = carbox.rubisco(find([carbox.rubisco.tag]==rubisco_tag)).x; % find its x position
        y_position = carbox.rubisco(find([carbox.rubisco.tag]==rubisco_tag)).y; % find its y position
        z_position = carbox.rubisco(find([carbox.rubisco.tag]==rubisco_tag)).z; % find its z position
        rubisco_positions(rubisco_tags == rubisco_tag, :) = [x_position y_position z_position]; % store them in array
    end
end

function chain_segments = generate_segments(rubisco_positions, carbox, rubisco_tags, rubisco_diameter)
% Creates the segments that make up a chain. For a chain with N rubisco
% there are N+1 segments. We make N-1 interior segments between rubisco and
% two exterior segments which extend from an end-of-chain rubisco center to
% its edge.

% Calculate interior segments between adjacent rubisco (not end-of-chain 
% extensions yet).
    chain_segments = zeros(size(rubisco_positions, 1) + 1, 7); % an array to hold chain segments
    
    for k = 2:length(rubisco_tags) % for each interior segment
        line_segment = rubisco_positions(k, :) - rubisco_positions(k - 1, :);
        chain_segments(k, 1:3) = rubisco_positions(k - 1, :); % segment start point
        chain_segments(k, 4) = 1; % lambda of segment
        chain_segments(k, 5:7) = line_segment; % direction of segment
    end

    % Determine the direction vectors of the end-of-chain extensions based
    % on the orientation vectors of the end-of-chain rubisco. The rubisco
    % orientations are in one of two possible directions. To get the one 
    % that extends the chain we check the dot product of the orientation 
    % vector with the next adjacent segment.

    % For the first rubisco in the chain, the vectors should point opposite
    % directions. If their dot product is positive we flip the rubisco
    % orientation vector and use that as the direction of the end-of-chain
    % extension. Else we take the orientation vector as is.
    if dot(carbox.rubisco(find([carbox.rubisco.tag]==rubisco_tags(1))).vector, chain_segments(2, 5:7)) > 0
        chain_segments(1, 5:7) = -1 * carbox.rubisco(find([carbox.rubisco.tag]==rubisco_tags(1))).vector;
    else
        chain_segments(1, 5:7) = carbox.rubisco(find([carbox.rubisco.tag]==rubisco_tags(1))).vector;
    end
    
    % For the last rubisco in the chain, the vectors should point in the 
    % same direction. If their dot product is negative we flip the rubisco
    % orientation vector and use that as the direction of the end-of-chain
    % extension. Else we take the orientation vector as is.
    if dot(carbox.rubisco(find([carbox.rubisco.tag]==rubisco_tags(end))).vector, chain_segments(end - 1, 5:7)) < 0
        chain_segments(end, 5:7) = -1 * carbox.rubisco(find([carbox.rubisco.tag]==rubisco_tags(end))).vector;
    else
        chain_segments(end, 5:7) = carbox.rubisco(find([carbox.rubisco.tag]==rubisco_tags(end))).vector;
    end
    
    % Set the end-of-chain segment xyz as the xyz of the end-of-chain 
    % rubisco and set their extension lengths to half of the diameter of a 
    % rubisco.
    chain_segments(1, 1:4) = [rubisco_positions(1, :) rubisco_diameter/2];
    chain_segments(end, 1:4) = [rubisco_positions(end, :) rubisco_diameter/2];
end

function closest_point = closest_point_on_segment_to_point(free_point, endpoints)
% Calculates the spot on a line segment which minimizes the distance to a
% free point. This is done by calculating t, the ratio of the scalar
% projection (dot product) of the point on the segment divided by the
% length of the segment.

    segment = endpoints(2, :) - endpoints(1, :); % a vector representing the segment
    segment_to_point = free_point - endpoints(1, :); % a vector from one segment endpoint to the free point

    t = dot(segment, segment_to_point) / dot(segment, segment); % how far along the segment the closest point is

    % for t<0 and t>1 the closest point is out of the bounds of the
    % segment. In this case we take the nearest endpoint (represented when
    % t=0 or t=1) as the closest point.
    if t < 0
        t = 0;
    elseif t > 1
        t = 1;
    end

    % Calculates the xyz coordinates of the closest point
    closest_point = endpoints(1, :) + t * segment;
end