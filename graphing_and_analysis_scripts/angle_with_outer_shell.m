function [] = angle_with_outer_shell(carboxysome_data, bin_width)
% This function plots a histogram of the angles between rubiscos and the
% carboxysome shell. It fetches from the input data the orientation and
% average normal vector of each rubisco and calculates the angle between
% them. These angles are plotted for all outer rubiscos in the dataset. In
% addition, random angles are calculated and overlaid on the plot for
% comparison.
%
% Inputs
% carboxysome_data - an array of carboxysome objects containing data
%                    through at least local_global_alignment.m
% bin_width - the width to make the bins in the histogram (in degrees)
%
% angle_with_outer_shell.m © 2025 is licensed under CC BY-NC-SA 4.0
    
    %% get rubisco average normals and orientations from data
    orientations = []; % holds the orientations of the rubiscos
    ave_normals = []; % holds the average normal vectors of facets surrounding rubiscos
    concentrations = []; % holds the parent carboxysome concentrations

    for carb = carboxysome_data % for every carboxysome
        for rubisco = carb.rubisco % for every rubisco
            if ~rubisco.inside % if it is an outer rubisco
                orientations = [orientations; rubisco.vector]; % get rubisco orientation
                ave_normals = [ave_normals; rubisco.ave_normal]; % get average normal vector of surrounding facets
                concentrations(end+1) = carb.concentration; % get parent carboxysome concentration
            end
        end
    end

    %% get random angles on sphere
    % get as many random vectors as you have real data points
    random_vecs = compute_random_vectors_marsaglia_72(length(concentrations));

    %% calculate real and random angles
    vec_wall = [1 0 0]; % a vector to compare the random vectors with
    data_angles = []; % the angles made by the real data
    random_angles = []; % the angles made by the random vectors

    % calculate the angles for the real data and random data
    for i = 1:length(concentrations)
        data_angles(end+1) = calc_angle(orientations(i,:), ave_normals(i,:)); % real data
        random_angles(end+1) = calc_angle(random_vecs(i,:), vec_wall); % random angles
    end

    %% plot the data
    custom_bins = 0:bin_width:90; % the edges of the custom bins
    [~, ~, data_bin_relation] = histcounts(data_angles, custom_bins); % finds to which bin each angle data point went
    [random_counts, ~, ~] = histcounts(random_angles, custom_bins); % finds to which bin each random angle went
    plotz = NaN(length(custom_bins) - 1, length(data_angles)); % this array will be used to make the colorbar
    
    % Group the values in concentrations by the vertical bar their
    % corresponding angle values belong to. Each row in plotz 
    % represents a vertical bar. Each column in plotz represents a layer in
    % the bar.
    for i = 1:length(data_angles)
        plotz(data_bin_relation(i), i) = concentrations(i);
    end

    % Reorder the data in plotz so each column holds only one value that is
    % not NaN and each row is sorted from largest to smallest
    bookmark = 1; % keep track of the end of the data from the last row
    for i = 1:size(plotz, 1) % for each row in plotz
        data_to_move = sort(plotz(i, ~isnan(plotz(i, :))), 'descend')'; % extract all the data from the row
        plotz(i, :) = nan; % clear the row
        plotz(i, bookmark:bookmark + length(data_to_move) - 1) = data_to_move; % paste in the extracted data starting at the bookmark
        bookmark = bookmark + length(data_to_move); % update the bookmark
    end

    z_values = unique(sort(concentrations, 'descend'), 'stable'); % all the unique z values from large to small
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
    b = bar(1:length(custom_bins) - 1, plotdata, 'stacked', 'FaceColor', 'flat'); % make a stacked bar graph
    hold on;
    scatter(1:length(custom_bins) - 1, random_counts, 150, 'filled', 'MarkerFaceColor', [1 0.4 0.3], 'MarkerFaceAlpha', 0.7); % plot random data

    % Link the concentration data to a colormap if there are multiple
    % carboxysomes
    if length(carboxysome_data) > 1
        cmap = colormap('winter');
        zmap = linspace(min(concentrations), max(concentrations), length(cmap));

        % Color each data point based on where it is between the min and max
        for i = 1:length(b)
            b(i).CData = interp1(zmap, cmap, z_values(i)); % set the proportional color to the bar's color
            b(i).EdgeColor = 'none'; % remove the edges of the bars
        end
    end

    % Make some labels for the plot
    x_axis_label = {'Angle between the Rubisco orientation and the average normal';...
        'vector of neighboring facets from convex hull (deg)'};
    title('Angles Between Outer Rubisco and Shell');
    xlabel(x_axis_label);
    ylabel('Number of Rubiscos');

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

    if length(carboxysome_data) > 1
        % Edit the colorbar
        c = colorbar('Location', 'southoutside'); % colorbar location
        c.Label.String = 'Rubisco Concentration(\muM)'; % colorbar title
        c.Ticks = [0, 0.25, 0.5, 0.75, 1];
        c.TickLabels = round(linspace(min(concentrations), max(concentrations), 5), 0); % colorbar labels
    end
end

%% helper functions
function angle = calc_angle(vec1, vec2)
% calculate the angle in degrees between two vectors vec1 and vec2 in R^3
    angle = acos(abs((dot(vec1, vec2))/(norm(vec1)*norm(vec2))))*(180/pi);
end

function mat = compute_random_vectors_marsaglia_72(N)
% compute N random vectors in R^3 for comparison with rubisco data. mat is
% an N by 3 matrix where each row is a complete random vector.
    mat = zeros(N, 3);
    for i = 1:N
        init = 2;
        v1 = 0;
        v2 = 0;
        while init >= 1
            v1 = 2*(rand() - 0.5);
            v2 = 2*(rand() - 0.5);
            init = v1*v1 + v2*v2;
        end
        vec = [2*v1*sqrt(1 - init) 2*v2*sqrt(1 - init) 1 - 2*init];
        normed = vec/norm(vec);
        mat(i, :) = normed; % return the normalized random vectors
    end
end
