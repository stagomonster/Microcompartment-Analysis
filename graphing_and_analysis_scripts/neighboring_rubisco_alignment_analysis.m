function [] = neighboring_rubisco_alignment_analysis(carboxysome_data, max_distance, plot_only_inner, min_rubiscos_inner)
% This function creates a plot of the probability density of the angle 
% between nearby rubiscos, grouped by parent carboxysome. The user can 
% select the maximum distance between rubiscos to be considered in this 
% plot, whether to plot only inner rubisco pairs, and the minimum number of
% rubiscos a carboxysome needs to have to be plotted. The density is a 
% normal Kernel probability density function with a bandwidth of 5. The 
% plot's data are colored based on the total concentration of the rubiscos'
% parent carboxysome. For reference, a plot of random angles is stacked
% on top of the real data in red.
%
% Inputs
% carboxysome_data - an array of carboxysome objects filled with data at
%                    least through local_global_alignment.m
% max_distance - the maximum distance two rubiscos can be from each other
%                and still be plotted in this analysis (pixels)
% plot_only_inner - a boolean which filters out rubisco pairs involving
%                   outer rubiscos if set to true
% min_rubiscos_inner - the minimum number of inner rubiscos a carboxysome 
%                      needs to have to be plotted in this analysis
%
% neighboring_rubisco_alignment_analysis.m © 2025 is licensed under CC BY-NC-SA 4.0

    %% Get Angles Between Rubiscos Data
    % initialize arrays to hold data to plot
    angles = [];
    carb_ids = [];
    concentrations = [];

    for carb = carboxysome_data % for each carboxysome
        % filter out rubisco pairs further away than max_distance
        valid_pairs = carb.rubisco_pairs([carb.rubisco_pairs.distance] < max_distance);

        % filter out rubisco pairs from carboxysomes containing fewer than
        % min_rubiscos_inner interior rubiscos
        valid_pairs = valid_pairs([valid_pairs.num_rubisco_inner] >= min_rubiscos_inner);

        if plot_only_inner % if the user wants to plot only inner rubiscos
            % filter out pairs involving outer rubiscos
            valid_pairs = valid_pairs([valid_pairs.inner]);
        end

        % add relevant data from remaining valid pairs to data to plot
        angles = [angles, [valid_pairs.angle]];
        carb_ids = [carb_ids, [valid_pairs.reg]];
        concentrations = [concentrations, [valid_pairs.concentration]];
    end
            
    %% Angles Between Rubiscos Plot
    % get bandwidth from user
    bandwidth = input('Please enter the bandwidth for the plot (recommended value is 5): ');

    % create a kernel probability density function for each carboxysome
    [pdf_array, ~, group_names] = fitdist(angles', 'Kernel', 'By', carb_ids, 'Kernel', 'Normal', 'Width', bandwidth);

    % start a new plot
    figure;
    hold on;

    % Link the concentration data to a colormap if there are multiple
    % carboxysomes
    if length(carboxysome_data) > 1
        cmap = colormap('winter');
        zmap = linspace(min(concentrations), max(concentrations), length(cmap));
    end

    % for each carboxysome, make a plot
    for i = 1:length(group_names)
        % filter out angles and concentrations in other carboxysomes
        x_data = angles(carb_ids == str2double(group_names{i}));
        z_data = concentrations(carb_ids == str2double(group_names{i}));

        % create 200 evenly spaced data points to sample the density plot
        density_x = linspace(min(x_data), max(x_data), 200);

        this_pdf = pdf(pdf_array{i}, density_x); % calculate the y values from the pdf

        % set the color of the line
        if length(carboxysome_data) > 1
            plot_color = interp1(zmap, cmap, unique(z_data));
        else
            plot_color = 'blue';
        end

        plot(density_x, this_pdf, 'LineWidth', 2, 'color', plot_color); % plot the density
    end

    %% Generate Random Angles Data
    % calculate a number of random vectors to generate
    rubisco_per_carb_average = sum([carboxysome_data.num_rubisco])/length(carboxysome_data);
    N = ceil(rubisco_per_carb_average);

    % generate random vectors
    rand_vectors = compute_random_vectors_marsaglia_72(N);
    random_angles = []; % initialize an array to hold random angles

    % for each possible combination of random vectors
    for k = 1:N
        for j = k+1:N
            % calculate the angle between them
            random_angles(end+1) = calc_angle(rand_vectors(k, :), rand_vectors(j, :));
        end
    end

    %% Plot Random Angles Density
    % create a kernel probability density function for the random angles
    random_pdf = fitdist(random_angles', 'Kernel', 'Kernel', 'Normal', 'Width', bandwidth);

    % create 200 evenly spaced data points to sample the density plot
    density_x = linspace(min(random_angles), max(random_angles), 200);

    random_y = pdf(random_pdf, density_x); % calculate the y values from the pdf
    plot(density_x, random_y, 'LineWidth', 5, 'color', 'r'); % plot the random density in red

    % Make some labels for the plot
    title('Probability Density of Angles Between Nearby RuBisCOs');
    xlabel('Angle Between a Reference RuBisCO and its Nearest Neighbors (deg)');
    ylabel('Probability Density');

    % create the colorbar if needed
    if length(carboxysome_data) > 1
        caxis([min(concentrations), max(concentrations)]); % set colorbar tick labels
        c = colorbar('Location', 'southoutside');
        c.Label.String = 'Rubisco Concentration (\muM)'; % colorbar title
    end
end

    %% Helper Functions
function mat = compute_random_vectors_marsaglia_72(N)
% compute N random vectors in R^3 for comparison with rubisco data. mat is
% an N by 3 matrix where each row is a complete random vector.
    mat = zeros(N, 3);
    for l = 1:N
        init = 2;
        v1 = 0;
        v2 = 0;
        while init >= 1
            v1 = 2*(rand() - 0.5);
            v2 = 2*(rand() - 0.5);
            init = v1*v1 + v2*v2;
        end
        vec = [2*v1*sqrt(1 - init), 2*v2*sqrt(1 - init), 1 - 2*init];
        normed = vec/norm(vec);
        mat(l, :) = normed;
    end
end

function angle = calc_angle(vec1, vec2)
% calculate the angle in degrees between two vectors vec1 and vec2 in R^3
    angle = acos(abs(dot(vec1, vec2)/(norm(vec1)*norm(vec2)))) * (180/pi);
end
