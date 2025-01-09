function [] = rubisco_alignment_analysis(carboxysome_data, max_angle, max_distance, min_inner_conc)
% This function generates a graph showing the relationship between the
% number of aligned rubiscos in a carboxysome and the parameters of the
% carboxysome. Specifically, carboxysome volume, total concentration, and
% inner concentration are relevant values of the data points. A pair of
% rubiscos are "aligned" if they pass all filters specified by the user.
%
% Inputs
% carboxysome_data - an array of carboxysome objects filled with data
%                    through at least local_global_alignment.m
% max_angle - the maximum angle to allow between aligned rubiscos (degrees)
% max_distance - the maximum distance to allow between aligned rubiscos
%                (pixels, recommended to use the same value as in 
%                local_global_alignment.m)
% min_inner_conc - the minimum inner concentration of a carboxysome to
%                  allow on the plot (micro molar)
%
% rubisco_alignment_analysis.m © 2025 is licensed under CC BY-NC-SA 4.0

    % Initialize arrays to hold data to plot
    valid_rubiscos = [];
    volumes = [];
    inner_concentrations = [];
    total_concentrations = [];

    for carb = carboxysome_data % loop through every carboxysome
        num_valid_pairs = 0; % number of rubisco pairs that satisfy every condition

        are_inner = [carb.rubisco_pairs.inner]; % both rubiscos are inner rubiscos
        within_angle = [carb.rubisco_pairs.angle] <= max_angle; % angle between rubiscos is less than user-defined max
        within_dist = [carb.rubisco_pairs.distance] <= max_distance; % distance between rubiscos is less than user-defined max

        if carb.inner_concentration >= min_inner_conc % if carboxysome inner concentration meets user-defined threshold
            % calculate how many rubiscos satisfy every condition
            num_valid_pairs = length(carb.rubisco_pairs(are_inner & within_angle & within_dist));
        end

        if num_valid_pairs > 0 % if carboxysome has any valid rubisco pairs
            valid_rubiscos(end+1) = num_valid_pairs; % store number of valid pairs
            volumes(end+1) = carb.volume * (10^18); % store carboxysome volume
            inner_concentrations(end+1) = carb.inner_concentration; % store carboxysome inner concentration
            total_concentrations(end+1) = carb.concentration; % store carboxysome total concentration
        end
    end

    % Number of aligned rubiscos vs. Inner Concentration Plot
    bubblechart(inner_concentrations, valid_rubiscos, volumes, total_concentrations); % make bubble chart
    set(gca, 'YScale', 'log'); % log scale y axis

    title('Rubisco Alignment as a function of Carboxysome Parameters'); % graph title
    ylabel(['Count of Neighboring RuBisCO pairs aligned within ', num2str(max_angle), ' degrees']); % y axis title
    xlabel('Inner RuBisCO Concentration (\muM)'); % x axis title

    colormap('winter');
    c = colorbar('Location', 'westoutside'); % color bubbles by total concentration
    c.Label.String = 'Total RuBisCO Concentration (\muM)'; % colorbar title
    bubblelegend('Carboxysome Volume (\mum^3)', 'Location', 'northeastoutside'); % size points by carboxysome volume
end