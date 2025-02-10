function [] = rubisco_spacing_analysis(carboxysome_data, min_chain_length)
% This function creates a plot of the probability density of the distance 
% between adjacent rubiscos in a chain, grouped by parent carboxysome. The 
% user can select the minimum chain length to be considered in this plot. 
% The density is a normal Kernel probability density function with a 
% bandwidth of 0.3.The plot's data are colored based on the inner 
% concentration of the rubiscos' parent carboxysome.
%
% Inputs
% carboxysome_data - an array of carboxysome objects filled with data at
%                    least through chain_maker.m
% min_chain_length - the minimum length a chain can be and still be
%                    included in this analysis
%
% rubisco_spacing_analysis.m © 2025 is licensed under CC BY-NC-SA 4.0

    %% important constants
    % load constants from external file
    CONSTANTS = constants();

    %% get relevant data from carboxysome array
    % initalize arrays to hold relevant data
    distances = [];
    inner_concentrations = [];
    carb_ids = [];
    last_index = 0;

    for carb = carboxysome_data % for each carboxysome
        for chain = carb.chains([carb.chains.length] >= min_chain_length) % for each chain over the minimum length
            iterations_per_chain = ((length(chain.indices) * (length(chain.indices) - 1 )) / 2);
            distances = [distances, zeros(iterations_per_chain)];
            inner_concentrations = [inner_concentrations, zeros(iterations_per_chain)];
            carb_ids = [carb_ids, zeros(iterations_per_chain)];
            for i = 1:length(chain.indices) - 1 % for each rubisco in the chain
                for j = i+1:length(chain.indices) % for each other rubisco in the chain
                    % get the indices of the two rubisco to find a linkage
                    % for
                    rubisco_i = chain.indices(i);
                    rubisco_j = chain.indices(j);
    
                    % find any existing linkage where rubisco_i is the 
                    % I_index and rubisco_j is the J_index
                    linkage = carb.links([carb.links.I_index] == rubisco_i & [carb.links.J_index] == rubisco_j);

                    if isempty(linkage) % if there is no such linkage
                        % find any existing linkage where rubisco_i is the
                        % J_index and rubisco_j is the I_index
                        linkage = carb.links([carb.links.J_index] == rubisco_i & [carb.links.I_index] == rubisco_j);
                    end
                    
                    if ~isempty(linkage) % if a linkage exists between the rubiscos
                        % fetch the distance, inner concentration, and
                        % parent carboxysome
                        distances(last_index+1) = linkage.distance * 1e9 * CONSTANTS.PIXEL_SIZE; % distance in nm
                        inner_concentrations(last_index+1) = carb.inner_concentration;
                        carb_ids(last_index+1) = carb.carb_index;
                    end
                end
            end
        end
    end
    if last_index == 0
        distances = [];
        inner_concentrations = [];
        carb_ids = [];
    else
        distances = distances(1:last_index);
        inner_concentrations = inner_concentrations(1:last_index);
        carb_ids = carb_ids(1:last_index);
    end

    %% Distances Between Adjacent Rubiscos Plot
    % get bandwidth from user
    bandwidth = input('Please enter the bandwidth for the plot (recommended value is 0.3): ');

    % create a kernel probability density function for each carboxysome
    [pdf_array, ~, group_names] = fitdist(distances', 'Kernel', 'By', carb_ids, 'Kernel', 'Normal', 'Width', bandwidth);

    % start a new plot
    figure;
    hold on;

    % Link the inner concentration data to a colormap
    cmap = colormap('winter');
    zmap = linspace(min(inner_concentrations), max(inner_concentrations), length(cmap));

    % for each carboxysome, make a plot
    for i = 1:length(group_names)
        % filter out distances and concentrations in other carboxysomes
        x_data = distances(carb_ids == str2double(group_names{i}));
        z_data = inner_concentrations(carb_ids == str2double(group_names{i}));

        % create 200 evenly spaced data points to sample the density plot
        density_x = linspace(min(x_data), max(x_data), 200);

        this_pdf = pdf(pdf_array{i}, density_x); % calculate the y values from the pdf
        plot_color = interp1(zmap, cmap, unique(z_data)); % set the color of the line
        plot(density_x, this_pdf, 'LineWidth', 2, 'color', plot_color); % plot the density
    end

    % Make some labels for the plot
    title('Probability Density of Distance Between RuBisCOs in Chains');
    xlabel('Distance Between Rubiscos (nm)');
    ylabel('Probability Density');

    % create the colorbar
    caxis([min(inner_concentrations), max(inner_concentrations)]); % set colorbar tick labels
    c = colorbar('Location', 'southoutside');
    c.Label.String = 'Inner Rubisco Concentration (\muM)'; % colorbar title
end