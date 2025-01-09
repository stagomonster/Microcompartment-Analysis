function [] = chain_lengths_analysis(carboxysome_data, min_chain_length, min_linkages)
% This function plots a histogram of chain lengths present in the dataset.
% It allows you to filter out shorter chains by setting a minimum chain
% length that a chain must meet to be plotted. It also allows you to filter
% out chains that aren't in lattices by setting a minimum number of
% linkages the chain must participate in to be plotted. The data points are
% colored according to the total rubisco concentration of their parent
% carboxysomes.
%
% Inputs
% carboxysome_data - an array of carboxysome objects filled with data at
%                    least through chain_maker.m
% min_chain_length - the minimum length of a chain that should still be
%                    included in this analysis
% min_linkages - the minimum number of linkages a chain must have to other
%                chains to be included in this analysis. Set this to zero 
%                if you want every chain in the dataset. Carboxysome
%                objects must be filled with data through at least
%                chain_linkages.m to use this filter with any value other
%                than zero.
%
% chain_lengths_analysis.m © 2025 is licensed under CC BY-NC-SA 4.0

    % Initialize arrays to hold the data
    chain_lengths = [];
    total_concentrations = [];

    % find the data in the carboxysome objects
    for carb = carboxysome_data % for each caroxysome
        for chain = carb.chains % for each chain
            if chain.length >= min_chain_length % filter out short chains
                %filter out chains not in lattices if user chooses
                if sum([carb.chain_links.I_index] == chain.index | [carb.chain_links.J_index] == chain.index) >= min_linkages
                    chain_lengths(end+1) = chain.length; % save the length of the chain
                    total_concentrations(end+1) = carb.concentration; % save the concentration of its parent carboxysome
                end
            end
        end
    end

    x = min_chain_length:max(chain_lengths); % x is what will be used to make the bar chart x axis
    plotz = NaN(length(x), length(chain_lengths)); % this array will be used to make the colorbar
    
    % Group the values in total_concentrations by the vertical bar their
    % corresponding chain_length values belong to. Each row in plotz 
    % represents a vertical bar. Each column in plotz represents a layer in
    % the bar.
    for i = 1:length(chain_lengths)
        plotz(x == chain_lengths(i), i) = total_concentrations(i);
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

    z_values = unique(sort(total_concentrations, 'descend'), 'stable'); % all the unique z values from large to small
    plotdata = zeros(length(x), length(z_values)); % the data that will be bar heights

    % calculates how many identical values exist in each row of plotz and
    % condenses them into one column. If there are n copies of a value, 
    % then num_repeats will have value n for that row and store n in some
    % position in plotdata
    for z_value = z_values
        num_repeats = sum(plotz == z_value, 2); % the number of times z_value appears in each row of plotz
        plotdata(:, z_values == z_value) = num_repeats; % each column of plotdata contains instances of the same z_value
    end

    figure;
    b = bar(x, plotdata, 'stacked', 'FaceColor', 'flat'); % make a stacked bar graph
    cmap = colormap('winter');
    zmap = linspace(min(z_values), max(z_values), length(cmap));

    % Color each data point based on where it is between the min and max
    for i = 1:length(b)
        % make the bar's color proportional to its z value's distance between z_min and z_max
        b(i).CData = interp1(zmap, cmap, z_values(i));
        b(i).EdgeColor = 'none'; % remove the edges of the bars
    end

    % Make some labels for the plot
    if min_linkages == 0
        title(['Chains of Length >=', num2str(min_chain_length), '']);
    else
        title(['Chains of Length >=', num2str(min_chain_length), ' with >=', num2str(min_linkages), ' Linkages']);
    end
    xlabel('Chain Length (number of RuBisCos)');
    ylabel('Number of Chains');

    % Create and edit the colorbar
    c = colorbar('Location', 'southoutside'); % colorbar location
    c.Label.String = 'Rubisco Concentration(\muM)'; % colorbar title
    c.Ticks = [0, 0.25, 0.5, 0.75, 1];
    c.TickLabels = round(linspace(min(total_concentrations), max(total_concentrations), 5), 0); % colorbar labels
end