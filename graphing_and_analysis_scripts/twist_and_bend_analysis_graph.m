function [] = twist_and_bend_analysis_graph(carboxysome_data, min_chain_length)
% This function plots the twist angle between two rubiscos linked in a 
% chain as a function of the bend angle between those rubiscos. The data
% points are colored according to the inner rubisco concentration of their
% parent carboxysome.
%
% Inputs
% carbosyxome_data - an array of carboxysome objects run through the main
%                    pipeline until at least chain_maker.m
% min_chain_length - the minimum length of a chain whose data you want
%                    plotted

    % create arrays to hold the data that will be plotted
    bends = [];
    twists = [];
    inner_concentrations = [];

    % for each carboxysome in the dataset
    for carb = carboxysome_data 
        % for each chain with length >= min_chain_length
        for chain = carb.chains([carb.chains.length] >= min_chain_length)
            % for each rubisco linkage in the chain
            for i = 1:length(chain.indices) - 1

                % get the two adjacent rubisco objects
                rubisco_i = carb.rubisco(chain.indices(i));
                rubisco_j = carb.rubisco(chain.indices(i+1));

                % calculate their bend, twist, and inner concentration and
                % save it to the arrays to be plotted
                bends(end+1) = calc_bend(rubisco_i, rubisco_j);
                twists(end+1) = calc_twist(rubisco_i, rubisco_j, chain.average_vector);
                inner_concentrations(end+1) = carb.inner_concentration;
            end
        end
    end

    % Plot the data in a scatter plot
    figure;
    scatter(bends, twists, [], inner_concentrations, "filled", 'MarkerFaceAlpha', 0.75);
    % plot title
    title(['Relative Orientations of Consecutive Rubiscos in Chains of Length >= ', num2str(min_chain_length), '']);
    xlabel('Bending Angle to Next RuBisCo in the Chain (deg)', 'fontsize', 14); % x axis label
    ylabel('Twist (deg)', 'fontsize', 14); % y axis label
    xlim([-inf inf]); % set x and y limits
    ylim([-inf inf]);
    colormap('winter'); % set the colormap
    c = colorbar('fontsize', 14); % make a colorbar
    c.Label.String = 'Inner RuBisCo Concentration (\muM)'; % colorbar title
end