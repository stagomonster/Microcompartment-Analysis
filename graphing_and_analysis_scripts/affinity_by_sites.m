function affinity_by_sites(filename, carboxysome_data, label, min_chain_length)
% This function generates two graphs of the ratio of the number of bound to 
% unbound rubisco binding sites in the carboxysomes vs. the carboxysome 
% concentrations. One plot has total concentration on the x axis and the
% other has inner concentration on the x axis. It fits an exponential to
% the total concentration graph data and creates a .csv table with relevant
% carboxysome properties.
%
% Inputs
% filename - the name of the file that holds the raw rubisco data
% carboxysome_data - an array of carboxysome objects filled with data
%                    through at least the chain_maker.m script
% label - a string to put in the titles of the plots
% min_chain_length - for the inner concentration plot, only carboxysomes
%                    with at least one chain of length min_chain_length 
%                    will be plotted

    % this script requires the chain length to be 1
    addpath('../main_pipeline/');
    carboxysome_data = chain_maker(filename, carboxysome_data, 1);

    affinity_data = zeros(length(carboxysome_data), 7); % data to output in .csv file
    xData = []; % total concentration
    xData_inner = []; % inner concentration
    participation_ratios = []; % the fraction of rubiscos participating in a chain
    volumes = []; % carboxysome volume
    yData = []; % ratio of bound to unbound sites
    yData_inner = []; % ratio of bound to unbound sites

    % goes through each carboxysome, by number in the dataset and not index
    for i = 1:length(carboxysome_data)
        carb = carboxysome_data(i);
        rubiscos = carb.rubisco;
        % fills the first three columns with preliminary data
        affinity_data(i, 1) = carb.carb_index;
        affinity_data(i, 2) = carb.num_rubisco;
        affinity_data(i, 3) = carb.concentration;

        % used to plot the affinity figure later
        xData(end+1) = carb.concentration;
        xData_inner(end+1) = carb.inner_concentration;
        volumes(end+1) = carb.volume * 10^18; % convert from m^3 to micrometers^3

        bound = 0; % number of bound sites
        unbound = 0; % number of unbound sites
        longestChain = 0; % length of longest chain in carboxysome
        interior_unbound = 0; % number of unbound inner rubiscos
        num_singles = 0; % number of rubiscos not in a chain

        % uses the fact that each chain has one unbound site on each end
        % and the other sites are all bound, no matter the length of the
        % chain
        for chain = carb.chains
            indices = chain.indices;
            % check to get the end rubiscos of the chain
            if chain.length > 1
                ends = [rubiscos(indices(1)) rubiscos(indices(end))];
            else
                num_singles = num_singles + 1;
                ends = rubiscos(indices);
            end
            % knows that the two ending rubiscos are the only ones that
            % will have unbound sites
            for rubisco = ends
                if rubisco.inside
                    interior_unbound = interior_unbound+1;
                end
            end

            % the number of bound rubisco sites will be 2x - 2, where x is
            % the length of the chain
            bound = bound + ((chain.length * 2) - 2);

            % there is always only 2 unbound sites per chain
            unbound = unbound + 2;

            % update longest chain
            if chain.length > longestChain
                longestChain = chain.length;
            end
        end
        % used to plot affinity data later
        participation_ratios(end+1) = 1 - num_singles / carb.num_rubisco;
        yData(end+1) = bound / unbound;
        yData_inner(end+1) = bound / unbound;

        % fills in the rest of the data
        affinity_data(i, 4) = unbound;
        affinity_data(i, 5) = bound;
        total = unbound + bound;
        affinity_data(i, 6) = unbound / total;
        affinity_data(i, 7) = bound / total;
        affinity_data(i, 8) = longestChain;

        % for inner concentration plot, remove data from carboxysomes with
        % no chains of length min_chain_length
        if longestChain < min_chain_length
            xData_inner = xData_inner(1:end-1);
            yData_inner = yData_inner(1:end-1);
            participation_ratios = participation_ratios(1:end-1);
            volumes = volumes(1:end-1);
        end
    end

    % creates a .csv table of Carboxysome data
    affinity_table = array2table(affinity_data);
    affinity_table.Properties.VariableNames = {'reg (Carboxysome No.)', 'RubiscoTotal', 'Concentration', 'Unbound', 'Bound', 'Rel. Unbound', 'Rel. Bound', 'Longest Chain'};

    writetable(affinity_table, ['polymerizationTable' label '.csv'])

    ft = fittype('a*exp(b*x)');

    % first the total concentration plot data to an exponential
    opts = fitoptions('Method', 'NonlinearLeastSquares');
    opts.Lower = [0, 0];  % Lower bounds for a and b
    opts.Upper = [0.1, 0.1];  % Upper bounds for a and b

    [fitresult, gof] = fit(xData', yData', ft, opts);

    if ~isempty(label)
        label = ['(' label ')'];
    end

    % Bound/Unbound vs. Total Concentration Plot
    figure;
    subplot(2, 1, 1)
    scatter(xData, yData, "filled");
    set(gca, 'YScale', 'log'); % log scale y axis
    title(['Binding Affinity of Rubiscos ' label]);
    ylabel('Bound/Unbound');
    xlabel('Concentration (\muM)');
    hold on;

    % Bound/Unbound vs. Inner Concentration Plot
    subplot(2, 1, 2);
    bubblechart(xData_inner, yData_inner, volumes, participation_ratios);
    set(gca, 'YScale', 'log'); % log scale y axis
    title(['Binding Affinity of Rubiscos ' label]);
    ylabel('Bound/Unbound');
    xlabel('Inner Concentration (\muM)');
    c = colorbar('Location', 'westoutside'); % color points by fraction of rubiscos participating in chains
    c.Label.String = 'Ratio of RuBisCos Participating In Chain'; % colorbar title
    colormap('winter');
    bubblelegend('Carboxysome Volume \mum^3', 'Location', 'northeastoutside'); % size points by carboxysome volume
    
    hold off;

    % display fit results to the command window
    disp(gof);
    disp(fitresult);
end
