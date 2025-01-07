function concentration_volume_graphs(carboxysome_data)
% This funciton creates two plots. First it plots the number of rubiscos on
% the outer layer of a carboxysome as a function of carboxysome volume. The
% data in this first plot are colored according to the rubisco 
% concentration of the carboxysome. Second, it plots the number of rubiscos
% on the outer layer of a carboxysome as a function of rubisco
% concentration.
%
% Inputs
% carboxysome_data - an array of caroxysome objects filled with data at
%                    least through convex_hull_and_volume.m

    % get the data to be plotted from the array of carboxysomes
    vol = ([carboxysome_data.volume]')*(10^18); % in cubic micrometers
    num_outer = [carboxysome_data.num_rubisco_outer]';
    conc = [carboxysome_data.concentration]'; % in micromolar

    % print plots of concentration vs. numbers
    figure;
    sgtitle('Population of the Shell-Adjacent Rubisco Layer');
   
    %% First Subplot
    % first subplot plots # of outer rubisco vs. carboxysome volume
    subplot(2, 1, 1);
    scatter(vol, num_outer, [], conc, 'filled');
    hold on;

    % give first subplot a title, axis labels, and a colorbar
    title('Carboxysome Volume vs. Count of RuBisCOs in Outer Layer');
    xlabel('Carboxysome Volume ($\mu m^3$)', 'fontsize',14,'interpreter','latex');
    ylabel('Count of RuBisCOs in Outer Layer');
    colormap('winter');
    bar = colorbar;
    set(get(bar,'Title'),'String','Concentration ($\mu M$)', 'interpreter', 'latex');
    hold off;
    
    %% Second Subplot
    % second subplot plots # of outer rubisco vs. rubisco concentration
    subplot(2, 1, 2);
    scatter(conc, num_outer, 'filled');
    hold on;

    % give second subplot a title and axis labels
    title('RuBisCO Concentration vs. Count of RuBisCOs in Outer Layer');
    xlabel('RuBisCO Concentration ($\mu M$)', 'fontsize',14,'interpreter','latex');
    ylabel('Count of RuBisCOs in Outer Layer');
    hold off;
end