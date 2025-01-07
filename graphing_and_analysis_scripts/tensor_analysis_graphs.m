function tensor_analysis_graphs(carboxysomes, min_num_rubiscos)
% This function plots the scalar S orientation parameter vs 
% concentration/inner concentration for an array of Carboxysomes containing
% the necessary data.
%
% Inputs
% caroxysomes - an array of carboxysome objects that have been processed
%               through at least local_global_alignment.m
% min_num_rubiscos - the minimum number of rubisco a carboxysome can have 
%                    that you still want to include in this plot
    
    % Filter out carboxysomes with less than min_num_rubiscos rubiscos
    carbs_to_include = carboxysomes([carboxysomes.num_rubisco] >= min_num_rubiscos);

    figure;
    sgtitle("Orientation Tensor Analysis of Carboxysomes With At Least " +num2str(min_num_rubiscos) + " RuBisCOs");
    
    % inner concentration plot
    subplot(2, 1, 1);
    hold on;
    title('Scalar Orientation Parameter S vs. Inner RuBisCO Concentration');
    s = [carbs_to_include.S_val];
    inner_concentration = [carbs_to_include.inner_concentration];
    scatter(inner_concentration, s, 'filled');
    xlabel('Inner RuBisCO Concentration ($\mu M$)', 'fontsize',14,'interpreter','latex');
    ylabel('Orientation Parameter S', 'fontsize',14,'interpreter','latex');
    xlim([-inf inf]);
    ylim([-inf inf]);
    
    % total concentration plot
    subplot(2, 1, 2);
    hold on;
    title('Scalar Orientation Parameter S vs. RuBisCO Concentration');
    s = [carbs_to_include.S_val];
    concentration = [carbs_to_include.concentration];
    scatter(concentration, s, 'filled');
    xlabel('RuBisCO Concentration ($\mu M$)', 'fontsize',14,'interpreter','latex');
    ylabel('Orientation Parameter S', 'fontsize',14,'interpreter','latex');
    xlim([-inf inf]);
    ylim([-inf inf]);
end
