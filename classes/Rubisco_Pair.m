classdef Rubisco_Pair < handle
% Model for a pair of Rubisco
%
% Rubisco_Pair.m © 2025 is licensed under CC BY-NC-SA 4.0

    properties
        reg = nan; % carb index
        num_rubisco = nan; % number of rubisco in carboxysome
        num_rubisco_inner = nan; % number of inner rubisco in carboxysome
        volume = nan; % volume of carboxysome
        inner = false; % true if both rubisco are inner rubisco
        I_index = nan; % index of i rubisco
        J_index = nan; % index of j rubisco
        I_tag = nan; % tag of i rubisco
        J_tag = nan; % tag of j rubsico
        angle = nan; % angle between rubisco (degrees)
        distance = nan; % distance between rubisco centers in pixels
        concentration = nan; % concentration of carboxysome
        inner_concentration = nan; % inner concentration of parent carboxysome
        projection = nan; % projection of vector from target rubisco to source rubisco onto source rubisco's vector
    end

    methods
        function pair = Rubisco_Pair(carb_index, num_rubisco, num_rubisco_inner, ...
                vol, in, I_index, J_index, I_tag, J_tag, ang, dist, conc, inner_conc)
            % Constructor: create an instance of Rubisco Pair
            if nargin == 0
                return;
            end
            pair.reg = carb_index;
            pair.num_rubisco = num_rubisco;
            pair.num_rubisco_inner = num_rubisco_inner;
            pair.volume = vol;
            pair.inner = in;
            pair.I_tag = I_tag;
            pair.J_tag = J_tag;
            pair.I_index = I_index;
            pair.J_index = J_index;
            pair.angle = ang;
            pair.distance = dist;
            pair.concentration = conc;
            pair.inner_concentration = inner_conc;
        end

        function disp(self)
            % Print information about rubisco pair to command window
            disp("Rubisco_Pair object with properties:")
            props = properties(self);
            for i = 1:length(props)
                name = props{i};
                value = self.(name);
                if isnan(value)
                    value = 'NaN';
                end
                disp("  " + name + " = " + value);
            end
        end

        function obj_copy = copy(obj)
            % Create a copy of a rubisco pair object
            obj_copy = Rubisco_Pair();
            props = properties(obj);
            for i = 1:length(props)
                obj_copy.(props{i}) = obj.(props{i});
            end
        end

    end
end