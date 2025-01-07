classdef Chain_Pair < handle
    % Model for a pair of Rubisco Chains
    properties
        reg = nan; % carb index
        I_index = nan; % index of i chain
        J_index = nan; % index of j chain
        I_tag = nan; % tag of i chain
        J_tag = nan; % tag of j rubsico
        angle = nan; % angle between chains
        distance = nan; % distance between chains in pixels
        concentration = nan; % concentration of carboxysome 
    end

    methods
        function pair = Chain_Pair(carb_index, I_index, J_index, ...
                I_tag, J_tag, ang, dist, conc)
            % Constructor: create an instance of Rubisco Chain Pair
            if nargin > 0
                pair.reg = carb_index;
                pair.I_tag = I_tag;
                pair.J_tag = J_tag;
                pair.I_index = I_index;
                pair.J_index = J_index;
                pair.angle = ang;
                pair.distance = dist;
                pair.concentration = conc;
            end
        end

        function disp(self)
            % Print information about a chian pair object to the command
            % window
            disp("Chain_Pair object with properties:")
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
            % make a copy of a chain pair object
            obj_copy = Chain_Pair();
            props = properties(obj);
            for i = 1:length(props)
                obj_copy.(props{i}) = obj.(props{i});
            end
        end
    end
end