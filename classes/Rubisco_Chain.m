classdef Rubisco_Chain < handle
% Rubisco Chain Model
%
% Rubisco_Chain.m © 2025 is licensed under CC BY-NC-SA 4.0

    properties
        length = nan; % length of rubisco chain
        reg = nan; % carboxysome id
        index = nan; % chain id
        tag = nan;
        tags uint32 = uint32([]); % tags of rubiscos in chain
        indices uint32 = uint32([]); % indices of rubiscos in chain
        % Primarily testing two ways - tags using casting initialization, indices using direct initialization in the constructor
        average_vector = []; % average vector of rubiscos in chain
        centroid = []; % centroid of chain
        eigenvalues = []; % eigenvalues of the S tensor
        S = []; % chain orientation tensor
        S_val = nan; % scalar orientation value; describes how uniform the rubisco directions are
        type = ChainType.free; % Default value using chain enumeration
        xy_angle = nan; % measures angle of vector's projection in xy plane
        z_angle = nan; % measures angle between the xy plane and the vector
    end

    methods
        % Constructor: builds a rubisco chain object
        function chain = Rubisco_Chain(reg, index, tag, tags, indices, type)
            if nargin > 0
                chain.reg = reg;
                chain.tags = tags;
                chain.length = length(tags);
                chain.indices = uint32(indices);
                chain.type = type;
                chain.index = index;
                chain.tag = tag;
            end
        end

        function disp(self)
            % Prints information about a rubisco chain to the command
            % window
            disp("Rubisco_Chain object with properties:")
            props = properties(self);
            for i = 1:length(props)
                name = props{i};
                value = self.(name);
                if name == "type"
                    value = to_str(value);
                elseif isnan(value)
                    value = 'NaN';
                elseif ~isscalar(value)
                    value = mat2str(value);
                end
                disp("  " + name + " = " + value);
            end
        end

        function obj_copy = copy(obj)
            % creates a copy of a rubisco chain object
            obj_copy = Rubisco_Chain();
            props = properties(obj);
            for i = 1:length(props)
                obj_copy.(props{i}) = obj.(props{i});
            end
        end
    end
end
