classdef ChainType
% Model for types of chains
%
% ChainType.m © 2025 is licensed under CC BY-NC-SA 4.0

    enumeration
        attached % chain has 1 or 2 outer rubiscos
        free % entirely made of inner rubiscos
        edge % has 3 or more outer rubiscos
    end

    methods
        % Prints information about chain's type
        function disp(obj)
            % Get the name of the enumeration value
            name = char(obj);
             
            % Display the class name and the specific enumeration value
            fprintf('ChainType: %s\n', name);
        end

        function str = to_str(obj)
            % converts an object to a string
            str = "ChainType: " + char(obj);
        end
    end
end
