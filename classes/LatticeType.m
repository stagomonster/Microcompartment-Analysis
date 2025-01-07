classdef LatticeType
    % Model for types of rubisco lattices
    enumeration
        full % rubisco chain is paired with 6 or more other chains
        incomplete % rubisco chain is paired with 4-5 other chains
        multiple % carboxysome has multiple lattices
        shell % most chains in carboxysome are edge chains
        none % carboxysome has no chains
    end

    methods
        % Prints information about lattice type to the command window
        function disp(obj)
            % Get the name of the enumeration value
            name = char(obj);
            
            % Display the class name and the specific enumeration value
            fprintf('LatticeType: %s\n', name);
        end

        function str = to_str(obj)
            % converts an object to a string
            str = "LatticeType: " + char(obj);
        end
    end
end