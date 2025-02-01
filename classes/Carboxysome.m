classdef Carboxysome < handle
% Model for Carboxysome
%
% Carboxysome.m © 2025 is licensed under CC BY-NC-SA 4.0

    properties
        carb_index uint16 = nan; % carboxysome tag
        rubisco = Rubisco.empty; % array of rubisco objects in the carboxysome
        num_rubisco uint16 = nan; % number of rubisco in carboxysome
        num_rubisco_inner uint16 = nan; % number of rubisco inside the boundary of the convex hull
        num_rubisco_outer uint16 = nan; % number of rubisco on the boundary of the convex hull
        tags_inside uint32 = uint32([]); % tags of the rubisco inside the boundary 
        convex_hull = []; % traingulation matrix for the convex hull
        vertices = []; % coordinates of the vertices of the convex hull. Each row is a row vector corresponding to the position of the vertex in R^3
        coordinates = []; % coordinates of the rubisco within the carboxysome.
        volume = nan; % volume of carboxysome
        concentration = nan; % concentration of rubisco in carboxysome (micro-molar)
        inner_concentration = nan; % concentration of rubisco within convex hull boundary (micro-molar)
        num_events uint16 = nan; % number of rubiscos which are neighbors and well aligned
        rubisco_pairs = Rubisco_Pair.empty; % pairs of Rubisco
        normals = []; % each row is a normal vector to the corresponding facet row in convex_hull
        centroids = []; % each row is the centroid of the facet in the corresponding row in convex_hull
        S = []; % S tensor describing global orientation
        eigenvalues = []; % eigenvalues of S tensor
        S_val = nan; % scalar orientation parameter; describes how similarly oriented the rubisco are
        ave_vector = []; % average vector of carboysome
        links = Rubisco_Pair.empty; %Rubsico Pairs that are actually linked together in a chain
        chains = Rubisco_Chain.empty; % list of rubisco chain objects in the carboxysome
        chain_links = Chain_Pair.empty; % list of chain pair objects in the carboxysome
        lattice = Rubisco_Chain.empty; % list of lattices in the carboxysome
        max_connections uint16 = 0; % max number of pairs to a single chain
        tomo uint16 = nan; % tomogram index this carboxysome belongs to
        lattice_type = LatticeType.none; % Type of carboxysome lattice
    end

    methods
        % Constructor: create an instance of a Carboxysome
        function carb = Carboxysome(index, rubiscos, num_rubiscos, num_rubiscos_inner, num_rubiscos_outer, ...
                tags_in, conv_hull, vol, verts, cords, conc, inner_conc, event, norms, cents, s, ave, pairs, links)
            if nargin == 0
                return;
            else
                carb.carb_index = index;
                carb.rubisco = rubiscos;
                carb.num_rubisco = num_rubiscos;
                carb.num_rubisco_inner = num_rubiscos_inner;
                carb.num_rubisco_outer = num_rubiscos_outer;
                carb.tags_inside = tags_in;
                carb.convex_hull = conv_hull;
                carb.volume = vol;
                carb.vertices = verts;
                carb.coordinates = cords;
                carb.concentration = conc;
                carb.inner_concentration = inner_conc;
                carb.num_events = event;
                carb.normals = norms;
                carb.centroids = cents;
                if nargin > 15
                    carb.S = s;
                    carb.eigenvalues = eigs(s);
                    carb.S_val = (3/2)*max(carb.eigenvalues);
                    carb.ave_vector = ave;
                    carb.rubisco_pairs = pairs;
                    carb.links = links;
                end
            end
        end

        % Displays information contained by the carboxysome
        function disp(self)
            disp("Carboxysome object with properties:")
            props = properties(self);
            for i = 1:length(props)
                name = props{i};
                value = self.(name);
                if name == "rubisco"
                    value = "Rubisco array of length " + string(length(value));
                elseif name == "rubisco_pairs" || name == "links"
                    value = "Rubisco_Pair array of length " + string(length(value));
                elseif name == "chains"
                    value = "Rubisco_Chain array of length " + string(length(value));
                elseif name == "chain_links" || name == "lattice"
                    value = "Chain_Pair array of length " + string(length(value));
                elseif name == "lattice_type"
                    value = to_str(value);
                elseif isnan(value)
                    value = "NaN";
                elseif ~isscalar(value)
                    value = mat2str(value);
                end
                disp("  " + name + " = " + value);
            end
        end

        function obj_copy = copy(obj)
            % Create a copy of the carboxysome so that we can modify new
            % objects without changing previously created ones.
            obj_copy = Carboxysome();
            props = properties(obj);
            obj_props = ["rubisco" "rubisco_pairs" "links" "chains" "chain_links" "lattice"];
            for i = 1:length(props)
                % Create deep copies of the properties that hold arrays of
                % other cutsom classes' objects.
                if ismember(props{i}, obj_props)
                    obj_copy.(props{i}) = arrayfun(@(x) x.copy(), obj.(props{i}));
                else
                    obj_copy.(props{i}) = obj.(props{i});
                end
            end
        end

        function coords = get_coords(self)
        % finds the mean of the vertices, which approximates the centroid
            coords = mean(self.vertices);
        end

        function visualize(self)
        % visualize the convex_hull of the carboxysome
            title_string = 'Convex Hull for Carboxysome with index = ';
            figure;
            axis equal;
            hold on;
            trisurf(self.convex_hull, self.vertices(:, 1), self.vertices(:, 2), self.vertices(:, 3), 'FaceAlpha', 0.1);
            plot3(self.coordinates(:, 1), self.coordinates(:, 2), self.coordinates(:, 3), '*', 'Color', 'k');
            title([title_string, num2str(self.carb_index)]);
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            hold off;
        end
       
        function visualize_normals(self)
        % visualize a carboxysome with all normal vectors on facets of the
        % convex hull.
        title_string = 'Convex Hull for Carboxysome with index = ';
        norms = self.normals;
        cents = self.centroids;
        norms = norms*500;
        figure;
        axis equal;
        hold on;
        trisurf(self.convex_hull, self.vertices(:, 1), self.vertices(:, 2), self.vertices(:, 3), 'FaceAlpha', 0.1);
        plot3(self.coordinates(:, 1), self.coordinates(:, 2), self.coordinates(:, 3), '*', 'Color', 'k');
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        for k = 1:size(norms, 1)
            hold on;
            plot3(cents(k, 1), cents(k, 2), cents(k, 3), 'o', 'Color', 'r');
            quiver3(cents(k, 1), cents(k, 2), cents(k, 3), norms(k, 1), norms(k, 2), norms(k, 3), 'Color', 'k');
        end
        title([title_string, num2str(self.carb_index)]);
        hold off;
        end
        
    end

end