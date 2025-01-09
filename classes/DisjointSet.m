classdef DisjointSet < handle
% data structure used to classify rubisco chains and later lattices
%
% DisjointSet.m © 2025 is licensed under CC BY-NC-SA 4.0

    properties
        % list of every element's parent in the disjoint set
        parent
        % list of how many children are in each root's sets
        rank
        % how many are in the total set
        size
        % actual set elements
        all
    end
    
    methods
        function set = DisjointSet(numElements)
        % generates a disjoint set of user specified size
            set.parent = zeros(1, numElements);
            set.all = zeros(1, numElements);
            % all elements have rank 1 at the beginning
            set.rank = ones(1, numElements);
            set.size = numElements;
            for i = 1:numElements
                % all elements' parent is themselves at initialization
                set.parent(i) = i;
                set.all(i) = i;
            end
        end

        function root = locate(set, x)
            % finds the most fundamental parent of an element
            if set.parent(x) ~= x
                set.parent(x) = set.locate(set.parent(x)); % concatonates the set to refer to its root as its parent
            end
            % when it finally finds an element that refers to itself as its own parent (the root)
            root = set.parent(x);
        end
        
        function unite(set, x, y)
        % used to put two sets in the disjoint set together
            rootX = set.locate(x);
            rootY = set.locate(y);
            if rootX == rootY
                return;
            end
            % makes the root of the smaller set's parent be the root of the
            % larger set and adds the rank of the smaller set to the larger set's
            if set.rank(rootX) < set.rank(rootY)
                set.parent(rootX) = rootY;
                set.rank(rootY) = set.rank(rootY) + set.rank(rootX);
                set.rank(rootX) = 1;
            else
                set.parent(rootY) = rootX;
                set.rank(rootX) = set.rank(rootX) + set.rank(rootY);
                set.rank(rootY) = 1;
            end
        end

        function children = find_children(set, parent)
        % used to obtain all members of a single set under one root
            children = [];
            for element = set.all
                root = set.locate(element);
                if parent == root
                    children(end+1) = element;
                end
            end
        end

        function chains = chain_builder(set)
            % builds a rubisco chain for each unique root
            chains = {};
            chain_map = containers.Map('KeyType', 'int32', 'ValueType', 'any');
    
            % puts the chains in a Map by the root of their set
            for element = 1:set.size
                root = set.locate(element);
                % makes the root into its own unique key and adds the root
                % as the first element in this chain
                if ~isKey(chain_map, root)
                    chain_map(root) = int32(root);
                    % and then adds the element that refers to the root to
                    % its map too
                    if element ~= root
                        chain_map(root) = [chain_map(root) element];
                    end
                % otherwise just adds the element to the root's map values
                else
                    if element ~= root
                        chain_map(root) = [chain_map(root) element];
                    end
                end
            end

            % combine all the chains into a single 2D cell array
            keys_list = keys(chain_map);
            for i = 1:length(keys_list)
                key = keys_list{i};
                chains{end+1} = chain_map(key);
            end
        end
    end
end
