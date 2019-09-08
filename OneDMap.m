classdef OneDMap < MapGeometry
    %OneDMap is descendant for piecewise linear map
    %   Constructor contains one argument which is number of nodes
    
    methods
        function map = OneDMap(N)
            % Create map
            map@MapGeometry(1);
            % Store size
            map.sizes = N;
            % Calculate internal coordinates of nodes
            map.internal = (1:N)';
            % Form array of links
            map.links = [map.internal(1:end-1), map.internal(2:end)];
            % Form array of ribs
            map.ribs = [map.internal(1:end-2), map.internal(2:end-1), map.internal(3:end)];
            % Set mapped coordinates to empty set
            map.mapped = [];
        end
        
        function newMap = extendPrim(map)
            % Get size of existing map
            n = map.sizes;
            % Create new map with greater size
            newMap = OneDMap(n + 2);
            % Copy known positions of nodes into new map
            newMap.mapped(2:n + 1, :) = map.mapped;
            % Calculate positions of the first and last nodes
            newMap.mapped(1, :) = 2 * newMap.mapped(2, :)...
                - newMap.mapped(3, :);
            newMap.mapped(n + 2, :) = 2 * newMap.mapped(n + 1, :)...
                - newMap.mapped(n, :);
        end
        
        function res = getBorder(map)
            res = [1, map.sizes];
        end
    end
    
end

