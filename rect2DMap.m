classdef rect2DMap < MapGeometry
    %OneDMap is descendant for piecewise linear map
    %   Constructor contains two argument which are number of rows and
    %   number of columns
    
    properties (SetAccess = protected)
        faces   %Array of faces to projections
    end
    
    methods
        function map = rect2DMap(rows, cols)
            if nargin < 2
                error('You MUST specify number of rows and columns for rect2DMap');
            end
            % Create map
            map@MapGeometry(2);
            % Store size
            map.sizes = [rows, cols];
            % Calculate internal coordinates of nodes
            N=rows*cols;
            a1=repmat(1:cols,rows,1);
            a2=repmat((1:rows)',1,cols);
            map.internal = [a1(:),a2(:)];
            % Form array of links
            A=reshape(1:N,rows,cols);
            B=A(1:end-1,:);
            C=A(2:end,:);
            D=A(:,1:end-1);
            E=A(:,2:end);
            map.links = [B(:), C(:); D(:), E(:)];
            % Form array of ribs
            B1=A(1:end-2,:);
            B2=A(2:end-1,:);
            B3=A(3:end,:);
            C1=A(:,1:end-2);
            C2=A(:,2:end-1);
            C3=A(:,3:end);
            map.ribs = [B1(:), B2(:), B3(:); C1(:), C2(:), C3(:)];
            % Form array of faces
            B1 = A(1:end-1,1:end-1);
            B2 = A(2:end,1:end-1);
            B3 = A(1:end-1,2:end);
            C1 = A(1:end-1,2:end);
            C2 = A(2:end,1:end-1);
            C3 = A(2:end,2:end);
            map.faces=[B1(:), B2(:), B3(:); C1(:), C2(:), C3(:)];
            % Set mapped coordinates to empty set
            map.mapped = [];
        end
        
        function face = getFaces(map)
            %Function to access to the faces of map.
            %face is k-by-3 matrix. Each row contains three numbers of
            %nodes which form one face.
            face=map.faces;
        end
        
        function newMap = extendPrim(map)
            % Get size of existing map
            n = map.sizes(1);
            m = map.sizes(2);
            N = n + 2;
            M = m + 2;
            NM = N * M;
            % Create new map with greater size
            newMap = rect2DMap(N, M);
            % Form list of old nodes
            ind = 1:NM;
            ind(newMap.getBorder) = [];
            % Copy known positions of nodes into new map
            newMap.mapped(ind, :) = map.mapped;
            % Calculate positions of the nodes added to left side
            ind = 2:n + 1;
            newMap.mapped(ind, :) = 2 * newMap.mapped(ind + N, :)...
                - newMap.mapped(ind + 2 * N, :);
            % Calculate positions of the nodes added to right side
            ind = NM - N + 2:NM - 1;
            newMap.mapped(ind, :) = 2 * newMap.mapped(ind - N, :)...
                - newMap.mapped(ind - 2 * N, :);
            % Calculate positions of the nodes added to bottom side
            ind = 1:N:NM;
            newMap.mapped(ind, :) = 2 * newMap.mapped(ind + 1, :)...
                - newMap.mapped(ind + 2, :);
            % Calculate positions of the nodes added to top side
            ind = N:N:NM;
            newMap.mapped(ind, :) = 2 * newMap.mapped(ind - 1, :)...
                - newMap.mapped(ind - 2, :);
        end
        
        function res = getBorder(map)
            % Get map sizes
            n = map.sizes(1);
            m = map.sizes(2);
            % Form list of border nodes
            res = [2:n - 1,...                  % left edge
                (m - 1) * n + 2:n * m - 1,...   % right edge
                (0:m - 1) * n + 1,...           % bottom edge
                (1:m) * n];                     % top edge
        end
    end
end

