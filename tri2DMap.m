classdef tri2DMap < MapGeometry
    %OneDMap is descendant for piecewise linear map
    %   Constructor contains two argument which are number of rows and
    %   number of columns. Each odd row contains cols nodes and each even
    %   row contains cols-1 nodes.
    
    properties (SetAccess = protected)
        faces   %Array of faces to projections
    end
    
    methods
        function map = tri2DMap(rows, cols)
            % Create map
            map@MapGeometry(2);
            % Store size
            map.sizes = [rows, cols];
            % Calculate internal coordinates of nodes
            vstep = sqrt(3)/2;
            % Calculate number of odd and even rows
            odd = round(rows/2 +0.001);
            even = rows - odd;
            N = odd*cols+even*(cols-1);
            % Calculate internal node coordinates.
            % X coordinates for two rows
            a1 = repmat([(1:cols)'; (1:cols-1)'+0.5],even,1);
            a2 = repmat([zeros(cols,1); repmat(vstep,cols-1,1)],1,even);
            a3 = (0:even-1)*(2*vstep);
            a2 = bsxfun(@plus,a2,a3);
            map.internal = [a1, a2(:)];
            if odd>even
                % Add last row
                map.internal = [map.internal; (1:cols)', repmat(even*2*vstep,cols,1)];
            end
            % Form array of links and ribs
            % Horizontal edges
            b = 2*cols-1;
            a1 = repmat((1:cols)',1,odd);
            a3 = (0:odd-1)*b;
            a1 = bsxfun(@plus,a1,a3);
            a2 = a1;
            B1 = a1(1:end-1,:);
            B2 = a1(2:end,:);
            map.links = [B1(:), B2(:)];
            B1 = a1(1:end-2,:);
            B2 = a1(2:end-1,:);
            B3 = a1(3:end,:);
            map.ribs = [B1(:), B2(:), B3(:)];
            a1=repmat((cols+1:2*cols-1)',1,even);
            a3=(0:even-1)*b;
            a1=bsxfun(@plus,a1,a3);
            B1=a1(1:end-1,:);
            B2=a1(2:end,:);
            map.links=[map.links; B1(:), B2(:)];
            B1 = a1(1:end-2,:);
            B2 = a1(2:end-1,:);
            B3 = a1(3:end,:);
            map.ribs = [map.ribs; B1(:), B2(:), B3(:)];
            % Form faces
            a3 = [reshape(a2(1:end-1,:),[],1);reshape(a1(1:end-1,:),[],1)];
            B1 = [a3, a3+1, a3+cols];
            ind = B1(:,2)<=N & B1(:,3)<=N;
            a3 = [reshape(a2(2:end-1,:),[],1);a1(:)];
            B2 = [a3,a3+cols,a3+(cols-1)];
            ind2 = B2(:,2)<=N & B2(:,3)<=N;
            map.faces = [B1(ind,:); B2(ind2,:)];
            % Left to right and up edges
            a3 = [reshape(a2(1:end-1,:),[],1);a1(:)];
            B1 = [a3,a3+cols];
            ind = B1(:,2)<=N;
            map.links = [map.links; B1(ind,:)];
            a3 = [reshape(a2(1:end-1,:),[],1);reshape(a1(1:end-1,:),[],1)];
            B1 = [a3,a3+cols,a3+2*cols];
            ind = B1(:,3)<=N;
            map.ribs = [map.ribs;B1(ind,:)];
            % Right to left up edges
            a3 = [reshape(a2(2:end,:),[],1);a1(:)];
            B1 = [a3,a3+(cols-1)];
            ind = B1(:,2)<=N;
            map.links = [map.links; B1(ind,:)];
            a3 = [reshape(a2(2:end,:),[],1);reshape(a1(2:end,:),[],1)];
            B1 = [a3,a3+(cols-1),a3+2*(cols-1)];
            ind = B1(:,3)<=N;
            map.ribs = [map.ribs;B1(ind,:)];
            % Set mapped coordinates to empty set
            map.mapped = [];
        end
        
        function face = getFaces(map)
            %Function to access to the feces of map
            %face is k-by-3 matrix. Each row contains three numbers of
            %nodes which form one face.
            face = map.faces;
        end
    
        function newMap = extendPrim(map)
            % Get size of existing map
            row = map.sizes(1);
            col = map.sizes(2);
            rowN = row + 2;
            colN = col + 1;
            % Calculate specific geometric information for new maps
            oddN = round(rowN / 2 + 0.001);
            evenN = rowN - oddN;
            NM = oddN * colN + evenN * (colN - 1);
            % Create new map with greater size
            newMap = tri2DMap(rowN, colN);
            % Creater array for mapped nodes
            newMap.mapped = zeros(size(newMap.internal, 1), size(map.mapped, 2));
            % Form list of old nodes
            ind = 1:NM;
            ind(newMap.getBorder) = [];
            % Copy known positions of nodes into new map
            newMap.mapped(ind, :) = map.mapped;
            
            % Expansion down left
            ind = 1:colN - 2;
            newMap.mapped(ind, :) = 2 * newMap.mapped(ind + colN, :)...
                - newMap.mapped(ind + 2 * colN, :);
            % Expansion down right
            ind = 3:colN;
            newMap.mapped(ind, :) = 2 * newMap.mapped(ind + colN - 1, :)...
                - newMap.mapped(ind + 2 * colN - 2, :);
            
            % Expansion up left
            if oddN > evenN
                ind = NM - colN + 1:NM - 2;
            else
                ind = NM - colN + 2:NM - 1;
            end
            newMap.mapped(ind, :) = 2 * newMap.mapped(ind - colN + 1, :)...
                - newMap.mapped(ind - 2 * colN + 2, :);
            % Expansion up right
            ind = NM - colN + 3:NM;
            newMap.mapped(ind, :) = 2 * newMap.mapped(ind - colN, :)...
                - newMap.mapped(ind - 2 * colN, :);
            
            % Expansion left
            od = oddN - 1 + evenN - oddN;
            ind = (1:od) * (2 * colN - 1) + 1;
            newMap.mapped(ind, :) = 2 * newMap.mapped(ind + 1, :)...
                - newMap.mapped(ind + 2, :);
            % Expansion right
            ind = (1:od) * (2 * colN - 1) + col + 1;
            newMap.mapped(ind, :) = 2 * newMap.mapped(ind - 1, :)...
                - newMap.mapped(ind - 2, :);
        end
        
        function res = getBorder(map)
            % Get map sizes
            row = map.sizes(1);
            col = map.sizes(2);
            
            % Calculate specific geometric information for old and new maps
            odd = round(row / 2 + 0.001);
            even = row - odd;
            N = odd * col + even * (col - 1);
            
            % Form list of border nodes
            od = odd - 1 + even - odd;
            res = [1:col,...                    % Bottom row
                N - col + 2 + even - odd:N,...  % Top row
                (1:od) * (2 * col - 1) + 1,...  % Left column
                (1:od) * (2 * col - 1) + col];  % Right column
        end
    end
end

