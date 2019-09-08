function drawMap(map, data, varargin)
%drawMap create figure and draw data and maps. If data dimension (number of
%columns) is less than 3 then original data coordinate system is used. If
%data dimension is greater than 3 then projection into space of the first
%three PCs is used for drawing.
%
%Usage
%   drawMap(map, data);
%   drawMap(__, Name, Value);
%
%Inputs:
%   map is object of class MapGeometry (descendant of this class).
%   data is n-by-dim matrix of data points. Each row is one point.
%   There are several possible 'Name', value pairs:
%       'axes' is array of three different nonzero integer. This argument
%           can be used for space with more than three coordinates only.
%           The first number is vector for x axis, the second number is
%           vector for y axis and the third number is vector for z axis.
%           Positive number means number of coordinate in original space.
%           Negative number mean specified principal component.
%       'classes' is n-by-1 vector of class labels for points. Each label
%           is positive integer number which is the index of cell array
%           with marker descriptions. Marker descriptions can be specified
%           by user in the arguments markColour, markShape, markSize.
%           Otherwise standard marker is used for all points.
%       'markColour' is K-by-1 vector of standard Matlab colours ('r', 'g',
%           'b', 'y', 'm', 'c', 'w', 'k'). Value K is number of defined
%           markers. If markColour is omitted then 'b' (blue) is used for
%           all markers. 
%       'markShape' is K-by-1 vector of standard Matlab marker shapes ('o',
%           '+', '*', '.', 'x', 's', 'd', '^', 'v', '<', '>', 'p', 'h').
%           Value K is number of defined markers. If markShape is omitted
%           then 's' (square) is used for all markers.
%       'markSize' is K-by-1 vector of positive numbers. Value K is number
%           of defined markers. If markSize is omitted then 6 is used for
%           all markers.
%       'nodeMarker' is one of the possible marker shape symbol ('o', '+',
%           '*', '.', 'x', 's', 'd', '^', 'v', '<', '>', 'p', 'h') or
%           'none'. Default value is 'none'.
%       'nodeMarkerSize' is positive number which is size of node marker.
%           Default value is 6;
%       'nodeColour' is one of the possible Matlab colours ('r', 'g', 'b',
%           'y', 'm', 'c', 'w', 'k'). Default value is 'r';
%       'lineWidth' is non-negative number for map line width. Zero means
%           absence of map at all (width of line is zero and 'nodeMarker'
%           is 'none'. Default value is 0.5.
%       'newFigure' is logical argument. Value true (default) causes
%           creation of new figure. Value false causes usage of current
%           active figure. This option can be used, for example, for
%           subplots.
%       'coloring' defines the type of data to colour. It can have
%               following values: 
%           [] (empty) means no colouring. It is default value.
%           'density' is density colouring. It is equivalent to vector of
%               ones.
%           fun is function handle of form function res = fun(X), where X
%               is n-by-d matrix with one data point per row and res is
%               n-by-1 vector with values to use. Coordinate of vector X
%               are defined in PREPROCESSED space.
%           k is positive integer number. In this case k is number of
%               coordinate to use.
%           k is negative integer number. In this case map is coloured by
%               value of projection on the k'th principal component.
%           vect is n-by-1 vector with data defined function. Each element
%               of vector corresponds to data point in matrix data.
%           matr is N-by-(d+1) matrix with one data point in the first d
%               columns of each row and value of function of this point in
%               the last column. For example to calculate density function
%               it is necessary to send data matrix with 1 in (d+1) column. 
%       'ColorMap' is colormap to use for current graph. This parameter has
%           meaning for maps with colouring only. Colormap is standard
%           Matlab colormap. It can be one of predefined colormaps or
%           user defined. For more details find Colormap in the Matlab
%           documentation. 
%       'ColoringInSpace' is logical attribute. True value means
%           calculation of density or another data defined colouring
%           without projection onto map. False (default) value means
%           firstly projection of data onto map and then calculation of
%           colouring.
%       'Smooth' is positive real number. It is parameter of smoothing of
%           function interpolation for densities and other functions
%           defined in the data points. Default value is 2. Data
%           interpolation is implemented in radial basis function like
%           fasion: for query point x we calculate squared distance to each
%           data point d(i). Value of function in the point x is sum of
%           exp(-d(i)/(smooth*s2)) where s2 is mean of variances of all
%           attributes.
%       'projType' is type of projection onto map. Projection is necessary
%           to create colouring from the function defined in the data
%           points including densities. projType can have following values: 
%               0 is projection to the nearest node,
%               1 is projection to the nearest edge,
%               2 is projection to the nearest face.
%           default value is 0.
%

    % Data preprocessing
    data = map.preprocessData(data);

    %Get data dimension
    [N, dim] = size(data);

    % Parse varargin
    classes = [];
    markColour = [];
    markShape = [];
    markSize = [];
    nodeMarker = 'none';
    nodeMarkerSize = 6;
    nodeColour = 'r';
    lineWidth = 0.5;
    newFigure = true;
    source = [];
    colMap = parula;
    ColouringInSpace = false;
    nodeMarkerSpecified = false;
    lineWidthSpecified = false;
    smooth = 0.15;
    projType = 0;
    axeS = [];
    
    for i=1:2:length(varargin)
        switch lower(varargin{i})
            case 'axes'
                axeS = varargin{i + 1};
            case 'classes'
                classes = varargin{i + 1};
            case 'markcolour'
                markColour = varargin{i + 1};
            case 'markshape'
                markShape = varargin{i + 1};
            case 'marksize'
                markSize = varargin{i + 1};
            case 'nodemarker'
                nodeMarker = varargin{i + 1};
                nodeMarkerSpecified = true;
            case 'nodemarkersize'
                nodeMarkerSize = varargin{i + 1};
                nodeMarkerSpecified = true;
            case 'nodecolour'
                nodeColour = varargin{i + 1};
            case 'linewidth'
                lineWidth = varargin{i + 1};
                lineWidthSpecified = true;
            case 'newfigure'
                newFigure = varargin{i + 1};
            case 'coloring'
                source = varargin{i + 1};
            case 'colormap'
                colMap = varargin{i + 1};
            case {'coloringinspace', 'colouringinspace'}
                ColouringInSpace = varargin{i + 1};
            case 'smooth'
                smooth = varargin{i + 1};
            case 'projtype'
                projType = varargin{i + 1};
            otherwise
                error(['Unknown argument "', varargin{i},...
                    '" at position ', num2str(i + 2)]);
        end
    end
    
    % Sanity check of arguments.
    if isempty(classes)
        classes = ones(N, 1);
        cls = 1;
    else
        cls = unique(classes);
    end
    
    nCls = length(cls);
    
    if isempty(markColour)
        markColour = repmat('b', nCls, 1);
    end

    if isempty(markShape)
        markShape = repmat('s', nCls, 1);
    end
    
    if isempty(markSize)
        markSize = repmat(6, nCls, 1);
    end
    
    if lineWidth == 0
        mapDraw = 0;
    else
        mapDraw = 1;
    end
    
    if ~isempty(axeS)
        if length(axeS) ~= 3 || ~isnumeric(axeS) || any(axeS == 0)...
                || length(unique(axeS)) ~= 3
            error(['"axes" must be array of three different nonzero',...
                ' integers. This argument can be used for space',...
                ' with more than three coordinates only. The',...
                ' first number is vector for x axis, the second',...
                ' number is vector for y axis and the third',...
                ' number is vector for z axis. Positive number',...
                ' means number of coordinate in original space.',...
                ' Negative number mean specified principal',...
                ' component.']);
        end
    end
                
    %Create figure
    if newFigure
        figure;
    end
    %Get map coordinates
    maps = map.getMappedCoordinates;

    % Check colouring parameters
    if ~isempty(source)
        % Is map appropriate for colouring
        if ~any(strcmp(methods(map), 'getFaces'))
            error('Map must implement method "getFaces" for colouring');
        end
        
        if map.getDimension() ~= 2
            error('Map coloring can be used for 2D maps only');
        end
        
        if ~isnumeric(smooth) || smooth <= 0
            error('Smooth must be positive number.');
        end
        
        % Get internal map nodes' coordinates
        intern = map.getInternalCoordinates;
        
        % Extract grid
        [gridReady, nodeMap, nodeInt] = formGrid(map.getFaces, maps, intern);
        
        % Check the correctness of source parameter
        % Identify type of source
        tmp = false;
        if isnumeric(source) 
            if isscalar(source)
                source = round(source);
                % It is number of coordinate or PCs
                source = round(source);
                if source > 0
                   % It is coordinate
                   % Get data in original space
                   temp = map.deprocessData(nodeMap);
                   if source > size(temp, 2)
                       error(['Number of coordinate (value of "coloring")',...
                           ' %d to draw must be positive and cannot be',...
                           ' greater than data space dimension %d'],...
                           source, size(temp, 2));
                   end
                   f = temp(:, source);
                elseif source == 0
                    error(['Number of coordinate (value of "coloring")',...
                        ' to draw must be positive and number of',...
                        ' principal component must be negative.',...
                        ' Value of "coloring" cannot be zero']);
                else
                    if -source > size(map.PCs, 2)
                       error(['Number of principal component (MINUS',...
                           ' value of "source") %d to draw',...
                           ' must be positive and cannot be greater',...
                           ' than %d which is the number of',...
                           ' principal component calculated at map',...
                           ' initialisation'],...
                           -source, size(map.PCs, 2));
                    end
                    if map.preproc
                        % Get required coordinate
                        f = nodeMap(:, -source);
                    else
                        % Calculate projection on required PC
                        f = bsxfun(@minus, nodeMap, map.means)...
                            * map.PCs(:, -source);
                    end
                end
            elseif isvector(source)
                source = source(:);
                if size(source, 1) ~= N
                    error(['Number of elements in vector "coloring"',...
                        ' %d must coincides with number of data points %d'],...
                        size(source, 1), N);
                end
                
                if ColouringInSpace
                    f = interpol(data, source, nodeMap, smooth);
                else
                    f = interpol(map.project(data, projType, 'internal'),...
                        source, nodeInt, smooth);
                end
            elseif ismatrix(source)
                temp = map.preprocessData(source(:, 1:d));
                if size(temp, 2) ~= dim
                    error(['Wrong dimension of matrix in source argument\n',...
                        'Matrix of data points with function to draw must be\n',...
                        'n-by-(d+1) matrix with one data point in the\n',...
                        'first d columns of each row and value of function\n',...
                        'of this point in the last column. For example\n',...
                        'to calculate density function it is necessary\n',...
                        'to send data matrix with 1 in (d+1) column.\n%s'],'');
                end
                % Calculate influence of data points in nodes
                if ColouringInSpace
                    f = interpol(temp, source(:, end),...
                        nodeMap, smooth);
                else
                    f = interpol(...
                        map.project(temp, projType, 'internal'),...
                        source(:, end), nodeInt, smooth);
                end
                clear temp;
            else
                tmp = true;
            end
        elseif ischar(source) && strcmp(source, 'density')
            if ColouringInSpace
                f = interpol(data, ones(N, 1), nodeMap, smooth);
            else
                f = interpol(map.project(data, projType, 'internal'),...
                    ones(N, 1), nodeInt, smooth);
            end
        elseif isa(source, 'function_handle')
            % Get values to draw
            f = source(nodeMap);
        else
            tmp = true;
        end
    
        % Throw error is necessary
        if tmp
            error(['Wrong type of source argument. Source must be',...
                ' either\n',...
                '[] (empty) means no colouring. It is default value.\n',...
                '''density'' is density colouring. It is equivalent to vector of\n',...
                '     ones.\n',...
                'fun is function handle of form function res = fun(X), where X\n',...
                '     is n-by-d matrix with one data point per row and res is\n',...
                '     n-by-1 vector with values to use. Coordinates of vector X\n',...
                '     are defined in original space even for map for preprocessed\n',...
                '     data.\n',...
                'k is positive integer number. In this case k is number of\n',...
                '     coordinate to use.\n',...
                'k is negative integer number. In this case map is coloured by\n',...
                '     value of projection on the k''th principal component.\n',...
                'vect is n-by-1 vector with data defined function. Each element\n',...
                '     of vector corresponds to data point in matrix data.\n',...
                'matr is N-by-(d+1) matrix with one data point in the first d\n',...
                '     columns of each row and value of function of this point in\n',...
                '     the last column. For example to calculate density function\n',...
                '     it is necessary to send data matrix with 1 in (d+1) column.\n',...
                ' For example to calculate density function\n',...
                '     it is necessary to send data matrix with 1 in',...
                ' (d+1)     column.\n%s'],'');
        end
        
        %Including of colouring by default removes drawing of nodes and
        %sets map edges line width to 0.5 if other options are not
        %specified by user.   
        %
        %3D map colouring automatically excludes data drawing, removes
        %nodes and sets map edges line width to 0.5 if other options are
        %not specified by user.  
        if ~nodeMarkerSpecified 
            nodeMarker = 'none';
        end
        if ~lineWidthSpecified
            lineWidth = 0.5;
        end
        if dim > 3 
            % Dimension is greater than 3 and we need to use PCs or
            % specified coordinates
            if isempty(axeS)
                % By default Project data onto three PCs if it is necessary
                if ~map.preproc
                    V = map.PCs(:, 1:3);
                    nodeMap = bsxfun(@minus, nodeMap, map.means) * V;
                end
            else
                % Create copy of node coordinates
                tmp = nodeMap;
                % Create new array for node coordinates
                nodeMap = zeros(size(nodeMap, 1), 3);
                % Check the necessity of original coordinates
                ind = axeS > 0;
                if sum(ind) > 0
                    % At least one coordinate will be used. We need to
                    % restore it
                    temp = map.deprocessData(tmp);
                    if any(axeS > size(temp, 2))
                        error(['Number of coordinate cannot be greater',...
                            ' than dimension of original space']);
                    end
                    nodeMap(:, ind) = temp(:, axeS(ind));
                end
                ind = ~ind;
                % Check the necessity of PCs
                if sum(ind) > 0
                    % At least on PCs will be used
                    if any(axeS < -size(map.PCs, 2))
                       error(['Number of principal component (MINUS',...
                           ' value of "axes") %d %d %d to use',...
                           ' must be positive and cannot be greater',...
                           ' than %d which is the number of',...
                           ' principal component calculated at map',...
                           ' initialisation'],...
                           -axeS(1), -axeS(2), -axeS(3), size(map.PCs, 2));
                    end
                    if map.preproc
                        nodeMap(:, ind) = tmp(:, axeS(ind));
                    else
                        % Calculate projection on required PC
                        nodeMap(:, ind) = bsxfun(@minus, tmp, map.means)...
                            * map.PCs(:, -axeS(ind));
                    end
                end
            end
        end
        if dim > 2
            % Now draw surface
            trisurf(gridReady, nodeMap(:, 1), nodeMap(:, 2),...
                nodeMap(:, 3), f, 'FaceColor', 'interp',...
                'EdgeColor', 'none');
        elseif dim == 2
            % Now draw surface in 2D graph
            trisurf(gridReady, nodeMap(:, 1), nodeMap(:, 2),...
                zeros(size(nodeMap, 1), 1), f, 'FaceColor', 'interp',...
                'EdgeColor', 'none');
        end
        hold on;
        colormap(gca, colMap);
    end
    
    %get map links
    links = map.getLinks;
    if dim > 3 
        % Dimension is greater than 3 and we need to use three PCs or user
        % specified axes
        if isempty(axeS)
            % Project data onto PCs if it is necessary
            if ~map.preproc
                data = bsxfun(@minus, data, map.means);
                V = map.PCs(:, 1:3);
                data = data * V;
                maps = bsxfun(@minus, maps, map.means) * V;
            end
        else
            % Create copy of node coordinates and data coordinates
            tmp = maps;
            dataTmp = data;
            % Create new arrays for node coordinates and data
            maps = zeros(size(maps, 1), 3);
            data = zeros(size(data, 1), 3);
            % Check the necessity of original coordinates
            ind = axeS > 0;
            if sum(ind) > 0
                % At least one coordinate will be used. We need to
                % restore it
                temp = map.deprocessData(tmp);
                if any(axeS > size(temp, 2))
                    error(['Number of coordinate cannot be greater',...
                        ' than dimension of original space']);
                end
                maps(:, ind) = temp(:, axeS(ind));
                temp = map.deprocessData(dataTmp);
                data(:, ind) = temp(:, axeS(ind));
            end
            ind = ~ind;
            % Check the necessity of PCs
            if sum(ind) > 0
                % At least on PCs will be used
                if any(axeS < -size(map.PCs, 2))
                    error(['Number of principal component (MINUS',...
                        ' value of "axes") %d %d %d to use',...
                        ' must be positive and cannot be greater',...
                        ' than %d which is the number of',...
                        ' principal component calculated at map',...
                        ' initialisation'],...
                        -axeS(1), -axeS(2), -axeS(3), size(map.PCs, 2));
                end
                if map.preproc
                    maps(:, ind) = tmp(:, axeS(ind));
                    data(:, ind) = dataTmp(:, axeS(ind));
                else
                    % Calculate projection on required PC
                    maps(:, ind) = bsxfun(@minus, tmp, map.means)...
                        * map.PCs(:, -axeS(ind));
                    data(:, ind) = bsxfun(@minus, dataTmp, map.means)...
                        * map.PCs(:, -axeS(ind));
                end
            end
        end
        %Draw data
        for k = 1:nCls
            ind = classes == cls(k);
            plot3(data(ind, 1), data(ind, 2), data(ind, 3),...
                [markColour(k), markShape(k)],...
                'MarkerFaceColor', markColour(k),...
                'MarkerSize', markSize(k));
            hold on
        end
        if mapDraw
            %Draw edges
            %Prepare arrays

            X=[maps(links(:,1), 1)';maps(links(:,2), 1)'];
            Y=[maps(links(:,1), 2)';maps(links(:,2), 2)'];
            Z=[maps(links(:,1), 3)';maps(links(:,2), 3)'];
            %Draw edges
            plot3(X, Y, Z, nodeColour, 'LineWidth', lineWidth);
            %Draw map nodes
            plot3(maps(:, 1), maps(:, 2), maps(:, 3), 'Marker', nodeMarker,...
                'MarkerFaceColor', nodeColour, 'MarkerEdgeColor', nodeColour,...
                'MarkerSize', nodeMarkerSize, 'LineStyle', 'none');
        end
        if isempty(axeS)
            xlabel('PC1');
            ylabel('PC2');
            zlabel('PC3');
        else
            if axeS(1) > 0 
                xlabel(['Attr', num2str(axeS(1))]);
            else
                xlabel(['PC', num2str(-axeS(1))]);
            end
            if axeS(2) > 0 
                ylabel(['Attr', num2str(axeS(2))]);
            else
                ylabel(['PC', num2str(-axeS(2))]);
            end
            if axeS(3) > 0 
                zlabel(['Attr', num2str(axeS(3))]);
            else
                zlabel(['PC', num2str(-axeS(3))]);
            end
        end
        axis equal;
    elseif dim == 3
        %3d data
        %Draw data
        for k = 1:nCls
            ind = classes == cls(k);
            plot3(data(ind, 1), data(ind, 2), data(ind, 3),...
                [markColour(k), markShape(k)],...
                'MarkerFaceColor', markColour(k),...
                'MarkerSize', markSize(k));
            hold on
        end
        if mapDraw
            %Draw edges
            %Prepare arrays
            size(maps)
            size(links)
            X=[maps(links(:,1), 1)';maps(links(:,2), 1)'];
            Y=[maps(links(:,1), 2)';maps(links(:,2), 2)'];
            Z=[maps(links(:,1), 3)';maps(links(:,2), 3)'];
            %Draw edges
            plot3(X, Y, Z, nodeColour, 'LineWidth', lineWidth);
            %Draw map nodes
            plot3(maps(:, 1), maps(:, 2), maps(:, 3), 'Marker', nodeMarker,...
                'MarkerFaceColor', nodeColour, 'MarkerEdgeColor', nodeColour,...
                'MarkerSize', nodeMarkerSize, 'LineStyle', 'none');
        end
        axis equal;
    elseif dim == 2
        %two dimensional data
        %Draw data
        for k = 1:nCls
            ind = classes == cls(k);
            plot(data(ind, 1), data(ind, 2),...
                [markColour(k), markShape(k)],...
                'MarkerFaceColor', markColour(k),...
                'MarkerSize', markSize(k));
            hold on
        end
        if mapDraw
            %Draw edges
            %Prepare arrays
            X=[maps(links(:,1),1)';maps(links(:,2),1)'];
            Y=[maps(links(:,1),2)';maps(links(:,2),2)'];
            %Draw edges
            plot(X, Y, nodeColour, 'LineWidth', lineWidth);
            %Draw map nodes
            plot(maps(:,1),maps(:,2), 'Marker', nodeMarker,...
                'MarkerFaceColor', nodeColour, 'MarkerEdgeColor', nodeColour,...
                'MarkerSize', nodeMarkerSize, 'LineStyle', 'none');
        end
        axis equal;
    else
        %one dimensional data
        %Draw data
        for k = 1:nCls
            ind = classes == cls(k);
            plot(data(ind, 1), 0,...
                [markColour(k), markShape(k)],...
                'MarkerFaceColor', markColour(k),...
                'MarkerSize', markSize(k));
            hold on
        end
        if mapDraw
            
            %Draw edges
            %Prepare arrays
            X=[maps(links(:,1),1)'; maps(links(:,2),1)'];
            Y=zeros(2,size(links,1));
            %Draw edges
            plot(X, Y, nodeColour, 'LineWidth', lineWidth);
            %Draw map nodes
            plot(maps(:,1),0, 'Marker', nodeMarker,...
                'MarkerFaceColor', nodeColour, 'MarkerEdgeColor', nodeColour,...
                'MarkerSize', nodeMarkerSize, 'LineStyle', 'none');
        end
    end
end

function res = interpol(X, y, nodes, r)
%interpol calculates value of function y defined in data points X in each
%node from nodes.
%
%Inputs:
%   X is n-by-d data matrix with one data point in each row.
%   y is n-by-1 vector with values of function. y(i) contains function
%       value for point X(i, :).
%   nodes is m-by-d matrix of nodes to calculate function values.
%   r is smoothing parameter for function calculation.
%
%   memSize is constant for quick 

    % Calculate variances for all attributes
    smooth = -1 / (r * mean(var(X)));
    % Calculate distances from each node to each data point
    dist = bsxfun(@plus, sum(X.^2,2), sum(nodes.^2, 2)') - 2 * (X * nodes');
    % Calclulate RBF in each point
    tmp = bsxfun(@times, exp(dist * smooth), y);
    % Calculate result
    res = sum(tmp)';
    % Normalise result
    mins = min(res);
    res = (res - mins) / (max(res) - mins);
end

function [grid, maps, inter] = formGrid(grid, maps, inter)
    % Step 1. Create list of all edges in grid
    % Unify description of triangles: sort nodes in ascend order
    grid = sort(grid, 2);
    % Form list of all edges in grid
    edges = [grid(:, 1:2); grid(:, 2:3); grid(:, [1, 3])];
    % Search unique values
    edges = unique(edges, 'rows');
    % Step 2. Form list of nodes in new node list
    nN = size(maps, 1);
    nE = size(edges, 1);
    maps = [maps; (maps(edges(:, 1), :) + maps(edges(:, 2), :)) / 2];
    inter = [inter; (inter(edges(:, 1), :) + inter(edges(:, 2), :)) / 2];
    % Form list of indexes for nodes
    ind = zeros(nE + nN);
    siz = size(ind);
    indL = sub2ind(siz, edges(:, 1), edges(:, 2));
    ind(indL) = nN + 1:nN + nE;
    % Step 3. Form list of six nodes for each face
    face = [grid, ind(sub2ind(siz, grid(:, 1), grid(:, 2))),...
        ind(sub2ind(siz, grid(:, 2), grid(:, 3))),...
        ind(sub2ind(siz, grid(:, 1), grid(:, 3)))];
    % Step 4. Form final list of triangles
    grid = [face(:, [2, 4, 5]); face(:, [1, 4, 6]);...
        face(:, [3, 5, 6]); face(:, [4, 5, 6])];
end