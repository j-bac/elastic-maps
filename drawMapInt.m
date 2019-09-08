function drawMapInt( map, data, projType, varargin )
%drawMapInt create figure (optional) colouring map (optional) and draw data
%and maps in the internal maps coordinates. 
%
%IMPORTANT! It is necessary to stress that colouring of map is possible for
%the 2D maps only. Including of colouring by default removes drawing of
%nodes and sets map edges line width to 0.5 if other options are not
%specified by user.  
%
%IMPORTANT! 3D map colouring automatically excludes data drawing, removes
%nodes and sets map edges line width to 0.5 if other options are not
%specified by user.  
%
%Usage
%   drawMapInt(map, data);
%   drawMapInt(map, data, projType);
%   drawMapInt(__, Name, Value);
%
%Inputs:
%   map is object of class MapGeometry (descendant of this class).
%   data is n-by-d matrix of data points. Each row is one point.
%   projType is the type of projection: 
%       0 is projection to the nearest node,
%       1 is projection to the nearest edge,
%       2 is projection to the nearest face.
%       default value is 0.
%   There are several additional parameters in form 'Name', value:
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
%           'none'. Default value is 'o'.
%       'nodeMarkerSize' is positive number which is size of node marker.
%           This parameter is ignored for the projection of the first type.
%           Default value is 6;
%       'nodeColour' is one of the possible Matlab colours ('r', 'g', 'b',
%           'y', 'm', 'c', 'w', 'k'). Default value is 'r';
%       'lineWidth' is non-negative number for map line width. Zero means
%           absence of map at all (width of line is zero and 'nodeMarker'
%           is 'none'. Zero value is reasonable for 2D maps with nonzero
%           project type. Default value is 2. 
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
%               n-by-1 vector with values to use. Coordinates of vector X
%               are defined in original space even for map for preprocessed
%               data. 
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
%       'flatColoring' is logical attribute. True value (default) assumes
%           2D graph. Value false assumes surface with intensity
%           corresponding to height of surface. Value false also prevent
%           drawing of data points, removed nodes markers and set
%           'lineWidth' to 0.5 if linewidth is not directly specified by
%           user. 
%       'ColorMap' is colormap to use for current graph. This parameter has
%           meaning for maps with colouring only. Colormap is standard
%           Matlab colormap. It can be one of predefined colormaps or
%           user defined one. For more details find Colormap in the Matlab
%           documentation. 
%       'ColoringInSpace' is logical attribute. True value means
%           calculation of density or another data defined colouring
%           without projection onto map. False (default) value means
%           firstly projection of data onto map and then calculation of
%           colouring.
%       'Smooth' is positive real number. It is parameter of smoothing of
%           function interpolation for densities and other functions
%           defined in the data points. Default value is 0.15. Data
%           interpolation is implemented in radial basis function like
%           fashion: for query point x we calculate squared distance to each
%           data point d(i). Value of function in the point x is sum of
%           exp(-d(i)/(smooth*s2)) where s2 is mean of variances of all
%           attributes.
%       'PseudoTimePlot' draw graph of specified function for 1D map. This
%           argument is similar to 'coloring' argument and can be
%           [] (empty) means no graphs. It is default value.
%           'density' is density graph. It is equivalent to vector of
%               ones.
%           fun is function handle of form function res = fun(X), where X
%               is n-by-d matrix with one data point per row and res is
%               n-by-1 vector with values to use. Coordinates of vector X
%               are defined in original space even for map for preprocessed
%               data. 
%           k is positive integer number. In this case k is number of
%               coordinate to use.
%           k is negative integer number. In this case graph presents 
%               value of projection on the k'th principal component.
%           vect is n-by-1 vector with data defined function. Each element
%               of vector corresponds to data point in matrix data.
%           matr is N-by-(d+1) matrix with one data point in the first d
%               columns of each row and value of function of this point in
%               the last column. For example to calculate density function
%               it is necessary to send data matrix with 1 in (d+1) column. 
%

    % Check presence of projType
    if ischar(projType)
        varargin = [{projType}, varargin];
        projType = 0;
    end
    
    % Data preprocessing
    data = map.preprocessData(data);
    % Get number of data points
    N = size(data, 1);

    % Default values of optional parameters
    classes = [];
    markColour = [];
    markShape = [];
    markSize = [];
    nodeMarker = 'o';
    nodeMarkerSize = 6;
    nodeColour = 'r';
    lineWidth = 2;
    newFigure = true;
    source = [];
    flatColoring = true;
    colMap = parula;
    ColouringInSpace = false;
    nodeMarkerSpecified = false;
    lineWidthSpecified = false;
    smooth = 0.15;
    graph = [];
    % Parse optional parameters
    for i=1:2:length(varargin)
        switch lower(varargin{i})
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
            case 'flatcoloring'
                flatColoring = varargin{i + 1};
            case {'colouringinspace', 'coloringinspace'}
                ColouringInSpace = varargin{i + 1};
            case 'smooth'
                smooth = varargin{i + 1};
            case 'pseudotimeplot'
                graph = varargin{i + 1};
            otherwise
                error(['Unknown argument at position "', num2str(i + 2)]);
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
    
    % Get map coordinates
    intern = map.getInternalCoordinates;
    
    % Project data onto map
    dataP = [];
    if ~isempty(data)
        dataP = map.project(data, projType, 'internal');
    end
    
    % Just in case value of function for colouring or graph is 
    f = 0;
    
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
        
        % Get data space dimension and mapped coordinates
        mapped = map.getMappedCoordinates;
        d = size(mapped, 2);
        
        % Extract grid
        [gridReady, nodeMap, nodeInt] = formGrid(map.getFaces, mapped, intern);
        
        % Check the correctness of source parameter
        % Identify type of source
        tmp = false;
        if isnumeric(source) 
            if isscalar(source)
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
                    f = interpol(dataP, source, nodeInt, smooth);
                end
            elseif ismatrix(source)
                temp = map.preprocessData(source(:, 1:d));
                if size(temp, 2) ~= d
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
                f = interpol(dataP, ones(N, 1), nodeInt, smooth);
            end
        elseif isa(source, 'function_handle')
            % Get values to draw
            f = source(map.deprocess(nodeMap));
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
                'matr is N-by-(d+1) matrix with one data point in the nfirst d\n',...
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
    end
    
    % Check graph parameters
    if ~isempty(graph)
        
        if map.getDimension() ~= 1
            error('Pseudo time plot can be used for 1D maps only');
        end
        
        if ~isnumeric(smooth) || smooth <= 0
            error('Smooth must be positive number.');
        end
        
        % Get data space dimension and mapped coordinates
        mapped = map.getMappedCoordinates;
        d = size(mapped, 2);
        
        % Extract detailed grid for graph
        [nodeMap, nodeInt] = formGraphGrid(mapped, intern);
        
        % Check the correctness of graph parameter
        % Identify type of source
        tmp = false;
        if isnumeric(graph) 
            if isscalar(graph)
                % It is number of coordinate or PCs
                graph = round(graph);
                if graph > 0
                   % It is coordinate
                   % Get data in original space
                   temp = map.deprocessData(nodeMap);
                   if graph > size(temp, 2)
                       error(['Number of coordinate (value of "PseudoTimePlot")',...
                           ' %d to draw must be positive and cannot be',...
                           ' greater than data space dimension %d'],...
                           graph, size(temp, 2));
                   end
                   f = temp(:, graph);
                elseif graph == 0
                    error(['Number of coordinate (value of "PseudoTimePlot")',...
                        ' to draw must be positive and number of',...
                        ' principal component must be negative.',...
                        ' Value of "PseudoTimePlot" cannot be zero']);
                else
                    % -graph is number of PC
                    if -graph > size(map.PCs, 2)
                       error(['Number of principal component (MINUS',...
                           ' value of "PseudoTimePlot") %d to draw',...
                           ' must be positive and cannot be greater',...
                           ' than %d which is the number of',...
                           ' principal component calculated at map',...
                           ' initialisation'],...
                           -graph, size(map.PCs, 2));
                    end
                    if map.preproc
                        % Get required coordinate
                        f = nodeMap(:, -graph);
                    else
                        % Calculate projection on required PC
                        f = bsxfun(@minus, nodeMap, map.means)...
                            * map.PCs(:, -graph);
                    end
                end
            elseif isvector(graph)
                source = source(:);
                if size(source, 1) ~= size(data, 1)
                    error(['Number of elements in vector "graph"',...
                        ' %d must coincides with number of data points %d'],...
                        size(source, 1), size(data, 1));
                end
                
                if ColouringInSpace
                    f = interpolGraph(data, source, nodeMap, smooth);
                else
                    f = interpolGraph(dataP, source, nodeInt, smooth);
                end
            elseif ismatrix(graph)
                temp = map.preprocessData(source(:, 1:d));
                if size(temp, 2) ~= d
                    error(['Wrong dimension of matrix in graph argument\n',...
                        'Matrix of data points with function to draw must be\n',...
                        'n-by-(d+1) matrix with one data point in the\n',...
                        'first d columns of each row and value of function\n',...
                        'of this point in the last column. For example\n',...
                        'to calculate density function it is necessary\n',...
                        'to send data matrix with 1 in (d+1) column.\n%s'],'');
                end
                % Calculate influence of data points in nodes
                if ColouringInSpace
                    f = interpolGraph(temp, source(:, end),...
                        nodeMap, smooth);
                else
                    f = interpolGraph(...
                        map.project(temp, projType, 'internal'),...
                        source(:, end), nodeInt, smooth);
                end
                clear temp;
            else
                tmp = true;
            end
        elseif ischar(graph) && strcmp(graph, 'density')
            if ColouringInSpace
                f = interpolGraph(data, ones(N, 1), nodeMap, smooth);
            else
                f = interpolGraph(dataP, ones(N, 1), nodeInt, smooth);
            end
        elseif isa(source, 'function_handle')
            % Get values to draw
            f = source(nodeMap);
        else
            tmp = true;
        end
    
        % Throw error is necessary
        if tmp
            error(['Wrong type of PseudoTimePlot argument. PseudoTimePlot',...
                ' must be either\n',...
                '[] (empty) means no PseudoTimePlot. It is default value.\n',...
                '''density'' is density PseudoTimePlot. It is equivalent to vector of\n',...
                '     ones.\n',...
                'fun is function handle of form function res = fun(X), where X\n',...
                '     is n-by-d matrix with one data point per row and res is\n',...
                '     n-by-1 vector with values to use. Coordinates of vector X\n',...
                '     are defined in original space even for map for preprocessed\n',...
                '     data.\n',...
                'k is positive integer number. In this case k is number of\n',...
                '     coordinate to use.\n',...
                'k is negative integer number. In this case graph is drawn by\n',...
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
        
        %Including of graph by default removes drawing of nodes and
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
    end
    
    % Remove data since it is not necessary
    data = dataP;
    
    % Create figure
    if newFigure
        figure;
    end
    % Get limits of map
    lims = zeros(1, 4);
    lims([1, 3]) = min(intern) - 1;
    lims([2, 4]) = max(intern) + 1;
    % Get map links
    links = map.getLinks;
    
    % Get map dimension
    dim = size(intern,2);
    if dim > 3
        
    elseif dim == 3
        % 3d map
        
    elseif dim == 2
        % Two dimensional map
        hold off;
        %Form edges
        %Prepare arrays
        X=[intern(links(:,1),1)';intern(links(:,2),1)'];
        Y=[intern(links(:,1),2)';intern(links(:,2),2)'];
        
        % Draw colouring
        if ~isempty(source)
            if flatColoring
                % Create 2D graph and fix it.
                plot(intern(:, 1), intern(:, 2), '.');
                hold on;
                axis square;
                minn = min(intern);
                maxx = max(intern);
                axis([minn(1), maxx(1), minn(2), maxx(2)]);
                % Now draw surface
                trisurf(gridReady, nodeInt(:, 1), nodeInt(:, 2),...
                    zeros(size(nodeInt, 1), 1), f, 'FaceColor', 'interp',...
                    'EdgeColor', 'none');
                %Draw edges
                if mapDraw
                    plot(X, Y, nodeColour, 'LineWidth', lineWidth);
                end
            else
                % Draw surface
                trisurf(gridReady, nodeInt(:, 1), nodeInt(:, 2),...
                    f, 'FaceColor', 'interp', 'EdgeColor', 'none');
                hold on;
                %Draw edges
                if mapDraw
                    plot3(X, Y,...
                        [f(links(:,1),1)';f(links(:,2),1)'],...
                        nodeColour, 'LineWidth', lineWidth);
                end
                data = [];
            end
            colormap(gca, colMap);
        else
            %Draw edges
            if mapDraw
                plot(X, Y, nodeColour, 'LineWidth', lineWidth);
            end
        end
        hold on;
        %Project data to map
        if projType == 0 && mapDraw
            if ~isempty(data)
                %Search unique points
                [dat, ~, ic] = unique(data,'rows');
                count = accumarray(ic, 1);
                ma = max(count);
                %Draw map nodes
                if ~strcmpi(nodeMarker, 'none')
                    if nCls == 1 % No classes
                        scatter(dat(:, 1), dat(:, 2), count/ma * 400,...
                            nodeColour, 'filled', nodeMarker);
                    else
                        % We have classes We need to draw each marker by it's
                        % own colour or several colours
                        % First of all calculate number of points of each class
                        % in each node.
                        nDat = size(dat, 1);
                        props = zeros(nDat, nCls);
                        for k = 1:nCls
                            ind = classes == cls(k);
                            props(:, k) = accumarray(ic(ind), 1, [nDat, 1]);
                        end
                        % Now we put pie Chart in each necessary node.
                        for k = 1:nDat
                            drawPieChart(dat(k, 1), dat(k, 2),...
                                0.5 * count(k) / ma, props(k, :), markColour);
                        end
                    end
                end
            end
            axis(lims);
        else
            if mapDraw && ~strcmpi(nodeMarker, 'none')
                %Draw maps nodes
                plot(intern(:,1), intern(:,2), 'Marker', nodeMarker,...
                    'MarkerFaceColor', nodeColour, 'MarkerEdgeColor',...
                    nodeColour, 'MarkerSize', nodeMarkerSize,...
                    'LineStyle', 'none');
            end

            %Draw data points
            if ~isempty(data)
                for k = 1:nCls
                    ind = classes == cls(k);
                    plot(data(ind, 1), data(ind, 2),...
                        [markColour(k), markShape(k)],...
                        'MarkerFaceColor', markColour(k),...
                        'MarkerSize', markSize(k));
                    hold on
                end
            end
        end
    else
        %one dimensional data
        hold on
        %Draw edges
        %Prepare arrays
        X=[intern(links(:,1),1)';intern(links(:,2),1)'];
        Y=zeros(2,size(links,1));
        %Draw edges
        plot(X, Y, nodeColour, 'LineWidth', lineWidth);
        if projType == 0
            if ~isempty(data)
                %Search unique points
                [dat, ~, ic] = unique(data,'rows');
                count = accumarray(ic,1);
                ma=max(count);
                %Draw map nodes
                if nCls == 1 % No classes
                    scatter(dat(:,1), zeros(size(dat,1),1),...
                        count / ma * 400, nodeColour, 'filled', nodeMarker);
                else
                    % We have classes We need to draw each marker by it's
                    % own colour or several colours
                    % First of all calculate number of points of each class
                    % in each node.
                    nDat = size(dat, 1);
                    props = zeros(nDat, nCls);
                    for k = 1:nCls
                        ind = classes == cls(k);
                        props(:, k) = accumarray(ic(ind), 1, [nDat, 1]);
                    end
                    % Now we put pie Chart in each necesssary node.
                    for k = 1:nDat
                        drawPieChart(dat(k, 1), 0,...
                            0.5 * count(k) / ma, props(k, :), markColour);
                    end
                    lims(3:4) = lims(1:2) - sum(lims(1:2)) / 2;
                    axis(lims);
                end
                % Renormalise function to use
                mins = min(f);
                maxs = max(f);
                if mins < 0 
                    mins = lims(3) / mins;
                else
                    mins = 0;
                end
                if maxs > 0 
                    maxs = lims(4) / maxs;
                else
                    maxs = 0;
                end
                if mins > maxs || mins == 0
                    mins = maxs;
                end
                f = f * mins;
            end
        else
            %Draw maps nodes
            plot(intern(:,1), 0, 'Marker', nodeMarker,...
                'MarkerFaceColor', nodeColour, 'MarkerEdgeColor',...
                nodeColour, 'MarkerSize', nodeMarkerSize,...
                'LineStyle', 'none');
            %Draw data points
            if ~isempty(data)
                for k = 1:nCls
                    ind = classes == cls(k);
                    plot(data(ind, 1), 0,...
                        [markColour(k), markShape(k)],...
                        'MarkerFaceColor', markColour(k),...
                        'MarkerSize', markSize(k));
                end
            end
        end
    end
    % Draw graph if required
    if ~isempty(graph)
        plot(nodeInt, f, 'k-');
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
    % Calculate RBF in each point
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

function [nodeMap, nodeInt] = formGraphGrid(mapped, intern)
    % This is standard 1D map - polyline with standard numeration from 1 to
    % N where N is the last node.
    % We use divider nDiv for map
    nDiv = 10;
    step = 1 / nDiv;
    arr = 0:step:1;
    arr = arr(1:nDiv);
    % Form internal
    nodeInt = bsxfun(@plus, intern, arr)';
    nodeInt = nodeInt(:);
    nodeInt = nodeInt(1:end - nDiv + 1);
    % Form external
    [n, m] = size(mapped);
    nodeMap = zeros((n - 1) * nDiv + 1, m);
    % Fill extern
    arr = arr(2:end)';
    for k = 1:n - 1
        nodeMap((k - 1) * nDiv + 1, :) = mapped(k, :);
        nodeMap((k - 1) * nDiv + 2:k * nDiv, :) = ...
            bsxfun(@plus, mapped(k, :),...
            bsxfun(@times, mapped(k + 1, :) - mapped(k, :), arr));
    end
    nodeMap(end, :) = mapped(end, :);
end

function res = interpolGraph(X, y, nodes, r)
%interpolGraph calculates value of function y defined in data points X in
%each node from nodes.
%
%Inputs:
%   X is n-by-d data matrix with one data point in each row.
%   y is n-by-1 vector with values of function. y(i) contains function
%       value for point X(i, :).
%   nodes is m-by-d matrix of nodes to calculate function values.
%   r is smoothing parameter for function calculation.
%

    % Calculate distances from each node to each data point
    dist = bsxfun(@plus, sum(X.^2,2), sum(nodes.^2, 2)') - 2 * (X * nodes');
    % Calculate RBF in each point
    tmp = bsxfun(@times, exp(- dist * r), y);
    % Calculate result
    res = sum(tmp)';
    % Normalise result
    mins = min(res);
    res = (res - mins) / (max(res) - mins);
end
