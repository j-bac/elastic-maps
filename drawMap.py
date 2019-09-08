from matlab_like_funcs import accumarray, isnumeric, histc
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objs as go
from plotly.offline import init_notebook_mode, iplot
init_notebook_mode(connected=True)



def drawMap(_map, 
            data,
            classes = [], 
            markColour = [], 
            markShape = [], 
            markSize = [], 
            nodeMarker = 'none', 
            nodeMarkerSize = 2, 
            nodeColour = 'r', 
            lineWidth = 1, 
            newFigure = True, 
            source = [], 
            colMap = 'Viridis', 
            ColouringInSpace = False, 
            nodeMarkerSpecified = False, 
            lineWidthSpecified = False, 
            smooth = 0.15, 
            projType = 0, 
            axeS = [],
           ):

    #drawMap create figure and draw data and maps. If data dimension (number of
    #columns) is less than 3 then original data coordinate system is used. If
    #data dimension is greater than 3 then projection into space of the first
    #three PCs is used for drawing.
    #
    #Usage
    #   drawMap(map, data) 
    #   drawMap(__, Name, Value) 
    #
    #Inputs:
    #   map is object of class MapGeometry (descendant of this class).
    #   data is n-by-dim matrix of data points. Each row is one point.
    #   There are several possible 'Name', value pairs:
    #       'axes' is array of three different nonzero integer. This argument
    #           can be used for space with more than three coordinates only.
    #           The first number is vector for x axis, the second number is
    #           vector for y axis and the third number is vector for z axis.
    #           Positive number means number of coordinate in original space.
    #           Negative number mean specified principal component.
    #       'classes' is n-by-1 vector of class labels for points. Each label
    #           is positive integer number which is the index of cell array
    #           with marker descriptions. Marker descriptions can be specified
    #           by user in the arguments markColour, markShape, markSize.
    #           Otherwise standard marker is used for all points.
    #       'markColour' is K-by-1 vector of standard Matlab colours ('r', 'g',
    #           'b', 'y', 'm', 'c', 'w', 'k'). Value K is number of defined
    #           markers. If markColour is omitted then 'b' (blue) is used for
    #           all markers. 
    #       'markShape' is K-by-1 vector of standard Matlab marker shapes ('o',
    #           '+', '*', '.', 'x', 's', 'd', '^', 'v', '<', '>', 'p', 'h').
    #           Value K is number of defined markers. If markShape is omitted
    #           then 's' (square) is used for all markers.
    #       'markSize' is K-by-1 vector of positive numbers. Value K is number
    #           of defined markers. If markSize is omitted then 6 is used for
    #           all markers.
    #       'nodeMarker' is one of the possible marker shape symbol ('o', '+',
    #           '*', '.', 'x', 's', 'd', '^', 'v', '<', '>', 'p', 'h') or
    #           'none'. Default value is 'none'.
    #       'nodeMarkerSize' is positive number which is size of node marker.
    #           Default value is 6 
    #       'nodeColour' is one of the possible Matlab colours ('r', 'g', 'b',
    #           'y', 'm', 'c', 'w', 'k'). Default value is 'r' 
    #       'lineWidth' is non-negative number for map line width. Zero means
    #           absence of map at all (width of line is zero and 'nodeMarker'
    #           is 'none'. Default value is 0.5.
    #       'newFigure' is logical argument. Value true (default) causes
    #           creation of new figure. Value false causes usage of current
    #           active figure. This option can be used, for example, for
    #           subplots.
    #       'coloring' defines the type of data to colour. It can have
    #               following values: 
    #           [] (empty) means no colouring. It is default value.
    #           'density' is density colouring. It is equivalent to vector of
    #               ones.
    #           fun is def handle of form def res = fun(X), where X
    #               is n-by-d matrix with one data point per row and res is
    #               n-by-1 vector with values to use. Coordinate of vector X
    #               are defined in PREPROCESSED space.
    #           k is positive integer number. In this case k is number of
    #               coordinate to use.
    #           k is negative integer number. In this case map is coloured by
    #               value of projection on the k'th principal component.
    #           vect is n-by-1 vector with data defined def. Each element
    #               of vector corresponds to data point in matrix data.
    #           matr is N-by-(d+1) matrix with one data point in the first d
    #               columns of each row and value of def of this point in
    #               the last column. For example to calculate density def
    #               it is necessary to send data matrix with 1 in (d+1) column. 
    #       'ColorMap' is colormap to use for current graph. This parameter has
    #           meaning for maps with colouring only. Colormap is standard
    #           Matlab colormap. It can be one of predefined colormaps or
    #           user defined. For more details find Colormap in the Matlab
    #           documentation. 
    #       'ColoringInSpace' is logical attribute. True value means
    #           calculation of density or another data defined colouring
    #           without projection onto map. False (default) value means
    #           firstly projection of data onto map and then calculation of
    #           colouring.
    #       'Smooth' is positive real number. It is parameter of smoothing of
    #           def interpolation for densities and other defs
    #           defined in the data points. Default value is 2. Data
    #           interpolation is implemented in radial basis def like
    #           fasion: for query point x we calculate squared distance to each
    #           data point d(i). Value of def in the point x is sum of
    #           exp(-d(i)/(smooth*s2)) where s2 is mean of variances of all
    #           attributes.
    #       'projType' is type of projection onto map. Projection is necessary
    #           to create colouring from the def defined in the data
    #           points including densities. projType can have following values: 
    #               0 is projection to the nearest node,
    #               1 is projection to the nearest edge,
    #               2 is projection to the nearest face.
    #           default value is 0.

    # Data preprocessing
    data = _map.preprocessData(data) 

    #Get data dimension
    [N, dim] = data.shape 
    
    # Sanity check of arguments.
    if not classes.size > 0:
        classes = np.ones((N, 1)) 
        _cls = 1 
        nCls = 1
    else:
        _cls = np.unique(classes) 
        nCls = len(_cls)
        
    
    if markColour == []:
        markColour = np.tile('b', (nCls, 1)) 

    if markShape == []:
        markShape = np.tile('s', (nCls, 1)) 
    
    if markSize == []:
        markSize = [nodeMarkerSize]*nCls
    
    if lineWidth == 0:
        mapDraw = 0 
    else:
        mapDraw = 1 
    
    if not(axeS == []):
        if len(axeS) != 3 or not(isnumeric(axeS)) or any(axeS == 0) or len(np.unique(axeS)) != 3:
            raise ValueError(['"axes" must be array of three different nonzero'+
                ' integers. This argument can be used for space'+
                ' with more than three coordinates only. The'+
                ' first number is vector for x axis, the second'+
                ' number is vector for y axis and the third'+
                ' number is vector for z axis. Positive number'+
                ' means number of coordinate in original space.'+
                ' Negative number mean specified principal'+
                ' component.']) 

                
    #Create figure
    if newFigure:
        layout=go.Layout(
            height=800,
            width=800, 
            title="Elastic Map & Data",
            scene={"aspectmode": "cube",
                   "xaxis": {"title": "PC1", "showbackground":True},
                   "yaxis": {"title": "PC2", "showbackground":True},
                   "zaxis": {"title": "PC3", "showbackground":True}
                  }
        )
        
        fig = go.Figure(layout=layout)
        
    #Get map coordinates
    maps = _map.getMappedCoordinates() 

    # Check colouring parameters
    if not source == []:
        # Is map appropriate for colouring
        if not hasattr(_map, 'getFaces'):
            raise ValueError('Map must implement method "getFaces" for colouring') 

        if _map.getDimension() != 2:
            raise ValueError('Map coloring can be used for 2D maps only') 

        
        if not(isnumeric(smooth)) or smooth <= 0:
            raise ValueError('Smooth must be a positive number') 
        
        # Get internal map nodes' coordinates
        intern = _map.getInternalCoordinates() 
        
        # Extract grid
        [gridReady, nodeMap, nodeInt] = formGrid(_map.getFaces(), maps, intern) 
        
        # Check the correctness of source parameter
        # Identify type of source
        tmp = False 
        if isnumeric(source): 
            if np.isscalar(source):
                # It is number of coordinate or PCs
                source = np.round(source) 
                
                if source > 0:
                   # It is coordinate
                   # Get data in original space
                    temp = _map.deprocessData(nodeMap) 
                    if source > temp.shape[1]:
                        raise ValueError(['Number of coordinate (value of "coloring")'+
                       ' #d to draw must be positive and cannot be'+
                       ' greater than data space dimension #d']+
                       source, temp.shape[1]) 
                    f = temp[:, source]
                    
                elif source == 0:
                    raise ValueError(['Number of coordinate (value of "coloring")'+
                        ' to draw must be positive and number of'+
                        ' principal component must be negative.'+
                        ' Value of "coloring" cannot be zero']) 
                else:
                    if -source > _map.PCs.shape[1]:
                        raise ValueError(['Number of principal component (MINUS'+
                           ' value of "source") #d to draw'+
                           ' must be positive and cannot be greater'+
                           ' than #d which is the number of'+
                           ' principal component calculated at map'+
                           ' initialisation']+
                           -source, _map.PCs.shape[1])
                        
                    elif _map.preproc:
                        # Get required coordinate
                        f = nodeMap[:, -source] 
                    else:
                        # Calculate projection on required PC
                        f = (nodeMap- _map.means) @ _map.PCs[:, -source] 

            elif len(source.shape) == 1 or (1 in source.shape):
                if any(source.shape) == 1:
                    source = source.reshape(-1,1)
                else: 
                    source = source[:,None]
                    
                if source.shape[0] != N:
                    raise ValueError(['Number of elements in vector "coloring"'+
                        ' #d must coincides with number of data points #d']+
                        source.shape[0], N) 
                
                if ColouringInSpace:
                    f = interpol(data, source, nodeMap, smooth) 
                else:
                    f = interpol(_map.project(data, projType, 'internal'),
                        source, nodeInt, smooth) 
                    
            elif len(source.shape) == 2:
                temp = _map.preprocessData(source[:,:d]) 
                if temp.shape[1] != dim:
                    raise ValueError(['Wrong dimension of matrix in source argument\n'+
                        'Matrix of data points with def to draw must be\n'+
                        'n-by-(d+1) matrix with one data point in the\n'+
                        'first d columns of each row and value of def\n'+
                        'of this point in the last column. For example\n'+
                        'to calculate density def it is necessary\n'+
                        'to send data matrix with 1 in (d+1) column.\n#s']) 

                # Calculate influence of data points in nodes
                if ColouringInSpace:
                    f = interpol(temp, source[:,-1]+
                        nodeMap, smooth) 
                else:
                    f = interpol(_map.project(temp, projType, 'internal'),
                        source[:, end], nodeInt, smooth) 

            else:
                tmp = true 

        elif isinstance(source,str) and source.lower() == 'density':
            if ColouringInSpace:
                f = interpol(data, np.ones(N, 1), nodeMap, smooth) 
            else:
                f = interpol(_map.project(data, projType, 'internal'),
                    np.ones(N, 1), nodeInt, smooth) 

        else:
            tmp = True 
        # Throw raise ValueError is necessary
        if tmp:
            raise ValueError(['Wrong type of source argument. Source must be'+
                ' either\n'+
                '[] (empty) means no colouring. It is default value.\n'+
                "'density' is density colouring. It is equivalent to vector of\n"+
                '     ones.\n'+
                'fun is def handle of form def res = fun(X), where X\n'+
                '     is n-by-d matrix with one data point per row and res is\n'+
                '     n-by-1 vector with values to use. Coordinates of vector X\n'+
                '     are defined in original space even for map for preprocessed\n'+
                '     data.\n'+
                'k is positive integer number. In this case k is number of\n'+
                '     coordinate to use.\n'+
                'k is negative integer number. In this case map is coloured by\n'+
                '     value of projection on the k''th principal component.\n'+
                'vect is n-by-1 vector with data defined def. Each element\n'+
                '     of vector corresponds to data point in matrix data.\n'+
                'matr is N-by-(d+1) matrix with one data point in the first d\n'+
                '     columns of each row and value of def of this point in\n'+
                '     the last column. For example to calculate density def\n'+
                '     it is necessary to send data matrix with 1 in (d+1) column.\n'+
                ' For example to calculate density def\n'+
                '     it is necessary to send data matrix with 1 in'+
                ' (d+1)     column.\n#s']) 

        
        #Including of colouring by default removes drawing of nodes and
        #sets map edges line width to 0.5 if other options are not
        #specified by user.   
        #
        #3D map colouring automatically excludes data drawing, removes
        #nodes and sets map edges line width to 0.5 if other options are
        #not specified by user.  
        if not nodeMarkerSpecified: 
            nodeMarker = 'none' 

        if not lineWidthSpecified:
            lineWidth = 0.5 

        if dim > 3: 
            # Dimension is greater than 3 and we need to use PCs or
            # specified coordinates
            if axeS == []:
                # By default Project data onto three PCs if it is necessary
                if _map.preproc:
                    V = _map.PCs[:, :3] 
                    nodeMap = (nodeMap - _map.means) * V 

            else:
                # Create copy of node coordinates
                tmp = nodeMap 
                # Create new array for node coordinates
                nodeMap = np.zeros((nodeMap.shape[0], 3)) 
                # Check the necessity of original coordinates
                ind = axeS > 0 
                if np.sum(ind,axis=0) > 0:
                    # At least one coordinate will be used. We need to
                    # restore it
                    temp = _map.deprocessData(tmp) 
                    if any(axeS > temp.shape[1]):
                        raise ValueError(['Number of coordinate cannot be greater'+
                            ' than dimension of original space']) 
                        
                    nodeMap[:, ind] = temp[:, axeS[ind]]

                ind =  not ind 
                # Check the necessity of PCs
                if np.sum(ind,axis=0) > 0:
                    # At least on PCs will be used
                    if any(axeS < -_map.PCs.shape[1]):
                           raise ValueError(['Number of principal component (MINUS'+
                               ' value of "axes") #d #d #d to use'+
                               ' must be positive and cannot be greater'+
                               ' than #d which is the number of'+
                               ' principal component calculated at map'+
                               ' initialisation']+
                               -axeS[0], -axeS[1], -axeS[2],_map.PCs.shape[1]) 

                    if _map.preproc:
                        nodeMap[:, ind] = tmp[:, axeS[ind]] 
                    else:
                        # Calculate projection on required PC
                        nodeMap[:, ind] =(tmp- _map.means)* _map.PCs[:, -axeS[ind]]

        if dim > 2:
            # Now draw surface
            trisurf(gridReady, nodeMap[:, 0], nodeMap[:, 1]+
                nodeMap[:, 2], f, 'FaceColor', 'interp'+
                'EdgeColor', 'none') 
        elif dim == 2:
            # Now draw surface in 2D graph
            trisurf(gridReady, nodeMap[:, 0], nodeMap[:, 1]+
                np.zeros((nodeMap.shape[0], 1)), f, 'FaceColor', 'interp'+
                'EdgeColor', 'none') 
    
    #get map links
    links = _map.getLinks() 
    if dim > 3: 
        # Dimension is greater than 3 and we need to use three PCs or user
        # specified axes
        if axeS == []:
            # Project data onto PCs if it is necessary
            if not _map.preproc:
                data = data - _map.means 
                V = _map.PCs[:,:3] 
                data = data * V 
                maps = (maps - _map.means) * V 
                
        else:
            # Create copy of node coordinates and data coordinates
            tmp = maps 
            dataTmp = data 
            # Create new arrays for node coordinates and data
            maps = np.zeros((maps.shape[0], 3)) 
            data = np.zeros((data.shape[0], 3)) 
            # Check the necessity of original coordinates
            ind = axeS > 0 
            if np.sum(ind,axis=0) > 0:
                # At least one coordinate will be used. We need to
                # restore it
                temp = _map.deprocessData(tmp) 
                if any(axeS > temp.shape[1]):
                    raise ValueError(['Number of coordinate cannot be greater'+
                        ' than dimension of original space']) 

                maps[:, ind] = temp[:, axeS[ind]] 
                temp = _map.deprocessData(dataTmp) 
                data[:, ind] = temp[:, axeS[ind]] 

            ind = ~ind 
            # Check the necessity of PCs
            if np.sum(ind,axis=0) > 0:
                # At least on PCs will be used
                if any(axeS < -_map.PCs.shape[1]):
                    raise ValueError(['Number of principal component (MINUS'+
                        ' value of "axes") #d #d #d to use'+
                        ' must be positive and cannot be greater'+
                        ' than #d which is the number of'+
                        ' principal component calculated at map'+
                        ' initialisation']+
                        -axeS[0], -axeS[1], -axeS[2], _map.PCs.shape[1]) 

                if _map.preproc:
                    maps[:, ind] = tmp[:,axeS[ind]]
                    data[:, ind] = dataTmp[:, axeS[ind]]
                else:
                    # Calculate projection on required PC
                    maps[:, ind] = (tmp- _map.means)* _map.PCs[:, -axeS[ind]] 
                    data[:, ind] = (dataTmp- _map.means)* _map.PCs[:, -axeS[ind]]

        #Draw data
        if nCls > 1:
            for k in range(nCls):
                ind = (classes == _cls[k]).squeeze()
               # plot3(data(ind, 1), data(ind, 2), data(ind, 3)+
                #    [markColour(k), markShape(k)]+
                #  'MarkerFaceColor', markColour(k)+
                 #   'MarkerSize', markSize(k))
                fig.add_trace(go.Scatter3d(x=data[ind, 0],
                                           y=data[ind, 1],
                                           z=data[ind, 2],
                                           mode='markers',
                                           marker=dict(
                                               size=markSize[k],
                                               color=markColour[k]
                                                 )
                               ))

        else:
            k = 0
            ind = (classes == _cls).squeeze()
            fig.add_trace(go.Scatter3d(x=data[ind, 0],
                                       y=data[ind, 1],
                                       z=data[ind, 2],
                                       mode='markers',
                                       marker=dict(
                                           size=markSize[k],
                                           color=markColour[k]
                                             )
                           ))
        if mapDraw:
            
            #Draw map nodes
            fig.add_trace(go.Scatter3d(x=maps[:, 0],
                                       y=maps[:, 1],
                                       z=maps[:, 2],
                                       mode='markers',
                                       marker=dict(
                                           size=1,
                                           color='coral'
                                       )
                           ))
            #Draw edges
            #Prepare arrays
            #pairs_XYZ = []
            #for j in range(3):
            #    a=np.array([maps[:,j][i] for i in links[:,0]])[None]
            #    b=np.array([maps[:,j][i] for i in links[:,1]])[None]
            #    pairs_XYZ.append(np.concatenate([a,b]).T)
            
            x_lines = list()
            y_lines = list()
            z_lines = list()

            #create the coordinate list for the lines
            for l in links:
                for i in range(2):
                    x_lines.append(maps[l[i],0])
                    y_lines.append(maps[l[i],1])
                    z_lines.append(maps[l[i],2])
                x_lines.append(None)
                y_lines.append(None)
                z_lines.append(None)

            fig.add_trace(go.Scatter3d(
                x=x_lines,
                y=y_lines,
                z=z_lines,
                mode='lines',
                name='map edges',
                marker=dict(color='coral')
            ))
            


        if axeS == []:
            pass
            
            
        
        else:
            if axeS(1) > 0: 
                xlabel(['Attr', num2str(axeS(1))]) 
            else:
                xlabel(['PC', num2str(-axeS(1))]) 
            
            if axeS(2) > 0: 
                ylabel(['Attr', num2str(axeS(2))]) 
            else:
                ylabel(['PC', num2str(-axeS(2))]) 
            
            if axeS(3) > 0: 
                zlabel(['Attr', num2str(axeS(3))]) 
            else:
                zlabel(['PC', num2str(-axeS(3))]) 
        
        fig.show()


    elif dim == 3:
        #3d data
        #Draw data
        for k in range(nCls):
            ind = classes == cls(k) 
            plot3(data(ind, 1), data(ind, 2), data(ind, 3)+
                [markColour(k), markShape(k)]+
                'MarkerFaceColor', markColour(k)+
                'MarkerSize', markSize(k)) 

        if mapDraw:
            #Draw edges
            #Prepare arrays
            X=np.concatenate([maps[links[: 0], 0].T, maps[links[: 1], 0].T])
            Y=np.concatenate([maps[links[: 0], 1].T, maps[links[: 1], 1].T]) 
            Z=np.concatenate([maps[links[: 0], 2].T, maps[links[: 1], 2].T]) 

            #Draw edges
            plot3(X, Y, Z, nodeColour, 'LineWidth', lineWidth) 
            #Draw map nodes
            plot3(maps[:, 1], maps[:, 2], maps[:, 3], 'Marker', nodeMarker+
                'MarkerFaceColor', nodeColour, 'MarkerEdgeColor', nodeColour+
                'MarkerSize', nodeMarkerSize, 'LineStyle', 'none') 

    elif dim == 2:
        #two dimensional data
        #Draw data
        for k in range(nCls):
            ind = classes == _cls(k) 
            plot(data(ind, 1), data(ind, 2)+
                [markColour(k), markShape(k)]+
                'MarkerFaceColor', markColour(k)+
                'MarkerSize', markSize(k)) 

        if mapDraw:
            #Draw edges
            #Prepare arrays
            X=np.concatenate([maps[links[: 0], 0].T, maps[links[: 1], 0].T])
            Y=np.concatenate([maps[links[: 0], 1].T, maps[links[: 1], 1].T]) 
            #Draw edges
            plot(X, Y, nodeColour, 'LineWidth', lineWidth) 
            #Draw map nodes
            plot(maps[:, 1], maps[:, 2], 'Marker', nodeMarker+
                'MarkerFaceColor', nodeColour, 'MarkerEdgeColor', nodeColour+
                'MarkerSize', nodeMarkerSize, 'LineStyle', 'none') 

    else:
        #one dimensional data
        #Draw data
        for k in range(nCls):
            ind = classes == _cls(k) 
            plot(data(ind, 1), 0+
                [markColour(k), markShape(k)]+
                'MarkerFaceColor', markColour(k)+
                'MarkerSize', markSize(k)) 

        if mapDraw:
            #Draw edges
            #Prepare arrays
            X=np.concatenate([maps[links[: 0], 0].T, maps[links[: 1], 0].T])
            Y=np.zeros((2,links.shape[1])) 
            #Draw edges
            plot(X, Y, nodeColour, 'LineWidth', lineWidth) 
            #Draw map nodes
            plot(maps[:,1],0, 'Marker', nodeMarker+
                'MarkerFaceColor', nodeColour, 'MarkerEdgeColor', nodeColour+
                'MarkerSize', nodeMarkerSize, 'LineStyle', 'none') 


def interpol(X, y, nodes, r):
    #interpol calculates value of def y defined in data points X in each
    #node from nodes.
    #
    #Inputs:
    #   X is n-by-d data matrix with one data point in each row.
    #   y is n-by-1 vector with values of def. y(i) contains def
    #       value for point X(i, :).
    #   nodes is m-by-d matrix of nodes to calculate def values.
    #   r is smoothing parameter for def calculation.
    #
    #   memSize is constant for quick 

    # Calculate variances for all attributes
    smooth = -1 / (r * np.mean(np.var(X,axis=0),axis=0)) 
    # Calculate distances from each node to each data point
    dist = (np.sum(X**2,axis=1) + np.sum(nodes**2,axis=1).T) - 2 * (X @ nodes.T) 
    # Calclulate RBF in each point
    tmp = np.exp(dist * smooth) * y 
    # Calculate result
    res = np.sum(tmp,axis=0).T
    # Normalise result
    mins = np.min(res,axis=0) 
    res = (res - mins) / (np.max(res,axis=0) - mins) 
    return res

def formGrid(grid, maps, inter):
    # Step 1. Create list of all edges in grid
    # Unify description of triangles: sort nodes in ascend order
    grid = np.sort(grid, axis=1) 
    # Form list of all edges in grid
    edges = np.concatenate([grid[:,:2],  grid[:, 1:3],  grid[:, [0, 2]]]) 
    # Search unique values
    edges = np.unique(edges, axis=0) 
    # Step 2. Form list of nodes in new node list
    nN = maps.shape[0] 
    nE = edges.shape[0] 
    maps = np.concatenate([maps,  (maps[edges[:, 0], :] + maps[edges[:, 1], :]) / 2]) 
    inter = np.concatenate([inter,  (inter[edges[:, 0], :] + inter[edges[:, 1], :]) / 2]) 
    # Form list of indexes for nodes
    ind = np.zeros((nE + nN,nE + nN))
    #siz = ind.shape
    #indL = np.ravel_multi_index((edges[:, 0], edges[:, 1]), (siz))
    ind[edges[:, 0], edges[:, 1]] = nN + np.array(range(nN)) + nE 
    # Step 3. Form list of six nodes for each face
    face = np.concatenate([grid, ind[np.ravel_multi_index((grid[:, 0], grid[:, 1]),siz)]+
        ind[np.ravel_multi_index((grid[:, 1], grid[:, 2]),siz)]+
        ind[np.ravel_multi_index((grid[:, 0], grid[:, 2]),siz)]]) 
    # Step 4. Form final list of triangles
    grid = np.concatenate([face[:, [1, 3, 4]],  face[:, [0, 3, 5]], 
        face[:, [2, 4, 5]],  face[:, [3, 4, 5]]]) 
    return grid, maps, inter