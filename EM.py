from rect2DMap import rect2DMap
from MapGeometry import MapGeometry
from matlab_like_funcs import accumarray, isnumeric, histc
import numpy as np

def constStretch(*args):
    stretch = constStretching 
    return stretch

def constBend(*args):
    bend = constBending 
    return bend

def EM(
    _map, 
    data, 
    strFun = constStretch,
    constStretching = 0.7,
    bendFun = constBend,
    constBending = 0.7,
    weights = [],
    func = [],
    intervals = [],
    nInt = 5,
    delta = 1,
    _type = 'none',
    **kwargs
      ):
    
    #EM is def to fit map to data
    #   Syntax
    #   EM( map, data )
    #   EM( __, Name, Value )
    #   
    #Inputs:
    #   map is an object of MapGeometry class or subclass of this class.
    #   data is n-by-m matrix which is the data set to fit. Each row contains
    #       m coordinates of one data point.
    #   Name can be one of the following:
    #       'type' is one of the following strings:
    #          'hard' is hard map with stretch = 1 and bend = 1
    #           'medium' is more flexible map with stretch = 0.7 and bend = 0.7
    #           'soft' is soft map with stretch = 0.5 and bend = 0.5
    #           If 'type', 'stretch' and 'bending' are omitted then 'medium' is
    #           used. 
    #       'stretch' is a positive numeric value which is the value of
    #           stretching modulo or a def with syntax
    #               val = stretch( epoch )
    #           where epoch is number of epoch (see epoch definition below) and
    #           val is the nonnegative stretching modulo to use on specified
    #           epoch. Epochs are numerated from 1.
    #           Default value corresponds to type 'medium'
    #       'bend' is a positive numeric value which is the value of the
    #           bending modulo or a def with syntax 
    #               val = bend( epoch )
    #           where epoch is number of epoch (see epoch definition below) and
    #           val is the bending modulo to use on specified epoch. Epochs are
    #           numerated from 1. 
    #           Default value corresponds to type 'medium'
    #       'weights' is n-by-1 vector of weights for data points. Weights must
    #           be nonnegative.
    #       'intervals', intervals serves to specify user defined intervals.
    #           intervals is row vector. The first element must be zero. By
    #           default is created by usage of 'number_of_intervals' and
    #           'intshrinkage'. Maximal value M is calculated as maximum among
    #           data points of distance from data points to the nearest node of
    #           map after initiation. Then this value is multiplied by
    #           'intshrinkage'. All other borders are calculated as r(i) =
    #           M*i^2/p^2, where p is number_of_intervals'. Ignored if
    #           'potential' is not specified.
    #       'Number_of_intervals' specifies the number of intervals to
    #           automatic interval calculation. Default value is 5. Ignored if
    #           'potential' is not specified. 
    #       'intshrinkage' is fraction of maximal distance from data points to
    #           original map which is used for intervals shrinkage (see
    #           argument delta in defineIntervals). Default value is 1 (no
    #           shrinkage). Ignored if 'potential' is not specified.
    #       'potential' is majorant def for PQSQ. L2 distance without
    #           shrinkage is used if 'potential' is not specified.
    #
    # One epoch is fitting of map with fixed values of stretching and bending
    # modulo. This process can include several iterations of two step
    # algorithm:
    #   1. associate each data point with nearest node.
    #   2. recalculate node position.
    # Process of map fitting is stopped if new values of stretching and bending
    # modulo are the same as on previous epoch OR if both stretching and
    # bending modulo are zero.

    # Check the number of input attributes and types of the two first
    # attributes.
    if not issubclass(type(_map),MapGeometry):
        raise ValueError('Incorrect type of the "map" argument, it must be MapGeometry') 
    if not len(data.shape) == 2:
        raise ValueError('Incorrect type of the "data" argument, data must be a matrix') 
    
    # Data preprocessing
    data = _map.preprocessData(data)
    
    # Get sizes of data
    [n, dim] = data.shape

    # Default values of customisable variables
    
    # Decode varargin            
    if _type.casefold() == 'hard'.casefold():
        strFun = constStretch 
        constStretching = 1 
        bendFun = constBend 
        constBending = 1 
    elif _type.casefold() == 'medium'.casefold():
        strFun = constStretch 
        constStretching = 0.7 
        bendFun = constBend 
        constBending = 0.7 
    elif _type.casefold() == 'soft'.casefold():
        strFun = constStretch 
        constStretching = 0.5 
        bendFun = constBend 
        constBending = 0.5 
    
    # Check _type and length of weights
    if weights == []:
        weights = np.ones((n, 1)) 

    # Define total weights
    TotalWeight = np.sum(weights,axis=0) 
    weigh = weights 
    pFunc = [] 

    # Analyse PQSQ request
    if not(func==[]):
        if intervals==[]:
            #def has to create intervals by automatic way
            #nInt must be positive integer scalar
            if not(isinstance(nInt,int)) or not(np.isfinite(nInt)) or nInt < 1:
                raise ValueError(['Incorrect value of "number_of_intervals" argument'+
                    'It must be positive integer scalar']) 
            else:
                nInt = np.floor(nInt) 
                
            #delta has to be positive real scalar
            if isinstance(delta,complex) or not(np.isfinite(delta)) or delta < 0:
                raise ValueError(['Incorrect value of "intshrinkage" argument' +
                    'It must be positive real scalar']) 
            pFunc = definePotentialdef(_map.getDisp(), nInt, func, delta) 
        else:
            #intervals must contains non negative values in ascending order.
            #The first value must be zero.
            if intervals[0]!=0 or not(all(np.isfinite(intervals))) or any((intervals[1:-1]-intervals[0:-2])<=0):
                raise ValueError(['Incorrect values in argument intervals: intervals must'+
                    ' contains finite non negative values in ascending order.'+
                    ' The first value must be zero.']) 
                
            pFunc.intervals = [intervals.T, np.inf] 
            [pFunc.A, pFunc.B] = computeABcoefficients(intervals, func) 
    
    #Get initial state of nodes
    nodes = _map.getMappedCoordinates() 
    if nodes.shape[1] != dim:
        raise ValueError('Dimensions of mapped nodes and data must be the same') 
                                      
    N = nodes.shape[0] 
    #Form matrices B and C
    tmp = _map.getLinks() 
    B = np.diag(np.bincount(tmp.reshape(-1,order='F')))        
    tmp = accumarray(tmp, 1, size=[N, N]) 
    B = B - tmp - tmp.T 
    
    tmp = _map.getRibs() 
    C = np.diag(np.bincount(np.concatenate([tmp[:, [0]],  tmp[:, [2]]]).squeeze(),minlength=N)
        +np.bincount(tmp[:, 1],minlength=N)*4) 
    w = accumarray(tmp[:, [0,2]], 1,[N, N]) 
    tmp = accumarray(np.concatenate([tmp[:,[0,1]],tmp[:,[1,2]]]), 2, [N, N]) 
    C = C + w + w.T - tmp - tmp.T 

    # Start iterative process
    epoch = 1  # Number of iteration
    ass = np.zeros((n, 1))  # Initial associations. It is impossible combination
    qInd = np.zeros((n, 1)) 
    # Get initial modulo
    #stretch = strFun(epoch)
    #bend = bendFun(epoch) 
    stretch = constStretching
    bend = constBending
    while True:
        # Save old associations and q indices.
        oldAss = ass 
        oldQInd = qInd 
        # Find new associations
        [dist, ass] = _map.associate(nodes, data) 
        # Find indices for PQSQ if required
        if not(pFunc==[]):
            qInd = np.digitize(dist,pFunc.sqint)-1
            weigh =  weights * pFunc.A[qInd][:,None]

        # If nothing has changed then we have end of epoch
        if np.all(oldAss == ass) and np.all(oldQInd == qInd):
            epoch = epoch + 1 
            #tmp = strFun(epoch) 
            #tmp1 = bendFun(epoch) 
            tmp = stretch
            tmp1 = bend
            if tmp == 0 and tmp1 == 0:
                break 
            if abs(tmp - stretch) + abs(tmp1 - bend) == 0:
                break 
            stretch = tmp 
            bend = tmp1 
        
        # Form matrix A
        # For further robustness and so on we consider possibility of zeros in
        # ass and create dummy element
        ass = ass + 1 
        # Calculate number of points for each node
        tmp = np.bincount(ass.squeeze(), weigh.squeeze(), minlength=N + 1)
        # Normalise and remove dummy element
        NodeClusterRelativeSize = tmp[1:] / TotalWeight 
        # Create centroids
        NodeClusterCenters = np.zeros((N + 1, dim)) 
        for k in range(dim):
            NodeClusterCenters[:, k] = np.bincount(ass.squeeze(), (data[:, [k]] * weigh).squeeze(), minlength=N + 1) / TotalWeight 

        # Remove dummy element
        NodeClusterCenters = NodeClusterCenters[1:,:]
        
        # form SLAE
        SLAUMatrix = np.diag(NodeClusterRelativeSize) + stretch * B + bend * C 
        nodes = np.linalg.lstsq(SLAUMatrix, NodeClusterCenters,rcond=None)[0] 
    
        # Restore ass
        ass = ass - 1 
    
    # Put new nodes into map
    _map.putMapped(nodes) 

def definePotentialdef( 
    x,
    number_of_intervals, 
    potential_def_handle, 
    delta=1
    ):
    #definePotentialdef defines "uniform in square" intervals for trimming
    #threshold x and specified number_of_intervals.
    #   x is upper boundary of the interval last but one.
    #   number_of_intervals is required number of intervals.
    #   potential_def_handle is def handler for coefficients
    #       calculation.
    #   delta is coefficient of shrinkage which is greater than 0 ang not
    #       greater than 1.
    #Output argument potentialdef is structure with three fields:
    #   intervals is matrix m-by-number_of_intervals. Each row contains
    #       number_of_intervals values of thresholds for intervals and one
    #       additional value Inf
    #   A and B are the m-by-number_of_intervals matrices with quadratic
    #       defs coefficients
    class empty_class:
        pass

    potentialdef = empty_class()
    p = int(number_of_intervals - 1) 
    
    #intervals is the product of row and maximal coefficient multiplied by delta:
    intervals = (x * delta) * (np.array(range(p+1)) / p) ** 2 
    
    potentialdef.intervals = np.concatenate([intervals, np.array([np.inf])],axis=None) 
    potentialdef.sqint = potentialdef.intervals ** 2 
    [potentialdef.A,potentialdef.B] = computeABcoefficients(intervals, potential_def_handle) 
    return potentialdef

def computeABcoefficients(intervals, potential_def_handle):
    #PQSQR_computeABcoefficients calculates the coefficients a and b for
    #quadratic fragments of potential def.
    #   intervals is the 1-by-K matrix of intervals' boundaries without final
    #       infinity boundary.
    #   potential_def_handle is a handle of majorant def.
    #Get dimensions of intervals
    p = len(intervals)

    #Preallocate memory
    A = np.zeros(p) 
    B = np.zeros(p) 

    #Calculate value of def all boundaries
    pxk = potential_def_handle(intervals) 
    sxk = intervals**2 

    A[:p-1] = (pxk[:p-1]-pxk[1:p])/(sxk[:p-1]-sxk[1:p]) 
    B[:p-1] = (pxk[1:p]*sxk[:p-1]-pxk[:p-1]*sxk[1:p]) / (sxk[:p-1]-sxk[1:p]) 
    B[p-1] = pxk[p-1] 
    return A,B