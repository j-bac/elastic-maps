import numpy as np
from copy import deepcopy

def accumarray(idx_arr,val,size):
    ''' only for 2D index'''
    out = np.zeros(size)
    if hasattr(val,'shape'):
        for i in range(len(idx_arr)):
            out[idx_arr[i,0],idx_arr[i,1]]+= val[i]
        return out
    else:
        for i in idx_arr:
            out[i[0],i[1]]+=val
        return out

class MapGeometry:
    def __init__(self,dim):
        self.dimension = dim
        self.preproc = False
        self.PCs = []
    
    def getDimension(self):
        dim = self.dimension
        return dim
    
    def getInternalCoordinates(self):
        coord = self.internal
        return coord
    
    def getMappedCoordinates(self):
        coord = self.mapped
        return coord
    
    def getLinks(self):
        # def to access edges of map. 
        # links is k-by-2 matrix. Each row contains two numbers of
        # nodes which form one edge.
        links = self.links
        return links
    
    def getRibs(self):
        # Function to access ribs of map
        # ribs is k-by-3 matrix. Each row contains three numbers of 
        # nodes which form this rib.
        ribs = self.ribs
        return ribs
    
    def getDisp(self):
        # def to access disp of map
        # disp is non-negative number which presents the maximal
        # distance from data point to nearest original node.
        disp = self.disp
        return disp
    
    def init(self,data,_type='pci',_reduce=0):
        
        data = self.preprocessDataInit(data, _reduce)
        
        if 'random'.casefold() == _type.casefold():
            #Random generation
            #Calculate intervals for each coordinates
            mini = np.min(data,axis=0)
            maxi = np.max(data,axis=0)-mini
            #Generate random coordinates
            data = np.random.uniform(size=(self.internal.shape[0],data.shape[1]))
            #Scale coordinates
            data = data*maxi
            #shift coordinates
            self.mapped = data+mini
        elif 'randomSelection'.casefold() == _type.casefold():
            #Random selection
            #Generate vector of prototype numbers
            data = np.random.uniform(high = data.shape[0],size=(self.internal.shape[0],1))
            #Get selected points and put them as mapped coordinates
            self.mapped = data[data,:]
        elif 'pci'.casefold() == _type.casefold():
            # Principal component initialization
            # Get mean and PCs
            if self.preproc:
                # Data were preprocessed
                V = np.eye(data.shape[1], self.dimension)
                tmp = data[:,:self.dimension]
                meanDat = np.zeros((1, data.shape[1]))
            else:
                # Get required number of PCs:
                meanDat = self.means 
                V = self.PCs[:, :self.dimension+1]
                tmp = data @ V
            #Calculate mean and dispersion along each PCs
            mini = np.min(tmp,axis=0)
            maxi = np.max(tmp,axis=0)
            meant = np.mean(tmp,axis=0)
            disper = np.min(np.array([[meant - mini],[maxi - meant]]),axis=0)
            #Calculate mean and dispersion along internal coordinates
            minI = np.min(self.internal,axis=0)
            maxI = np.max(self.internal,axis=0)
            meanI = np.sum(self.internal,axis=0) / self.internal.shape[0]
            disP = np.min(np.array([[meanI - minI],[maxI - meanI]]),axis=0)
            #auxiliary calculations
            V = V * (disper / disP)
            #final values
            self.mapped=((self.internal- meanI) @ V.T) + meanDat
        else:
            raise ValueError('type "' +_type+ '" is not recognized as valid type of initialization')
        data, _  = self.associate(self.mapped, data)
        self.disp = np.sqrt(np.max(data,axis=0))
    
    def preprocessDataInit(self, data, reduce):

        #Perform data preprocessing if necessary (see description of
        #'reduce')
        #
        #Inputs:
        #   map is MapGeometry object to initialise map.
        #   data is n-by-m matrix with n data points and m coordinates for
        #       each point (each row is one data point)
        #   reduce is integer. If 'reduce' is positive and is
        #       less than n then specified number of the first principal
        #       components are used. If 'reduce' is zero and m>n then the
        #       first n-1 principal components is used. If 'reduce' is
        #       positive and is greater or equal to n or 'reduce' is zero
        #       and n>m then dimensionality reduction is not performed. If
        #       reduce is negative then -reduce PCs are calculated but
        #       dimensionality reduction is not performed.

        # Get sizes
        [n, m] = data.shape
        reduce = np.floor(reduce)
        if reduce >= m or (reduce == 0 and n > m):
            # Calculate 3 PCs and mean but do not apply preprocessing
            k = 3
            if k > m:
                k = m
            
            self.preproc = False
        elif reduce < 0:
            k = -reduce
            if k > m:
                k = m
            
            self.preproc = False
        else:
            # Define required number of PCs
            k = n - 1
            if reduce > 0 and reduce > k:
                k = reduce
            
            self.preproc = True
        # Search required number of PCs
        self.means = np.mean(data,axis=0)
        [_, D, VH] = np.linalg.svd(data-self.means,full_matrices=False)
        V = VH.T   # for consistency with matlab version
        D = D[:k]
        ind = D[::-1].argsort()
        V = V[:,ind]
        # Standardise direction of PCs
        ind = np.diag(V) < 0
        V[:, ind] = -V[:, ind]
        # Store results
        self.PCs = V
        
        # Preprocess data if it is required
        if self.preproc:
            data = self.preprocessData(data)
            
        return data
    
    def preprocessData(self, data):
        if not self.preproc:
            return
        data = (data - self.means) @ self.PCs
        return data

    def deprocessData(self, data):
        if not self.preproc:
            return

        data = (data  @ self.PCs.T) + self.means
        return data
    
    def project(self, points, _type, kind):
        #Project is the def to calculate projection of data point
        #(points) into map. There are d+1 types of projection for d
        #dimensional map: 0 means projection into nearest node of map, 1
        #means projection onto nearest edge of map, 2 means projection onto
        #nearest face of map. Projection can be calculated in the internal
        #or mapped coordinates. There are three input arguments for this
        #method: set of point to project, type of projection (integer
        #number) and coordinates space for projection: internal or
        #mapped.
        #
        #Inputs:
        #   map is MapGeometry object to use
        #   points is n-by-m matrix where m is number of mapped coordinates
        #       and n is number of points to project.
        #   type is type of projection: 0 or 1 for 1D maps, 0,1 or 2 for 2D
        #       maps.
        #   kind is one of words 'internal' for internal coordinates and
        #       'mapped'for mapped coordinates.
        coord, _ = projectPrime(self, self.mapped, points, _type, kind)
        return coord
        
    def projectPrime(self, nodes, points, _type, kind):
        #projectPrime is the def to calculate projection of data point
        #(points) into map. There are d+1 types of projection for d
        #dimensional map: 0 means projection into nearest node of map, 1
        #means projection onto nearest edge of map, 2 means projection onto
        #nearest face of map. Projection can be calculated in the internal
        #or mapped coordinates. There are three input arguments for this
        #method: set of point to project, type of projection (integer
        #number) and coordinates space for projection: internal or
        #mapped.
        #
        #Inputs:
        #   map is MapGeometry object to use
        #   nodes is the current state of the mapped nodes. It can be
        #       diffed from map.mapped. It is useful for estimation of
        #       calculated mapped nodes without fixing it into map object
        #   points is n-by-m matrix where m is number of mapped coordinates
        #       and n is number of points to project.
        #   type is type of projection: 0 or 1 for 1D maps, 0,1 or 2 for 2D
        #       maps.
        #   kind is one of words 'internal' for internal coordinates and
        #       'mapped'for mapped coordinates.
        #
        #Outputs:
        #   coord is the set of requested projections.
        #   dist is the vector of distances from points to map.
        # Check which type of coordinates is necessary to return
        cType = 'mapped'.casefold() == kind.casefold()
        N = points.shape[0]

        if _type == 0:
            #Projection to the nearest node
            #Search the nearest node
            [dist, num] = self.associate(nodes, points)
            if cType:
                coord = nodes[num,:]
            else:
                coord = self.internal[num,:]
        elif _type == 1:
            #projection onto nearest edge
            #Get array of edges end
            V2 = nodes[self.links[:,2],:].T
            #Form matrix of edges directions
            _dir = nodes[self.links[:,1],:].T -V2
            #Calculate squared length of edge directions
            _len = np.sum(_dir**2,axis=0)
            #Calculate projections length (l in documentation, matrix
            #analogue of (2))
            pr = (points*_dir)-np.sum(V2*_dir,axis=0)
            #Copy projections to normalize (l* in documentation)
            prn = pr/_len
            #Non negativity
            prn[prn<0] = 0
            #Cut too long projections (it is the same as step 3 of
            #algorithm in documentation)
            prn[prn>1] = 1
            #Calculate distances:
            dist = (np.sum(points**2,axis=1))+(np.sum(V2**2,axis=0)
                -2*points*V2+prn*(prn*_len-2*pr))
            #Select the nearest edge
            [dist, edge] = np.min(dist,axis=1)
            #form index to find length of projections
            ind = np.ravel_multi_index(np.concatenate((list(range(N+1)),edge)),([N,self.links.shape[0]]))
            if cType:
                coord = ((1-prn(ind))*nodes[self.links[edge,1],:])+(prn(ind)*nodes[self.links[edge,0],:])
            else:
                coord = (1-prn[ind])*self.internal[self.links[edge,1],:]+prn[ind]*self.internal[self.links[edge,0],:]
        elif _type == 2:
            #Projections onto face
            if not any(methods(self).casefold() ==  'getFaces'.casefold()):
                raise ValueError('request of the projection onto face for MapGeometry without methods "getFaces"')
            #Get faces
            face = self.getFaces
            #form auxiliary vectors
            Y2 = (nodes[face[:,2],:]).T
            Y20 = Y2-(nodes[face[:,0],:]).T
            Y21 = Y2-(nodes[face[:,1],:]).T
            Y10 = (nodes[face[:,1],:]-nodes[face[:,0],:]).T
            Y20Y20 = np.sum(Y20**2,axis=0)
            Y21Y20 = np.sum(Y20*Y21,axis=0)
            Y21Y21 = np.sum(Y21**2,axis=0)
            #Calculate projections
            A20 = np.sum(Y2*Y20,axis=0)-points*Y20
            A21 = np.sum(Y2*Y21)-points*Y21
            A10 = (A20-A21-np.sum(Y21*Y10,axis=0))/np.sum(Y10**2,axis=0)
            tmp = Y20Y20*Y21Y21-Y21Y20**2
            A0 = (A20*Y21Y21-A21*Y21Y20)/tmp
            A1 = (A21*Y20Y20-A20*Y21Y20)/tmp
            A20 = A20/Y20Y20
            A21 = A21/Y21Y21
            #Normalize projections
            A20N = A20
            A20N[A20N<0] = 0
            A20N[A20N>1] = 1
            A21N = A21
            A21N[A21N<0] = 0
            A21N[A21N>1] = 1
            A10[A10<0] = 0
            A10[A10>1] = 1
            tmp = A0<0
            A0[tmp] = 0
            A1[tmp] = A21N[tmp]
            tmp = A1<0
            A0[tmp] = A20N[tmp]
            A1[tmp] = 0
            tmp = (1-(A0+A1))<0
            A0[tmp] = A10[tmp]
            A1[tmp] = 1-A10[tmp]
            #Calculate distances
            dist = (np.sum(points**2,axis=1)+np.sum(Y2**2,axis=0)-2*points*Y2
                +(A0*Y20Y20)*(A0-2*A20)
                +(A1*Y21Y21)*(A1-2*A21)
                +2*(A0*A1)*Y21Y20)
            #Select the nearest face
            [dist, tmp] = np.min(dist,axis=1)
            #form index to find length of projections
            ind = np.ravel_multi_index(np.concatenate((list(range(N+1)),tmp)),[N,face.shape[0]])

            if cType:
                coord = (A0[ind]*nodes[face[tmp,0],:]
                    +A1[ind]*nodes[face[tmp,1],:]
                    +(1-A0[ind]-A1[ind])*nodes[face[tmp,2],:])
            else:
                coord = (A0[ind]*self.internal[face[tmp,0],:]
                    +A1[ind]*self.internal[face[tmp,1],:]
                    +(1-A0[ind]-A1[ind])*self.internal[face[tmp,2],:])
        else:
            raise ValueError('unacceptable type or projections')
        return coord, dist

    def associate(self,node, data):
        #associate identify the nearest node for each data point and
        #return the squared distance between selected node and data
        #point and number of nearest node.
        #
        #Inputs:
        #   node is n-by-k matrix of mapped coordinates for tested state
        #       of map, where n is number of nodes and m is dimension of
        #       data space.
        #   data is m-by-k data points to test, where m is number of
        #       points and k is dimension of data space.
        #
        #Outputs:
        #   dist is m-by-1 matrix of squared distances from data point to
        #       nearest node
        #   klass is m-by-1 vector which contains number of nearest node
        #       for each data point.

        dist = (np.sum(data ** 2, axis=1,keepdims=1)+np.sum(node**2, axis=1,keepdims=1).T) - 2 * (data @ node.T) 
        [dist, klas] = np.min(dist, axis=1), np.argmin(dist, axis=1) 
        return dist, klas

    def extend(self, data, val=1):
        #extend create extended version of map to reduce/prevent border
        #effect.
        #
        #Inputs:
        #   map is MapGeometry object to extend
        #   val is optional parameter to customise process:
        #       is greater of equal 1 is used to add val ribbons to each
        #           side of map.
        #       is positive number between 0 and 1 means maximal acceptable
        #           fraction of points which are projected onto map border.
        #       default value is 1
        #   data is n-by-m data points to test, where n is number of
        #       points and m is dimension of data space.

        # Check the val value
        if val < 0:
            raise ValueError(['Value of val attribute must be positive value'+
                ' between 0 and 1 for fraction restriction or'+
                ' positive integer to add val ribbons to each'+
                ' side of map']) 
        if val < 1:
            # Restriction for fraction of border cases
            dat = self.preprocessData(data) 
            newMap = deepcopy(self) 
            while newMap.borderCases(dat, newMap.getBorder()) > val:
                newMap = newMap.extendPrim() 
        else:
            val = np.floor(val) 
            newMap = self.extendPrim() 
            for k in range(2,val+1):
                newMap = newMap.extendPrim() 

        return newMap

    def FVU(self, data, node=None, _type=1):
        #Calculate fraction of variance unexplained for specified data and
        #nodes.
        #
        #Inputs:
        #   map is MapGeometry object to use
        #   data is set of data points
        #   node is the set of considered mapped nodes. If t is omitted or
        #       empty then the map.mapped is used.
        #   type is the type of projection: 0 means projection into nearest
        #       node of map, 1 means projection onto nearest edge of map, 2
        #       means projection onto nearest face of map. If this argument
        #       is omitted then 1 is used.
        if node==None:
            node = self.mapped
        #Calculate base variance
        N = data.shape[0] 
        meanS = np.sum(data,axis=0)/N 
        base = np.sum(np.reshape(data,(-1),order='F')**2)-N*np.sum(meanS**2) 

        #Get distances to map
        [_, dist] = self.projectPrime(node, data, _type, 'mapped') 

        #Calculate FVU
        fvu = np.sum(dist,axis=0)/base 

        return fvu

    def putMapped(self, newMapped):
        #This def is used for the putting the fitted mapped
        #coordinates of map.
        #
        #Inputs:
        #   newMapped is new matrix of mapped coordinates. It must have
        #       the same size as previously defined matrix
        if not self.mapped.shape==newMapped.shape:
            raise ValueError('Matrix newMapped must have the same size as matrix mapped') 
        self.mapped = newMapped 


    def borderCases(self, data, _list):
        #borderCases calculates fraction of border cases among all data
        #points.
        #
        #Inputs:
        #   map is MapGeometry object to use
        #   data is n-by-m data points to test, where n is number of
        #       points and m is dimension of data space.
        #   list is list of indices of border nodes.

        [_, ass] = self.associate(self.getMappedCoordinates, data) 
        # Calculate number of points for each node
        N = self.getMappedCoordinates.shape[0]
        tmp = accumarray(ass + 1, 1, [N + 1, 1]) 
        # Normalise and remove dummy element
        tmp = tmp[1:] 
        frac = np.sum(tmp[_list],axis=0) / data.shape[0] 
        return frac