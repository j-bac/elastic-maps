from MapGeometry import MapGeometry
from matlab_like_funcs import accumarray, isnumeric, histc
import numpy as np


class rect2DMap(MapGeometry):
    #OneDMap is descendant for piecewise linear map
    #   Constructor contains two argument which are number of rows and
    #   number of columns
    def __init__(self,rows,cols):
        # Create map
        super().__init__(2)
        # Store size
        self.sizes = [rows, cols] 
        # Calculate internal coordinates of nodes
        N=rows*cols 
        a1=np.tile(np.array(range(cols)),(rows,1)) 
        a2=np.tile(np.array(range(rows))[:,None],(1,cols)) 
        self.internal = np.array([a2.reshape(-1),a1.reshape(-1)]).T        
        # Form array of links
        A=np.reshape(np.array(range(N)),(rows,cols),order='F') 
        B=A[:-1,:] 
        C=A[1:,:] 
        D=A[:,:-1] 
        E=A[:,1:]
        self.links = np.concatenate([
            np.concatenate([B.reshape((-1,1),order='F'), C.reshape((-1,1),order='F')],axis=1),
            np.concatenate([D.reshape((-1,1),order='F'), E.reshape((-1,1),order='F')],axis=1)
        ])
        # Form array of ribs
        B1=A[0:-2,:] 
        B2=A[1:-1,:] 
        B3=A[2:,:]
        C1=A[:,:-2]
        C2=A[:,1:-1] 
        C3=A[:,2:] 
        self.ribs = np.concatenate([
            np.concatenate([B1.reshape((-1,1),order='F'), B2.reshape((-1,1),order='F'), B3.reshape((-1,1),order='F')],axis=1),
            np.concatenate([C1.reshape((-1,1),order='F'), C2.reshape((-1,1),order='F'), C3.reshape((-1,1),order='F')],axis=1)
        ]) 
        # Form array of faces
        B1 = A[:-1,:-1]
        B2 = A[1:,:-1]
        B3 = A[:-1,1:]
        C1 = A[:-1,1:]
        C2 = A[1:,:-1]
        C3 = A[1:,1:]
        self.faces = np.concatenate([
            np.concatenate([B1.reshape((-1,1),order='F'), B2.reshape((-1,1),order='F'), B3.reshape((-1,1),order='F')],axis=1),
            np.concatenate([C1.reshape((-1,1),order='F'), C2.reshape((-1,1),order='F'), C3.reshape((-1,1),order='F')],axis=1)
        ]) 
        # Set mapped coordinates to empty set
        self.mapped = [] 

    def getFaces(self):
        #def to access to the faces of map.
        #face is k-by-3 matrix. Each row contains three numbers of
        #nodes which form one face.
        face=self.faces
        return face

    def extendPrim(self):
        # Get size of existing map
        n = self.sizes(1) 
        m = self.sizes(2) 
        N = n + 2 
        M = m + 2 
        NM = N * M 
        # Create new map with greater size
        newMap = rect2DMap(N, M) 
        # Form list of old nodes
        ind = np.array(range(1,NM)) 
        ind[newMap.getBorder] = []
        # Copy known positions of nodes into new map
        newMap.mapped[ind, :] = self.mapped 
        # Calculate positions of the nodes added to left side
        ind = np.array(range(1,n)) + 1 
        newMap.mapped[ind, :] = (2 * newMap.mapped[ind + N, :]
            - newMap.mapped[ind + 2 * N, :] )
        # Calculate positions of the nodes added to right side
        ind = NM - N + np.array(range(2,NM+1)) - 1 
        newMap.mapped[ind, :] = (2 * newMap.mapped[ind - N, :]
            - newMap.mapped[ind - 2 * N, :] )
        # Calculate positions of the nodes added to bottom side
        ind = np.array(range(0,NM,N)) 
        newMap.mapped[ind, :] = (2 * newMap.mapped[ind + 1, :]
            - newMap.mapped[ind + 2, :] )
        # Calculate positions of the nodes added to top side
        ind = np.array(range(N-1,NM,N)) 
        newMap.mapped[ind, :] = (2 * newMap.mapped[ind - 1, :]
            - newMap.mapped[ind - 2, :] )
        
        return newMap

    def getBorder(self):
        # Get map sizes
        n = self.sizes(1) 
        m = self.sizes(2) 
        # Form list of border nodes
        res = [np.array(range(2,n)) - 1,                  # left edge
            (m - 1) * n + np.array(range(2,n)) * m - 1,   # right edge
            (np.array(range(m)) - 1) * n + 1,             # bottom edge
            np.array(range(1,m)) * n]                     # top edge
        return res 