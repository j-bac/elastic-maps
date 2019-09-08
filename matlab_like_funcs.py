import numpy as np

### implementing matlab functions missing in python ###

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

def isnumeric(obj):
    try:
        obj + 0
        return True
    except:
        return False
    
def histc(X, bins):
    map_to_bins = np.digitize(X,bins)
    r = np.zeros((len(X[0,:]),len(bins)))
    for j in range(len(map_to_bins[0,:])):
        for i in map_to_bins[:,j]:
            r[j,i-1] += 1
    return [r, map_to_bins]