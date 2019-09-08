import numpy as np

def LSqrt(x):
    #Lsqrt is the sqrt energy potential
    return np.sqrt(x)

def L1(x):
    #L1 is the L1 norm energy potential
    return np.abs(x)

def L1_5(x):
    #L1_5 is the L1.5 norm energy potential
    return np.abs(x)**(3/2)

def L2(x):
    #L2 is the L2 norm energy potential
    return x**2

def LLog(x):
    #L1 is the logarithmic energy potential
    return np.log1p(x)