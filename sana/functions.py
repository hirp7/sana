import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

pi = np.pi
sin = np.sin
cos = np.cos
log = np.log
inv = np.linalg.inv


def FDM(x0,xn,n,a,b,c,f,alpha,beta):
    """
    Finite Differential Method in 1D coordination system
    ay'' + by' + c = f(x)
    
    alpha and beta are boundary condition at the begin and the end x0 and xn, respectively
    f(x0) = alpha, f(xn) = beta
    
    a(x),b(x),c(x) can be either a function of x or constants. (Scalar or Vector with the same length of x
    """
        
    inv = np.linalg.inv
    x,h = np.linspace(x0,xn,n,retstep = True)
    def sc2vec(a):
        if np.isscalar(a)==True:
            a = np.ones(len(x))*a
        return a
    a,b,c = list(map(sc2vec,[a,b,c]))
    #define A
    diagonal = h**2*c - 2*a
    diagonal_side = a[1:] + h*b[1:]/2
    A = np.diag(diagonal) + np.diag(diagonal_side,1) + np.diag(diagonal_side,-1)   
    A = 1/h**2*A

    #define F
    F = f(x)
    F[0] = F[0] - (a[0]/h**2 - b[0]/2/h)*alpha
    F[n-1] = F[n-1] - (a[n-1]/h**2 + b[n-1]/2/h)*beta
    Y = inv(A)@F
    return Y,A,F,x
