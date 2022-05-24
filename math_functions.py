import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from numpy import sin,cos,tan,pi, sqrt, log, log10
from scipy import integrate

import scipy as sp
import numpy as np
import skrf as rf

#from scipy.optimize import fsolve,brentq,minimize, differential_evolution
#from scipy.optimize import newton

"""
def quad_complex(fun,a,b):
    Re = integrate.quad(lambda x:fun(x).real,a,b)[0]
    Im = integrate.quad(lambda x:fun(x).imag,a,b)[0]
    return Re + 1j*Im
"""
"""
Fourier series is an expansion of periodic signal as a linear combination of sines and cosines, while Fourier transform is the process or function used to convert signals from the time domain to frequency domain.
"""

"""
def fourier(fun,a,n):
    tol = 1e-6
    x = np.linspace(-a,a,100)
    L = a
    cn = np.array([1/(L)*(integrate.quad(lambda x:fun(x)*cos(i*x*pi/L),-L,L)[0] - 1j* \
                            integrate.quad(lambda x:fun(x)*sin(i*x*pi/L),-L,L)[0]) \
                   for i in np.arange(n)+1 ])
    approx = (fun(0+tol)+fun(0-tol))/2 + np.array([cn[i]*np.exp(1j*pi*j*x/L) for i,j in enumerate(np.arange(n)+1) ]).sum(axis = 0)
    #approx =  np.array([cn[i]*np.exp(1j*pi*j*x/L) for i,j in enumerate(np.arange(n)+1) ]).sum(axis = 0)

    #Im = integrate.quad(lambda x:fun(x).imag,a,b)[0]
    return approx,cn
"""

def fourier_expansion(fun,L,n):
    """
    fun:a periodic function
    L:interval length
    n:the number of harmonics
    """
    #tol = 1e-6
    x = np.linspace(-L,L,100)
    

    
    a0 = 1/L/2 * integrate.quad(fun,-L,L)[0]
    an = 1/L * np.array([integrate.quad(lambda x:fun(x)*cos(i*x*pi/L),-L,L)[0] for i in np.arange(n)+1])
    bn = 1/L * np.array([integrate.quad(lambda x:fun(x)*sin(i*x*pi/L),-L,L)[0] for i in np.arange(n)+1])
    approx = a0 + np.array([an[i]*cos(pi*j*x/L) + bn[i]*sin(pi*j*x/L) for i,j in enumerate(np.arange(n)+1)]).sum(axis = 0)

    """
    c0 = 1/L/2 * integrate.quad(fun,-a,a)[0]

    cn = np.array([1/(L)*(integrate.quad(lambda x:fun(x)*cos(i*x*pi/L),-L,L)[0] - 1j* \
                            integrate.quad(lambda x:fun(x)*sin(i*x*pi/L),-L,L)[0]) \
                   for i in np.arange(n)+1 ])
    approx = c0 + np.array([cn[i]*np.exp(1j*pi*j*x/L) for i,j in enumerate(np.arange(n)+1) ]).sum(axis = 0)
    """
    return approx,np.array([an,bn])



def plot_fourier(fun,a,n):
    """
    fun: the vectorized function
    a(float): the interval length
    n(list): the list of the required harmonics number, e.g. [1,5,10].
    """
    x = np.linspace(-a,a,100)
    fig,ax = plt.subplots()
    ax.plot(x,fun(x),label = 'exact')
    approx,cn = np.array([fourier_expansion(fun,a,i) for i in n]).T
    [ax.plot(x,j,'--', label = f'n = {n[i]}') for i,j in enumerate(approx)]
    [ax.set_xlabel('x'),ax.set_ylabel('f(x)')]
    ax.legend()
    """
    VisibleDeprecationWarning: Creating an ndarray from ragged nested sequences
    (which is a list-or-tuple of lists-or-tuples-or ndarrays with different lengths
     or shapes) is deprecated. If you meant to do this, you must specify 
    'dtype=object' when creating the ndarray.
    """
    return fig,ax



"""
Solving differential equation using fourier series expansion
"""

#Q2.1 f''(x) +3f'(x) = cost x(0) = 0, x'(0) = 0



"""
    
There is convplution function in numpy as well. how to use is shown below
np.convolve([1, 2, 3], [0, 1, 0.5])
Out[22]: array([0. , 1. , 2.5, 4. , 1.5])
mode{‘full’, ‘valid’, ‘same’}, optional

    ‘full’:

        By default, mode is ‘full’. This returns the convolution at each point of overlap, with an output shape of (N+M-1,). At the end-points of the convolution, the signals do not overlap completely, and boundary effects may be seen.
    ‘same’:

        Mode ‘same’ returns output of length max(M, N). Boundary effects are still visible.
    ‘valid’:

        Mode ‘valid’ returns output of length max(M, N) - min(M, N) + 1. The convolution product is only given for points where the signals overlap completely. Values outside the signal boundary have no effect.




"""
x = np.linspace(0,10,100)
laplace = lambda x: integrate.quad(lambda x:x**2*np.exp(-x),0,np.inf)[0]



y = np.array([laplace(i) for i in x])
fig,ax = plt.subplots(2,1)
ax[0].plot(x)
ax[1].plot(y)

integrate.quad(lambda x:np.exp(-2*x),0,np.inf)

freq = 1
f = lambda x,freq:0 if (x >= -3 and x <= 0) else np.sin(2*pi*freq*x)  # else 1 if x >= 0 and x < 1  
f = np.vectorize(f)
a = 0.3
n = [2,5,10,30]
fig,ax = plot_fourier(lambda x:f(x,freq),a,n)
fig.suptitle('f(x) = 0 (x >= -3 and x <= 1) else 10 ', fontsize=16)



f = lambda x,a:a*x# else 1 if x >= 0 and x < 1  
f = np.vectorize(f)

L =2
approx,[an,bn] = fourier_expansion(lambda x:f(x,2),L,20)

dev = f(np.linspace(-L,L,100),2) - approx

fig,ax = plt.subplots(2,1)
ax[0].plot(dev)
ax[1].plot(np.abs(bn))

a = pi
n = [2,5,20]
fig,ax = plot_fourier(lambda x:f(x,2),a,n)
fig.suptitle(' f(x) = |sin(x)| ', fontsize=16)

"""
Q1.1
f = lambda x:0 if (x >= -3 and x <= 1) else 10 # else 1 if x >= 0 and x < 1  
f = np.vectorize(f)
a = 3
n = [1,5,10]
fig,ax = plot_fourier(f,a,n)
fig.suptitle('f(x) = 0 (x >= -3 and x <= 1) else 10 ', fontsize=16)


Q1.2
f = lambda x:np.abs(sin(x))# else 1 if x >= 0 and x < 1  
f = np.vectorize(f)
a = pi
n = [2,5,10]
fig,ax = plot_fourier(f,a,n)
fig.suptitle(' f(x) = |sin(x)| ', fontsize=16)

Q1.3
f = lambda x:np.abs(x)# else 1 if x >= 0 and x < 1  
f = np.vectorize(f)
a = 2
n = [1,5,10]
fig,ax = plot_fourier(f,a,n)
fig.suptitle('f(x) = |x|', fontsize=16)
"""

"""
to fourier transform
"""
L = 2*pi
f = lambda x:sin(x) + sin(3*x) + sin(10*x)
test = fourier_expansion(f,L,20)
x = np.linspace(-L,L,100)

fig,ax = plt.subplots(2,1)
ax[0].plot(x,f(x))
[ax[0].set_xlabel('t'),ax[0].set_ylabel('magnitude')]
ax[1].plot(test[1][1],'o')

[ax[1].set_xlabel('frequency'),ax[1].set_ylabel('magnitude')]
fig.tight_layout(pad = 1.2)

L = 6.35e-3
d = 11.43e-3
t = 1.27e-3
h = 24.13e-3
w = 0.635e-3
Er = 8.875

#basis function?

Jz_d_t = lambda x: 1/sqrt(w**2 - x**2)
Jx_d_t = lambda x: x*sqrt(w**2 - x**2)
Ez = lambda x: sqrt(a**2 - x**2)
Ex = lambda x: x/sqrt(a**2 - x**2)

#We need to transform these functions into that in frequency domain
#using bessel functions

#space -> spectral domain by fourier transform

#AssertionErrorhttps://homepages.inf.ed.ac.uk/rbf/HIPR2/fourier.htm










