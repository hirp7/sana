# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 15:43:00 2022

@author: hiro7
"""
from scipy import integrate
import matplotlib.pyplot as plt
import numpy as np
from numpy import sin,cos,tan


import os
import skrf as rf
import numpy as np
import scipy as sp
from numpy import real, log10, sum, absolute, pi, sqrt
import matplotlib.pyplot as plt
from scipy.optimize import fsolve,brentq,minimize, differential_evolution
from skrf.media import DefinedGammaZ0,RectangularWaveguide,CircularWaveguide
from scipy.optimize import newton


#integrate.quad(np.exp(1j*i*x*pi/L) for i in (np.arange(n)+1)

               

integrate.quad(lambda x: fun(x).re)


#path = r"C:\Users\hiro7\Desktop\Phideliti"
#os.chdir(path)

pi,sqrt,sin,cos,exp,jv,kv,yv,iv = [np.pi,np.sqrt,np.sin,np.cos,np.exp,sp.special.jv,sp.special.kv,sp.special.yv,sp.special.iv] 

inv = np.linalg.inv
abs = np.abs

kn,jvp,kvp,yvp,ivp = [sp.special.kn,sp.special.jvp,sp.special.kvp,sp.special.yvp,sp.special.ivp]


hankel2 = sp.special.hankel2
fminbound = sp.optimize.fminbound #This is used to determine u in dispersion relations to find x-value giving minimum y
gradient = np.gradient
e0 = 8.85e-12
m0 = 4*pi*1e-7
c0 = 2.9979e8
det = np.linalg.det

atan = np.arctan


root = sp.optimize.root_scalar


coth = lambda x:1/tanh(x)




def plot_3D(x,y,Z,plot_type = 'heat'):
    if plot_type == 'heat':
        fig, ax = plt.subplots()
        ax.imshow(Z, cmap='hot', interpolation='nearest')
    else:
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        X,Y = np.meshgrid(x,y)
        ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    return fig,ax


#definition of domain
"""
Solving Helmholtz Equation Over Rectangular with homogeneuous B.C.
using Green function
"""
#25.11
x = np.linspace(0,1,21)
y = np.linspace(0,1,21)

m = np.arange(3)+1
n = np.arange(3)+1
k0 = 10
a,b = [1,1]
def calc_x(x,y,eta,a,b,m,n,k0):
    fun_green = lambda x,xi,y,eta,a,b,m,n,k0:\
        sin(n*pi*x/a)*sin(n*pi*xi/a)*sin(m*pi*y/b)*sin(m*pi*eta/b)\
            /((m*pi/a)**2 + (n*pi/b)**2 - k0**2)
    return np.array([[fun_green(i,y,j,eta,a,b,m,n,k0) for i in x] for j in x])
#Homogeneus B.C. G_b.c. = 0
lis = []
[[lis.append(np.array([[calc_x(x,k,l,1,1,mi,nj,k0) for k in y] for l in y]).reshape(len(x)*len(y),len(x)*len(y))) for mi in m] for nj in n]

stacked = np.dstack(lis).transpose(2,0,1) 
summed = np.sum(stacked,axis = 0)*4/(a*b)

source = lambda xi,eta: 1 if np.round(xi,2) == np.round(eta,2) == 0.35 else 0
F = np.array([[source(i,j) for i in x]for j in y]).reshape(len(x)*len(y),1)

U = (summed@F).reshape(len(x),len(y))

fig,ax = plot_3D(x/a,y/b,U,plot_type = 'heat')
ax.set_xlabel("x/a")
ax.set_ylabel("y/b")
ax.set_title('Helmholz Equation with H.B.C.')



x = np.linspace(0,1,20)
y = np.linspace(0,1,10)




x = np.linspace(0,5,100)
#f_x = lambda x: 5 + integrate.quad()
A = lambda x:-15/7*x**2 + 20/7

pi = 3.141592
t = np.linspace(-pi,pi,100)

binary_generate = lambda x : 1 if(x >= 0.5) else 0
f1 = lambda x:-1 if (x >= -pi and x < 0) else 1 # else 1 if x >= 0 and x < 1  
#f1 = lambda x:-1 if (x <= 0) 
#plt.plot(t,np.array([f1(i) for i in t]))
# else 1 if x >= 0 and x < 1  



def test(x):
    return x
    #return np.exp(1j*x)

def quad_complex(fun,a,b):
    Re = integrate.quad(lambda x:fun(x).real,a,b)[0]
    Im = integrate.quad(lambda x:fun(x).imag,a,b)[0]
    return Re + 1j*Im

def fourier(fun,a,n):
    x = np.linspace(-a,a,100)
    L = a
    cn = np.array([1/(L)*(integrate.quad(lambda x:fun(x)*cos(i*x*pi/L),-L,L)[0] - 1j* \
                            integrate.quad(lambda x:fun(x)*sin(i*x*pi/L),-L,L)[0]) \
                   for i in np.arange(n)+1 ])
    approx = fun(0)/2 + np.array([cn[i]*np.exp(1j*pi*j*x/L) for i,j in enumerate(np.arange(n)+1) ]).sum(axis = 0)
    #Im = integrate.quad(lambda x:fun(x).imag,a,b)[0]
    return approx,cn

def run_fourier(fun,a,n): #fun is already vectorized n is a list, e.g. [1,5,10].
    x = np.linspace(-a,a,100)
    fig,ax = plt.subplots()
    ax.plot(x,f(x))
    approx,cn = np.array([fourier(f,a,i) for i in n]).transpose(1,0)
    [ax.plot(x,j,'--', label = f'n = {n[i]}') for i,j in enumerate(approx)]
    [ax.set_xlabel('x'),ax.set_ylabel('f(x)')]
    ax.legend()
    return 0
"""
Q1.1
f(x) = 0 (-5<x<0) 3 (0<x<5) 

"""
f = lambda x:0 if (x >= -5 and x < 0) else 3 # else 1 if x >= 0 and x < 1  
f = np.vectorize(f)
a = 5
n = [1,5,10]
"""
x = np.linspace(-a,a,100)
fig,ax = plt.subplots()
ax.plot(x,f(x))
approx,cn = np.array([fourier(f,a,i) for i in n]).transpose(1,0)
[ax.plot(x,j,'--', label = f'n = {n[i]}') for i,j in enumerate(approx)]
ax.legend()
"""
run_fourier(f,a,n)



"""
Q1.2 
f(x) = x [-2,2]
"""
f = lambda x:x# else 1 if x >= 0 and x < 1  
f = np.vectorize(f)
a = 2
n = [1,5,10]
run_fourier(f,a,n)


"""
Q1.3
f = |sin(x)| [-pi,pi]
"""
f = lambda x:np.abs(sin(x))# else 1 if x >= 0 and x < 1  
f = np.vectorize(f)
a = pi
n = [2,5,10]
run_fourier(f,a,n)






x = np.linspace(-pi,pi,100)
f_test,cn_test = fourier(test,pi,50)
fig,ax = plt.subplots()
ax.plot(x,test(x))
ax.plot(x,fourier(test,pi,50)[0])

plt.plot(cn_test.imag,'o')


"""
4/pi * sum_(n=1)^(inf) sin(2k-1)x/(2k-1)
"""
def fourier(fun,x,n):
    list0 = np.array([fun(x,i+1) for i in np.arange(n)])
    return list0.sum(axis = 0)

rect = lambda x,n: 4/pi*sin((2*n-1)*x)/(2*n-1) #
tri = lambda x,n: -(1/n/pi)**2*cos(n*x)  + (1/n/pi)*sin(n*x) + 2/3 




#a2 = lambda x: 4/pi*sin(x) + 4/3/pi*sin(3*x) + 4/5/pi*sin(5*x)
n_list = np.array([1,3,5,10])
fig,ax = plt.subplots(2,1)
ax[0].plot(t,np.array([f1(i) for i in t]),'--',label = 'original')
[ax[0].plot(t,fourier(sin0,t,i),label = f'n = {i}') for i in n_list]
#ax.plot(t,np.array([a2(i) for i in t]),label = 'analytical2')
ax[0].legend()

x1 = np.linspace(0,1,100)
#ax[1].plot(x1,np.array([lambda x: , ]))
ax[1].plot(x1,list(map(lambda x:x**2, x1)))
[ax[1].plot(x1,fourier(tri,x1,i),label = f'n = {i}') for i in n_list]


test_fun = lambda x,y: np.exp(x+y)
a = 3

test1 = integrate.quad(lambda x:np.exp(x),-a,a)
test2 = integrate.quad(lambda y:np.exp(y),-a,a)

test12 = integrate.dblquad(test_fun,-a,a,lambda x:-a,lambda x:a)

print(f'double integral = {test12[0]}, separated = {test1[0]*test2[0]}')
#Prada_test = integrate.dblquad(lambda theta,phi:U_test(theta,phi)*sin(theta),0,2*pi,lambda theta:0,lambda theta:pi)[0]

