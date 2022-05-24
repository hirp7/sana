# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 09:35:04 2021

@author: hiro7
"""


"""
Abstruct:
DWG propagation characteristics such as alpha and beta are analytically calculated.
Attributes and methods are summarized in a class 'CircularDWG' or 'RectangularDWG'
One can define each DWG class with given parameters: a, b, Er1, Er2, f1, f2, and npoints,
where a and b is longline and shortline, respectively, in a rectangular structure, or the diameter of 
circular structure. Er1 and Er2 are dielectric constant of core and cladding material. tanD is tangent delta of core material.
f1 and f2 are start and stop frequecy analyzed, and finally npoints is number of points.
e.g.) 
    CDWG1 = CircularDWG(a=1.58e-3/2, Er1 = 2.25, Er2 = 1.0,tanD = 0.0003,f1 = 70e9,f2 = 150e9,n=1000)
    RDWG2 = RectangularDWG(a=3.1e-3,b = 1.55e-3, Er1 = 2.25, Er2=1,tanD = 0.0003,f1=70e9,f2=150e9,n = 1000)

Once you define DWG class, you can get its attributes(analytically calculated results) and run some functions (methods) such as plotting.

e.g.) -Attributes
    CDWG1.beta, CDWG1.alpha, CDWG1,gamma #give us phase, attenuation, and propagation constant in complex 
    CDWG1.GD #giving us the group delay per meter
    CDWG1.media #give us a media class of scikit-rf
      -Methods(function)
    CDWG1.plot_summary() #give us the plot summary showing Omega-beta diagram, normalized V-b diagram, Group delay, and loss 
    CDWG1.S(length,z0) #give us a Scikit-RF Network class for the defined DWG with a given length and characteristic impedance

Scikit-RF Network class has a wide variety functions and attributes. look at the original website https://scikit-rf.readthedocs.io/en/latest/index.html
Export as a touchstone file is also available by write_touchstone() command #source'https://scikit-rf.readthedocs.io/en/latest/api/generated/skrf.network.Network.write_touchstone.html'
    
e.g.) CDWG1.S(length,z0).write_touchstone('test.s2p')    
    Note: characteristic impedance analysis is currently unavailable. Instead of exact value, a typical constant value 110 Ohm is used.
"""

import os
import skrf as rf
import numpy as np
import scipy as sp
from numpy import real, log10, sum, absolute, pi, sqrt
import matplotlib.pyplot as plt
from scipy.optimize import fsolve,brentq,minimize, differential_evolution
from skrf.media import DefinedGammaZ0,RectangularWaveguide,CircularWaveguide
from scipy.optimize import newton


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

a = [1,2,3,4,5]
#[abs(x) for x in a] is equal to
map(abs,a)


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


coth = lambda x:1/tanh(x)



sin = np.sin
coth = lambda x:1/np.tanh(x)
#Parseval
W = 1e-3
Er = 3.6
e0 = 8.85e-12
h = 0.5e-3


def perseval(alpha):
    roh = lambda alpha:2*sin(alpha*W/2)/(alpha*W/2) - (sin(alpha*W/4)/(alpha*W/4))**2 
    g = lambda alpha: 1/(np.abs(alpha)*(1 + Er*coth(np.abs(alpha)*h)))

    return 1/e0*roh(alpha)*roh(alpha)*g(alpha)

alpha = np.linspace(1,5e4,1000)
fig,ax = plt.subplots()
ax.plot(alpha,perseval(alpha))

#    a0 = 1/L/2 * integrate.quad(fun,-L,L)[0]
lim = np.array([1e3,5e3,1e4,5e4,1e5,5e5,1e6])
tol = 1e-1
per_int = lambda lim: quad(perseval,-lim,tol)[0] + quad(perseval,tol,lim)[0]
ROH = np.array([per_int(i) for i in lim])

C = 2*pi/ROH
plt.plot(C*1e12)
#We can see it will converge to around 70 pF

W = 0.5
er = 3.6
e0 = 8.85e-12
u0 = 4*pi*1e-7
ur = 1
k1 = sqrt(omega**2*er*ur*e0*u0)
k2 = sqrt(omega**2*e0*u0)

k0 = k2

b11 = lambda alpha,beta: 1j*alpha*((k0**2 - beta**2)/(k1**2 - beta**2) - 1)
b12 = lambda alpha,beta: omega*u0*gamma1/beta*(gamma2/gamma1  + (k0**2 - beta**2)/(k1**2 - beta**2)*tanh(gamma1*h))
b21 = lambda alpha,beta: omega*e0*gamma1/beta*(gamma2/gamma1 + er*(k0**2 - beta**2)/(k1**2 - beta**2)*coth(gamma1*h))
b22 = lambda alpha,beta: -b11(alpha,beta)
det = lambda alpha,beta: b11(alpha,beta)*b22(alpha,beta) - b12(alpha,beta)*b21(alpha,beta)
F1 = lambda beta:omega*u0*gamma1*tanh(gamma1*h)/(1j*(k0**2 - beta**2))

G11 = lambda alpha,beta: 1/det(alpha,beta)*(F1(beta)*b22(alpha,beta) \
                        + alpha*beta/(k1**2 - beta**2)*b12(alpha,beta))
G12 = lambda alpha,beta: b12(alpha,beta)/det(alpha,beta)
G21 = lambda alpha,beta: gamma2/det(alpha,beta)*(F1(beta)*b21(alpha,beta)\
                        + alpha*beta/(k1**2 - beta**2)*b11(alpha,beta))
G22 = lambda alpha,beta: gamma2*b11(alpha,beta)/det(alpha,beta)


Jz1 = lambda x: 1 if (abs(x)<=W/2) else 0
Jz2 = lambda x: 2/W*x if (x >= 0 and x < W/2) else -2/W*x if (x < 0 and x > -W/2) else 0
Jx1 = lambda x: 1 if (x >= 0 and x < W/2) else -1 if (x < 0 and x > -W/2) else 0
#Jx2 = lambda x: 4/W*x if (x >= -W/4 and x < W/4) else -4/W*x if ((x < -W/4 and x > -W/2) or (x < W/2 and x > W/4)) else 0
Jx2 = lambda x: 2/W*x if (x >= -W/4 and x < W/4) else -2/W*x - 4*W/2 if (x < -W/4 and x > -W/2) else -2/W*x + 4*W/2  if (x < W/2 and x > W/4) else 0


def test(alpha,beta,f):
    omega  = 2*pi*f
    k1 = sqrt(omega**2*er*ur*e0*u0)
    k2 = sqrt(omega**2*e0*u0)
    k0 = k2

    gamma1 = alpha**2 + beta**2 - k1**2
    gamma2 = alpha**2 + beta**2 - k2**2

    b11 = lambda alpha,beta: 1j*alpha*((k0**2 - beta**2)/(k1**2 - beta**2) - 1)
    b12 = lambda alpha,beta: omega*u0*gamma1/beta*(gamma2/gamma1  + (k0**2 - beta**2)/(k1**2 - beta**2)*tanh(gamma1*h))
    b21 = lambda alpha,beta: omega*e0*gamma1/beta*(gamma2/gamma1 + er*(k0**2 - beta**2)/(k1**2 - beta**2)*coth(gamma1*h))
    b22 = lambda alpha,beta: -b11(alpha,beta)
    det = lambda alpha,beta: b11(alpha,beta)*b22(alpha,beta) - b12(alpha,beta)*b21(alpha,beta)
    F1 = lambda beta:omega*u0*gamma1*tanh(gamma1*h)/(1j*(k0**2 - beta**2))
    
    G11 = lambda alpha,beta: 1/det(alpha,beta)*(F1(beta)*b22(alpha,beta) \
                            + alpha*beta/(k1**2 - beta**2)*b12(alpha,beta))
    G12 = lambda alpha,beta: b12(alpha,beta)/det(alpha,beta)
    G21 = lambda alpha,beta: gamma2/det(alpha,beta)*(F1(beta)*b21(alpha,beta)\
                            + alpha*beta/(k1**2 - beta**2)*b11(alpha,beta))
    G22 = lambda alpha,beta: gamma2*b11(alpha,beta)/det(alpha,beta)





    return G11(alpha,beta),det(alpha,beta),F1(beta)     
W = 5
Jz1 = lambda x: 1 if (abs(x)<=W/2) else 0
Jz2 = lambda x: 2/W*x if (x >= 0 and x < W/2) else -2/W*x if (x < 0 and x > -W/2) else 0
Jx1 = lambda x: 1 if (x >= 0 and x < W/2) else -1 if (x < 0 and x > -W/2) else 0
#Jx2 = lambda x: 4/W*x if (x >= -W/4 and x < W/4) else -4/W*x if ((x < -W/4 and x > -W/2) or (x < W/2 and x > W/4)) else 0
Jx2 = lambda x: 2/W*x if (x >= -W/4 and x < W/4) else -2/W*x - 4*W/2 if (x < -W/4 and x > -W/2) else -2/W*x + 4*W/2  if (x < W/2 and x > W/4) else 0



def fourier(fun,f):
    operator_e = lambda f,t: cos(2*pi*f*t)
    operator_o = lambda f,t: sin(2*pi*f*t)

    even = integrate.quad(lambda t:fun(t)*operator_e(f,t),-np.inf,np.inf)[0]
    odd = integrate.quad(lambda t:fun(t)*operator_o(f,t),-np.inf,np.inf)[0]
    return even + 1j*odd

Jz1_FT,Jz2_FT,Jx1_FT,Jx2_FT = [lambda f: fourier(i,f) for i in [Jz1,Jz2,Jx1,Jx2]]
W = 8

fig,ax = plt.subplots(2,2)
ax = ax.reshape(4)
x = np.linspace(-10,10,100)
[i.plot(x,list(map(j,x))) for i,j in zip(ax,[Jz1,Jz2,Jx1,Jx2])]
k = np.linspace(-1,1,100)

fig,ax = plt.subplots(2,2)
ax = ax.reshape(4)
[i.plot(k,np.array(list(map(j,k))).imag) for i,j in zip(ax,[Jz1_FT,Jz2_FT,Jx1_FT,Jx2_FT])]


plt.plot(k,np.array(list(map(Jz2_FT,k))).imag)

w = 3.0
l = w/2

Exm_x = lambda x,r:cos(2*(r-1)*pi*x/w)/sqrt((w/2)**2 - x**2) if abs(x)<w/2 else 0 
Exm_z = lambda z,s:cos((2*s-1)*pi*z/l)/sqrt((l/2)**2 - z**2) if abs(z)<l/2 else 0

Ezn_x = lambda x,r:sin(2*r*pi*x/w)/sqrt((w/2)**2 - x**2) if abs(x)<w/2 else 0
Ezn_z = lambda z,s:sin((2*s-1)*pi*z/l)/sqrt((l/2)**2 - z**2) if abs(z)<l/2 else 0 
Exm = lambda x,z,r,s: Exm_x(x,r)*Exm_z(z,s)
Ezn = lambda x,z,r,s: Ezn_x(x,r)*Ezn_z(z,s)

EXM_X = lambda f:fourier(lambda t:Exm_x(t,1),f)
EZN_X = lambda f:fourier(lambda t:Ezn_x(t,1),f)

EXM_Z = lambda f:fourier(lambda t:Exm_z(t,1),f)
EZN_Z = lambda f:fourier(lambda t:Ezn_z(t,1),f)


x = np.linspace(-5,5,100)
alpha = np.linspace(-1,1,100)


fig,ax = plt.subplots(2,1)
ax[0].plot(x,list(map(lambda x:Exm_x(x,1),x)),label = "Exm_x")
ax[0].plot(x,list(map(lambda x:Ezn_x(x,1),x)),label = 'Ezn_x')
ax[0].legend()
ax[1].plot(x,list(map(lambda x:EXM_X(x),alpha)),label = "F(Exm_x)")
ax[1].plot(x,np.array(list(map(lambda x:EZN_X(x),alpha))).imag,label = "F(Ezn_x)")
ax[1].legend()
fig.suptitle("Basis function and fourier transform")

fig,ax = plt.subplots(2,1)
ax[0].plot(x,list(map(lambda x:Exm_z(x,1),x)),label = "Exm_z")
ax[0].plot(x,list(map(lambda x:Ezn_z(x,1),x)),label = 'Ezn_z')
ax[0].legend()
ax[1].plot(x,list(map(lambda x:EXM_Z(x),alpha)),label = "F(Exm_z)")
ax[1].plot(x,np.array(list(map(lambda x:EZN_Z(x),alpha))).imag,label = "F(Ezn_z)")
ax[1].legend()
fig.suptitle("Basis function and fourier transform")


#ax.plot(x,list(map(Ex3,x)))


test_fun = lambda t,a: sin(a*t)/t if abs(t)<20 else 0
Test_fun = lambda f:fourier(lambda t:test_fun(t,2),f)

t = np.linspace(-2*pi,2*pi,100)
f = np.linspace(-1,1,100)

step = lambda x: 0 if x < 0 else 1
fig,ax = plt.subplots(2,1)
ax[0].plot(t,np.array(list(map(lambda t:test_fun(t,2),t))).real)
ax[1].plot(f,np.array(list(map(Test_fun,f))).real)
ax[1].plot(f,np.array([step(2-i*2*pi)*sqrt(pi/2) for i in f]))


fig,ax = plt.subplots()
x = np.linspace(-5,5,100)
ax.plot(x,np.array([step(i) for i in x]))
ax.plot(x,np.array([step(2-i) for i in x]))



fig,ax = plt.subplots()
ax.plot(alpha,np.array(list(map(Ex1_FT,alpha))).imag)

"""
fig,ax = plt.subplots(2,2)
ax = ax.reshape(4)
k = np.linspace(-10,10,100)
[i.plot(k,list(map(j,k))) for i,j in zip(ax,[Jz1_FT,Jz2_FT,Jx1_FT,Jx2_FT])]
"""



from scipy.fft import fft, ifft, fftfreq, fftshift, ifftshift

#signal = np.array([-2, 8, 6, 4, 1, 0, 3, 5], dtype=float)
s1 = lambda t: sin(30*t) + sin(50*t)
n = 30
t,d = np.linspace(-1,1,n,retstep = True)
fourier = fft(s1(t))



freqs = fftfreq(n, d=d)
f = fftshift(freqs)
fig,ax = plt.subplots(2,1)
ax[0].plot(t,s1(t))
ax[1].plot(f,fourier)
#freqs = np.fft.fftfreq(10, 0.1)
#f = np.fft.fftshift(freqs)



a= lambda t:fft(sin(t))

t = np.linspace(-pi,pi,100)
fig,ax = plt.subplots(2,1)
A = a(t)
ax[0].plot(t,sin(t))
ax[1].plot(A.imag)

Jzm = lambda x,z,w,l,r,s:cos((r-1)*pi*2*x/w)/sqrt((w/2)**2 - x**2) \
    * cos((2*s-1)*pi*z/l)/sqrt((l/2)**2 - z**2) if (abs(x)<w/2 and abs(z)<l/2) else 0 #a function of x and z
Jxm = lambda x,z,w,l,r,s:sin(2*r*pi*x/w)/sqrt((w/2)**2 - x**2)\
    *sin((2*s-1)*pi*z/l)/sqrt((l/z)**2 - z**2) if (abs(x)<w/2 and abs(z)<l/2) else 0
w = 3.6e-3
l = 1.8e-3
Jz1_patch = lambda x,z:Jzm(x,z,w,l,1,1)
Jz2_patch = lambda x,z:Jzm(x,z,w,l,1,2)

Jx1_patch = lambda x,z:Jxm(x,z,w,l,1,1)
Jx2_patch = lambda x,z:Jxm(x,z,w,l,2,1)

def calc_2D(x,y,f): #(f(x,y))
    return np.array([[f(i,j) for i in x] for j in y])

x = np.linspace(-4e-3,4e-3,100)
z = np.linspace(-3e-3,3e-3,100)
fig,ax = plt.subplots(2,1)
#ax.plot
from matplotlib import cm
plot_3D(x,z,calc_2D(x,z,Jz1_patch),plot_type = '3D')
plot_3D(x,z,calc_2D(x,z,Jz2_patch),plot_type = '3D')
plot_3D(x,z,calc_2D(x,z,Jx1_patch),plot_type = '3D')
plot_3D(x,z,calc_2D(x,z,Jx2_patch),plot_type = '3D')


#fourier lambda alpha,beta,fun: integrate(lambda:quad*exp(1j*(alpha*x+beta*y)),)


def fourier2D(alpha,beta,fun):
    return integrate.dblquad(lambda x, y: fun(x,y)*np.exp(1j*(alpha*x+beta*y)), -np.inf, np.inf, lambda x: -np.inf, lambda x: np.inf)[0]

rect = lambda x,y:1 if(abs(x)<1 and abs(y)<1) else 0
x = np.linspace(-2,2,20)
y = np.linspace(-2,2,20)
Z = np.array([[rect(i,j) for i in x] for j in y])
plot_3D(x,y,Z)
alpha = 1/x
beta = 1/y

Z_f = np.array([[fourier2D(i,j,rect) for i in alpha] for j in beta])
plot_3D(alpha,beta,Z_f.real,plot_type = '3D')
#fourier2D()


def fft_2D(x,y,fun):
    
    Y = np.meshgrid(x,y)
    fft2()

#fullwave analysis microstrip resonator and radiation using spectral domain immitance matrix approach 
def test(alpha,beta,k):
    
    
    
    omega = 2*pi*f0
    er = 3.1 

    e1 = er
    e2 = 1
    k = omega/c0
    gamma1 = sqrt(alpha**2 + beta**2 - e1*k**2)
    gamma2 = sqrt(alpha**2 + beta**2 - e2*k**2)

    Ytm1 = 1j*omega*e0*e1/gamma1 #may be a function of gamma
    Ytm2 = 1j*omega*e0*e2/gamma2
    Yte1 = gamma1/(1j*omega*u0)
    Yte2 = gamma2/(1j*omega*u0)
    
    
    Ye_p = Ytm1*coth(gamma1*d) #Ye+
    Yh_p = Yte1*coth(gamma1*d) #Yh+
    
    Ye_m = Ytm2
    Yh_m = Yte2


    Z0e = 1/(Ye_p + Ye_m)
    Z0h = 1/(Yh_p + Yh_m)
    
    Nx = alpha/sqrt(alpha**2 + beta**2)
    Nz = beta/sqrt(alpha**2 + beta**2)
    Zxx = Nx**2*Z0e +Nz**2*Z0h
    Zxz = Nx*Nz*(-Z0e + Z0h)
    Zzz = Nz**2*Z0e + Nx**2*Z0h
    
    Jy
w = np.linspace(-0.5,0.5,100)
fig,ax = plt.subplots(2,2)
ax[0,0].plot(w,np.array(list(map(Jz1,w))))
ax[0,1].plot(w,np.array(list(map(Jz2,w))))
ax[1,0].plot(w,np.array(list(map(Jx1,w))))
ax[1,1].plot(w,np.array(list(map(Jx2,w))))




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



"""
Fourier transform for even and odd function
"""
def fourier(fun,f):
    operator_e = lambda f,t: cos(2*pi*f*t)
    operator_o = lambda f,t: sin(2*pi*f*t)

    even = integrate.quad(lambda t:fun(t)*operator_e(f,t),-np.inf,np.inf)[0]
    odd = integrate.quad(lambda t:fun(t)*operator_o(f,t),-np.inf,np.inf)[0]
    return even + 1j*odd

f1 = lambda x,b: 1 if (abs(x)<b) else 0
x = np.linspace(-5,5,100)
b = 1
ans = np.array([f1(i,b) for i in x])

freq = np.linspace(-3,3,100)
ans2 = np.array([fourier(lambda t:f1(t,b),i) for i in freq])
f1a = lambda w,b: sqrt(2/pi)*sin(b*w)/w


fig,ax = plt.subplots(2,1)
ax[0].plot(x,ans)
ax[1].plot(freq,ans2.real,label = 'real')
ax[1].plot(freq,ans2.imag,label = 'imag')
ax[1].plot(freq,list(map(lambda t:f1a(t,b),freq*2*pi)),label = 'analytical',linestyle = 'dashed')
ax[1].legend()

f2 = lambda x,b: x if(abs(x)<b) else 0
x = np.linspace(-5,5,100)
b = 2

F2 = np.array([fourier(lambda t:f2(t,b),i) for i in freq])
fig,ax = plt.subplots(2,1)
ax[0].plot(x,list(map(lambda x:f2(x,b),x)))
ax[1].plot(freq,F2.real,label = 'real')
ax[1].plot(freq,F2.imag,label = 'imag')
ax[1].legend()





"""
Slot antenna analysis using spectral domain method of immitance approach
"""

Ex1 = lambda x,z: cos(2*(r-1)*pi*x/w)/sqrt((w/2)**2 - x**2) *\
    cos(2*(s-1)*pi*z/l)/sqrt((l/2)**2 - z**2) if (abs(x)<w/2 and abs(z)<l/2) else 0 
Ez1 = lambda x,z: sin(2*r*pi*x/w)/sqrt((w/2)**2 - x**2) *\
    sin((2*x-1)*pi*z/l)/sqrt((l/2)**2 - z**2)
    



    
"""
Laplace transformation
"""
from scipy import integrate
def laplace(fun,s):
    operator = lambda s,t:np.exp(-s*t)
    calc = list(map(lambda s:integrate.quad(lambda t:fun(t)*\
                            operator(s,t),0,np.inf)[0],s))   
    return np.array(calc)

f1 = lambda t:1
s1 = np.linspace(1,10,10)
F1 = laplace(f1,s1)


f2 = lambda t,a:np.exp(a*t)
a = 0.1
s2 = np.linspace(a+1,10+a+1,10)
F2 = laplace(lambda t:f2(a,t),s2)



#F2 = laplace(lambda t:f2(t,a),s2)


F1 = integrate.quad(lambda t:f1(t)*operator(s,t),0,np.inf)[0]   