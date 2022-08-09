# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 14:39:40 2022

@author: hiro7
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from numpy import sqrt,sin,cos,pi
from matplotlib import cm
dB = lambda x:10*np.log10(x)
from scipy import integrate

import mpl_toolkits.mplot3d.axes3d as axes3d

c0 = 2.9979e8

E0 = 1
f = 30e9
k0 = 2*pi*f/c0

def find_nearest(a, a0):
    "Element in nd array `a` closest to the scalar value `a0`"
    idx = np.abs(a - a0).argmin()
    return idx,a.flat[idx]

a_test = np.random.randn(10)
#D(theta,phi)
def plot_antenna(D,var,hold = 'phi',ax = None):
    fig,ax = plt.subplots(2,1)
    if hold == 'phi':
        ax[0].plot(theta*180/pi,D[:,find_nearest(phi,var[0])[0]])
        ax[1].plot(theta*180/pi,D[:,find_nearest(phi,var[1])[0]])
    else: #hold is theta
        ax[0].plot(phi*180/pi,D[find_nearest(phi,var[0])[0],:])
        ax[1].plot(phi*180/pi,D[find_nearest(phi,var[1])[0],:])

    #ax[0].plot(theta,D[:,find_nearest(phi,var)[0]])
    return fig,ax 
    #ax[0,0].plot
"""
def plot_polar(D,var,hold = 'phi',ax = None):
    fig,ax = plt.subplots(2,1)
    if hold == 'phi':
        ax[0].plot(theta*180/pi,D[:,find_nearest(phi,var[0])[0]])
        ax[1].plot(theta*180/pi,D[:,find_nearest(phi,var[1])[0]])
    else: #hold is theta
        ax[0].plot(phi*180/pi,D[find_nearest(phi,var[0])[0],:])
        ax[1].plot(phi*180/pi,D[find_nearest(phi,var[1])[0],:])

    #ax[0].plot(theta,D[:,find_nearest(phi,var)[0]])
    return fig,ax 
    #ax[0,0].plot
"""



def plot_3D(x,y,Z,plot_type = 'heat'):
    if plot_type == 'heat':
        fig, ax = plt.subplots()
        ax.imshow(Z, cmap='hot', interpolation='nearest')
    elif plot_type == 'Polar': #x = theta, y = phi,

        fig = plt.figure()
        THETA, PHI = np.meshgrid(theta, phi)
        #R = np.cos(PHI**2)
        X =  np.sin(PHI) * np.cos(THETA)# * R
        Y =  np.sin(PHI) * np.sin(THETA)# * R
        #Z = R * np.cos(PHI)

        ax = fig.add_subplot(1,1,1, projection='3d')
        plot = ax.plot_surface(
            X, Y, Z, rstride=1, cstride=1, cmap=plt.get_cmap('jet'),
            linewidth=0, antialiased=False, alpha=0.5)
            
    else:
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        X,Y = np.meshgrid(x,y)
        ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    return fig,ax





#matrice = lambda x:[x]
xi = 377

#E_theta_aperture = lambda theta,phi,a,b:1j*a*b*k0*E0/(2*pi)*(sin(phi)*(sin(X)/X)*(sin(Y)/Y))#*exp(-1j*k*r)/r
#E_phi_aperture = lambda theta,phi,a,b:1j*a*b*k0*E0/(2*pi)*(cos(theta)*cos(phi)*(sin(X)/X)*(sin(Y)/Y))#*exp(-1j*k*r)/r

def test_fun(a,b):
    return a+b,a*b

def aperture(theta,phi,a,b): #analytical solution
    a = a*c0/f #q wavelength
    b = b*c0/f
    xi = 377
    X = k0*a/2*sin(theta)*cos(phi)
    Y = k0*b/2*sin(theta)*sin(phi)
    E_theta = 1j*a*b*k0*E0/(2*pi)*(sin(phi)*(sin(X)/X)*(sin(Y)/Y))#*exp(-1j*k*r)/r
    E_phi = 1j*a*b*k0*E0/(2*pi)*(cos(theta)*cos(phi)*(sin(X)/X)*(sin(Y)/Y))#*exp(-1j*k*r)/r
    H_theta = -E_phi/xi
    H_phi = E_theta/xi
    U = 1/2/xi*(np.abs(E_theta)**2 + np.abs(E_phi)**2)
    return E_theta,E_phi,H_theta,H_phi,U
#class Antenna(object):
#    def __init__(self,phi,theta):
 
#Aperture
theta = np.linspace(-pi/2,pi/2,100)
phi = np.linspace(0.01,pi/2,100)

       
#E_theta,E_phi,H_theta,H_phi = np.array([[aperture(i,j) for i in theta] for j in phi]).transpose(2,0,1)
a = 3
b = 2
E_thetaa,E_phia,H_thetaa,H_phia,U_a = np.array([[aperture(i,j,a,b) for i in theta] for j in phi]).transpose(2,0,1)
Prada = integrate.dblquad(lambda theta,phi:aperture(theta,phi,a,b)[4]*sin(theta),0,pi,lambda theta:0.00,lambda theta:pi)[0]
D = U_a/Prada
D0 = dB(4*pi*np.max(U_a)/Prada) 
#D = np.array([[directivity(E_theta[i[0],j[0]],E_phi[i[0],j[0]],H_theta[i[0],j[0]],H_phi[i[0],j[0]]) for i in enumerate(theta)] for j in enumerate(phi)])
D_normalize = D/np.max(D)

print(f'D0 = {D0}')

plot_3D(theta,phi,dB(np.abs(D_normalize)),plot_type = 'heat')
fig,ax = plot_antenna(dB(D_normalize),var=[0,pi/2],hold = 'phi')
ax[0].set_xlim(-90,90)
ax[1].set_xlim(-90,90)
#ax[0].set_ylim(-40,0)
#ax[1].set_ylim(-40,0)





"""
#Balanis p.43 correct.
U_test = lambda theta,phi: sin(theta)
U = np.array([[U_test(i,j) for i in np.linspace(0,pi,100)] for j in np.linspace(0,2*pi,100)])
Prada_test = integrate.dblquad(lambda theta,phi:U_test(theta,phi)*sin(theta),0,2*pi,lambda theta:0,lambda theta:pi)[0]

D0_test = dB(4*pi*np.max(U)/Prada_test)
"""
#TypeError: can't convert complex to float
def patch(theta,phi,W,L,h,f0):
    k0 = 2*pi*f0/c0
    lambda0 = c0/f0
    #W = W * f/c0 #a wavelength
    #L = L * f/c0
    #h = h * f/c0
    xi = 120*pi
    X = k0*(h*lambda0)/2*sin(theta)*cos(phi)
    Z = k0*(W*lambda0)/2*cos(theta)
    E_theta = 0
    #E_theta = 1j*a*b*k*E0/(2*pi)*(sin(phi)*(sin(X)/X)*(sin(Y)/Y))#*exp(-1j*k*r)/r
    E_phi = 1j*k0*(h*lambda0)*(W*lambda0)*E0/pi*(sin(theta)*sin(X)/X*sin(Z)/Z)# AF->*cos(k0*(L*lambda0)/2*sin(theta)*sin(phi))
    #E_phi = 1j*a*b*k*E0/(2*pi)*(cos(theta)*cos(phi)*(sin(X)/X)*(sin(Y)/Y))#*exp(-1j*k*r)/r
    H_theta = -E_phi/xi
    H_phi = E_theta/xi
    U = 1/2/xi*(np.abs(E_theta)**2 + np.abs(E_phi)**2)
    return E_theta,E_phi,H_theta,H_phi,U
#patch
theta = np.linspace(0.001,pi,100)
phi = np.linspace(0.001,pi/2,100)
f0 = 10e9
k0 = 2*pi*f0/c0
lambda0 = c0/f0
E0 = 1
W = 11.86e-3/lambda0 #
L = 9.06e-3/lambda0
h = 1.588e-3/lambda0
Er = 2.2
k1 = k0
E_theta,E_phi,H_theta,H_phi,U_patch = np.array([[patch(i,j,W,L,h,f0) for i in theta] for j in phi]).transpose(2,0,1)
Prad = lambda theta: (sin(k1*W*lambda0/2*cos(theta))/cos(theta))**2*sin(theta)**3
Prad_patch_tb = integrate.quad(Prad,0,pi)
#G1 = Prad_patch_tb[0]/(pi*377)
G1 = Prad_patch_tb[0]/(120*pi**2) # approx 1/120*(W*lambda0/lambda0)



Prad_patch = integrate.dblquad(lambda theta,phi:patch(theta,phi,W,L,h,f0)[4]*sin(theta),0,2*pi,lambda theta:0.00,lambda theta:pi)[0]
D = U_patch/Prad_patch
D0 = dB(4*pi*np.max(U_patch)/Prad_patch) 
#D = np.array([[directivity(E_theta[i[0],j[0]],E_phi[i[0],j[0]],H_theta[i[0],j[0]],H_phi[i[0],j[0]]) for i in enumerate(theta)] for j in enumerate(phi)])
D_normalize = D/np.max(D)
print(f'D0 = {D0}')



plot_3D(theta,phi,dB(np.abs(D_normalize)),plot_type = 'heat')
plot_antenna(D_normalize,[pi/2,pi/4] ,hold = 'phi')





fig, ax = plt.subplots(2,1,subplot_kw={'projection': 'polar'})
#fig, ax = plt.subplots(2,1)

ax[0].plot(phi,dB(D_normalize[50,:]),label = 'E-plane')
ax[1].plot(theta, dB(D_normalize[:,50]),label = 'H-plane')

#ax.plot(theta, dB(D_normalize[:,0]),'b-')
#ax.plot(theta, dB(D_normalize[:,-1]),'r-')
ax[0].set_rmax(0)
ax[0].set_rmin(-40.0)
#ax.set_theta_zero_location("W")
ax[0].set_theta_offset(pi/2)
ax[1].set_rmax(0)
ax[1].set_rmin(-40.0)
#ax.set_theta_zero_location("W")
#ax[1].set_theta_offset(pi/2)

#Exact definition of directivity. However, directivity in far-field is already implemented inside a function
"""
def directivity(E_theta,E_phi,H_theta,H_phi):
    a_theta = np.array([1,0])
    a_phi = np.array([0,1])
    U = 1/2*np.cross(a_theta*E_theta + a_phi*E_phi,(a_theta*H_theta + a_phi*H_phi).conj().T)
    return U
"""


#example of integration
"""
from scipy import integrate

f = lambda y, x, c : x*y**2 + c

t = integrate.dblquad(f, 0, 2, lambda x: 0, lambda x: 1,args = [1])
#    (0.6666666666666667, 7.401486830834377e-15)
invexp = lambda x: np.exp(-x)

integrate.quad(invexp, 0, np.inf)
#(1.0, 5.842605999138044e-11)



integrate.quad(lambda x,c: test_fun(x,c)[0],1,5,args =(1))
f_ = lambda theta,phi: sin(theta)
test1 = integrate.dblquad(f_,0,2*pi,lambda theta: 0, lambda theta: pi/6)
test2 = integrate.quad(f_,0,pi/6,args = [pi/2])
print(f'test1 = {test1[0]}, test2 = {test2[0]}')


"""



#axisb: axis -1 is out of bounds for array of dimension 0
D = np.array([[directivity(E_theta[i[0],j[0]],E_phi[i[0],j[0]],H_theta[i[0],j[0]],H_phi[i[0],j[0]]) for i in enumerate(theta)] for j in enumerate(phi)])
D_normalize = D/np.max(D)
E_mag = np.array([[sqrt(np.abs(np.array([1,0])*E_theta[i[0],j[0]])**2 + np.abs(np.array([0,1])*E_phi[i[0],j[0]])**2) for i in enumerate(theta)] for j in enumerate(phi)])
"""
#x = fig.add_subplot(projection='3d')

# Create the mesh in polar coordinates and compute corresponding Z.
theta, phi = np.linspace(0, 2 * np.pi, 40), np.linspace(0, np.pi, 40)
THETA, PHI = np.meshgrid(theta, phi)
R = np.cos(PHI**2)
X = R * np.sin(PHI) * np.cos(THETA)
Y = R * np.sin(PHI) * np.sin(THETA)
Z = R * np.cos(PHI)
fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection='3d')
plot = ax.plot_surface(
    X, Y, Z, rstride=1, cstride=1, cmap=plt.get_cmap('jet'),
    linewidth=0, antialiased=False, alpha=0.5)
"""


"""
fig,ax = plt.subplots()
ax.plot(theta*180/pi,dB(D_normalize[:,0]),label = r'$\phi$ = 0')
#ax.plot(theta,D_normalize[:,-1])
#ax.plot(theta*180/pi,dB(D_normalize[:,50]),label = r'$\phi$ = $\pi$/4')
#ax.plot(theta*180/pi,dB(D_normalize[:,25]),label = r'$\phi$ = $\pi$/6')

ax.plot(theta*180/pi,dB(D_normalize[:,-1]),label = r'$\phi$ = $\pi$/2')

ax.set_xlabel('theta')
ax.set_ylabel('Directivity')
ax.set_ylim(-40,0)
ax.legend()
"""


"""
Green function circular disk capacitance
"""
from scipy.integrate import quad
e0 = 8.85e-12
d = 1e-3
er = 3
Green = lambda alpha:1/(e0*alpha*(1 + er*coth(alpha*d)))
Rho1 = lambda alpha,a: a*jv(1,alpha*a)/alpha
Phi = lambda alpha,a: a/alpha*jv(1,alpha*a)

a = 1e-3
k11 = lambda alpha,a: Rho1(alpha,a)*Green(alpha)*Rho1(alpha,a)
K11 = quad(lambda k:k11(k,a),0,np.inf)
a1 = quad(lambda k:Rho1(k,a)*Phi(k,a),0,np.inf)

C = 2*pi*a1[0]**2/K11[0]


"""
Via barrel-plate capacitance from TMz0 calculation
"""

from scipy.special import gamma,hankel1,hankel2,jv


def Cb(a,b,R,h,er,f,N): #N is list
    k0 = 2*pi*f/2.9979e8
    kn = lambda n:sqrt(0j+(k0*er)**2 - (n*pi/h)**2)
    e = er*e0

    #gamma_Rn = -hankel2(0,kn*R)/jv(0,kn*R) #if PEC
    #gamma_Rn = -hankel2(1,kn*R)/jv(1,kn*R) #if PMC
    gamma_Rn = 0
    
    #gamma_an = -jv(0,kn*a)/hankel2(0,kn*a) #if PEC
    #gamma_an = -jv(1,kn*a)/hankel2(1,kn*a) #if PMC
    gamma_an = 0

    def Cb_b(n):
        
        cb_n = (1 - gamma_Rn*gamma_an)**(-1)/(kn(n)**2*hankel2(0,kn(n)*a))*((hankel2(0,kn(n)*b) \
                                                    - hankel2(0,kn(n)*a)) + gamma_Rn*(jv(0,kn(n)*b) - jv(0,kn(n)*a)))
        return cb_n    
    cb = 8*pi*e/h/log(b/a)*np.array([Cb_b(i) for i in N])
    cb_sum = cb.sum()
    return cb_sum,cb
R =  1e-3
a,b = [0.1e-3,0.35e-3]
er = 3.84
e0 = 8.85e-12
h = 0.228e-3
f = 1e9
#Cb = 8*pi*e/h/log(b/a)*(1 - gamma_Rn*gamma_an)**(-1)/(kn**2*hankel2(0,kn*a))*((hankel2(0,kn*b) - hankel2(0,kn*a)) + gamma_Rn*(jv(0,kn*b) - jv(0,kn*a)))

Cb(a,b,R,h,er,1e9,[1,3,5,7,9])

print(f'Cb,complex:{Cb},abs:{np.abs(Cb)}')
l = 0.025e-3
Ca = 2*pi*e*l/log(b/a)
