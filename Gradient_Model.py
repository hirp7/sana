import numpy as np
from scipy.special import hyp2f1 #arguments: a,b,c, and z
from mpmath import hyp3f2,hyp2f1
#mpmath.hyp3f2(a1, a2, a3, b1, b2, z)
import matplotlib.pyplot as plt
from skrf.media import MLine, DefinedGammaZ0
import skrf as rf

import os
path = r"D:\ATS\takahashi\Desktop\Analysis\Roughness"
os.chdir(path)


inv = np.linalg.inv
pi = np.pi 
sqrt = np.sqrt
exp = np.exp
log10 = np.log10
log = np.log
tanh = np.tanh
atanh = np.arctanh
u0 = 4*pi*1e-7
e0 = 8.85*1e-12
gradient = np.gradient


Sigma0 = 5.8e7
"""
Gradient
Analytical Model using closed form
"""

#Gradient Model function input Rq. Frequency

#from mpmath import hyp2f1
from scipy.special import hyp2f1 #arguments: a,b,c, and z


Mag_noise,Phase_noise = [0.001,0.01]
#Standard means that the RMS of the copper roughness is 1 um,while smooth is smooth 
Waveguide10mm = rf.Network('StandardWaveguide_10mm.s2p')
Waveguide10mm = rf.Network('SmoothWaveguide_10mm.s2p')

Waveguide10mm.add_noise_polar(Mag_noise,Phase_noise)
Waveguide5mm = rf.Network('StandardWaveguide_5mm.s2p')
Waveguide5mm = rf.Network('SmoothWaveguide_5mm.s2p')

Waveguide5mm.add_noise_polar(Mag_noise,Phase_noise)

Waveguide_Loss = np.abs(Waveguide10mm.s21.s[:,0,0])- \
                np.abs(Waveguide5mm.s21.s[:,0,0])
                
#curve fittingh
fig,ax = plt.subplots()
Waveguide10mm.s21.plot_s_db(ax = ax)
Waveguide5mm.s21.plot_s_db(ax = ax)

#hyp3f2(a1, a2, a3, b1, b2, z)
calc_F = lambda i,alpha,beta,ksi:hyp3f2(1+alpha-beta,2+alpha+beta,i+alpha,1+2*alpha,1+i+alpha,ksi) 
class GradientG():
    """
    Reference: Grujuc 2018 'Closed Form Solution of Rough Conductor Surface Impedance'
    """
    def __init__(self,Rq,freq):
        psi = 1/2
        sigma0 = 5.8e7
        u0 = 4*pi*1e-7
        self.sigma0 = sigma0
        self.f = freq
        self.Rq = Rq
        self.omega = 2*pi*self.f
        self.x,self.x_step = np.linspace(-5*Rq,5*Rq,100,retstep = True)
        self.chi = sqrt(2)*Rq
        self.sigma = sigma0/4*(1 + tanh(self.x/self.chi + psi))**2
        """
        self.B = self.ksi**self.alpha*hyp2f1(self.alpha.real + \
            self.beta.real, self.alpha.real - self.beta.real - \
                1, 1 + 2*self.alpha.real,self.ksi)
        """
        self.B = np.array([calc_By(i,self.Rq,self.f) for i in self.x])
        #self.Jz = gradient(self.B,self.x_step)/u0
        #self.P_diss = np.abs(self.Jz)**2/(2*self.sigma)
        self.Z_surf = calc_Zsurf(self.Rq,self.f)
        self.sigma_eff = u0*self.omega/(2*self.Z_surf.real**2)
        self.u_eff = 2*sigma0*self.Z_surf.imag**2/self.omega
        
def calc_By(x,Rq,f):
    omega = f
    chi = sqrt(2)*Rq
    psi = 1/2

    alpha = (1+1j)/2*Rq*sqrt(u0*omega*Sigma0)
    beta = 1/2*(sqrt(1 + 4*alpha**2) - 1)
    ksi = (1 + np.e**(2*(x/chi + psi)))**-1
    B = ksi**alpha*hyp2f1(alpha.real + \
        beta.real, alpha.real - beta.real - \
            1, 1 + 2*alpha.real,ksi)
    
    return np.complex(B)

def calc_Zsurf(Rq,f):
    omega = 2*pi*f
    chi = sqrt(2)*Rq
    psi = 1/2
    x0 = -5*Rq
    
    alpha = (1+1j)/2*Rq*sqrt(u0*omega*Sigma0)
    beta = 1/2*(sqrt(1 + 4*alpha**2) - 1)
    ksi = lambda x:(1 + np.e**(2*(x/chi + psi)))**-1
    Bg = lambda x:chi/2*ksi(x)**alpha*(ksi(x)/(1 + alpha)*calc_F(1,alpha,beta,ksi(x)) - 1/alpha*calc_F(0,alpha,beta,ksi(x)))    
    #calc_By(x0,Rq,f)
    Zsurf = -1j*omega*u0*Bg(x0)/calc_By(x0,Rq,f)
    
    #Zsurf = 1j*u0*omega*integrate.quad(lambda x:calc_By(x,Rq,f).imag,x0,np.inf)[0]/calc_By(x0,Rq,f)
    return np.complex(Zsurf)
            
def calc_RoughParamf(Rq,Frequency): #Z_surf, sigma_eff, ur_eff
    Zsurf = []
    sigma_eff = []
    ur_eff = []
    for i in Frequency.f:
        x = GradientG(Rq,i)
        Zsurf.append(x.Z_surf)
        sigma_eff.append(x.sigma_eff)
        ur_eff.append(x.u_eff)
    
    Zsurf,sigma_eff,ur_eff = np.array([np.array(i) for i in [Zsurf,sigma_eff,ur_eff]])
    return Zsurf, np.abs(sigma_eff),np.abs(ur_eff)


"""
for discrete frequency analysis
"""
u0 = 4*pi*1e-7
Rq = 1e-6
depth = np.linspace(-5*Rq,5*Rq,100)

frequencies = [1e9,10e9,100e9]
GradientG1,GradientG10,GradientG100 = [GradientG(Rq,i) for i in frequencies]
Gradient_set = [GradientG1,GradientG10,GradientG100]
#fig_B = plt.figure(figsize = (13,7))
fig_B,[[ax_B,ax_Jz],[ax_Breal,ax_Bimag]] = plt.subplots(2,2) #ax_P = ax_Jz.twinx()
ax_set = [ax_B,ax_Jz,ax_Breal,ax_Bimag]
set_xlabels = 'x [\u03BC m]'
set_ylabels = ['|B|','Current Dens. Jz','Re(B)','Im(B)']#'Dissipated Power']
color = ['orange','blue','green']
for j,i in enumerate(Gradient_set):
    ax_B.plot(depth*1e6,np.abs(i.B)/np.abs(i.B[0]),color = color[j],label = str(frequencies[j]/1e9) + ' [GHz]')
    ax_Breal.plot(depth*1e6,(i.B/i.B[0]).real,color = color[j],label = str(frequencies[j]/1e9) + ' [GHz]')
    ax_Bimag.plot(depth*1e6,(i.B/i.B[0]).imag,color = color[j],label = str(frequencies[j]/1e9) + ' [GHz]')
    ax_Jz.plot(depth*1e6,np.abs(i.Jz)/np.abs(np.max(i.Jz)),color = color[j],label = str(frequencies[j]/1e9) + ' [GHz]')
    #ax_P.plot(depth*1e6,np.abs(i.P_diss),color = color[j],linestyle = '--')
[i.legend(loc = 'upper right',fontsize = 6) for i in ax_set]
[[i.set_xlabel(set_xlabels),i.set_ylabel(set_ylabels[j])] for j,i in enumerate(ax_set) ]
fig_B.suptitle(f'Conductor Surface Roughness with RMS = {Rq*1e6} [\u03BCm]',fontsize = 12,y = 0.99)
fig_B.tight_layout(pad = 1.2)



"""
Surface Impedance with discrete Rq
"""

#calc_RoughParamf(Rq,Frequency) #-> Zsurf,sigma_eff,ur_eff


Rq = [1e-6,0.5e-6,0.25e-6]
u0 = 4*pi*1e-7
Frequency = rf.Frequency(0.1,100,100)
Zsurf,Sigma_eff,ur_eff = np.array([calc_RoughParamf(i,Frequency) for i in Rq]).transpose(1,0,2)
fig_Zsurf = plt.figure(figsize = (6,7))
fig,ax = plt.subplots(2,1,sharex = True)
ax02 = ax[0].twinx()
ax12 = ax[1].twinx()

ax_sets = [ax[0],ax02,ax[1],ax12]
ylabels = ['Re(Zsurf) ','Im(Zsurf)','\u03C3_eff','\u03BCr_eff']
#[ax[0].plot(Frequency.f/1e9,i,label = f'{np.round(Rq[j]*1e6,2)} \u03BC m') for j,i in enumerate(Zsurf.real)]
[ax[0].plot(Frequency.f/1e9,i) for j,i in enumerate(np.array(Zsurf.real))]

#[ax[0].plot(Frequency.f/1e9,i) for j,i in Zsurf.real]

[ax02.plot(Frequency.f/1e9,i,linestyle = '--') for i in Zsurf.imag]

[ax[1].plot(Frequency.f/1e9,i) for i in Sigma_eff]
[ax12.plot(Frequency.f/1e9,i/u0,linestyle = '--') for i in ur_eff]
ax[1].set_xlabel('Frequency [GHz]')
[i.set_ylabel(j) for i,j in zip(ax_sets,ylabels)]
#ax[0].legend(['1 \u03BC m','0.5 \u03BC m','0.25 \u03BC'])
ax[0].legend()




def f(x):
    return 0*x



def FDM(x0,xn,n,a,b,c,f,alpha,beta): #Dirichlet 
    inv = np.linalg.inv
    x = np.linspace(x0,xn,n,retstep = True)
    h = x[1]
    x = x[0]
    #define A
    diagonal = h**2*c - 2*a
    diagonal_side = a[1:] + h*b[1:]/2
    A = np.diag(diagonal) + np.diag(diagonal_side,1) + np.diag(diagonal_side,-1)   
    A = 1/h**2*A

    #define F
    F = f(x)
    F[0] = F[0] - (a[0]/h**2 - b[0]/2/h)*alpha
    F[n-1] = F[n-1] - (a[n-1]/h**2 + b[n-1]/2/h)*beta
    U = inv(A)@F
    return U,A,F,x

from scipy.stats import norm
"""
fig,ax = plt.subplots()
Rq = 1e-6
x = np.linspace(-5*Rq,5*Rq,100)
ax.plot(x,[norm.cdf(i,scale = Rq) for i in x])
"""
CDF_norm = lambda x,Rq:norm.cdf(x,scale = Rq)
def Gradient_Model(frequency,Rq):
    Sigma0 = 5.8e7
    Omega = 2*pi*frequency
    n = 1000
    xmin =-Rq*5 #should be determined by frequency and Rq
    #x0
    xn = 10e-6
    x = np.linspace(xmin,xn,n)
    Sigma_x = Sigma0*CDF_norm(x,Rq)
    a = np.ones(n)
    b = -gradient(Sigma_x)
    c = -1j*Omega*u0*Sigma_x**2
    #c = 0
    U,A,F,x = FDM(xmin,xn,n,a,b,c,f,1,0)
    return U,x,Sigma_x


U,x,Sigma_x = Gradient_Model(1e9,1e-6)

from scipy.special import jv

FDM(0,10,)



"""
#np.insert(A,index,value)
def FDM_N(x0,xn,n,a,b,c,f,alpha_dash,beta): #Neuman
    inv = np.linalg.inv
    x = np.linspace(x0,xn,n,retstep = True)
    h = x[1]
    x = x[0]
    #define A [mx2]^2
    diagonal = h**2*c - 2*a
    diagonal = np.insert(diagonal,0,-h)
    diagonal = np.append(diagonal,h**2)
    
    diagonal_left = a[:] + h*b[:]/2
    #diagonal_left = np.insert(diagonal_side,0,a[0] - h*b[0]/2)
    diagonal_left = np.append(diagonal_left,0)
    #diagonal_left = np.insert(diagonal_left,0,a[0] - h*b[0]/2)

    diagonal_right = a[:] + h*b[:]/2
    diagonal_right = np.insert(diagonal_right,0,h)
    #diagonal_left = np.insert(diagonal_side,0,a[0] - h*b[0]/2)
    #diagonal_right = np.append(diagonal_side,0
    A = np.diag(diagonal) + np.diag(diagonal_right,1) + np.diag(diagonal_left,-1)   
    A = 1/h**2*A

    #define F
    F = f(x)
    F = np.insert(F,0,alpha_dash)
    #F[0] = F[0] - (a[0]/h**2 - b[0]/2/h)*alpha
    #F[n-1] = F[n-1] - (a[n-1]/h**2 + b[n-1]/2/h)*beta
    F = np.append(F,beta)
    U = inv(A)@F
    U = U[1:-1]
    return U,A,F,x

"""
"""
Gradient Model with surface roughness copper (analytical solution unavailable)
"""
"""
def Gradient_Model(frequency,Rq):
    Sigma0 = 5.8e7
    
    Omega = 2*pi*frequency
    n = 1000
    xmin =-Rq*5 #should be determined by frequency and Rq
    
    #x0
    xn = 10e-6
    x = np.linspace(xmin,xn,n)
    Sigma_x = Sigma0*CDF_norm(x,Rq)
    a = np.ones(n)
    b = -gradient(Sigma_x)*0
    #b = np.ones(n)*0
    c = -1j*Omega*u0*Sigma_x**2*0
    #c = 0
    U,A,F,x = FDM(xmin,xn,n,a,b,c,f,1,0)
    return U,x,Sigma_x

def Gradient_Impedance_Model(frequency,Rq):
    Sigma0 = 5.8e7
    
    Omega = 2*pi*frequency
    n = 1000
    xmin =-Rq*5 #should be determined by frequency and Rq
    xn = 10e-6
    x = np.linspace(xmin,xn,n)
    Sigma_x = Sigma0*CDF_norm(x,Rq)
    Sigma_x = Sigma0*CDF_norm(x,Rq)
    a = np.zeros(n)
    b = np.ones(n)
    c = Sigma_x
    return

def Gradient_Model_FS(frequency,Rq): #free space
    Sigma0 = 5.8e7/5.8e7
    
    Omega = 2*pi*frequency
    n = 1000
    xmin = -10e-6 #should be determined by frequency and Rq
 
    #x0
    xn = 0
    x = np.linspace(xmin,xn,n)
    Sigma_x = Sigma0*CDF_norm(x,Rq)
    k = Omega*sqrt(u0*e0)*sqrt(1-1j*Sigma_x/Omega/e0)
    a = np.ones(n)*1
    b = np.ones(n)*0
    #b = -gradient(log(Sigma_x),10e-6)
    #c = 1j*Omega*u0*Sigma_x
    c = np.ones(n)*k**2
    #c = 0
    U,A,F,x = FDM_N2(xmin,xn,n,a,b,c,f,0,1)
    return U,x,Sigma_x


def Gradient_Model_N(frequency,Rq):
    Sigma0 = 5.8e7/5.8e7
    
    Omega = 2*pi*frequency
    n = 1000
    xmin = 0 #should be determined by frequency and Rq
    
    #x0
    xn = 10e-6
    x = np.linspace(xmin,xn,n)
    Sigma_x = Sigma0*CDF_norm(x,Rq)
    a = np.ones(n)*1
    b = -gradient(log(Sigma_x),0.1e-6)
    c = -1j*Omega*u0*Sigma_x
    
    U,A,F,x = FDM_N2(xmin,xn,n,a,b,c,f,1,0)
    return U,x,Sigma_x

"""

"""
    def plot_sigma(self,ax):
        ax.plot(self.x/self.Rq,self.sigma/self.sigma0)
        ax.set_xlabel('x')
        ax.set_ylabel(r'Conductivity $\sigma(x)$ [S]')
        
    def plot_B(self,ax):
        ax.plot(self.x/self.Rq,self.B/self.B[0])
        ax.set_xlabel('x')
        ax.set_ylabel('Normalized B')
        
    def plot_Jz(self,ax):
        ax.plot(self.x/self.Rq,np.abs(self.Jz)/max(np.abs(self.Jz)))
        ax.set_xlabel('x')
        ax.set_ylabel('Normalized Jz')
            
    def plot_summary(self):
        self.fig = plt.figure(figsize = (12,9))
        self.ax1 = self.fig.add_subplot(221)
        self.ax2 = self.fig.add_subplot(222)
        self.ax3 = self.fig.add_subplot(223)
        
        self.plot_sigma(self.ax1)
        self.plot_B(self.ax2)
        self.plot_Jz(self.ax3)
        #self.ax3 = self.fig.add_subplot(223)
        #self.ax4 = self.fig.add_subplot(224)
        #self.ax5 = self.fig.add_subplot(235)
        return self.fig, [self.ax1,self.ax2,self.ax3]
"""
"""
    def plot_Z_surf(self,ax1,ax2):
        ax1.plot(self.f/1e9,self.Z_surf.real)
        ax2.plot(self.f/1e9,self.Z_surf.imag)
        ax1.set_xlabel('Frequency [GHz]')
        ax1.set_ylabel(r'Re{$Z_{surf}$}')
        ax2.set_ylabel(r'Im{$Z_{surf}$}')
        
    def plot_prop(self,ax1,ax2):
        ax1.plot(self.f/1e9,self.sigma_eff)
        ax2.plot(self.f/1e9,self.u_eff)
        ax1.set_xlabel('Frequency [GHz]')
        ax1.set_ylabel(r'Effective Conductivity')
        ax2.set_ylabel(r'Effective permeability')


        
        
    def plot_summary(self):
        self.fig = plt.figure()
        self.ax1 = self.fig.add_subplot(211)
        self.ax2 = self.ax1.twinx()
        self.plot_Z_surf(self.ax1,self.ax2)
        
        self.ax3 = self.fig.add_subplot(212)
        self.ax4 = self.ax3.twinx()
        self.plot_prop(self.ax3,self.ax4)
        
        return self.fig,[self.ax1,self.ax2,self.ax3,self.ax4]
    
    
    def tab_write(self,name,data):
        csv = np.array([self.f/1e9,data]).T
        np.savetxt(name + '.tab',csv, delimiter = '\t')
    
"""
