# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 18:34:23 2022

@author: hiro7
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy.integrate import quad
from matplotlib.animation import FuncAnimation
import os


pi,sin,cos,exp,sqrt = [np.pi,np.sin,np.cos,np.exp,np.sqrt]


path = r"C:\Users\hiro7\github\sana"
os.chdir(path)

def main():
    
    #f = lambda x:1 if 0.4<=x<=0.6 else 0
    #Amp = 0.1
    
    """
    f(x) = k(sinx - 1/2 sin2x)
    u_ana(x,t) = k(costsinx - 1/2 cos2tsin2x)
    """
    L = pi
    k = 1
    f1 = lambda x: k*(sin(x) - 0.5*sin(2*x))
    g1 = lambda x: 0
    modes = 10
    analytical1 = lambda x,t: k*(cos(t)*sin(x) - 0.5*cos(2*t)*sin(2*x))
        
    string1 = Wave_1D(0,L,100,0,5,20,f1,g1,1,modes,1,analytical1)
    
    #fig,ax = plt.subplots()
    #ax.plot(string2.x,string2.f_x)
    string1.animation("gif/string1")

    L = 1
    k = L/2
    c = 1
    modes = 20
    f2 = lambda x:2*k*x/L if 0<=x<=L/2 else 2*k/L*(L-x) if L/2 <=x<=L else 0  
    g2 = lambda x: 0
    #analytical2 = lambda x,t: 8*k/pi**2*(sin(pi*x/L)*cos(pi*c*t/L) - 1/3**2*sin(3*pi*x/L)*cos(3*pi*c*t/L))
    analytical2 = lambda x,t: 8*k/pi**2*(sin(pi*x/L)*cos(pi*c*t/L) - 1/3**2*sin(3*pi/L*x)*cos(3*pi*c/L*t) )

    string2 = Wave_1D(0,L,100,0,pi,30,f2,g2,c,modes,1,analytical2)
    string2.animation("gif/string2")
    
    
    
    f = lambda x:100*sin(pi*x/80)
    analytical2 = lambda x,t:100*sin(pi*x/80)*exp(-0.001785*t)
    rod_thermal = Diff_1D(0,80,100,0,700,20,f,sqrt(1.158),3,1,analytical2)
    rod_thermal.animation("gif/test2")
    
    return 0

def frame_update(frame,ax,x,y,analytical_flag = None, analysis = None):
    ax.cla() # ax をクリア
    ax.plot(x,y[frame,:],'red')
    #ax.plot(x,test_ana[frame],'blue',linestyle = '--')
    ax.set_ylim(-y.max(),y.max())
    #ax.text(0.5,0.5,f'{t[frame]}')

    if analytical_flag != None:
        ax.plot(x,analysis[frame],'blue',linestyle = '--')
        ax.legend(['numerical','analytical'],loc = 'upper right')



class Wave_1D(object):
    def __init__(self,x0,x1,nx,t0,t1,nt,f,g,c,modes,analytical_flag = None,analysis = None):
        self.x = np.linspace(x0,x1,nx)
        self.t = np.linspace(t0,t1,nt)
        self.f_x = np.array([f(i) for i in self.x])
        self.g_x = np.array([g(i) for i in self.x])
        self.modes = modes
        self.sol_all = np.array([[wave_1D(i_x,i_t,f,g,x1,1,self.modes) for i_x in self.x] for i_t in self.t]) #[time frame*x*modes]
        self.sol = self.sol_all.sum(axis=2)
        self.analytical_flag = analytical_flag
        if analytical_flag != None:
            self.analysis = np.array([[analysis(i,j) for i in self.x] for j in self.t])
        else:
            self.analysis = None
        
    #def plot(self,mode):
        
    def plot(self,mode):
        return 0
        
    def animation(self,name):        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(self.x,self.f_x)
        anim = FuncAnimation(fig, lambda x:frame_update(x,ax,self.x,self.sol,self.analytical_flag,self.analysis), frames=np.arange(self.t.shape[0]), interval=200)
        
        anim.save(name+".gif", writer="imagemagick")
        plt.close()
        return 0


#https://www.yutaka-note.com/entry/matplotlib_artist_anim#Matplotlib%E3%81%A7%E3%82%A2%E3%83%8B%E3%83%A1%E3%83%BC%E3%82%B7%E3%83%A7%E3%83%B3%E3%82%92%E4%BD%9C%E6%88%90%E3%81%99%E3%82%8B%EF%BC%92%E3%81%A4%E3%81%AE%E6%96%B9%E6%B3%95
def wave_1D(x,t,f,g,L,c,n): 
    """
    f(x),g(x) is a function
    c is a constant c**2 = T/roh T:power, roh: mass per length of the string
    """
    b1 = lambda n: 2/L*quad(lambda x:f(x)*sin(n*pi*x/L),0,L)[0]
    b2 = lambda n: 2/(c*n*pi)*quad(lambda x:g(x)*sin(n*pi*x/L),0,L)[0]
    #lambdan = c*n*pi/L
    u_sum = np.array([(b1(i)*cos(c*i*pi/L*t) + \
                       b2(i)*sin(c*i*pi/L*t))*sin(i*pi*x/L)\
                      for i in (np.arange(n)+1) ])
    return u_sum

class Diff_1D(object):
    def __init__(self,x0,x1,nx,t0,t1,nt,f,c,modes,analytical_flag = None,analysis = None):
        self.x = np.linspace(x0,x1,nx)
        self.t = np.linspace(t0,t1,nt)
        self.f_x = np.array([f(i) for i in self.x])
        self.modes = modes
        self.c = c
        self.sol_all = np.array([[diff_1D(i_x,i_t,f,x1,self.c,self.modes) for i_x in self.x] for i_t in self.t]) #[time frame*x*modes]
        self.sol = self.sol_all.sum(axis=2)
        self.analytical_flag = analytical_flag
        if analytical_flag != None:
            self.analysis = np.array([[analysis(i,j) for i in self.x] for j in self.t])
        else:
            self.analysis = None
        
    #def plot(self,mode):
        
    def plot(self,mode):
        return 0
        
    def animation(self,name):        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        anim = FuncAnimation(fig, lambda x:frame_update(x,ax,self.x,self.sol,self.analytical_flag,self.analysis), frames=np.arange(self.t.shape[0]), interval=200)
        
        anim.save(name+".gif", writer="imagemagick")
        plt.close()
        return 0
def diff_1D(x,t,f,L,c,modes): 
    """
    f(x),g(x) is a function
    c is a constant c**2 = T/roh T:power, roh: mass per length of the string
    """
    b1 = lambda n: 2/L*quad(lambda x:f(x)*sin(n*pi*x/L),0,L)[0]
    #b2 = lambda n: 2/(c*n*pi)*quad(lambda x:g(x)*sin(n*pi*x/L),0,L)[0]
    #lambdan = c*n*pi/L
    u_sum = np.array([b1(i)*sin(i*pi/L*x)*exp(-(c*i*pi/L)**2*t)\
                      for i in (np.arange(modes)+1) ])
    return u_sum


    


if __name__ == '__main__':    
    main()



    
def u_2D(x,y,t,f,L,c,modes): 
    """
    f(x),g(x) is a function
    c is a constant c**2 = T/roh T:power, roh: mass per length of the string
    """
    b1 = lambda n: 2/L*quad(lambda x:f(x)*sin(n*pi*x/L),0,L)[0]
    #b2 = lambda n: 2/(c*n*pi)*quad(lambda x:g(x)*sin(n*pi*x/L),0,L)[0]
    #lambdan = c*n*pi/L
    u_sum = np.array([b1(i)*sin(i*pi/L*x)*exp(-(c*i*pi/L)**2*t)\
                      for i in (np.arange(modes)+1) ])
    return u_sum
